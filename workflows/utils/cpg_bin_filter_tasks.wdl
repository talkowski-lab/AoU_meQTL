version 1.0

# Shared tasks for BuildFilteredCpGMatrix: given N per-sample BinCpGs-shaped
# beds (chrom, start, end, name, num_cpgs, total_coverage,
# weighted_mean_methylation, unweighted_mean_methylation), compute per-bin
# summary statistics across samples -- mean/variance/min/max of a value
# column (plus delta, their difference), mean number of CpGs backing that
# value, and presence (the fraction of samples the bin was even called in)
# -- the same running-total approach as CpGSummaryStats (count/sum/sum-of-
# squares/min/max, folded in two passes so no stage ever joins samples into
# a wide matrix), then label each bin PASS/FAIL against QC thresholds on
# those statistics. Tabular reductions are delegated to single-purpose Polars
# scripts in the project data-manipulation image.

task ComputeBinStatsBatch {
    input {
        Array[File] BinnedBeds
        # 1-based column index into each bed for the value to summarize --
        # BuildFilteredCpGMatrix defaults this to 8 (unweighted_mean_methylation).
        Int ValueCol
        String OutputName = "batch_bin_stats.tsv"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
        Int MemoryGB = 8
        Int CPU = 2
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BinnedBeds, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        binned_beds=(~{sep=' ' BinnedBeds})

        : > binned_manifest.tsv
        for f in "${binned_beds[@]}"; do
            printf '%s\n' "$f" >> binned_manifest.tsv
        done

        python3 /scripts/compute_bin_batch_stats.py \
            --manifest binned_manifest.tsv \
            --value-col ~{ValueCol} \
            --output "~{OutputName}"
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File BatchStats = OutputName
    }
}

task CombineBinStats {
    input {
        Array[File] BatchStats
        # Total number of samples fed into the workflow -- the denominator
        # for presence (n / TotalSamples), not just however many samples
        # happened to have a given bin.
        Int TotalSamples
        String OutputName = "bin_stats.bed"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
        Int MemoryGB = 8
        Int CPU = 2
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BatchStats, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        batch_stats=(~{sep=' ' BatchStats})
        : > stats_manifest.tsv
        for f in "${batch_stats[@]}"; do
            printf '%s\n' "$f" >> stats_manifest.tsv
        done

        python3 /scripts/combine_bin_stats.py \
            --manifest stats_manifest.tsv \
            --total-samples ~{TotalSamples} \
            --output "~{OutputName}"
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File StatsBed = OutputName
    }
}

task FilterBins {
    input {
        File StatsBed
        Float MinPresence
        Float MinMeanCpGs
        Float MinVariance
        Float MinDelta
        String OutputPrefix = "filtered_cpg_matrix"
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 1
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(StatsBed, "GB") * 2) + 5

    command <<<
        set -euo pipefail
        export LC_ALL=C

        # A bin PASSes only if it has enough variance and enough spread
        # (delta = max - min) to possibly show an association, enough CpGs
        # backing its value on average to be a reliable estimate, and is
        # present (called) in enough samples. Bins with fewer than 2 samples
        # have no variance (NA) and so automatically FAIL -- $7 != "NA" makes
        # that explicit, even though awk's numeric comparison on "NA"
        # (coerced to 0) would already fail the >= minvar check on its own.
        {
            head -n1 ~{StatsBed} | awk -F'\t' -v OFS='\t' '{ print $0, "filter" }'
            awk -F'\t' -v OFS='\t' -v minp=~{MinPresence} -v mincpg=~{MinMeanCpGs} -v minvar=~{MinVariance} -v mindelta=~{MinDelta} '
                !/^#/ {
                    pass = ($12 >= minp) && ($11 >= mincpg) && ($7 != "NA") && ($7 >= minvar) && ($10 >= mindelta)
                    print $0, (pass ? "PASS" : "FAIL")
                }
            ' ~{StatsBed}
        } > "~{OutputPrefix}_summary_stats.bed"

        awk -F'\t' -v OFS='\t' '!/^#/ && $13 == "PASS" { print $1, $2, $3 }' "~{OutputPrefix}_summary_stats.bed" > "~{OutputPrefix}_passing_bins.bed"
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File FilteredStatsBed = "~{OutputPrefix}_summary_stats.bed"
        File PassingBinsBed = "~{OutputPrefix}_passing_bins.bed"
    }
}

task FilterBedToPassingBins {
    input {
        File BinnedBed
        File PassingBinsBed
        String OutputName = "filtered_bins.bed"
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 1
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil((size(BinnedBed, "GB") + size(PassingBinsBed, "GB")) * 2) + 5

    command <<<
        set -euo pipefail
        export LC_ALL=C

        # Classic two-file awk filter: load passing bins' chrom:start:end
        # keys into memory first (NR==FNR, i.e. while reading the first
        # file), then keep only BinnedBed rows whose key is in that set.
        awk -F'\t' -v OFS='\t' '
            NR == FNR { keep[$1":"$2":"$3] = 1; next }
            !/^#/ && ($1":"$2":"$3) in keep
        ' ~{PassingBinsBed} ~{BinnedBed} > ~{OutputName}
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File FilteredBed = OutputName
    }
}
