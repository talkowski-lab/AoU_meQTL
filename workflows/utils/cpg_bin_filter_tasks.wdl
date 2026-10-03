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
# those statistics.

task ComputeBinStatsBatch {
    input {
        Array[File] BinnedBeds
        # 1-based column index into each bed for the value to summarize --
        # BuildFilteredCpGMatrix defaults this to 8 (unweighted_mean_methylation).
        Int ValueCol
        String OutputName = "batch_bin_stats.tsv"
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BinnedBeds, "GB") * 2) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        binned_beds=(~{sep=' ' BinnedBeds})

        # Reduce every sample's binned bed to chrom:start:end (a bin's
        # identity) and a per-sample (n=1, num_cpgs, value, value^2) triple
        # -- rows with a non-numeric value (e.g. weighted_mean_methylation
        # being "NA" for a bin with zero total coverage) are skipped rather
        # than corrupting the sums.
        for f in "${binned_beds[@]}"; do
            awk -F'\t' -v OFS='\t' -v val=~{ValueCol} '
                !/^#/ && $val != "NA" {
                    print $1":"$2":"$3, 1, $5, $val, $val * $val
                }
            ' "$f"
        done | sort -k1,1 > combined.tsv

        # Collapse to one row per bin by summing n/num_cpgs/value/value^2 (and
        # tracking min/max value) across whatever samples in this batch had
        # that bin -- a bin absent from a given sample's binned bed (e.g.
        # dropped by BinCpGs for having zero coverage there) just doesn't
        # contribute.
        awk -F'\t' -v OFS='\t' '
            function emit() {
                if (key != "") print key, n, sum_cpgs, sum_val, sumsq_val, min_val, max_val
            }
            {
                if ($1 != key) {
                    emit()
                    key = $1; n = 0; sum_cpgs = 0; sum_val = 0; sumsq_val = 0; min_val = ""; max_val = ""
                }
                n += $2; sum_cpgs += $3; sum_val += $4; sumsq_val += $5
                if (min_val == "" || $4 < min_val) min_val = $4
                if (max_val == "" || $4 > max_val) max_val = $4
            }
            END { emit() }
        ' combined.tsv > ~{OutputName}
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
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BatchStats, "GB") * 2) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        cat ~{sep=' ' BatchStats} | sort -k1,1 > combined.tsv

        # Same per-bin reduction as ComputeBinStatsBatch, just summing
        # already-batched (n, sum_cpgs, sum_val, sumsq_val) quadruples (and
        # taking the min/max of already-batched min/max values) across
        # batches instead of per-sample ones -- then converts the final
        # totals to mean/variance of value (Bessel's correction, n-1),
        # min/max/delta (max - min) of value, mean number of CpGs, and
        # presence. Bins seen in only one sample (n=1) get a mean/min/max/
        # delta but no variance (NA, division by zero).
        awk -F'\t' -v OFS='\t' -v total=~{TotalSamples} '
            function emit() {
                if (key == "") return
                split(key, coord, ":")
                name = coord[1]":"(coord[2] + 1)"-"coord[3]
                mean_numcpgs = sum_cpgs / n
                presence = n / total
                delta = max_val - min_val
                if (n > 1) {
                    mean_val = sum_val / n
                    variance_val = (sumsq_val - (sum_val * sum_val) / n) / (n - 1)
                    printf "%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", coord[1], coord[2], coord[3], name, n, mean_val, variance_val, min_val, max_val, delta, mean_numcpgs, presence
                } else {
                    mean_val = (n == 1) ? sum_val / n : "NA"
                    printf "%s\t%s\t%s\t%s\t%d\t%s\tNA\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", coord[1], coord[2], coord[3], name, n, mean_val, min_val, max_val, delta, mean_numcpgs, presence
                }
            }
            {
                if ($1 != key) {
                    emit()
                    key = $1; n = 0; sum_cpgs = 0; sum_val = 0; sumsq_val = 0; min_val = ""; max_val = ""
                }
                n += $2; sum_cpgs += $3; sum_val += $4; sumsq_val += $5
                if (min_val == "" || $6 < min_val) min_val = $6
                if (max_val == "" || $7 > max_val) max_val = $7
            }
            END { emit() }
        ' combined.tsv | sort -k1,1 -k2,2n -k3,3n > body.tsv

        {
            printf '#chrom\tstart\tend\tname\tn\tmean_value\tvariance_value\tmin_value\tmax_value\tdelta\tmean_num_cpgs\tpresence\n'
            cat body.tsv
        } > ~{OutputName}
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
