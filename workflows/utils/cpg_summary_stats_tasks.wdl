version 1.0

# Shared tasks for CpGSummaryStats: fold N per-sample bed files (the
# standard pb-CpG-tools combined pileup bed -- chrom, start, end,
# modification_probability, haplotype, coverage) into per-CpG mean and
# variance across samples, without ever joining samples into a wide matrix.
# The workflow splits the sample list into chunks (in pure WDL -- see its
# own comment on why); ComputeBatchStats reduces each chunk to one row per
# CpG (count/sum/sum-of-squares, the additive sufficient statistics for
# variance); CombineStats sums those same per-CpG statistics across batches
# and converts the totals to mean and sample variance. Every intermediate
# file has at most one row per CpG, regardless of how many samples feed in.
# A CpG's identity is just its chrom:start:end coordinate -- there's no
# separate feature-name column to track.

task ComputeBatchStats {
    input {
        Array[File] BedFiles
        # Defaults to column 4 (modification_probability) in the standard
        # pb-CpG-tools combined pileup bed; pass 6 (coverage) instead to get
        # variance of coverage rather than variance of methylation.
        Int ValueCol = 4
        String OutputName = "batch_stats.tsv"
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BedFiles, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        bed_files=(~{sep=' ' BedFiles})

        # Reduce every sample's bed to chrom:start:end (a single composite
        # sort key -- a CpG's identity, no separate feature-name column) and
        # a per-sample (n=1, sum=value, sum-of-squares=value^2) triple --
        # rows with a non-numeric value (e.g. "NA") are skipped rather than
        # corrupting the sums.
        for f in "${bed_files[@]}"; do
            awk -F'\t' -v OFS='\t' -v val=~{ValueCol} '
                !/^#/ && $val != "NA" {
                    print $1":"$2":"$3, 1, $val, $val * $val
                }
            ' "$f"
        done | sort -k1,1 > combined.tsv

        # Collapse to one row per CpG by summing n/sum/sumsq across whatever
        # samples in this batch had that CpG -- a CpG absent from a given
        # sample's bed just doesn't contribute, rather than requiring every
        # sample to share the same CpG set.
        awk -F'\t' -v OFS='\t' '
            function emit() {
                if (key != "") print key, n, sum, sumsq
            }
            {
                if ($1 != key) {
                    emit()
                    key = $1; n = 0; sum = 0; sumsq = 0
                }
                n += $2; sum += $3; sumsq += $4
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

task CombineStats {
    input {
        Array[File] BatchStats
        String OutputPrefix = "cpg_summary_stats"
        Boolean Bgzip = true
        # Only needs awk/sort/cat/bgzip, not a full bioinformatics image --
        # this is the legacy standalone tabix/bgzip package (not htslib),
        # same choice as BuildBedMatrix's JoinMatrices.
        String DockerImage = "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BatchStats, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        cat ~{sep=' ' BatchStats} | sort -k1,1 > combined.tsv

        # Same per-CpG reduction as ComputeBatchStats, just summing
        # already-batched (n, sum, sumsq) triples across batches instead of
        # per-sample ones -- then converts the final totals to mean and
        # sample variance (Bessel's correction, n-1): mean = sum/n,
        # variance = (sumsq - sum^2/n) / (n-1). CpGs seen in only one sample
        # (n=1) get a mean but no variance (NA, division by zero). The name
        # column is just the CpG's own coordinate (1-based inclusive,
        # samtools/tabix style, same convention BinCpGs uses), not a
        # separately-tracked value.
        awk -F'\t' -v OFS='\t' '
            function emit() {
                if (key == "") return
                split(key, coord, ":")
                name = coord[1]":"(coord[2] + 1)"-"coord[3]
                if (n > 1) {
                    mean = sum / n
                    variance = (sumsq - (sum * sum) / n) / (n - 1)
                    printf "%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\n", coord[1], coord[2], coord[3], name, n, mean, variance
                } else {
                    mean = (n == 1) ? sum / n : "NA"
                    printf "%s\t%s\t%s\t%s\t%d\t%s\tNA\n", coord[1], coord[2], coord[3], name, n, mean
                }
            }
            {
                if ($1 != key) {
                    emit()
                    key = $1; n = 0; sum = 0; sumsq = 0
                }
                n += $2; sum += $3; sumsq += $4
            }
            END { emit() }
        ' combined.tsv | sort -k1,1 -k2,2n -k3,3n > body.tsv

        {
            printf '#chrom\tstart\tend\tname\tn\tmean\tvariance\n'
            cat body.tsv
        } > "~{OutputPrefix}.bed"

        if [[ "~{Bgzip}" == "true" ]]; then
            bgzip "~{OutputPrefix}.bed"
        fi
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
        File StatsBed = if Bgzip then "~{OutputPrefix}.bed.gz" else "~{OutputPrefix}.bed"
    }
}
