version 1.0

# Shared tasks for CpGSummaryStats: fold N per-sample bed files (the
# standard pb-CpG-tools combined pileup bed -- chrom, start, end,
# modification_probability, haplotype, coverage) into per-CpG mean,
# variance, min, and max across samples, without ever joining samples into a
# wide matrix.
# The workflow splits the sample list into chunks (in pure WDL -- see its
# own comment on why); ComputeBatchStats reduces each chunk to one row per
# CpG (count/sum/sum-of-squares/min/max); CombineStats sums those same
# per-CpG statistics across batches and converts the totals to mean and
# sample variance while carrying min/max through. Every intermediate file
# has at most one row per CpG, regardless of how many samples feed in.
# A CpG's identity is its chrom/start/end BED coordinate. Tabular reductions
# are delegated to single-purpose Polars scripts in the project data-
# manipulation image.

task ComputeBatchStats {
    input {
        Array[File] BedFiles
        # Defaults to column 4 (modification_probability) in the standard
        # pb-CpG-tools combined pileup bed; pass 6 (coverage) instead to get
        # variance of coverage rather than variance of methylation.
        Int ValueCol = 4
        String OutputName = "batch_stats.tsv"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BedFiles, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        bed_files=(~{sep=' ' BedFiles})

        : > bed_manifest.tsv
        for f in "${bed_files[@]}"; do
            printf '%s\n' "$f" >> bed_manifest.tsv
        done

        python3 /scripts/compute_cpg_batch_stats.py \
            --manifest bed_manifest.tsv \
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

task CombineStats {
    input {
        Array[File] BatchStats
        String OutputPrefix = "cpg_summary_stats"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
        Int MemoryGB = 2
        Int CPU = 1
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

        python3 /scripts/combine_cpg_stats.py \
            --manifest stats_manifest.tsv \
            --output-prefix "~{OutputPrefix}"
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
        File StatsBed = "~{OutputPrefix}.bed"
    }
}

task Bgzip {
    input {
        File InputFile
        String DockerImage = "quay.io/biocontainers/htslib:1.22--h566b1c6_0"
    }

    Int disk_size = ceil(size(InputFile, "GB") * 2) + 10

    command <<<
        set -euo pipefail

        bgzip -c ~{InputFile} > ~{basename(InputFile)}.gz
    >>>

    runtime {
        docker: DockerImage
        memory: "4G"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File Output = basename(InputFile) + ".gz"
    }
}
