version 1.0

# Shared tasks for BuildBedMatrix: turn N per-sample bed files into a single
# wide bed matrix (chrom, start, end, feature, then one value column per
# sample) suitable as QTL-mapping phenotype input. Batched so large sample
# counts don't all get joined in a single task: the workflow splits the
# sample list into chunks (in pure WDL -- see its own comment on why),
# BuildMatrixBatch builds one matrix per chunk, and JoinMatrices combines
# the per-batch matrices into the final one. Tabular merges are delegated to
# single-purpose Polars scripts in the project data-manipulation image.

task BuildMatrixBatch {
    input {
        Array[File] BedFiles
        Array[String] SampleIDs
        Int FeatureCol
        Int ValueCol
        String OutputName = "matrix_batch.bed"
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
        sample_ids=(~{sep=' ' SampleIDs})
        n=${#bed_files[@]}

        : > bed_manifest.tsv
        for ((i = 0; i < n; i++)); do
            printf '%s\t%s\n' "${sample_ids[$i]}" "${bed_files[$i]}" >> bed_manifest.tsv
        done

        python3 /scripts/build_matrix_batch.py \
            --manifest bed_manifest.tsv \
            --feature-col ~{FeatureCol} \
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
        File MatrixBed = OutputName
    }
}

task JoinMatrices {
    input {
        Array[File] BatchMatrices
        String OutputPrefix = "phenotype_matrix"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(BatchMatrices, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        export LC_ALL=C

        matrix_files=(~{sep=' ' BatchMatrices})
        n=${#matrix_files[@]}

        : > matrix_manifest.tsv
        for ((i = 0; i < n; i++)); do
            printf '%s\n' "${matrix_files[$i]}" >> matrix_manifest.tsv
        done

        python3 /scripts/join_matrices.py \
            --manifest matrix_manifest.tsv \
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
        File MatrixBed = "~{OutputPrefix}.bed"
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
