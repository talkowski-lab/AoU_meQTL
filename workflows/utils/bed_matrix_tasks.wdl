version 1.0

# Shared tasks for BuildBedMatrix: turn N per-sample bed files into a single
# wide bed matrix (chrom, start, end, feature, then one value column per
# sample) suitable as QTL-mapping phenotype input. Batched so large sample
# counts don't all get joined in a single task: the workflow splits the
# sample list into chunks (in pure WDL -- see its own comment on why),
# BuildMatrixBatch builds one matrix per chunk, and JoinMatrices combines
# the per-batch matrices into the final one.

task BuildMatrixBatch {
    input {
        Array[File] BedFiles
        Array[String] SampleIDs
        Int FeatureCol
        Int ValueCol
        String OutputName = "matrix_batch.bed"
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
        sample_ids=(~{sep=' ' SampleIDs})
        n=${#bed_files[@]}

        # Reduce each sample's bed to chrom:start:end (a single composite
        # join key, since `join` only supports a single field), feature, and
        # value -- the feature column is only kept from the first sample,
        # since `join` would otherwise just duplicate it for every sample.
        for ((i = 0; i < n; i++)); do
            if [[ "$i" -eq 0 ]]; then
                awk -F'\t' -v OFS='\t' -v f=~{FeatureCol} -v v=~{ValueCol} '
                    !/^#/ { print $1":"$2":"$3, $f, $v }
                ' "${bed_files[$i]}" | sort -k1,1 > "reduced_${i}.tsv"
            else
                awk -F'\t' -v OFS='\t' -v v=~{ValueCol} '
                    !/^#/ { print $1":"$2":"$3, $v }
                ' "${bed_files[$i]}" | sort -k1,1 > "reduced_${i}.tsv"
            fi
        done

        # Fold all samples together via repeated joins on the composite key,
        # accumulating one value column per sample as we go. This is an
        # inner join: a region missing from any one sample's bed (e.g.
        # dropped by BinCpGs for having zero coverage) is dropped from the
        # batch matrix entirely, rather than padded with NA.
        cp reduced_0.tsv acc.tsv
        for ((i = 1; i < n; i++)); do
            join -t $'\t' -j 1 acc.tsv "reduced_${i}.tsv" > acc_next.tsv
            mv acc_next.tsv acc.tsv
        done

        {
            printf '#chr\tstart\tend\tfeature'
            for sid in "${sample_ids[@]}"; do
                printf '\t%s' "$sid"
            done
            printf '\n'
            awk -F'\t' -v OFS='\t' '
                {
                    split($1, coord, ":")
                    line = coord[1] OFS coord[2] OFS coord[3]
                    for (c = 2; c <= NF; c++) line = line OFS $c
                    print line
                }
            ' acc.tsv | sort -k1,1 -k2,2n -k3,3n
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
        File MatrixBed = OutputName
    }
}

task JoinMatrices {
    input {
        Array[File] BatchMatrices
        String OutputPrefix = "phenotype_matrix"
        Boolean Bgzip = true
        # Only needs awk/sort/join/bgzip, not a full bioinformatics image --
        # this is the legacy standalone tabix/bgzip package (not htslib),
        # much lighter weight.
        String DockerImage = "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
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

        # Combine every batch matrix's header (its sample-ID columns, i.e.
        # everything after chr/start/end/feature) into one header, in file
        # order.
        {
            printf '#chr\tstart\tend\tfeature'
            for f in "${matrix_files[@]}"; do
                head -n1 "$f" | awk -F'\t' '{ for (i = 5; i <= NF; i++) printf "\t%s", $i }'
            done
            printf '\n'
        } > header.tsv

        # Same composite-key join as BuildMatrixBatch, just combining
        # already-built batch matrices (each with its own set of sample
        # columns) instead of raw per-sample beds -- again an inner join, so
        # a region missing from any batch is dropped from the final matrix.
        for ((i = 0; i < n; i++)); do
            if [[ "$i" -eq 0 ]]; then
                awk -F'\t' -v OFS='\t' '
                    !/^#/ {
                        line = $1":"$2":"$3 OFS $4
                        for (c = 5; c <= NF; c++) line = line OFS $c
                        print line
                    }
                ' "${matrix_files[$i]}" | sort -k1,1 > "reduced_${i}.tsv"
            else
                awk -F'\t' -v OFS='\t' '
                    !/^#/ {
                        line = $1":"$2":"$3
                        for (c = 5; c <= NF; c++) line = line OFS $c
                        print line
                    }
                ' "${matrix_files[$i]}" | sort -k1,1 > "reduced_${i}.tsv"
            fi
        done

        cp reduced_0.tsv acc.tsv
        for ((i = 1; i < n; i++)); do
            join -t $'\t' -j 1 acc.tsv "reduced_${i}.tsv" > acc_next.tsv
            mv acc_next.tsv acc.tsv
        done

        {
            cat header.tsv
            awk -F'\t' -v OFS='\t' '
                {
                    split($1, coord, ":")
                    line = coord[1] OFS coord[2] OFS coord[3]
                    for (c = 2; c <= NF; c++) line = line OFS $c
                    print line
                }
            ' acc.tsv | sort -k1,1 -k2,2n -k3,3n
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
        File MatrixBed = if Bgzip then "~{OutputPrefix}.bed.gz" else "~{OutputPrefix}.bed"
    }
}
