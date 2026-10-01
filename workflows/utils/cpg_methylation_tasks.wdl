version 1.0

# Shared tasks for the BinCpGMethylation and BinCpGMethylationCustomBins
# workflows: optionally restrict a pb-CpG-tools bed to a region, build a
# fixed-size window bed from a chrom.sizes file, and aggregate CpG calls
# into whatever bins bed is supplied.

task IntersectWithIntervals {
    input {
        File CpGBed
        File? IntervalBed
        String? IntervalString
        String OutputName = "filtered_cpgs.bed"
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(CpGBed, "GB") * 2) + 5

    command <<<
        set -euo pipefail

        INTERVAL_BED="~{IntervalBed}"
        if [[ -z "$INTERVAL_BED" ]]; then
            INTERVAL_STRING="~{IntervalString}"
            if [[ "$INTERVAL_STRING" == *:* ]]; then
                # A 1-based inclusive region (samtools/tabix style, e.g.
                # chr1:1000000-2000000) -- convert to 0-based BED coordinates.
                printf '%s\n' "$INTERVAL_STRING" \
                    | awk -F'[:-]' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3}' \
                    > interval_from_string.bed
            else
                # A bare chromosome name (e.g. chr18) -- the whole chromosome.
                printf '%s\t0\t2147483647\n' "$INTERVAL_STRING" > interval_from_string.bed
            fi
            INTERVAL_BED="interval_from_string.bed"
        fi

        bedtools intersect -u -a ~{CpGBed} -b "$INTERVAL_BED" > ~{OutputName}
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
        File FilteredCpGBed = OutputName
    }
}

task MakeWindows {
    input {
        File ChromSizes
        Int WindowSize
        String OutputName = "windows.bed"
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
        Int MemoryGB = 1
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(ChromSizes, "GB") * 2) + 5

    command <<<
        set -euo pipefail

        bedtools makewindows -g ~{ChromSizes} -w ~{WindowSize} > ~{OutputName}
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
        File WindowsBed = OutputName
    }
}

task BinCpGs {
    input {
        File CpGBed
        File BinsBed
        String OutputPrefix = "cpg_bins"
        Int MethCol = 4
        Int CovCol = 6
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
        Int MemoryGB = 4
        Int CPU = 1
        Int? DiskGB
    }

    # size(CpGBed) is the file's on-disk (possibly gzip-compressed) size, but
    # it gets decompressed to a plain copy below, so budget for an inflated
    # decompressed working set rather than just a small multiple of that.
    Int auto_disk_size = ceil(size(CpGBed, "GB") * 8) + 10

    command <<<
        set -euo pipefail

        # CpGBed may be gzip-compressed (e.g. straight from pb-CpG-tools,
        # which emits .bed.gz) -- zcat -f auto-detects gzip and decompresses
        # it, while passing already-plain-text input through unchanged, so
        # this works either way without needing to know up front.
        zcat -f ~{CpGBed} > cpg_bed_plain.bed

        # bedtools map has no weighted-mean operation, so append an
        # integer-coverage column and a coverage*methylation column to the
        # CpG bed; map then sums both per bin, and weighted_mean_methylation
        # is their ratio (computed below, after mapping down to per-bin rows
        # instead of per-CpG rows).
        n_orig_cols=$(head -n1 cpg_bed_plain.bed | awk -F'\t' '{print NF}')
        cov_int_col=$((n_orig_cols + 1))
        weight_col=$((n_orig_cols + 2))

        awk -F'\t' -v OFS='\t' -v m=~{MethCol} -v c=~{CovCol} '{
            cov_int = int($c)
            print $0, cov_int, $m * cov_int
        }' cpg_bed_plain.bed > cpg_with_weight.bed

        {
            printf '#chrom\tstart\tend\tnum_cpgs\ttotal_coverage\tweighted_mean_methylation\tunweighted_mean_methylation\n'
            bedtools map -a ~{BinsBed} -b cpg_with_weight.bed \
                    -c ~{MethCol},~{MethCol},"$cov_int_col","$weight_col" \
                    -o count,mean,sum,sum \
                | awk -F'\t' -v OFS='\t' '
                    $4 > 0 {
                        weighted_mean = ($6 > 0) ? sprintf("%.4f", $7 / $6) : "NA"
                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2, $3, $4, $6, weighted_mean, $5
                    }' \
                | sort -k1,1 -k2,2n
        } > ~{OutputPrefix}.bed
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
        File Bed = OutputPrefix + ".bed"
    }
}
