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
        String ImageTag = "latest"
        Int MemoryGB = 4
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(CpGBed, "GB") * 3) + 10

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
        docker: "ayenkin1871/aou_meqtl-bioinformatics:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 2
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
        String ImageTag = "latest"
        Int MemoryGB = 2
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(ChromSizes, "GB") * 2) + 5

    command <<<
        set -euo pipefail

        bedtools makewindows -g ~{ChromSizes} -w ~{WindowSize} > ~{OutputName}
    >>>

    runtime {
        docker: "ayenkin1871/aou_meqtl-bioinformatics:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 1
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
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(CpGBed, "GB") * 3) + 10

    command <<<
        set -euo pipefail

        # bedtools map has no weighted-mean operation, so append an
        # integer-coverage column and a coverage*methylation column to the
        # CpG bed; map then sums both per bin, and weighted_mean_methylation
        # is their ratio (computed below, after mapping down to per-bin rows
        # instead of per-CpG rows).
        n_orig_cols=$(head -n1 ~{CpGBed} | awk -F'\t' '{print NF}')
        cov_int_col=$((n_orig_cols + 1))
        weight_col=$((n_orig_cols + 2))

        awk -F'\t' -v OFS='\t' -v m=~{MethCol} -v c=~{CovCol} '{
            cov_int = int($c)
            print $0, cov_int, $m * cov_int
        }' ~{CpGBed} > cpg_with_weight.bed

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
        docker: "ayenkin1871/aou_meqtl-bioinformatics:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 4
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File Bed = OutputPrefix + ".bed"
    }
}
