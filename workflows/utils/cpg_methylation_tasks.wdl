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

    # size(CpGBed)/size(BinsBed) are the files' on-disk (possibly
    # gzip-compressed) sizes, but both get decompressed to plain copies
    # below, so budget for an inflated decompressed working set rather than
    # just a small multiple of that.
    Int auto_disk_size = ceil((size(CpGBed, "GB") + size(BinsBed, "GB")) * 8) + 10

    command <<<
        set -euo pipefail

        # CpGBed/BinsBed may each be gzip-compressed (e.g. CpGBed straight
        # from pb-CpG-tools, which emits .bed.gz; BinsBed from
        # ref/functional_regions.hg38.bed.gz) or plain text. The
        # image's zcat is BusyBox's, which -- unlike GNU gzip's -- refuses
        # to pass plain text through even with -f ("no gzip/bzip2/xz magic"),
        # so detect gzip ourselves via its magic bytes and only invoke zcat
        # when actually needed.
        maybe_decompress() {
            if [[ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" == "1f8b" ]]; then
                zcat "$1" > "$2"
            else
                cat "$1" > "$2"
            fi
        }
        maybe_decompress ~{CpGBed} cpg_bed_plain.bed
        maybe_decompress ~{BinsBed} bins_bed_plain.bed

        # Bins used for binning need a name in column 4 to carry an
        # identifier through to the output. If BinsBed is a plain 3-column
        # bed (chrom/start/end, e.g. bedtools makewindows output with no
        # name column), synthesize one from the interval itself, in the same
        # 1-based inclusive samtools/tabix style used elsewhere in these
        # workflows (e.g. IntervalString); otherwise keep whatever name is
        # already there (e.g. functional_regions.hg38.bed.gz's region_type
        # column).
        n_bins_cols=$(awk -F'\t' '!/^#/ {print NF; exit}' bins_bed_plain.bed)
        if [[ "$n_bins_cols" -eq 3 ]]; then
            awk -F'\t' -v OFS='\t' '!/^#/ {print $1, $2, $3, $1":"($2 + 1)"-"$3}' bins_bed_plain.bed > bins_named.bed
        else
            awk -F'\t' '!/^#/' bins_bed_plain.bed > bins_named.bed
        fi

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
            printf '#chrom\tstart\tend\tname\tnum_cpgs\ttotal_coverage\tweighted_mean_methylation\tunweighted_mean_methylation\n'
            bedtools map -a bins_named.bed -b cpg_with_weight.bed \
                    -c ~{MethCol},~{MethCol},"$cov_int_col","$weight_col" \
                    -o count,mean,sum,sum \
                | awk -F'\t' -v OFS='\t' '
                    $5 > 0 {
                        weighted_mean = ($7 > 0) ? sprintf("%.4f", $8 / $7) : "NA"
                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2, $3, $4, $5, $7, weighted_mean, $6
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
