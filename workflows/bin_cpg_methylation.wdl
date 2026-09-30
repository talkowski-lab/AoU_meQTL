version 1.0

workflow BinCpGMethylation {
    input {
        File CpGBed
        File ChromSizes
        Int WindowSize
        File? IntervalBed
        String? IntervalString
        String OutputPrefix = "cpg_bins"
        Int MethCol = 4
        Int CovCol = 6
        String ImageTag = "latest"
    }

    Boolean RestrictToIntervals = defined(IntervalBed) || defined(IntervalString)

    if (RestrictToIntervals) {
        call IntersectWithIntervals {
            input:
                CpGBed = CpGBed,
                IntervalBed = IntervalBed,
                IntervalString = IntervalString,
                ImageTag = ImageTag
        }
    }

    File CpGBedForBinning = select_first([IntersectWithIntervals.FilteredCpGBed, CpGBed])

    call BinCpGs {
        input:
            CpGBed = CpGBedForBinning,
            ChromSizes = ChromSizes,
            WindowSize = WindowSize,
            OutputPrefix = OutputPrefix,
            MethCol = MethCol,
            CovCol = CovCol,
            ImageTag = ImageTag
    }

    output {
        File BinnedMethylationBed = BinCpGs.Bed
    }
}

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
            # IntervalString is a 1-based inclusive region (samtools/tabix style,
            # e.g. chr1:1000000-2000000) -- convert to 0-based BED coordinates.
            printf '%s\n' "~{IntervalString}" \
                | awk -F'[:-]' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3}' \
                > interval_from_string.bed
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

task BinCpGs {
    input {
        File CpGBed
        File ChromSizes
        Int WindowSize
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

        bedtools makewindows -g ~{ChromSizes} -w ~{WindowSize} > windows.bed

        {
            printf '#chrom\tstart\tend\tnum_cpgs\ttotal_coverage\tweighted_mean_methylation\tunweighted_mean_methylation\n'
            bedtools intersect -wa -wb -a windows.bed -b ~{CpGBed} \
                | python3 /scripts/bin_cpg_methylation.py --meth-col ~{MethCol} --cov-col ~{CovCol} \
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
