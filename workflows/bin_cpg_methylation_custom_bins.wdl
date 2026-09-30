version 1.0

import "utils/cpg_methylation_tasks.wdl" as CpGTasks

workflow BinCpGMethylationCustomBins {
    input {
        File CpGBed
        File BinsBed
        File? IntervalBed
        String? IntervalString
        String OutputPrefix = "cpg_bins"
        Int MethCol = 4
        Int CovCol = 6
        String ImageTag = "latest"
    }

    Boolean RestrictToIntervals = defined(IntervalBed) || defined(IntervalString)

    if (RestrictToIntervals) {
        call CpGTasks.IntersectWithIntervals as IntersectWithIntervals {
            input:
                CpGBed = CpGBed,
                IntervalBed = IntervalBed,
                IntervalString = IntervalString,
                ImageTag = ImageTag
        }
    }

    File CpGBedForBinning = select_first([IntersectWithIntervals.FilteredCpGBed, CpGBed])

    call CpGTasks.BinCpGs as BinCpGs {
        input:
            CpGBed = CpGBedForBinning,
            BinsBed = BinsBed,
            OutputPrefix = OutputPrefix,
            MethCol = MethCol,
            CovCol = CovCol,
            ImageTag = ImageTag
    }

    output {
        File BinnedMethylationBed = BinCpGs.Bed
    }
}
