version 1.0

import "utils/cpg_methylation_tasks.wdl" as CpGTasks

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
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    }

    Boolean RestrictToIntervals = defined(IntervalBed) || defined(IntervalString)

    if (RestrictToIntervals) {
        call CpGTasks.IntersectWithIntervals as IntersectWithIntervals {
            input:
                CpGBed = CpGBed,
                IntervalBed = IntervalBed,
                IntervalString = IntervalString,
                DockerImage = DockerImage
        }
    }

    File CpGBedForBinning = select_first([IntersectWithIntervals.FilteredCpGBed, CpGBed])

    call CpGTasks.MakeWindows as MakeWindows {
        input:
            ChromSizes = ChromSizes,
            WindowSize = WindowSize,
            DockerImage = DockerImage
    }

    call CpGTasks.BinCpGs as BinCpGs {
        input:
            CpGBed = CpGBedForBinning,
            BinsBed = MakeWindows.WindowsBed,
            OutputPrefix = OutputPrefix,
            MethCol = MethCol,
            CovCol = CovCol,
            DockerImage = DockerImage
    }

    output {
        File BinnedMethylationBed = BinCpGs.Bed
    }
}
