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
    Boolean CpGBedIsGz = basename(CpGBed) != sub(basename(CpGBed), "\\.gz$", "")

    String FilteredSuffix = if RestrictToIntervals then "_filtered" else ""
    # Abbreviate round window sizes for naming (20000 -> 20k, 1000000 -> 1m);
    # anything not a whole multiple of 1000 falls back to the raw bp count.
    String WindowSizeLabel = if (WindowSize % 1000000 == 0) then "~{WindowSize / 1000000}m" else if (WindowSize % 1000 == 0) then "~{WindowSize / 1000}k" else "~{WindowSize}bp"
    String BinnedOutputPrefix = "~{OutputPrefix}_binned_methyl_~{WindowSizeLabel}~{FilteredSuffix}"

    if (CpGBedIsGz) {
        call CpGTasks.Unzip as UnzipCpGBed {
            input:
                InputFile = CpGBed
        }
    }

    File CpGBedPlain = select_first([UnzipCpGBed.DecompressedFile, CpGBed])

    if (RestrictToIntervals) {
        call CpGTasks.IntersectWithIntervals as IntersectWithIntervals {
            input:
                CpGBed = CpGBedPlain,
                IntervalBed = IntervalBed,
                IntervalString = IntervalString,
                DockerImage = DockerImage
        }
    }

    File CpGBedForBinning = select_first([IntersectWithIntervals.FilteredCpGBed, CpGBedPlain])

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
            OutputPrefix = BinnedOutputPrefix,
            MethCol = MethCol,
            CovCol = CovCol,
            DockerImage = DockerImage
    }

    output {
        File BinnedMethylationBed = BinCpGs.Bed
    }
}
