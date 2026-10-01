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
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    }

    Boolean RestrictToIntervals = defined(IntervalBed) || defined(IntervalString)
    Boolean CpGBedIsGz = basename(CpGBed) != sub(basename(CpGBed), "\\.gz$", "")
    Boolean BinsBedIsGz = basename(BinsBed) != sub(basename(BinsBed), "\\.gz$", "")

    String FilteredSuffix = if RestrictToIntervals then "_filtered" else ""
    String BinnedOutputPrefix = "~{OutputPrefix}_binned_methyl~{FilteredSuffix}"

    if (CpGBedIsGz) {
        call CpGTasks.Unzip as UnzipCpGBed {
            input:
                InputFile = CpGBed
        }
    }

    if (BinsBedIsGz) {
        call CpGTasks.Unzip as UnzipBinsBed {
            input:
                InputFile = BinsBed
        }
    }

    File CpGBedPlain = select_first([UnzipCpGBed.DecompressedFile, CpGBed])
    File BinsBedPlain = select_first([UnzipBinsBed.DecompressedFile, BinsBed])

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

    call CpGTasks.BinCpGs as BinCpGs {
        input:
            CpGBed = CpGBedForBinning,
            BinsBed = BinsBedPlain,
            OutputPrefix = BinnedOutputPrefix,
            MethCol = MethCol,
            CovCol = CovCol,
            DockerImage = DockerImage
    }

    output {
        File BinnedMethylationBed = BinCpGs.Bed
    }
}
