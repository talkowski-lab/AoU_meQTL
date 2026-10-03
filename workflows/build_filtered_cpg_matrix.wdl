version 1.0

import "utils/cpg_methylation_tasks.wdl" as CpGTasks
import "utils/cpg_bin_filter_tasks.wdl" as FilterTasks
import "utils/bed_matrix_tasks.wdl" as MatrixTasks

# End-to-end pipeline from raw per-sample pb-CpG-tools beds to a QC-filtered
# QTL-mapping phenotype matrix:
#   1. Bin every sample's CpGs (BinCpGs, reused from BinCpGMethylation*) into
#      either fixed-size windows (ChromSizes + WindowSize) or a custom bins
#      bed (BinsBed) -- BinsBed is preferred when both are given.
#   2. Fold the per-sample binned beds down to per-bin summary statistics
#      across samples -- mean/variance of a value column, mean number of
#      CpGs backing it, and presence (fraction of samples the bin was even
#      called in) -- via the same batched running-total reduction as
#      CpGSummaryStats, so this never joins samples into a wide matrix.
#   3. Label every bin PASS/FAIL against QC thresholds on those three
#      statistics, and restrict each sample's binned bed to only the PASS
#      bins before matrix-building, so BuildBedMatrix's per-sample joins
#      never have to consider bins that would be dropped anyway.
#   4. Build the final phenotype matrix (BuildMatrixBatch/JoinMatrices,
#      reused from BuildBedMatrix) from just those filtered beds.
workflow BuildFilteredCpGMatrix {
    input {
        Array[File] CpGBeds
        Array[String] SampleIDs

        # Bins: BinsBed is used if given (custom regions, e.g. gene bodies,
        # ref/functional_regions.hg38.bed.gz); otherwise ChromSizes and
        # WindowSize are both required to build fixed-size windows.
        File? BinsBed
        File? ChromSizes
        Int? WindowSize

        String? IntervalString 

        Int MethCol = 4
        Int CovCol = 6
        # 1-based column index into each BinCpGs output bed for the value
        # summarized/filtered on and carried into the final matrix: 8 is
        # unweighted_mean_methylation; pass 7 for weighted_mean_methylation,
        # or 5/6 to summarize num_cpgs/total_coverage instead.
        Int ValueCol = 8

        Int BatchSize = 50

        # QC filter thresholds -- starting points from common region-level
        # methylation/mQTL phenotype filtering practice, not hard standards:
        # tune per cohort/tissue. MinPresence/MinMeanCpGs guard against
        # unreliable bins (too few samples called, too few CpGs backing the
        # average); MinVariance and MinDelta both guard against bins with no
        # meaningful variability left to test for association -- MinVariance
        # on the overall spread (sensitive to a few outlier samples),
        # MinDelta (max - min across samples) on the raw range (a commonly
        # cited minimum "meaningful" methylation difference in the DNAm
        # literature is ~0.05-0.1 on a 0-1 scale). All four are plain
        # decimals, not scientific notation -- Dockstore's WDL parser
        # rejects "1e-4"-style float literals.
        Float MinPresence = 0.8
        Float MinMeanCpGs = 3.0
        Float MinVariance = 0.001
        Float MinDelta = 0.1

        Boolean Bgzip = true
        String OutputPrefix = "filtered_cpg_matrix"
        String BasicDockerImage = "debian:bookworm-slim"
        String BedtoolsDockerImage = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    }

    Int NumSamples = length(CpGBeds)
    Boolean UseCustomBins = defined(BinsBed)
    Boolean RestrictToIntervals = defined(IntervalString)

    # --- Resolve bins: prefer BinsBed if given, else build fixed windows ---
    if (UseCustomBins) {
        File BinsBedInput = select_first([BinsBed])
        Boolean BinsBedIsGz = basename(BinsBedInput) != sub(basename(BinsBedInput), "\\.gz$", "")

        if (BinsBedIsGz) {
            call CpGTasks.Unzip as UnzipBinsBed {
                input:
                    InputFile = BinsBedInput,
                    DockerImage = BasicDockerImage
            }
        }

        File CustomBinsBedPlain = select_first([UnzipBinsBed.DecompressedFile, BinsBedInput])
    }

    if (!UseCustomBins) {
        call CpGTasks.MakeWindows as MakeWindows {
            input:
                ChromSizes = select_first([ChromSizes]),
                WindowSize = select_first([WindowSize]),
                DockerImage = BedtoolsDockerImage
        }
    }

    File BinsBedFinal = select_first([CustomBinsBedPlain, MakeWindows.WindowsBed])

    # --- Per-sample binning ---
    scatter (i in range(NumSamples)) {
        File RawCpGBed = CpGBeds[i]
        Boolean CpGBedIsGz = basename(RawCpGBed) != sub(basename(RawCpGBed), "\\.gz$", "")

        if (CpGBedIsGz) {
            call CpGTasks.Unzip as UnzipCpGBed {
                input:
                    InputFile = RawCpGBed,
                    DockerImage = BasicDockerImage
            }
        }

        File CpGBedPlain = select_first([UnzipCpGBed.DecompressedFile, RawCpGBed])

        if (RestrictToIntervals) {
            call CpGTasks.IntersectWithIntervals as IntersectWithIntervals {
                input:
                    CpGBed = CpGBedPlain,
                    IntervalString = IntervalString,
                    DockerImage = BedtoolsDockerImage
            }
        }

        File CpGBedForBinning = select_first([IntersectWithIntervals.FilteredCpGBed, CpGBedPlain])

        call CpGTasks.BinCpGs as BinCpGs {
            input:
                CpGBed = CpGBedForBinning,
                BinsBed = BinsBedFinal,
                MethCol = MethCol,
                CovCol = CovCol,
                DockerImage = BedtoolsDockerImage
        }
    }

    Array[File] BinnedBeds = BinCpGs.Bed

    # --- Per-bin summary stats across samples, batched (see BuildBedMatrix's
    # own comment on why this is pure WDL rather than a manifest-writing
    # task: a task only sees each File's *localized* path, which breaks on
    # Cromwell/Terra's GCP backend if round-tripped through a file). ---
    Int NumStatsBatches = (NumSamples + BatchSize - 1) / BatchSize

    scatter (b in range(NumStatsBatches)) {
        Int StatsBatchStart = b * BatchSize
        Int StatsBatchEndRaw = StatsBatchStart + BatchSize
        Int StatsBatchEnd = if StatsBatchEndRaw < NumSamples then StatsBatchEndRaw else NumSamples

        scatter (i in range(NumSamples)) {
            if (i >= StatsBatchStart && i < StatsBatchEnd) {
                File StatsBatchBed = BinnedBeds[i]
            }
        }

        call FilterTasks.ComputeBinStatsBatch as ComputeBinStatsBatch {
            input:
                BinnedBeds = select_all(StatsBatchBed),
                ValueCol = ValueCol,
                DockerImage = BasicDockerImage
        }
    }

    call FilterTasks.CombineBinStats as CombineBinStats {
        input:
            BatchStats = ComputeBinStatsBatch.BatchStats,
            TotalSamples = NumSamples,
            DockerImage = BasicDockerImage
    }

    call FilterTasks.FilterBins as FilterBins {
        input:
            StatsBed = CombineBinStats.StatsBed,
            MinPresence = MinPresence,
            MinMeanCpGs = MinMeanCpGs,
            MinVariance = MinVariance,
            MinDelta = MinDelta,
            OutputPrefix = OutputPrefix,
            DockerImage = BasicDockerImage
    }

    # --- Restrict each sample's binned bed to passing bins before matrix
    # build, so the matrix join below never has to consider dropped bins. ---
    scatter (i in range(NumSamples)) {
        call FilterTasks.FilterBedToPassingBins as FilterBedToPassingBins {
            input:
                BinnedBed = BinnedBeds[i],
                PassingBinsBed = FilterBins.PassingBinsBed,
                DockerImage = BasicDockerImage
        }
    }

    Array[File] FilteredBinnedBeds = FilterBedToPassingBins.FilteredBed

    # --- Build the final matrix from the filtered per-sample beds, same
    # batched join as BuildBedMatrix. ---
    Int NumMatrixBatches = (NumSamples + BatchSize - 1) / BatchSize

    scatter (b in range(NumMatrixBatches)) {
        Int MatrixBatchStart = b * BatchSize
        Int MatrixBatchEndRaw = MatrixBatchStart + BatchSize
        Int MatrixBatchEnd = if MatrixBatchEndRaw < NumSamples then MatrixBatchEndRaw else NumSamples

        scatter (i in range(NumSamples)) {
            if (i >= MatrixBatchStart && i < MatrixBatchEnd) {
                File MatrixBatchBed = FilteredBinnedBeds[i]
                String MatrixBatchSampleID = SampleIDs[i]
            }
        }

        call MatrixTasks.BuildMatrixBatch as BuildMatrixBatch {
            input:
                BedFiles = select_all(MatrixBatchBed),
                SampleIDs = select_all(MatrixBatchSampleID),
                FeatureCol = 4,
                ValueCol = ValueCol,
                DockerImage = BasicDockerImage
        }
    }

    call MatrixTasks.JoinMatrices as JoinMatrices {
        input:
            BatchMatrices = BuildMatrixBatch.MatrixBed,
            OutputPrefix = OutputPrefix,
            DockerImage = BasicDockerImage
    }

    if (Bgzip) {
      call MatrixTasks.Bgzip as Bgzip_task {
          input:
              InputFile = JoinMatrices.MatrixBed
      }
    }

    output {
        File Matrix = select_first([Bgzip_task.Output, JoinMatrices.MatrixBed])
        File SummaryStats = FilterBins.FilteredStatsBed
    }
}
