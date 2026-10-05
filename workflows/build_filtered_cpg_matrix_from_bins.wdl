version 1.0

import "utils/cpg_bin_filter_tasks.wdl" as FilterTasks
import "utils/bed_matrix_tasks.wdl" as MatrixTasks

# Same QC/filter/matrix tail as BuildFilteredCpGMatrix, but starts from
# already-binned per-sample BinCpGs-shaped beds:
#   chrom, start, end, name, num_cpgs, total_coverage,
#   weighted_mean_methylation, unweighted_mean_methylation
workflow BuildFilteredCpGMatrix_FromBins {
    input {
        Array[File] BinnedBeds
        Array[String] SampleIDs

        # 1-based column index into each binned bed for the value
        # summarized/filtered on and carried into the final matrix: 8 is
        # unweighted_mean_methylation; pass 7 for weighted_mean_methylation,
        # or 5/6 to summarize num_cpgs/total_coverage instead.
        Int ValueCol = 8

        Int BatchSize = 50

        Float MinPresence = 0.8
        Float MinMeanCpGs = 3.0
        Float MinVariance = 0.001
        Float MinDelta = 0.1

        Boolean Bgzip = true
        String OutputPrefix = "filtered_binned_cpg_matrix"
        String BasicDockerImage = "debian:bookworm-slim"
        String DataDockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
    }

    Int NumSamples = length(BinnedBeds)

    # --- Per-bin summary stats across samples, batched. Batching is pure WDL
    # indexing to preserve File identity across Cromwell/Terra task
    # localization; don't round-trip localized paths through a manifest task.
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
                DockerImage = DataDockerImage
        }
    }

    call FilterTasks.CombineBinStats as CombineBinStats {
        input:
            BatchStats = ComputeBinStatsBatch.BatchStats,
            TotalSamples = NumSamples,
            DockerImage = DataDockerImage
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

    # --- Build the final matrix from the filtered per-sample beds. ---
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
                DockerImage = DataDockerImage
        }
    }

    call MatrixTasks.JoinMatrices as JoinMatrices {
        input:
            BatchMatrices = BuildMatrixBatch.MatrixBed,
            OutputPrefix = OutputPrefix,
            DockerImage = DataDockerImage
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
