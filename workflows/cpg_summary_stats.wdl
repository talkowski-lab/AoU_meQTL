version 1.0

import "utils/cpg_summary_stats_tasks.wdl" as StatsTasks

# Computes, per CpG, the mean and variance of a value column across N
# per-sample pb-CpG-tools combined pileup beds (chrom, start, end,
# modification_probability, haplotype, coverage) -- without ever joining the
# samples into a wide matrix (unlike BuildBedMatrix). Instead, samples are
# folded down to the additive sufficient statistics for variance (count,
# sum, sum-of-squares) in two reduction passes: batched (ComputeBatchStats)
# then across batches (CombineStats). Every intermediate file has one row
# per CpG, not one column per sample, so memory/disk use doesn't grow with
# sample count the way a full matrix would.
workflow CpGSummaryStats {
    input {
        Array[File] BedFiles
        # Defaults to column 4 (modification_probability); pass 6 (coverage)
        # instead to get variance of coverage rather than methylation.
        Int ValueCol = 4
        Int BatchSize = 50
        Boolean Bgzip = true
        String OutputPrefix = "cpg_summary_stats"
        String DockerImage = "debian:bookworm-slim"
    }

    Int NumSamples = length(BedFiles)
    Int NumBatches = (NumSamples + BatchSize - 1) / BatchSize

    # Batching is done here in pure WDL expressions, not via a task that
    # writes BedFiles' paths into a manifest for a later task to read back
    # -- see BuildBedMatrix's build_bed_matrix.wdl for why that breaks on
    # Cromwell/Terra's GCP backend (a task only sees each File's *localized*
    # path inside its own sandbox, not a portable reference). Indexing
    # BedFiles[i] directly, as below, always keeps the original tracked
    # File reference.
    scatter (b in range(NumBatches)) {
        Int BatchStart = b * BatchSize
        Int BatchEndRaw = BatchStart + BatchSize
        Int BatchEnd = if BatchEndRaw < NumSamples then BatchEndRaw else NumSamples

        scatter (i in range(NumSamples)) {
            if (i >= BatchStart && i < BatchEnd) {
                File BatchBedFile = BedFiles[i]
            }
        }

        call StatsTasks.ComputeBatchStats as ComputeBatchStats {
            input:
                BedFiles = select_all(BatchBedFile),
                ValueCol = ValueCol,
                DockerImage = DockerImage
        }
    }

    call StatsTasks.CombineStats as CombineStats {
        input:
            BatchStats = ComputeBatchStats.BatchStats,
            OutputPrefix = OutputPrefix,
            Bgzip = Bgzip
    }

    output {
        File SummaryStatsBed = CombineStats.StatsBed
    }
}
