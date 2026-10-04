version 1.0

import "utils/bed_matrix_tasks.wdl" as MatrixTasks

# Builds a wide bed matrix (chrom, start, end, feature, then one value
# column per sample) from N per-sample bed files -- e.g. turning per-sample
# BinCpGMethylation output into a phenotype bed suitable as QTL-mapping
# input. Batched: samples are split into chunks of BatchSize, a matrix is
# built per chunk, then the chunks are joined into the final matrix.
workflow BuildBedMatrix {
    input {
        Array[File] BedFiles
        Array[String] SampleIDs
        Int FeatureCol
        Int ValueCol
        Int BatchSize = 50
        Boolean Bgzip = true
        String OutputPrefix = "phenotype_matrix"
        String DockerImage = "ayenkin1871/aou_meqtl-data-manipulation:latest"
    }

    Int NumSamples = length(BedFiles)
    Int NumBatches = (NumSamples + BatchSize - 1) / BatchSize

    # Batching is done here in pure WDL expressions -- not via a task that
    # writes BedFiles' paths into a manifest for a later task to read back
    # -- because a task only sees each File's *localized* path inside its
    # own sandbox, not a portable reference. Round-tripping that path
    # through a file and re-coercing it to File downstream works on a local
    # backend (shared filesystem) but breaks on Cromwell/Terra's GCP
    # backend, where each task localizes inputs into its own container: the
    # re-coerced "File" is just a bare local-looking path from a different
    # task's filesystem, which Cromwell can't resolve to the original GCS
    # object. Indexing BedFiles[i]/SampleIDs[i] directly, as below, always
    # keeps the original tracked File reference.
    scatter (b in range(NumBatches)) {
        Int BatchStart = b * BatchSize
        Int BatchEndRaw = BatchStart + BatchSize
        Int BatchEnd = if BatchEndRaw < NumSamples then BatchEndRaw else NumSamples

        scatter (i in range(NumSamples)) {
            if (i >= BatchStart && i < BatchEnd) {
                File BatchBedFile = BedFiles[i]
                String BatchSampleID = SampleIDs[i]
            }
        }

        call MatrixTasks.BuildMatrixBatch as BuildMatrixBatch {
            input:
                BedFiles = select_all(BatchBedFile),
                SampleIDs = select_all(BatchSampleID),
                FeatureCol = FeatureCol,
                ValueCol = ValueCol,
                DockerImage = DockerImage
        }
    }

    call MatrixTasks.JoinMatrices as JoinMatrices {
        input:
            BatchMatrices = BuildMatrixBatch.MatrixBed,
            OutputPrefix = OutputPrefix,
            DockerImage = DockerImage
    }

    if (Bgzip) {
        call MatrixTasks.Bgzip as Bgzip_task {
            input:
                InputFile = JoinMatrices.MatrixBed
        }
    }

    output {
        File PhenotypeMatrix = select_first([Bgzip_task.Output, JoinMatrices.MatrixBed])
    }
}
