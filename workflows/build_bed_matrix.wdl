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
        String DockerImage = "debian:bookworm-slim"
    }

    call MatrixTasks.MakeBatches as MakeBatches {
        input:
            BedFiles = BedFiles,
            SampleIDs = SampleIDs,
            BatchSize = BatchSize,
            DockerImage = DockerImage
    }

    scatter (manifest in MakeBatches.BatchManifests) {
        Array[Array[String]] ManifestRows = read_tsv(manifest)
        Array[Array[String]] ManifestCols = transpose(ManifestRows)
        Array[File] BatchBedFiles = ManifestCols[0]
        Array[String] BatchSampleIDs = ManifestCols[1]

        call MatrixTasks.BuildMatrixBatch as BuildMatrixBatch {
            input:
                BedFiles = BatchBedFiles,
                SampleIDs = BatchSampleIDs,
                FeatureCol = FeatureCol,
                ValueCol = ValueCol,
                DockerImage = DockerImage
        }
    }

    call MatrixTasks.JoinMatrices as JoinMatrices {
        input:
            BatchMatrices = BuildMatrixBatch.MatrixBed,
            OutputPrefix = OutputPrefix,
            Bgzip = Bgzip
    }

    output {
        File PhenotypeMatrix = JoinMatrices.MatrixBed
    }
}
