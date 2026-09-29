version 1.0

workflow VCFToPlink {
    input {
        File VCF
        Float MinAF = 0.01
        Float HWEPvalThreshold = 1e-6
        String OutputPrefix = "plink_out"
        String ImageTag = "latest"
    }

    call Plink2MakeBed {
        input:
            VCF = VCF,
            MinAF = MinAF,
            HWEPvalThreshold = HWEPvalThreshold,
            OutputPrefix = OutputPrefix,
            ImageTag = ImageTag
    }

    output {
        File Bed = Plink2MakeBed.Bed
        File Bim = Plink2MakeBed.Bim
        File Fam = Plink2MakeBed.Fam
        File Log = Plink2MakeBed.Log
    }
}

task Plink2MakeBed {
    input {
        File VCF
        Float MinAF = 0.01
        Float HWEPvalThreshold = 1e-6
        String OutputPrefix = "plink_out"
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(VCF, "GB") * 4) + 20

    command <<<
        set -euo pipefail

        plink2 \
            --vcf ~{VCF} \
            --double-id \
            --allow-extra-chr \
            --maf ~{MinAF} \
            --hwe ~{HWEPvalThreshold} \
            --make-bed \
            --out ~{OutputPrefix}
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
        File Bim = OutputPrefix + ".bim"
        File Fam = OutputPrefix + ".fam"
        File Log = OutputPrefix + ".log"
    }
}
