version 1.0

workflow VCFToPlink {
    input {
        File VCF
        String InputFormat = "vcf"
        Float MinAF = 0.01
        Float HWEPvalThreshold = 0.000001
        String OutputPrefix = "plink_out"
        String ImageTag = "latest"
    }

    call Plink2MakePgen {
        input:
            VCF = VCF,
            InputFormat = InputFormat,
            MinAF = MinAF,
            HWEPvalThreshold = HWEPvalThreshold,
            OutputPrefix = OutputPrefix,
            ImageTag = ImageTag
    }

    output {
        File Pgen = Plink2MakePgen.Pgen
        File Pvar = Plink2MakePgen.Pvar
        File Psam = Plink2MakePgen.Psam
        File Log = Plink2MakePgen.Log
    }
}

task Plink2MakePgen {
    input {
        File VCF
        String InputFormat = "vcf"
        Float MinAF = 0.01
        Float HWEPvalThreshold = 0.000001
        String OutputPrefix = "plink_out"
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(VCF, "GB") * 4) + 20

    command <<<
        set -euo pipefail

        case "~{InputFormat}" in
            vcf|bcf) ;;
            *)
                echo "InputFormat must be 'vcf' or 'bcf', got '~{InputFormat}'" >&2
                exit 1
                ;;
        esac

        plink2 \
            --~{InputFormat} ~{VCF} \
            --double-id \
            --allow-extra-chr \
            --maf ~{MinAF} \
            --hwe ~{HWEPvalThreshold} \
            --make-pgen \
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
        File Pgen = OutputPrefix + ".pgen"
        File Pvar = OutputPrefix + ".pvar"
        File Psam = OutputPrefix + ".psam"
        File Log = OutputPrefix + ".log"
    }
}
