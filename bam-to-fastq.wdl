version 1.0

workflow BamToFastq {
    input {
        File bam
        String sample_name
        String samtools_docker = "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
        Int cpu = 4
    }

    call SamtoolsFastq {
        input:
            bam = bam,
            sample_name = sample_name,
            docker = samtools_docker,
            cpu = cpu
    }

    output {
        File fastq_r1 = SamtoolsFastq.fq1
        File fastq_r2 = SamtoolsFastq.fq2
    }
}

task SamtoolsFastq {
    input {
        File bam
        String sample_name
        String docker
        Int cpu = 4
        Float mem_per_cpu_gb = 4
        Int boot_disk_gb = 10
    }

    Float total_mem_gb = mem_per_cpu_gb * cpu
    Int estimated_disk_gb = ceil(4 * size(bam, "GB") + 10)

    command <<<
        set -euo pipefail
        samtools sort -n -@ ~{cpu} ~{bam} \
            | samtools fastq \
                -@ ~{cpu} \
                -1 ~{sample_name}.R1.fastq.gz \
                -2 ~{sample_name}.R2.fastq.gz \
                -0 /dev/null -s /dev/null -n -c 7 \
        >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{total_mem_gb} GB"
        disks: "local-disk ~{estimated_disk_gb} HDD"
        bootDiskSizeGb: boot_disk_gb
    }

    output {
        File fq1 = "~{sample_name}.R1.fastq.gz"
        File fq2 = "~{sample_name}.R2.fastq.gz"
    }
}
