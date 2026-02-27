version 1.0

workflow deplete_host_reads_bowtie {
    input {
        File read1
        File read2
        File host_genome
        String sample_name
        String docker = "staphb/bowtie2:2.5.4"
        Int cpu = 8
    }

    call bowtie_depletion {
        input:
            read1 = read1,
            read2 = read2,
            host_genome = host_genome,
            sample_name = sample_name,
            cpu = cpu,
            docker = docker
    }

    output {
        File unmapped_r1 = bowtie_depletion.unmapped_read1
        File unmapped_r2 = bowtie_depletion.unmapped_read2
        Int read_counts_pre = bowtie_depletion.read_counts_pre
        Int read_counts_post = bowtie_depletion.read_counts_post
    }
}

task bowtie_depletion {
    input {
        File read1
        File read2
        File host_genome
        String sample_name
        String docker
        Int cpu = 8
        Float mem_per_cpu_gb = 2
        Int boot_disk_gb = 20
    }
    String template = "~{sample_name}_unmapped_%.fastq.gz"
    Float total_mem_gb = mem_per_cpu_gb * cpu
    Int estimated_disk_gb = ceil(2 * (size(host_genome, "GB") + size(read1, "GB") + size(read2, "GB")) + 10)

    command <<<
        set -euo pipefail
        zcat ~{read1} | wc -l | tee INITIAL_READS
        bowtie2-build --threads ~{cpu} ~{host_genome} host
        bowtie2 -x host \
            -1 ~{read1} \
            -2 ~{read2} \
            --un-conc-gz ~{template} \
            --very-sensitive-local \
            -p ~{cpu} \
            -S /dev/null
        zcat ~{sample_name}_unmapped_1.fastq.gz | wc -l | tee DEPLETED_READS
        >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{total_mem_gb} GB"
        disks: "local-disk ~{estimated_disk_gb} HDD"
        bootDiskSizeGb: boot_disk_gb
    }

    output {
        File unmapped_read1 = "~{sample_name}_unmapped_1.fastq.gz"
        File unmapped_read2 = "~{sample_name}_unmapped_2.fastq.gz"
        Int read_counts_pre = read_int("INITIAL_READS") / 2
        Int read_counts_post = read_int( "DEPLETED_READS" ) / 2
    }
}
