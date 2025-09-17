version 1.0

workflow downsample_fastqs {
    input {
        File fastq1
        File? fastq2
        String sample_name
        Float sampling_fraction

        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 4
    }
    call downsample_fastqs {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            sample_name = sample_name,
            sampling_fraction = sampling_fraction,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size
    }

    output {
        File    downsampled_read1 = downsample_fastqs.read1
        File?   downsampled_read2 = downsample_fastqs.read2
    }
}

task downsample_fastqs {
    input {
        File fastq1
        File? fastq2
        String sample_name
        Float sampling_fraction
        Int cpu = 4
        Int memory = 16
        String docker = "staphb/seqtk:1.5"
        Int disk_size = 16
    }

    Boolean is_paired_end = defined(fastq2)

    command <<<
        echo "Subsampling reads for ~{sample_name}."
        seqtk sample -s100 ~{fastq1} ~{sampling_fraction} | gzip > ~{sample_name}.R1.fastq.gz
        echo "Created ~{sample_name}.R1.fastq.gz."

        if [ ~{is_paired_end} = true ] ; then
          echo "Subsampling second reads for ~{sample_name}."
          seqtk sample -s100 ~{fastq2} ~{sampling_fraction} | gzip > ~{sample_name}.R2.fastq.gz
          echo "Created ~{sample_name}.R2.fastq.gz."
        fi
    >>>
    output {
        File downsampled_read1 = "~{sample_name}.R1.fastq.gz"
        File? downsampled_read2 = "~{sample_name}.R2.fastq.gz"
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
    }
}

