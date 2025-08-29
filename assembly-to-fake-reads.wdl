version 1.0

workflow ncbi_to_fake_reads {
    input {
        String ncbi_accession
        Int read_length = 150
        Int required_coverage = 30
        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 4
        Int cpu = 1
    }
    call ncbi_datasets_download_genome_accession {
        input:
            ncbi_accession = ncbi_accession,
            disk_size = disk_size,
            memory = memory,
            cpu = cpu,
    }
    call bbmap_randomreads {
        input:
            ncbi_accession = ncbi_accession,
            assembly = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta,
            read_length = read_length,
            required_coverage = required_coverage,
            disk_size = disk_size,
            memory = memory,
            cpu = cpu,
    }

    output {
        File    assembly_read1 = bbmap_randomreads.read1
        File    assembly_read2 = bbmap_randomreads.read2
    }
}

task ncbi_datasets_download_genome_accession {
    input {
        String ncbi_accession
        Int cpu = 1
        Int memory = 4
        String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" # not the latest version, but it's hard to keep up w/ the frequent releases
        Int disk_size = 50
    }
    meta {
        # added so that call caching is always turned off
        volatile: true
    }
    command <<<
        set -euo pipefail

        date | tee DATE
        datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

        #### download FASTA file using ncbi_accession ####
        datasets download genome accession \
          ~{ncbi_accession} \
          --filename ~{ncbi_accession}.zip \
          --include genome

        # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable
        unzip ~{ncbi_accession}.zip
        cp -v ncbi_dataset/data/~{ncbi_accession}*/~{ncbi_accession}*.fna ./~{ncbi_accession}.fasta
    >>>
    output {
        File ncbi_datasets_assembly_fasta = "~{ncbi_accession}.fasta"
        String ncbi_datasets_version = read_string("DATASETS_VERSION")
        String ncbi_datasets_docker = docker
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
        maxRetries: 3
    }
}

task bbmap_randomreads {
    input {
        String ncbi_accession
        File assembly
        Int read_length = 150
        Int required_coverage = 30
        Int cpu = 1
        Int memory = 4
        String docker = "evolbioinfo/bbmap:v39.01"
        Int disk_size = 16
    }

    Int overlap = read_length / (0.5 * required_coverage )

    command <<<
    randomreads.sh \
        ref=~{assembly} \
        out1=~{ncbi_accession}.R1.fastq.gz \
        out2=~{ncbi_accession}.R2.fastq.gz \
        length=~{read_length} \
        coverage=~{required_coverage} \
        overlap=~{overlap} \
        illuminanames=t \
        paired=t
    >>>
    output {
        File read1 = "~{ncbi_accession}.R1.fastq.gz"
        File read2 = "~{ncbi_accession}.R2.fastq.gz"
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
    }
}

