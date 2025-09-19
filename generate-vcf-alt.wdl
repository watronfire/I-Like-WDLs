version 1.0

workflow generate_vcf_alt {
    input {
        Array[File] consensus_sequences
        File reference = "gs://bacpage-resources/vc_reference.fasta"
        Array[String] sample_names
        Array[String] replacement_names

        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 4
    }

    Boolean rename_samples = defined( replacement_names )

    call generate_vcf {
        input:
            consensus_sequences = consensus_sequences,
            reference = reference,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size
    }

    if (rename_samples) {
        call rename_vcf {
            input:
                vcf = generate_vcf.vcf,
                original_names = sample_names,
                replacement_names = replacement_names
        }
    }


    output {
        File    vcf = generate_vcf.vcf
        File?   renamed_vcf = rename_vcf.renamed_vcf
    }
}

task generate_vcf {
    input {
        Array[File] consensus_sequences
        File reference
        Int cpu = 4
        Int memory = 16
        String docker = "quay.io/biocontainers/ucsc-fatovcf:482--hdc0a859_1"
        Int disk_size = 16
    }

    command <<<
        # concatenate reference
        echo "Concatenating reference"
        sed '1h;/>/d;H;$!d;x;s/\\n/@/;s/\\n//g;s/@/\\n/' ~{reference} | sed -e '$a\\' > concat_reference.fasta
        REFERENCE=$(head -n1 ~{reference} | cut -f2 -d \> | cut -f1 -d" ")

        # Generate fasta alignment
        echo "Generating fasta alignment"
        cat concat_reference.fasta ~{sep=" " consensus_sequences} > temp_alignment.fasta

        # Convert to VCF
        echo "Generating VCF alignment with snp-sites"
        faToVcf -ambiguousToN -vcfChrom=${REFERENCE} temp_alignment.fasta alignment.vcf
    >>>
    output {
        File vcf = "alignment.vcf"
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
    }
}

task rename_vcf {
    input {
        File vcf
        Array[String] original_names
        Array[String] replacement_names
        Int cpu = 1
        Int memory = 4
        String docker = "staphb/bcftools:1.22"
        Int disk_size = 4
    }

    Array[Array[String]] names = transpose( [original_names, replacement_names] )
    File names_f = write_tsv( names )

    command <<<
        bcftools view -Ou ~{vcf} | bcftools reheader --samples ~{names_f} --output renamed_alignment.vcf.gz -
    >>>
    output {
        File renamed_vcf = "renamed_alignment.vcf.gz"
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
    }
}

