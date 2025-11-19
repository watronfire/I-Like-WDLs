version 1.0

workflow generate_tree_from_vcf {
    input {
        File alignment
        File reference = "gs://bacpage-resources/vc_reference.union.fasta"
        String tree_prefix
        String? outgroup
        String docker = "staphb/bcftools:1.22"
        Int disk_size = 16
        Int memory = 16
    }

    call vcf_to_fasta {
        input:
            alignment = alignment,
            reference = reference,
            docker = docker,
            disk_size = disk_size,
            memory = memory
    }
    call gubbins {
        input:
            alignment = vcf_to_fasta.expanded_alignment,
            cluster_name = tree_prefix,
            outgroup = outgroup
    }

    output {
        String  bcftools_version = vcf_to_fasta. bcftools_version
        String  gubbins_version = gubbins.version
        File    gubbins_final_tree = gubbins.gubbins_final_tree
        File    gubbins_final_labelled_tree = gubbins.gubbins_final_labelled_tree
        File    gubbins_polymorphic_fasta = gubbins.gubbins_polymorphic_fasta
        File    gubbins_recombination_gff = gubbins.gubbins_recombination_gff
        File    gubbins_branch_stats = gubbins.gubbins_branch_stats
    }
}

task vcf_to_fasta {
    input {
        File alignment
        File reference = "gs://bacpage-resources/vc_reference.union.fasta"
        String docker = "staphb/bcftools:1.22"
        Int disk_size = 16
        Int memory = 16
    }

    command <<<
        bcftools --version | head -n1 | cut -f2 -d' ' | tee BCFTOOLS_VERSION

        bcftools view -Ob -o compressed_alignment.bcf.gz ~{alignment}
        bcftools index compressed_alignment.bcf.gz

        for sample in $(bcftools query -l compressed_alignment.bcf.gz); do
            bcftools consensus --mark-del N -f ~{reference} -s ${sample} compressed_alignment.bcf.gz |\
            sed -E "s,>.+$,>${sample},g" >> expanded_alignment.fasta
        done
    >>>

    output {
        String bcftools_version = read_string( "BCFTOOLS_VERSION" )
        File expanded_alignment = "expanded_alignment.fasta"
    }
    runtime {
        memory: "~{memory} GB"
        docker: docker
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
    }
}

task gubbins {
    input {
        File alignment
        String cluster_name
        String docker = "sangerpathogens/gubbins"
        Int? filter_percent = 25 #default is 25%
        Int? iterations = 5
        String? tree_builder = "raxml"
        String? tree_args
        String? nuc_subst_model = "GTRGAMMA"
        Int? bootstrap = 0
        String? outgroup
        File? dates_file
    }

    command <<<
        # date and version control
        date | tee DATE
        run_gubbins.py --version | tee VERSION

        run_gubbins.py \
            ~{alignment} \
            --prefix ~{cluster_name} \
            --filter-percentage ~{filter_percent} \
            --iterations ~{iterations} \
            --tree-builder ~{tree_builder} \
            ~{'--tree-args ' + tree_args} \
            ~{'--model ' + nuc_subst_model} \
            --bootstrap ~{bootstrap} \
            ~{'--outgroup ' + outgroup} \
            ~{'--date ' + dates_file} \
            --threads 4
    >>>
    output {
        String date = read_string("DATE")
        String version = read_string("VERSION")
        File gubbins_final_tree = "~{cluster_name}.final_tree.tre"
        File gubbins_final_labelled_tree = "~{cluster_name}.node_labelled.final_tree.tre"
        File gubbins_polymorphic_fasta = "~{cluster_name}.filtered_polymorphic_sites.fasta"
        File gubbins_recombination_gff = "~{cluster_name}.recombination_predictions.gff"
        File gubbins_branch_stats = "~{cluster_name}.per_branch_statistics.csv"
        File? gubbins_timetree = "~{cluster_name}.final_tree.timetree.tre"
        File? gubbins_timetree_stats = "~{cluster_name}.lsd.out"
    }
    runtime {
        docker: "~{docker}"
        memory: "32 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
        preemptible: 0
        maxRetries: 1
    }
}