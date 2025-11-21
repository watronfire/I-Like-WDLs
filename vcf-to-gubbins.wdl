version 1.0

workflow generate_tree_from_vcf {
    input {
        File alignment
        File reference = "gs://bacpage-resources/vc_reference.union.fasta"
        Array[String] alignment_names
        Array[String] dates
        String tree_prefix
        String? outgroup
        Int iqr_clock_filter = 3
        String docker = "staphb/bcftools:1.22"
        String iqtree_model = "GTR+G"
        Int bootstraps = 1000
        String iqtree_opts = ""
        String iqtree_docker = "staphb/gubbins:3.4.1"
        Int disk_size = 16
        Int memory = 16
        Int gubbins_disk_size = 64
        Int gubbins_memory = 16
        Int gubbins_cpu = 16
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
            outgroup = outgroup,
            disk_size = gubbins_disk_size,
            memory = gubbins_memory,
            cpu = gubbins_cpu
    }
    call generate_masked_vcf {
        input:
            alignment = alignment,
            recombinant_sites = gubbins.gubbins_recombination_gff
    }
    call sparsify_alignment {
        input:
            alignment = gubbins.gubbins_polymorphic_fasta
    }

    call generate_tree {
        input:
            sparse_alignment = sparsify_alignment.sparse_alignment,
            reference = reference,
            starting_tree = gubbins.gubbins_final_tree,
            cluster_name = tree_prefix,
            outgroup = outgroup,
            iqtree_model = iqtree_model,
            iqtree_bootstraps = bootstraps,
            iqtree_opts = iqtree_opts,
            docker = iqtree_docker,
    }

    call clock_rate_filter {
        input:
            ml_tree = generate_tree.ml_tree,
            sample_names = alignment_names,
            sample_dates = dates,
            iqr = iqr_clock_filter
    }

    call build_usher_tree {
        input:
            rooted_tree = clock_rate_filter.rooted_tree,
            vcf_alignment = generate_masked_vcf.masked_alignment
    }

    output {
        String  bcftools_version = vcf_to_fasta. bcftools_version
        String  gubbins_version = gubbins.version
        File    gubbins_final_tree = gubbins.gubbins_final_tree
        File    gubbins_final_labelled_tree = gubbins.gubbins_final_labelled_tree
        File    gubbins_polymorphic_fasta = gubbins.gubbins_polymorphic_fasta
        File    gubbins_recombination_gff = gubbins.gubbins_recombination_gff
        File    gubbins_branch_stats = gubbins.gubbins_branch_stats
        File    masked_vcf = generate_masked_vcf.masked_alignment
        File    ml_tree = generate_tree.ml_tree
        File    rooted_tree = clock_rate_filter.rooted_tree
        File    rtt_distances = clock_rate_filter.rtt_distances
        File    rtt_plot = clock_rate_filter.rtt_plot
        File    protobuf_file = build_usher_tree.protobuf_file
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
        String docker = "staphb/gubbins:3.4.1"
        Int? filter_percent = 25 #default is 25%
        Int? iterations = 5
        String? tree_builder = "hybrid"
        String? tree_args
        String? nuc_subst_model = "GTRGAMMA"
        Int? bootstrap = 0
        String? outgroup
        File? dates_file
        Int disk_size = 64 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }

    command <<<
        # date and version control
        date | tee DATE
        run_gubbins.py --version | tee VERSION

        run_gubbins.py \
            ~{alignment} \
            --prefix ~{cluster_name} \
            --iterations ~{iterations} \
            --tree-builder ~{tree_builder} \
            ~{'--tree-args ' + tree_args} \
            ~{'--model ' + nuc_subst_model} \
            --bootstrap ~{bootstrap} \
            ~{'--outgroup ' + outgroup} \
            ~{'--date ' + dates_file} \
            --threads ~{cpu}
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
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 50
    }
}

task generate_masked_vcf {
    input {
        File alignment
        File recombinant_sites
        String docker = "watronfire/vibecheck:2025.07.30"
    }
    command <<<
        bcftools view -Oz -o input.vcf.gz ~{alignment}
        bcftools view -h ~{alignment} > header.txt

        python << CODE
        import pandas as pd
        import io
        import gzip

        def read_vcf(path):
            with gzip.open(path, "rt") as f:
                lines = [l for l in f if not l.startswith('##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t'
            )

        gff_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        gff = pd.read_csv( "~{recombinant_sites}", sep="\t", header=None, comment="#", names=gff_columns )
        gff["taxa"] = gff["attribute"].str.extract( r'taxa="([^"]*)"' )

        vcf = read_vcf( "input.vcf.gz" )

        for _, entry in gff.iterrows():
            vcf.loc[vcf["POS"].between( entry["start"], entry["end"]), entry["taxa"].split()] = "."

        with open( "input.masked.vcf", "wt" ) as f:
            f.write("##fileformat=VCFv4.2\n")
            vcf.to_csv( f, sep="\t", index=False)
        CODE

        bcftools reheader -h header.txt input.masked.vcf | bcftools view -Oz -o masked_alignment.vcf.gz
    >>>
    output {
        File masked_alignment = "masked_alignment.vcf.gz"
    }
    runtime {
        docker: "~{docker}"
        cpu: 1
        memory: 8 + " GiB"
        disks: "local-disk 16 SSD"
        bootDiskSizeGb: 10
    }
}

task sparsify_alignment {
    input {
        File alignment
        String docker = "staphb/snp-sites"
    }
    command <<<
        snp-sites -o sparse_alignment.fasta ~{alignment}
    >>>
    output {
        File sparse_alignment = "sparse_alignment.fasta"
    }
    runtime {
        docker: "~{docker}"
        cpu: 1
        memory: 16 + " GiB"
        disks: "local-disk 24 SSD"
        bootDiskSizeGb: 50
    }
}

task generate_tree {
    input {
        File sparse_alignment
        String cluster_name
        File? starting_tree
        File reference = "gs://bacpage-resources/vc_reference.union.fasta"
        String iqtree_model = "GTR+G" # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
        String iqtree_bootstraps = 1000 #  Ultrafast bootstrap replicates
        String? iqtree_opts = ""
        String? outgroup
        String docker = "staphb/gubbins:3.4.1"
    }
    command <<<
        python << CODE
        counts = {
            "A" : 0,
            "C" : 0,
            "G" : 0,
            "T" : 0
        }

        with open( "~{reference}", "rt" ) as f:
            for line in f:
                if not line.startswith(">"):
                    upper_line = line.upper()
                    for char in upper_line:
                        if char in counts:
                            counts[char] += 1

        with open( "sites.txt", "wt" ) as f:
            f.write( ",".join([str( counts[c] ) for c in ["A", "C", "G", "T"]]) )
        CODE

        cp ~{sparse_alignment} msa.fasta
        iqtree \
            -nt AUTO \
            -s msa.fasta \
            -m ~{iqtree_model} \
            -bb ~{iqtree_bootstraps} \
            -fconst $(cat sites.txt) \
            ~{'-t ' + starting_tree} \
            ~{'-o ' + outgroup} \
            ~{iqtree_opts}

        cat msa.fasta.treefile | sed  's/|_/|?/g' > ~{cluster_name}_iqtree.tree
    >>>
    output {
        File ml_tree = "~{cluster_name}_iqtree.tree"
    }
    runtime {
        docker: "~{docker}"
        memory: "16 GB"
        cpu: 16
        disks: "local-disk 24 SSD"
    }
}

task clock_rate_filter {
    input {
        File ml_tree
        Array[String] sample_names
        Array[String] sample_dates
        Int iqr = 3
        String docker = "nextstrain/base:build-20251119T000157Z"
        Int memory = 24
        Int cpu = 4
        Int disk_space = 24
    }
    Array[Array[String]] dates = transpose( [sample_names, sample_dates] )
    File dates_f = write_tsv( dates )
    command <<<

        echo "node_name\tdate\n" > dates.tsv
        cat ~{dates_f} >> dates.tsv

        treetime clock \
            --tree ~{ml_tree} \
            --dates dates.tsv \
            --clock-filter ~{iqr} \
            --keep-root \
            --prune-outliers \
            --outdir clock_result
    >>>
    output {
        File rooted_tree = "clock_results/rerooted.newick"
        File rtt_distances = "clock_results/rtt.csv"
        File rtt_plot = "clock_results/root_to_tip_regression.pdf"
    }
runtime {
        docker: docker
        memory: memory + " GB"
        cpu: cpu
        disks:  "local-disk " + disk_space + " HDD"
        preemptible: 0
    }
}

task build_usher_tree {
    input {
        File rooted_tree
        File vcf_alignment
        String docker = "us-docker.pkg.dev/general-theiagen/pathogengenomics/usher:0.6.2"
        Int memory = 16
        Int cpu = 4
        Int disk_size = 24
    }
    command <<<
        usher -t ~{rooted_tree} -v ~{vcf_alignment} -o global_tree.pb -c -d output/
    >>>
    output {
        File protobuf_file = "global_tree.pb"
    }
    runtime {
        docker: docker
        memory: memory + " GB"
        cpu :  cpu
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        dx_instance_type: "mem3_ssd1_v2_x4"
    }
}