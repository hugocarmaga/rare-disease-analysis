configfile: "config/config.yaml"

intersects = [i for i in config["intersections"].keys()]
graph = config["aug_graph"]
ref = config["reference"]
callsets = config["callsets"]

rule extract_exclusive_vcf:
    input:
        inter = lambda wildcards: config["intersections"][wildcards.intersect]['table'],
        vcf = lambda wildcards: config["intersections"][wildcards.intersect]['vcf']
    output:
        "results/vcfs/{intersect}/{intersect}_exclusive.vcf.gz"
    conda:
        "../envs/extraction.yaml"
    wildcard_constraints:
        intersect = "|".join(intersects)
    params:
        fields_excl = " && ".join(['${}=="True"'.format(i + 10) for i in range(callsets[0])]),
        fields_miss = " && ".join(['${}=="False"'.format(i + 10) for i in range(callsets[0], callsets[1])]),
        id_rep = 10 + callsets[1]
    shell:
        """
        awk -F "\\t" '{params.fields_excl} && {params.fields_miss} {{print${params.id_rep}}}' {input.inter} | bcftools view --include ID=@/dev/stdin -Oz -o {output} {input.vcf}
        """

rule extract_alleles_not_found:
    input:
        inter = lambda wildcards: config["intersections"][wildcards.intersect]['table'],
        vcf = lambda wildcards: config["intersections"][wildcards.intersect]['vcf']
    output:
        vcf = "results/vcfs/{intersect}/{intersect}_not-found.vcf.gz",
        fasta = "results/fasta/{intersect}_alleles_not-found.fa"
    conda:
        "../envs/extraction.yaml"
    wildcard_constraints:
        intersect = "|".join(intersects)
    params:
        fields_excl = " && ".join(['${}=="True"'.format(i + 10) for i in range(callsets[1]) if i != callsets[0]]),
        fields_miss = '${}=="False"'.format(callsets[0] + 10),
        id_rep = 10 + callsets[1]
    shell:
        """
        awk -F "\\t" '{params.fields_excl} && {params.fields_miss} {{print${params.id_rep}}}' {input.inter} | bcftools view --include ID=@/dev/stdin -Oz -o {output.vcf} {input.vcf}
        python workflow/scripts/extract_fasta_alleles_fastx.py -f {ref} -v {output.vcf} -o {output.fasta}
        """

rule realign_to_graph:
    input:
        "results/fasta/{intersect}_alleles_not-found.fa"
    output:
        initial = temp("results/realignment/gaf/{intersect}_tmp.gaf"),
        final = "results/realignment/gaf/{intersect}_not-found_wholeLen.gaf"
    conda:
        "../envs/extraction.yaml"
    threads: 4
    shell:
        """
        GraphAligner -g {graph} -f {input} -a {output.initial} -t {threads} --multimap-score-fraction 1 --seeds-mxm-windowsize 5000 --seeds-mxm-length 30 --seeds-mem-count 10000 -C -1 --bandwidth 15 --discard-cigar
        awk -F "\\t" '$4-$3==$2{{a[$1]++; b[$1]=$0}}END{{for(i in a) if(a[i]==1) print b[i]}}' {output.initial} > {output.final}
        """

rule plot_seqIDs:
    input:
        "results/realignment/gaf/{intersect}_not-found_wholeLen.gaf"
    output:
        "results/realignment/plots/seqID_distribution-log_{intersect}.svg"
    conda:
        "../envs/plotting.yaml"
    shell:
        "python workflow/scripts/seqIdentity_plots.py -g {input} -o {output} -s {wildcards.intersect}"

rule plot_af_distribution:
    input:
        table = lambda wildcards: config["intersections"][wildcards.intersect]['table'],
        af_ids = "workflow/resources/SAGA_ID_AF.tsv"
    output:
        "results/af_comparison/af-distribution-saga-{intersect}.svg"
    conda:
        "../envs/plotting.yaml"
    wildcard_constraints:
        intersect = "|".join(intersects)
    params:
        calls = ",".join(str(x) for x in callsets)
    shell:
        """python workflow/scripts/af-distribution.py -i {input.af_ids} -t {input.table} -o {output} -s {wildcards.intersect} -c {params.calls}"""
