rule all:
    input:
        "auspice/RSV_A_tree.json",
        "auspice/RSV_A_meta.json",
        "auspice/RSV_B_tree.json",
        "auspice/RSV_B_meta.json"


rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = "data/RSV-{lineage}_ref_vp1_align_genotype.fas"
    output:
        sequences = "results/sequences_{lineage}.fasta",
        metadata = "results/metadata_{lineage}.tsv"
    params:
        fasta_fields = ["strain","virus", "id", "host", "country", "region", "date"]
    run:
    	from Bio import SeqIO
    	from datetime import datetime

    	meta = {}
    	seqs = list(SeqIO.parse(input.sequences, 'fasta'))
    	for seq in seqs:
    		m = {k:v for k,v in zip(params.fasta_fields, seq.name.split('|'))}
    		seq.name = m['strain']
    		seq.id= m['strain']
    		seq.description=''
    		print(m)
    		if len(m['date'])==4:
    			m['date'] = str(m['date'])+'-XX-XX'
    		else:
    			m['date'] = datetime.strptime(m['date'], '%d/%m/%Y').strftime('%Y-%m-%d')
    		meta[m['strain']] = m
    	SeqIO.write(seqs, output.sequences, 'fasta')
    	with open(output.metadata, 'w') as fh:
    		fh.write('\t'.join(params.fasta_fields)+'\n')
    		for m in meta.values():
	    		fh.write('\t'.join([m[x] for x in params.fasta_fields])+'\n')


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.parse.output.sequences,
        reference = "data/reference_{lineage}.gb"
    output:
        alignment = "results/aligned_{lineage}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw_{lineage}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{lineage}.nwk",
        node_data = "results/branch_lengths_{lineage}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts_{lineage}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "data/reference_{lineage}.gb"
    output:
        node_data = "results/aa_muts_{lineage}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{lineage}.json",
    params:
        columns = "region country",
        sampling_bias_correction = 2
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """

rule lbi:
	input:
		tree = rules.refine.output.tree,
		bl = rules.refine.output.node_data
	output:
		node_data = "results/LBI_{lineage}.json"
	params:
		tau = 0.7,
		window = 50,
		attr = 'lbi'
	shell:
		"""
		 augur lbi  --tree {input.tree} \
		 			--branch-lengths {input.bl}\
		 			--output {output.node_data} \
		 			--attribute-names {params.attr} \
		 			--tau {params.tau} --window {params.window}
		"""


rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        lbi = rules.lbi.output.node_data,
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"
    output:
        auspice_tree = "auspice/RSV-{lineage}_tree.json",
        auspice_meta = "auspice/RSV-{lineage}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.lbi} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta} --minify-json
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "data "
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
