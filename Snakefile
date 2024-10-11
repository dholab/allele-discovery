from datetime import datetime

# because clustering of shared allele is nondeterminisitc, add timestamp to output FASTA
timestamp=datetime.now()
current_time=timestamp.strftime("%Y%m%d%H%M%S")

## CONFIG ##

configfile: "config/config.yaml"

## FILES ##

input_dir = config["params"]["input_dir"]
ccs_fastq, = glob_wildcards(f'{input_dir}/{{id}}.fastq.gz')

## RULES ##

rule all:
	input:
		expand('results/04-trimmed/{id}.fastq.gz', id=ccs_fastq),
		'results/15-genotyping/genotyping.xlsx',
		'results/13-muscle-cdna-identical/cdna_matches.fasta',
		'results/13-muscle-cdna-identical/cdna.gff',
		'results/14-novel/novel.fasta',
		'results/14-novel/novel.gff'
	run:
		# remove intermediate files
		shell('rm -rf results/06-pbaa/*')

## pbaa clustering ##

rule map_ccs_to_reference:
	'''
	Map CCS HiFi to a set of reference sequences using minimap2 and output in SAM format.
	'''
	input:
		f'{input_dir}/{{id}}.fastq.gz',
		config['pbaa']['mapping_reference_fasta']
	output:
		temp('results/03-mapped/{id}.sam'),
		'logs/01_map_ccs_to_reference/map_ccs_to_reference_{id}.log'
	threads: 1
	run:
		import time

		# wait 3 seconds because there are I/O errors otherwise when too many map at once
		time.sleep(3)

		# run minimap2 in spliced mode, return SAM format
		# this seems to give more comprehensive results than bbmap
		# do not retain secondary alignments
		shell('minimap2 \
		-t {threads} \
		{input[1]} \
		{input[0]} \
		-ax map-hifi --secondary=no -A 1 -B 1 --end-bonus 30 \
		> {output[0]} \
		2> {output[1]}')

rule extract_soft_clipped_contigs:
	'''
	I wrote a custom parser to extract aligned sequences that have soft trims both upstream and downstream of the aligned region.
	This should extract sequences fully spanning the reference sequence.
	'''
	input:
		'results/03-mapped/{id}.sam'
	output:
		temp('results/04-trimmed/{id}.fastq')
	threads: 1
	script:
		'scripts/extract_soft_clipped_contigs.py'

rule gzip_compress_trimmed_fastq:
	'''
	compress trimmed FASTQ so these can be included in final results
	'''
	input:
		'results/04-trimmed/{id}.fastq'
	output:
		'results/04-trimmed/{id}.fastq.gz',
		'logs/03_gzip_compress_trimmed_fastq/gzip_compress_trimmed_fastq_{id}.log'
	threads: 1
	run:
		shell('gzip -c {input[0]} > {output[0]} 2> {output[1]}')

rule index_fastq:
	'''
	make index from each FASTQ file
	indexing is necessary for pbaa
	'''
	input:
		'results/04-trimmed/{id}.fastq'
	output:
		temp('results/05-index-fastq/{id}.fastq'),
		temp('results/05-index-fastq/{id}.fastq.fai'),
		'logs/04_index_fastq/index_fastq_{id}.log'
	threads: 1
	run:
		# copy FASTQ and index
		shell('cp \
		{input[0]} {output[0]} \
		&& samtools faidx {output[0]} 2> {output[2]}')

rule run_pbaa:
	'''
	run pbaa on each FASTQ that has been soft-clipped to full-length sequences
	'''
	input:
		config['pbaa']['guide_fasta'],
		'results/05-index-fastq/{id}.fastq',
		'results/05-index-fastq/{id}.fastq.fai'
	output:
		'results/06-pbaa/{id}_passed_cluster_sequences.fasta',
		'logs/05_run_pbaa/run_pbaa_{id}.log'
	params:
		pbaa_output_basename = '{id}'
	threads: 4
	run:
		shell('/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads {threads} \
		--skip-chimera-detection \
		{input[0]} \
		{input[1]} \
		results/06-pbaa/{params.pbaa_output_basename} \
		 2> {output[1]}')

rule rename_for_laa:
	'''
	pbaa adds extraneous sequence to the output fasta
	rename to remove
	'''
	input:
		'results/06-pbaa/{id}_passed_cluster_sequences.fasta'
	output:
		temp('results/07-pbaa-renamed-clusters/{id}.fasta')
	run:
		import shutil
		shutil.copy(input[0], output[0])

## find clusters shared between two or more animals ##

rule cluster_per_sample:
	'''
	run bbbmap dedupe.sh on each sample, finding exact matches and absorbing containments.
	use PacBio clustering parameters from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/ but change edit distance to 0.
	This should restrict clustering to identical sequences.
	outbest parameter is not documented in dedupe.sh but is used in the example from the PacBio clustering and saves one representative sequence per cluster.
	Since the sequences in the cluster are identical and have an edit distance of 0, this exemplar sequence should be sufficient.

	Some samples might not have output. Create empty output file and overwrite if data exists.
	'''
	input:
		'results/07-pbaa-renamed-clusters/{id}.fasta'
	output:
		temp('results/08-cluster_per_sample/{id}.fasta.gz'),
		'logs/07_cluster_per_sample/cluster_per_sample_{id}.log'
	threads: 1
	run:
		shell('gzip < /dev/null > {output[0]}')
		shell('dedupe.sh -Xmx1g ow=t \
		in={input[0]} \
		outbest={output[0]} \
		fo c \
		threads={threads} \
		2> {output[1]}')

rule rename_clusters:
	'''
	prepend sample name to cluster names. This makes it easier to track which sequences are associated with which clusters.
	'''
	input:
		'results/08-cluster_per_sample/{id}.fasta.gz'
	output:
		temp('results/09-renamed/{id}.fasta.gz'),
		'logs/08_rename_clusters/rename_clusters_{id}.log'
	params:
		sample_name = '{id}'
	threads: 1
	run:
		shell('rename.sh -Xmx1g \
		in={input[0]} \
		out={output[0]} \
		prefix={params.sample_name} \
		addprefix=t \
		threads={threads} \
		2> {output[1]}')

rule merge_per_animal_clusters:
	'''
	gymnastics with per-sample merged FASTA to merge them after all have been created - need to wait until rename_clusters completes before running this rule
	use of lambda function to join all samples into single string in command documented here: https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
	this may be broadly useful in cases when we need to merge together files but wait until all the input files exist
	'''
	input:
		expand('results/09-renamed/{per_sample}.fasta.gz', per_sample=ccs_fastq)
	output:
		temp('results/10-merged/merged_clusters.fasta.gz'),
		'logs/09_merge_per_animal_clusters/merge_per_animal_clusters.log'
	params:
		per_sample_files = lambda wildcards, input: " ".join(input)
	run:
		# cat of gzip-compressed FASTA files caused ASCII errors. So decompress FASTA and recompress
		shell('zcat {params.per_sample_files} | gzip > {output[0]} 2> {output[1]}')

rule shared_animals:
	'''
	run bbbmap dedupe.sh on each sample, finding exact matches, absorbing containments, and finding overlaps that are at least 3kb long.
	this requires a three step process:
	1. Run dedupe.sh and output a FASTA file containing all singletons and duplicated sequences
	2. Run dedupe.sh again and output a FASTA file containing only singletons.

	Now there are singleton sequences in both FASTA files, but sequences that are duplicates (and by definition found in more than one animal)
	in only one file. So if we find the sequences that are unique *between* these two files, we are left with only the duplicate sequences.
	So step 3:
	3. Remove singleton sequences found in both files

	A version of this is described in https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/
	'''
	input:
		'results/10-merged/merged_clusters.fasta.gz'
	output:
		temp('results/11-shared/shared.fasta'),
		temp('results/11-shared/all.fasta'),
		temp('results/11-shared/unique.fasta'),
		temp('results/11-shared/putative_alleles_temp.fasta'),
		'logs/10_shared_animals/shared_animals.log'
	threads: 1
	run:
		# create all singletons and unique sequences
		shell('dedupe.sh -Xmx1g \
		in={input[0]} \
		outbest={output[1]} \
		am=t ac=f arc=t fo c fcc nam=4 threads={threads} 2> {output[4]}')

		# find only unique sequences
		shell('dedupe.sh -Xmx1g \
		in={input[0]} \
		out={output[2]} \
		am=t ac=f arc=t fo fcc uniqueonly=t threads={threads} 2> {output[4]}')

		# find sequences that are not unique to steps 1 and 2
		# these are the sequences shared between two or more animals
		shell('dedupe.sh -Xmx1g \
		in={output[1]},{output[2]} \
		out={output[0]} \
		ac=f uniqueonly=t threads={threads} 2> {output[4]}')

		# run deduple on shared sequences to absorb near-identical matches
		shell('dedupe.sh -Xmx1g \
		in={output[0]} \
		out={output[3]} \
		ac=t threads={threads} 2> {output[4]}')

rule rename_putative_allele_clusters:
	'''
	add integer FASTA id to gdna_match FASTA name
	this simplifies downstream analyses, where the FASTA description can be complex and the FASTA id is simple
	since previous steps are nondeterminisitic, add timestamp to id too
	this makes it easier to clarify which files are used for genotyping
	'''
	input:
		'results/11-shared/putative_alleles_temp.fasta'
	output:
		temp('results/11-shared/putative_alleles.fasta')
	run:
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord
		from Bio.Seq import Seq

		# parse input FASTA
		with open(output[0], "w") as handle:
			for idx, record in enumerate(SeqIO.parse(input[0], "fasta")):
				record.id = str(current_time) + '-' + str(idx)
				SeqIO.write(record, handle, "fasta")

## classify MHC class I sequences into gdna matches, cdna extensions, and novel sequences ##

rule map_shared_clusters_to_full_length_gDNA:
	'''
	identify putative alleles whose sequences match full-length gDNA sequences already in IPD.
	save BAM file of reads that map to known IPD alleles and FASTA.gz file of reads that do not.
	'''
	input:
		config['classify']['gdna_reference_fasta'],
		'results/11-shared/putative_alleles.fasta'
	output:
		temp('results/12-gdna-identical/all_mappings.sam'),
		'logs/12_map_shared_clusters_to_full_length_gDNA/map_shared_clusters_to_full_length_gDNA.log'
	threads: 1
	run:
		# map to genomic DNA reference sequences
		# return SAM results
		# when returning secondary mappings, sequence only shows correctly for primary alignment
		# but since we only care about exact matches (NM=0), this isn't really a problem
		shell('minimap2 \
		-t {threads} \
		{input[0]} \
		{input[1]} \
		-ax splice -N 10000 \
		> {output[0]} \
		2> {output[1]}')

rule filter_exact_gdna_matches:
	'''
	filter mappings to only those that have NM:i:0 (no mismatches)
	use filterlines.sh tool
	'''
	input:
		'results/12-gdna-identical/all_mappings.sam',
		'results/11-shared/putative_alleles.fasta',
		config['classify']['gdna_reference_fasta']
	output:
		temp('results/12-gdna-identical/gdna_match.sam'),
		temp('results/13-muscle-cdna-identical/no-gdna_match.fasta'),
		temp('results/12-gdna-identical/' + os.path.basename(config['classify']['gdna_reference_fasta'])),
		'logs/12_map_shared_clusters_to_full_length_gDNA/map_shared_clusters_to_full_length_gDNA.log'
	run:
		# copy putative alleles FASTA file to 09-gdna-identical folder
		shell('cp {input[1]} {output[2]}')

		# filter SAM file
		# save mappings where NM=0
		shell('filterlines.sh \
		in={input[0]} \
		out={output[0]} \
		names=NM:i:0 \
		substring=t \
		include=t \
		2> {output[3]}')

		# create FASTA file of clusters that do not match to gDNA sequence
		shell('filterbyname.sh \
		in={input[1]} \
		names={output[0]} \
		out={output[1]} \
		2> {output[3]}')

		# copy gDNA reference to 12-gdna-identical
		shell('cp {input[2]} {output[2]}')

rule map_shared_clusters_to_cDNA_with_muscle:
	'''
	among sequences that don't match existing full-length gDNA sequences
	find ones that extend known cDNA sequences to gDNA versions

	rewritten in 23230 to exhaustively pairwise align gDNA sequences to cDNA library with MUSCLE
	and return alignments where the number of matched nucleotides exactly matches the length of the cDNA
	this is because the minimap2-based method I used previously failed to find cDNA matches that have a 5nt exon 8

	Split gDNA into separate files to process separately. That allows for pseudo-multiprocessing.
	'''
	input:
		'results/13-muscle-cdna-identical/no-gdna_match.fasta',
		config['classify']['cdna_reference_fasta']
	output:
		temp('results/13-muscle-cdna-identical/merged.aln'),
		'logs/13_map_shared_clusters_to_cDNA_with_muscle/map_shared_clusters_to_cDNA_with_muscle.log'
	run:
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord
		from Bio.Seq import Seq

		# create gDNA sequence object
		for gdna_record in SeqIO.parse(input[0], "fasta"):

			# create cDNA sequence object
			for cdna_record in SeqIO.parse(input[1], "fasta"):

				# write a tab delimited line containing the gdna_record.id, cdna_record.id, cdna_record length, and match length
				# use command substitution to get the match length
				# run MUSCLE on sequence objects, but do it in a stream to avoid a lot of file I/O
				# pipe output to CLUSTALW format
				# then count the number of '*' characters that denote matches between the two sequences
				# this works for class I
				# for class II, the cDNA can be longer than the gDNA so this doesn't work
				# if the count of natching characters equals the number of cDNA characters, write to file
				# add maxiters = 2 to accelerate processing per https://www.drive5.com/muscle/manual/compromise.html

				shell( 'echo "' + gdna_record.name  + '\t' + cdna_record.name + '\t' + str(len(gdna_record.seq)) + '\t' + str(len(cdna_record.seq)) + '\t' + ' \
					$(echo ">' + gdna_record.name + '\n' + str(gdna_record.seq) + '\n>' + cdna_record.name + '\n' + str(cdna_record.seq) + '" \
					| muscle -maxiters 2 -quiet -clwstrict 2> {output[1]} \
					| grep "^ " | grep "\*" -o | wc -l )" >> {output[0]}')

rule find_muscle_cdna_gdna_matches:
	'''
	run awk to find gdna sequences that fully match sequence of cdna
	'''
	input:
		'results/13-muscle-cdna-identical/merged.aln'
	output:
		temp('results/13-muscle-cdna-identical/matches.aln')
	run:
		shell('awk \'{{if( $4 == $5  ) print $0}}\' {input[0]} > {output[0]}')

rule rename_muscle_cdna_matched_fasta:
	'''
	input gDNA sequences
	change the name of the sequence header
	write output sequences that match cDNAs to new file
	'''
	input:
		'results/13-muscle-cdna-identical/no-gdna_match.fasta',
		'results/13-muscle-cdna-identical/matches.aln'
	output:
		'results/13-muscle-cdna-identical/cdna_matches.fasta'
	run:
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord
		from Bio.Seq import Seq
		import csv

		# create output file in case it is empty
		shell('touch {output[0]}')

		# read input FASTA line-by-line
		for record in SeqIO.parse(input[0], "fasta"):

			# parse file with gdna sequences that match cdna sequences
			with open(input[1]) as tsvfile:
				reader = csv.reader(tsvfile, delimiter='\t')
				for row in reader:

					# test if name of sequence in cdna match file matches gdna sequence name
					if row[0] == record.name:
						# update name of sequence in output file
						record.description = row[1] + '|' + row[0]

						# write to file
						with open(output[0], "a") as handle:
							SeqIO.write(record, handle, "fasta")

rule preliminary_exonerate:
	'''
	run exonerate as in 23188
	return GFF files with annotations

	this preliminary annotation relative to HLA-A is good but often has sequence before the starting methionine and after the stop codon

	therefore, parse the annotations with gffread to get protein sequence with correct stop codon.
	Use python regexp to trim protein sequence to starting methionine.

	then re-run exonerate in protein mode to get more accurate annotations
	'''
	input:
		config['classify']['hla_mrna_reference'],
		'results/13-muscle-cdna-identical/cdna_matches.fasta',
		config['classify']['hla_cds_annotation']
	output:
		temp('results/13-muscle-cdna-identical/mapped.gff'),
		'logs/14_preliminary_exonerate/preliminary_exonerate.log'
	run:
		# run exonerate
		shell('exonerate \
		--showtargetgff \
		--showalignment FALSE \
		--showvulgar FALSE \
		--model cdna2genome \
		--query {input[0]} \
		--target {input[1]} \
		--refine full \
		--annotation {input[2]} \
		> {output[0]} \
		2> {output[1]}')

rule preliminary_exonerate_process_gff:
	'''
	prepare exonerate GFF output for Geneious
	need to pass experiment and paths to needed tar.gz files
	script is in Docker container
	'''
	input:
		'results/13-muscle-cdna-identical/mapped.gff'
	output:
		temp('results/13-muscle-cdna-identical/processed.gff'),
		'logs/15_preliminary_exonerate_process_gff/preliminary_exonerate_process_gff.log'
	run:
		shell('bash /scripts/process-gff.sh \
			-e {input[0]}  \
			-p /scripts/21295-exonerate_gff_to_alignment_gff3.pl \
			-o {output[0]} \
			2> {output[1]}')

rule preliminary_exonerate_merge_cds:
	'''
	when Geneious-compatible GFF files are made,
	the CDS annotations are viewed as independent annotations.
	This step adds a name and ID annotation to CDS annotations
	so they view as a single annotation in Geneious and can be automatically translated.

	Since we only work with CDS annotations, only print these to final file.
	'''
	input:
		'results/13-muscle-cdna-identical/processed.gff'
	output:
		temp('results/13-muscle-cdna-identical/preliminary-annotations.gff')
	run:
		# extract cds annotations
		shell('awk \'{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' {input[0]} >> {output[0]}')

rule trim_annotations:
	'''
	I couldn't figure out how to trim the start and stop codons from exonerate annotations
	even after several days of trying.

	Evenutally I decided to write my own tool to do this.
	'''
	input:
		'results/13-muscle-cdna-identical/preliminary-annotations.gff',
		'results/13-muscle-cdna-identical/no-gdna_match.fasta',
	output:
		'results/13-muscle-cdna-identical/cdna.gff'
	script:
		'scripts/trim_annotations.py'

rule extract_novel_sequences:
	'''
	the previous steps processed clusters that perfectly match known IPD cDNA sequences.
	Because the sequences in this experiment are supported by PacBio cDNA and pbaa gDNA sequences,
	we can submit them to IPD even if they don't match a current IPD cDNA sequence.

	Map gDNA sequences to cDNAs and report reads that don't map

	'''
	input:
		'results/13-muscle-cdna-identical/no-gdna_match.fasta',
		'results/13-muscle-cdna-identical/cdna_matches.fasta',
	output:
		'results/14-novel/novel.fasta',
		'logs/17_extract_novel_sequences/extract_novel_sequences.log'
	run:
		shell('mapPacBio.sh in={input[0]} ref={input[1]} outu={output[0]} subfilter=0 2> {output[1]}')

rule novel_exonerate:
	'''
	run exonerate on novel sequences
	'''
	input:
		config['classify']['hla_mrna_reference'],
		'results/14-novel/novel.fasta',
		config['classify']['hla_cds_annotation']
	output:
		temp('results/14-novel/mapped.gff'),
		'logs/18_novel_exonerate/novel_exonerate.log'
	run:
		# test if novel file empty
		if os.stat(input[1]).st_size == 0:
			shell('touch {output[0]}')
		else:
			# run exonerate
			shell('exonerate \
			--showtargetgff \
			--showalignment FALSE \
			--showvulgar FALSE \
			--model cdna2genome \
			--query {input[0]} \
			--target {input[1]} \
			--refine full \
			--annotation {input[2]} \
			> {output[0]} \
			2> {output[1]}')

rule novel_exonerate_process_gff:
	'''
	prepare exonerate GFF output for Geneious
	need to pass experiment and paths to needed tar.gz files
	'''
	input:
		'results/14-novel/mapped.gff'
	output:
		temp('results/14-novel/processed.gff'),
		'logs/19_novel_exonerate_process_gff/novel_exonerate_process_gff.log'
	run:
		shell('bash /scripts/process-gff.sh \
			-e {input[0]}  \
			-p /scripts/21295-exonerate_gff_to_alignment_gff3.pl \
			-o {output[0]} \
			2> {output[1]}')

rule novel_exonerate_merge_cds:
	'''
	when Geneious-compatible GFF files are made,
	the CDS annotations are viewed as independent annotations.
	This step adds a name and ID annotation to CDS annotations
	so they view as a single annotation in Geneious and can be automatically translated.

	Since we only work with CDS annotations, only print these to final file.
	'''
	input:
		'results/14-novel/processed.gff'
	output:
		temp('results/14-novel/preliminary-annotations.gff')
	run:
		# extract cds annotations
		shell('awk \'{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' {input[0]} >> {output[0]}')

rule novel_trim_annotations:
	'''
	I couldn't figure out how to trim the start and stop codons from exonerate annotations
	even after several days of trying.

	Evenutally I decided to write my own tool to do this.
	'''
	input:
		'results/14-novel/preliminary-annotations.gff',
		'results/14-novel/novel.fasta',
	output:
		'results/14-novel/novel.gff'
	script:
		'scripts/trim_annotations.py'

rule merge_reads:
	'''
	concatenate files to align with clustal omega
	'''
	input:
		config['classify']['gdna_reference_fasta'],
		'results/14-novel/novel.fasta',
		'results/13-muscle-cdna-identical/cdna_matches.fasta'
	output:
		temp('results/14-novel/reads.fasta')
	run:
		shell('cat {input[0]} {input[1]} {input[2]} > {output[0]}')

rule clustal_align:
	'''
	align reads with clustal omega
	generate alignment in fasta format and distance matrix
	distance matrix can be used to parse closest matches to known alleles
	'''
	input:
		'results/14-novel/reads.fasta',
		'results/14-novel/novel.fasta',
	output:
		temp('results/14-novel/aligned.fasta'),
		temp('results/14-novel/distances.txt'),
		'logs/21_clustal_align/clustal_align.log'
	threads: workflow.cores
	run:
		# test if novel file empty
		if os.stat(input[1]).st_size == 0:
			shell('touch {output[0]} {output[1]}')
		else:
			shell('clustalo \
				--infile={input[0]} \
				--outfile={output[0]} \
				--distmat-out={output[1]} \
				--threads={threads} \
				--full \
				2> {output[2]}')

rule parse_distances:
	'''
	parse clustal omega distances to find alleles in reference database and cDNA extensions that are nearest neighbors to novel sequences
	'''
	input:
		'results/14-novel/distances.txt',
		'results/14-novel/novel.fasta',
		'results/13-muscle-cdna-identical/cdna_matches.fasta'
	output:
		temp('results/14-novel/novel_closest_matches.xlsx'),
		temp('results/14-novel/distances_tmp.txt')
	script:
		'scripts/parse_clustalo_distances.py'

## genotype against sequences in this dataset ##

rule create_genotyping_fasta:
	'''
	to check accuracy of novel sequences, want to genotype IPD gDNA matches, cDNA matches, and novel sequences
	need to create FASTA file that contains all of these putative_alleles along with their classification
	'''
	input:
		'results/11-shared/putative_alleles.fasta',
		config['classify']['gdna_reference_fasta'],
		'results/13-muscle-cdna-identical/matches.aln',
		'results/14-novel/novel_closest_matches.xlsx'
	output:
		temp('results/15-genotyping/classified.fasta')
	script:
		'scripts/create_genotyping_fasta.py'

rule genotype_ccs:
	'''
	idenitfy reference sequences that are fully contained in CCS reads
	use trimmed full-length reads for mapping
	use higher minratio and smaller maxindel per emails with Brian Bushnell
	'''
	input:
		'results/15-genotyping/classified.fasta',
		'results/04-trimmed/{id}.fastq'
	output:
		temp('results/15-genotyping/{id}.sam'),
		'logs/23_genotype_ccs/genotype_ccs_{id}.log'
	threads: 4
	run:
		# find reads with no substitutions relative to reference
		# return SAM file
		# minimap2 eqx flag is essential for reforma.sh to parse SAM file
		# per email with Brian Bushnell 14 July
		shell('minimap2 \
		{input[0]} \
		{input[1]} \
		-ax map-hifi --eqx -t 3 \
		2> {output[1]} \
		| reformat.sh \
		in=stdin.sam \
		out={output[0]} \
		ref={input[0]} \
		noheader=t \
		subfilter=0 \
		threads=1 \
		ow=t \
		-da \
		2> {output[1]}')

rule create_genotyping_csv:
	'''
	create CSV file containing sample name, allele name, and count of reads matching allele
	'''
	input:
		expand('results/15-genotyping/{per_sample}.sam', per_sample=ccs_fastq)
	output:
		temp('results/15-genotyping/genotyping.csv')
	params:
		# gymnastics with per-sample merged SAM to merge them after all have been created - need to wait until SAM files made before running this rule
		# use of lambda function to join all samples into single string in command documented here: https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
		per_sample_files = lambda wildcards, input: ",".join(input)
	script:
		'scripts/create_genotyping_csv.py'

rule create_genotyping_pivot:
	'''
	make Excel-formatted pivot table from genotyping CSV
	'''
	input:
		'results/15-genotyping/genotyping.csv'
	output:
		'results/15-genotyping/genotyping.xlsx'
	script:
		'scripts/genotyping.py'
