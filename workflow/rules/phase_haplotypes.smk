
checkpoint phase_all:
	'Phase all samples (parent-offspring duos or trios + unrelated samples) using SHAPEIT5 by chunks.'
        input:
                'results/geno/plink/all/all-chr{CHR}.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr{CHR}.b37.gmap.gz',
		'results/chunks/hg19_chunk_chr{CHR}.txt',
                'results/geno/plink/all/all-chr{CHR}.vcf.gz.tbi'
	output:
                directory('results/phased-haplotypes/temp/all/chunks/{CHR}/')
	params:
		'results/phased-haplotypes/temp/all/chunks/{CHR}/',
		'results/phased-haplotypes/temp/all/chunks/logs/{CHR}-',
		'results/phased-haplotypes/temp/all/chunks/logs/'
	threads: 10
	priority: 3
	run:
		chunk= pd.read_csv(input[3], sep= '\t', header= None, names= ['numi', 'CHR', 'pos1', 'pos2', 'x1', 'x2', 'x3', 'x4'])
		try:
                        os.makedirs(params[0])
		except FileExistsError:
                        # directory already exists
                        pass
		try:
			os.makedirs(params[2])
		except FileExistsError:
			# directory already exists
			pass
		for index, row in chunk.iterrows():
			pos1= row['pos1']
			print(pos1)
			outfile= params[0] + 'chunk' + str(row['numi']) + '.bcf'
			logs= params[1] + 'chunk' + str(row['numi']) + '.txt'
                	shell("/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region {pos1} --map {input[2]} --output-format bcf --output {outfile} --thread {threads} --log {logs}")

def aggregate_chunks(wildcards):
	checkpoint_output = checkpoints.phase_all.get(**wildcards).output[0]
	return expand("results/phased-haplotypes/temp/all/chunks/{{CHR}}/chunks{chunk_id}.bcf", chunk_id=glob_wildcards(os.path.join(checkpoint_output, "chunks{chunk_id}.bcf")).chunk_id)

rule list_chunk_files:
	'Write chunk file names in a new file.'
	input:
		aggregate_chunks
	output:
		'results/phased-haplotypes/temp/all/list_of_chunks/CHR{CHR}.txt'
	priority: 4
	run:
		d= input.sort()
		with open(output[0], 'w') as f:
			for line in lines:
				f.write(f'{line}\n')

rule ligate_all:
	'Merge bcf files from different chunks.'
	input:
		'results/phased-haplotypes/temp/all/list_of_chunks/CHR{CHR}.txt',
		'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
		aggregate_chunks
	output:
		'results/phased-haplotypes/temp/all/all-chr{CHR}.bcf'
	threads: 2
	log:
		'logs/phasing/shapeit5-ligate-all-chr{CHR}-log.txt'
	priority: 5
	shell:
		'~/soft/SHAPEIT5/ligate_static --input {input[0]} --pedigree {input[1]} --output {output[0]} --thread {threads} --log {log}'

rule convert_bcf:
	'Convert bcfs from related and unrelated samples into a vcf.'
	input:
		'results/phased-haplotypes/temp/all/all-chr{CHR}.bcf'
	output:
		temp('results/phased/haplotypes/delivery/temp/phased-MoBaPsychGen-chr{CHR}.vcf')
	threads: 2
	priority: 6
	shell:
		'bcftools convert -Ov -o {output[0]} {input[0]}'

rule bgzip_index_phased_vcf:
	'Block gzip and index phased vcf files.'
	input:
		'results/phased/haplotypes/delivery/temp/phased-MoBaPsychGen-chr{CHR}.vcf'
	output:
		'results/phased/haplotypes/delivery/phased-MoBaPsychGen-chr{CHR}.vcf.gz',
		'results/phased/haplotypes/delivery/phased-MoBaPsychGen-chr{CHR}.vcf.gz.tbi'
	priority: 7
	shell:
		'''
		bgzip -c {input[0]} > {output[0]}
		tabix -p vcf {output[0]}
		'''

