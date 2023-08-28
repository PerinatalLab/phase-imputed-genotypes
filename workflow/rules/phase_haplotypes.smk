
checkpoint phase_related:
	'Phase related samples (parten-offspring duos or trios) using SHAPEIT5 by chunks.'
        input:
                'results/geno/plink/related/related-chr{CHR}.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr{CHR}.b37.gmap.gz',
		'results/chunks/hg19_chunk_chr{CHR}.txt',
                'results/geno/plink/related/related-chr{CHR}.vcf.gz.tbi'
	output:
                directory('results/phased-haplotypes/temp/related/chunks/{CHR}/')
#		temp('results/phased-haplotypes/temp/related/related-chr{CHR}.bcf')
	params:
		'results/phased-haplotypes/temp/related/chunks/{CHR}/',
		'results/phased-haplotypes/temp/related/chunks/logs/{CHR}-'
	threads: 10
	priority: 2
	run:
		chunk= pd.read_csv(input[3], sep= '\t', header= None, names= ['num', 'CHR', 'pos1', 'pos2'])
		for index, row in chunk.iterrows():
			outfile= params[0] + 'chunk' + row['pos1'] + '.txt'
			logs= params[1] + 'chunk' + row['pos1'] + '.txt'
                	shell("/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region {row['pos1']} --map {input[2]} --output-format bcf --output {outfile} --thread {threads} --log {logs}")

def aggregate_chunks(wildcards):
	checkpoint_output = checkpoints.phase_related.get(**wildcards).output[0]
	return expand("results/phased-haplotypes/temp/related/chunks/{{CHR}}/chunks{chunk_id}.txt", chunk_id=glob_wildcards(os.path.join(checkpoint_output, "chunks{chunk_id}.txt")).chunk_id)

rule list_chunk_files:
	'Write chunk file names in a new file.'
	input:
		aggregate_chunks
	output:
		'results/phased-haplotypes/temp/related/list_of_chunks/CHR{CHR}.txt'
	run:
		d= input.sort()
		with open(output[0], 'w') as f:
			for line in lines:
				f.write(f'{line}\n')

rule ligate_relateds:
	'Merge bcf files from different chunks.'
	input:
		'results/phased-haplotypes/temp/related/list_of_chunks/CHR{CHR}.txt',
		'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
		aggregate_chunks
	output:
		'results/phased-haplotypes/temp/related/related-chr{CHR}.bcf'
	threads: 2
	log:
		'logs/phasing/shapeit5-ligate-related-chr{CHR}-log.txt'
	shell:
		'~/soft/SHAPEIT5/ligate_static --input {input[0]} --pedigree {input[1]} --output {output[0]} --thread {threads} --log {log}'

checkpoint phase_all_samples_chunks:
	'Phase all samples using related samples as a scaffold.'
	input:
                'results/geno/plink/all/all-chr{CHR}.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/genetic_maps.b37.tar.gz',
		'results/phased-haplotypes/temp/related/related-chr{CHR}.bcf',
		'results/chunks/hg19_chunk_chr{CHR}.txt',
                'results/geno/plink/all/all-chr{CHR}.vcf.gz.tbi'
	output:
                directory('results/phased-haplotypes/temp/all/chunks/{CHR}/')
	params:
                'results/phased-haplotypes/temp/all/chunks/{CHR}/',
                'results/phased-haplotypes/temp/all/chunks/logs/{CHR}-'
	threads: 10
	priority: 3
	run:
		chunk= pd.read_csv(input[4], sep= '\t', header= None, names= ['num', 'CHR', 'pos1', 'pos2'])
		for index, row in chunk.iterrows():
                	outfile= params[0] + 'chunk' + row['pos1'] + '.txt'
                	logs= params[1] + 'chunk' + row['pos1'] + '.txt'
                	shell("/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --region {row['pos1']} --map {input[2]} --scaffold {input[3]} --output-format bcf --output {outfile} --thread {threads} --log {logs}")


def aggregate_chunks_all():
        checkpoint_output = checkpoints.phase_all_samples_chunks.get(**wildcards).output[0]
        return expand("results/phased-haplotypes/temp/related/chunks/{{CHR}}/chunks{chunk_id}.txt", chunk_id=glob_wildcards(os.path.join(checkpoint_output, "chunks{chunk_id}.txt")).chunk_id)

rule list_chunk_files_all:
        'Write chunk file names in a new file.'
        input:
                aggregate_chunks_all
        output:
                'results/phased-haplotypes/temp/all/list_of_chunks/CHR{CHR}.txt'
        run:
                d= input.sort()
                with open(output[0], 'w') as f:
                        for line in lines:
                                f.write(f'{line}\n')

rule ligate_relateds_all:
        'Merge bcf files from different chunks.'
        input:
                'results/phased-haplotypes/temp/all/list_of_chunks/CHR{CHR}.txt',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                aggregate_chunks_all
        output:
                'results/phased-haplotypes/temp/all/related-chr{CHR}.bcf'
        threads: 2
        log:
                'logs/phasing/shapeit5-ligate-all-chr{CHR}-log.txt'
        shell:
                '~/soft/SHAPEIT5/ligate_static --input {input[0]} --pedigree {input[1]} --output {output[0]} --thread {threads} --log {log}'

rule convert_bcf:
	'Convert bcfs from related and unrelated samples into a vcf.'
	input:
		'results/phased-haplotypes/temp/all/all-chr{CHR}.bcf'
	output:
		temp('results/phased/haplotypes/delivery/temp/phased-MoBaPsychGen-chr{CHR}.vcf')
	threads: 2
	priority: 4
	shell:
		'bcftools convert -Ov -o {output[0]} {input[0]}'

rule bgzip_index_phased_vcf:
	'Block gzip and index phased vcf files.'
	input:
		'results/phased/haplotypes/delivery/temp/phased-MoBaPsychGen-chr{CHR}.vcf'
	output:
		'results/phased/haplotypes/delivery/phased-MoBaPsychGen-chr{CHR}.vcf.gz',
		'results/phased/haplotypes/delivery/phased-MoBaPsychGen-chr{CHR}.vcf.gz.tbi'
	shell:
		'''
		bgzip -c {input[0]} > {output[0]}
		tabix -p vcf {output[0]}
		'''

