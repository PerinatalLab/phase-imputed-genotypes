
rule phase_related_samples:
        'Phase related samples (parten-offspring duos or trios) using SHAPEIT5.'
        input:
                'results/geno/plink/related/related-chr{CHR}.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/genetic_maps.b37.tar.gz',
                'results/geno/plink/related/related-chr{CHR}.vcf.gz.tbi'
        output:
                temp('results/phased-haplotypes/temp/related/related-chr{CHR}.bcf')
        threads: 10
        log:
                'logs/phasing/shapeit5-phasing-related-chr{CHR}-log.txt'
        benchmark:
                'benchmarks/phasing/benchmark-shapeit5-phasing-related-chr{CHR}.txt'
	priority: 2
        shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region {wildcards.CHR} --map {input[2]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
                '''


rule phase_all_samples:
	'Phase all samples using related samples as a scaffold.'
	input:
                'results/geno/plink/all/all-chr{CHR}.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/genetic_maps.b37.tar.gz',
		'results/phased-haplotypes/temp/related/related-chr{CHR}.bcf',
                'results/geno/plink/all/all-chr{CHR}.vcf.gz.tbi'
	output:
                temp('results/phased-haplotypes/temp/all/all-chr{CHR}.bcf')
	threads: 10
	log:
                'logs/phasing/shapeit5-phasing-all-chr{CHR}-log.txt'
	benchmark:
                'benchmarks/phasing/benchmark-shapeit5-phasing-all-chr{CHR}.txt'
	priority: 3

	shell:
                '''
		/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --region {wildcards.CHR} --map {input[2]} --scaffold {input[3]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
                '''

rule convert_bcf:
	'Convert bcfs from related and unrelated samples into a vcf.'
	input:
		'results/phased-haplotypes/temp/all/all-chr{CHR}.bcf'
	output:
		temp('results/phased/haplotypes/delivery/temp/phased-MoBaPsychGen-chr{CHR}.vcf')
	threads: 2
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

