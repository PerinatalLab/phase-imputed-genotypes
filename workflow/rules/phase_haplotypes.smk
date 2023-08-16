
rule phase_related_samples:
        'Phase related samples (parten-offspring duos or trios) using SHAPEIT5.'
        input:
                'results/geno/plink/related/temp/related-chr{CHR}.vcf.gz',
                'results/aux/subsample/ids_SHAPEIT5/related-samples-IIDs.txt',
                'genetic_maps.b37.tar.gz',
                'results/geno/plink/related/temp/related-chr{CHR}.vcf.gz.tbi'
        output:
                temp('results/phased-haplotypes/temp/related/related-chr{CHR}.bcf'
        threads: 10
        log:
                'logs/phasing/shapeit5-phasing-related-chr{CHR}-log.txt'
        benchmark:
                'benchmarks/phasing/benchmark-shapeit5-phasing-related-chr{CHR}.txt'
        shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region {wildcards.CHR} --map {input[2]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
                '''


rule phase_unrelated_samples:
	'Phase unrelated samples using related samples as a scaffold.'
	input:
                'results/geno/plink/unrelated/temp/unrelated-chr{CHR}.vcf.gz',
                'results/aux/subsample/ids_SHAPEIT5/related-samples-IIDs.txt',
                'genetic_maps.b37.tar.gz',
		'results/phased-haplotypes/temp/related/related-chr{CHR}.bcf',
                'results/geno/plink/related/temp/unrelated-chr{CHR}.vcf.gz.tbi'
	output:
                temp('results/phased-haplotypes/temp/unrelated/unrelated-chr{CHR}.bcf')
	threads: 10
	log:
                'logs/phasing/shapeit5-phasing-unrelated-chr{CHR}-log.txt'
	benchmark:
                'benchmarks/phasing/benchmark-shapeit5-phasing-unrelated-chr{CHR}.txt'

	shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region {wildcards.CHR} --map {input[2]} --scaffold {input[3]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
                '''

rule merge_bcf:
	'Merge phased bcfs from related and unrelated samples into a vcf.'
	input:
		''
	output:
		''
	run:
		''

rule bgzip_index_phased_vcf:
	'Block gzip and index phased vcf files.'

