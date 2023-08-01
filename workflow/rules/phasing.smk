rule phase_related_samples:
	'Phase related samples (parten-offspring duos or trios) using SHAPEIT5.'
	input:
		'results/geno/subsample/related-chr22.vcf.gz',
		'results/aux/subsample/ids_SHAPEIT5/related-samples-IIDs.txt',
		'resources/chr22.b37.gmap.gz',
		'results/geno/subsample/related-chr22.vcf.gz.tbi'
	output:
		'results/geno/phased/related-chr22.bcf'
	threads: 10
	log:
		'logs/phasing/shapeit5-phasing-related-chr22-log.txt'
	benchmark:
		'benchmarks/phasing/shapeit5-phasing-related-chr22.txt'
	shell:
		'''
		/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region 22 --map {input[2]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
		'''
