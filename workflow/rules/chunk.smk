

rule chunk_imputed:
	''
	input:
		'results/geno/plink/related/related-chr{CHR}.vcf.gz',
		'resources/chr{CHR}.b37.gmap.gz'
	output:
		'results/chunks/hg19_chunk_chr{CHR}.txt'
	shell:
		'/home/pol.sole.navais/soft/GLIMPSE2/GLIMPSE2_chunk_static --input {input[0]} --map {input[1]} --region {wildcards.CHR} --sequential --output {output[0]}'
