rule filter_bfile:
	'Create one vcf file for each chromosome.'
	input:
		'results/aux/subsample/ids/related-samples-IIDs.txt',
                multiext('/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc', '.bed', '.bim', '.fam')
	output:
                temp('results/geno/plink/all/temp/all-chr{CHR}.vcf.gz')
	params:
                '/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc',
                'results/geno/plink/all/temp/all-chr{CHR}'
	log:
                'logs/geno/plink/all-filter-chr{CHR}.txt'
	threads: 4
	resources: 
		mem_mb=10000
	shell:
                '''
                /home/pol.sole.navais/soft/plink2 --bfile {params[0]} --max-alleles 2 --chr {wildcards.CHR} --export vcf-4.2 bgz --threads {threads} --memory {resources.mem_mb} --out {params[1]} > {log[0]}
                '''

rule add_AC:
	'Add allele count and index vcf - SHAPEIT5 requires both.'
	input:
                'results/geno/plink/all/temp/all-chr{CHR}.vcf.gz'
	output:
                temp('results/geno/plink/all/all-chr{CHR}.vcf.gz'),
                temp('results/geno/plink/all/all-chr{CHR}.vcf.gz.tbi')
	benchmark:
                'benchmarks/geno/plink/all-chr{CHR}-add-AC.txt'
	priority: 1
	shell:
                '''
                bcftools +fill-tags {input[0]} -Oz -o {output[0]}
                tabix -p vcf {output[0]}
                '''






