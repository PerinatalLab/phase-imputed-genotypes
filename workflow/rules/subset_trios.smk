rule select_subsample:
	'Select a subsample of duos/trios for testing phasing software in a small dataset.'
	input:
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam',
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-rel.kin'
	output:
		'results/aux/subsample/ids/related-samples-IIDs.txt',
		'results/aux/subsample/ids/sind-samples-IIDs.txt',
		'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt'
	script:
		'../scripts/select-samples.py'

rule split_bfile_related:
	'Extract samples and chr22 from the whole genotype file.'
	input:
		'results/aux/subsample/ids/related-samples-IIDs.txt',
		multiext('/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc', '.bed', '.bim', '.fam')
	output:
		temp('results/geno/subsample/temp/related-chr22.vcf.gz')
	params:
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc',
		'results/geno/subsample/temp/related-chr22'
	log:
		'logs/subsample/related-split.txt'
	shell:
		'''
		/home/pol.sole.navais/soft/plink2 --bfile {params[0]} --keep {input[0]} --max-alleles 2 --chr 22 --export vcf-4.2 bgz --out {params[1]} > {log}
		'''
rule add_allele_count:
	'Add allele count (AC) and index vcf.'
	input:
		'results/geno/subsample/temp/related-chr22.vcf.gz'
	output:
		'results/geno/subsample/related-chr22.vcf.gz',
		'results/geno/subsample/related-chr22.vcf.gz.tbi'
	shell:
		'''
		bcftools +fill-tags {input[0]} -Oz -o {output[0]}
		tabix -p vcf {output[0]}
		'''





