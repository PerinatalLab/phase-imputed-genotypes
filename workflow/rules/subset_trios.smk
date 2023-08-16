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

