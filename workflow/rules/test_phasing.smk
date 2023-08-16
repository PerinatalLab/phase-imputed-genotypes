rule test_split_bfile_related:
        'Extract samples and chr22 from the whole genotype file.'
        input:
                'results/aux/subsample/ids/related-samples-IIDs.txt',
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/500k_N_genotyped_1+_info_0.97785.snps',
                multiext('/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc', '.bed', '.bim', '.fam')
        output:
                temp('results/test/temp/all-related-chr22.vcf.gz'),
		temp('results/test/temp/good-related-chr22.vcf.gz')
        params:
                '/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc',
		'results/test/temp/all-related-chr22',
		'results/test/temp/good-related-chr22'
        log:
                'logs/test/all-related-split.txt',
		'logs/test/good-related-split.txt'
        shell:
                '''
		/home/pol.sole.navais/soft/plink2 --bfile {params[0]} --keep {input[0]} --max-alleles 2 --chr 22 --export vcf-4.2 bgz --out {params[1]} > {log[0]}
		/home/pol.sole.navais/soft/plink2 --bfile {params[0]} --keep {input[0]} --extract {input[1]} --max-alleles 2 --chr 22 --export vcf-4.2 bgz --out {params[2]} > {log[1]}
		'''

rule test_add_allele_count_related:
        'Add allele count (AC) and index vcf.'
        input:
                'results/test/temp/{quality}-related-chr22.vcf.gz'
        output:
                'results/test/subsample/{quality}-related-chr22.vcf.gz',
                'results/test/subsample/{quality}-related-chr22.vcf.gz.tbi'
	benchmark:
		'benchmarks/test/add_allele_count-{quality}.txt'
	shell:
                '''
                bcftools +fill-tags {input[0]} -Oz -o {output[0]}
                tabix -p vcf {output[0]}
                '''

rule test_phase_related_samples:
        'Phase related samples (parent-offspring duos or trios) using SHAPEIT5.'
        input:
                'results/test/subsample/{quality}-related-chr22.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr22.b37.gmap.gz',
                'results/test/subsample/{quality}-related-chr22.vcf.gz.tbi'
        output:
                'results/test/phased/{quality}-related-chr22.bcf'
        threads: 10
        log:
                'logs/test/shapeit5-phasing-{quality}-related-chr22-log.txt'
        benchmark:
                'benchmarks/test/phasing/shapeit5-phasing-{quality}-related-chr22.txt'
        shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region 22 --map {input[2]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
                '''

rule test_phase_with_scaffold:
	'Use scaffold generated above to phase all sites.'
	input:
		'results/test/subsample/all-related-chr22.vcf.gz',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr22.b37.gmap.gz',
		'results/test/phased/good-related-chr22.bcf',
                'results/test/subsample/all-related-chr22.vcf.gz.tbi'
	output:
                'results/test/phased/allsites-good/related-chr22.bcf'
	threads: 10
	log:
		'logs/test/allsites-good/shapeit5-phasing-related-chr22-log.txt'
	benchmark:
		'benchmarks/test/allsites-good/shapeit5-phasing-related-chr22.txt'
	shell:
		'''
		/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --pedigree {input[1]} --region 22 --map {input[2]} --scaffold {input[3]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
		'''

rule test_list_kids:
	'Get a list of child IIDs'
	input:
		'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt'
	output:
		'results/aux/subsample/fets_only/related-samples-IIDs.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= None, names= ['IID', 'dad', 'mom'])
		d.dropna(subset= ['dad', 'mom'], inplace= True, how= 'all')
		d.drop_duplicates('IID', keep= 'first', inplace= True)
		d.to_csv(output[0], sep= '\t', header= False, index= False, columns= ['IID'])

rule test_vcf_only_kids:
	'Extract child samples from vcf.'
	input:
		'results/test/subsample/all-related-chr22.vcf.gz',
		'results/aux/subsample/fets_only/related-samples-IIDs.txt'
	output:
		'results/test/subsample/fets_only/related-chr22.vcf.gz',
		'results/test/subsample/fets_only/related-chr22.vcf.gz.tbi'
	shell:
		'''
		bcftools view --force-samples -S {input[1]} -Oz -o {output[0]} {input[0]}
		tabix -p vcf {output[0]}
		'''

rule test_phase_only_kids:
	'Phase data using only child samples.'
        input:
                'results/test/subsample/fets_only/related-chr22.vcf.gz',
                'resources/chr22.b37.gmap.gz',
                'results/test/subsample/fets_only/related-chr22.vcf.gz.tbi'
	output:
		'results/test/phased/fets_only/related-chr22.bcf'
	threads: 10
	log:
		'logs/test/fets_only/shapeit5-phasing-related-chr22-log.txt'
	benchmark:
		'benchmarks/test/fets_only/shapeit5-phasing-related-chr22.txt'
	shell:
		'''
		/home/pol.sole.navais/soft/SHAPEIT5/phase_common_static --input {input[0]} --region 22 --map {input[1]} --output-format bcf --output {output[0]} --thread {threads} --log {log}
		'''


rule test_switch_error_good:
	'Calculate switch error rate for good markers.'
	input:
                'results/test/phased/allsites-good/related-chr22.bcf',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr22.b37.gmap.gz',
		'results/test/phased/fets_only/related-chr22.bcf'
	output:
		'results/test/switch/allsites-good/related-chr22.block.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.frequency.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.sample.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.sample.typing.txt.gz',
		'results/test/switch/allsites-good/related-chr22.variant.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.variant.typing.txt.gz'
	params:
		'results/test/switch/allsites-good/related-chr22'
	threads: 10
	log:
                'logs/test/switch/allsites-good/shapeit5-phasing-related-chr22-log.txt'
	benchmark:
                'benchmarks/test/switch/allsites-good/shapeit5-phasing-related-chr22.txt'
	shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/switch_static --validation {input[0]} --estimation {input[3]} --pedigree {input[1]} --region 22 --output {params[0]} --thread {threads} --log {log}
                '''

rule test_switch_error_all:
        'Calculate switch error rate for all markers.'
        input:
                'results/test/phased/all-related-chr22.bcf',
                'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
                'resources/chr22.b37.gmap.gz',
                'results/test/phased/fets_only/related-chr22.bcf'
        output:
                'results/test/switch/all/related-chr22.block.switch.txt.gz',
                'results/test/switch/all/related-chr22.frequency.switch.txt.gz',
                'results/test/switch/all/related-chr22.sample.switch.txt.gz',
                'results/test/switch/all/related-chr22.sample.typing.txt.gz',
                'results/test/switch/all/related-chr22.variant.switch.txt.gz',
                'results/test/switch/all/related-chr22.variant.typing.txt.gz'
	params:
		'results/test/switch/all/related-chr22'
	threads: 10
	log:
                'logs/test/switch/all/shapeit5-phasing-related-chr22-log.txt'
	benchmark:
                'benchmarks/test/switch/all/shapeit5-phasing-related-chr22.txt'
	shell:
                '''
                /home/pol.sole.navais/soft/SHAPEIT5/switch_static --validation {input[0]} --estimation {input[3]} --pedigree {input[1]} --region 22 --output {params[0]} --thread {threads} --log {log}
                '''

rule test_report_ser:
	'Report on switch error rate.'
	input:
		'results/aux/subsample/pedigree_SHAPEIT5/pedigree-parent-offspring-IIDs.txt',
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bim',
		'/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/500k_N_genotyped_1+_info_0.97785.snps',
		'results/test/switch/all/related-chr22.sample.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.sample.switch.txt.gz',
		'results/test/switch/allsites-good/related-chr22.frequency.switch.txt.gz',
		'results/test/switch/all/related-chr22.frequency.switch.txt.gz',
		'benchmarks/test/add_allele_count-all.txt',
		'benchmarks/test/phasing/shapeit5-phasing-all-related-chr22.txt',
		'benchmarks/test/add_allele_count-good.txt',
		'benchmarks/test/phasing/shapeit5-phasing-good-related-chr22.txt',
		'benchmarks/test/allsites-good/shapeit5-phasing-related-chr22.txt'
	output:
		'report/test-phase-pipeline.html'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/report-select-phasing-pipeline.Rmd'
