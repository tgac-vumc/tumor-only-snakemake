configfile: "config.yaml"
(Samples,) = glob_wildcards("../bam/{sample}_indelq_sorted.bam")

rule all:
	input:
		expand("../vcf/somatic/non_functional/{sample}_non_functional_somatics.csv",sample=Samples)
		
#TL je hoeft alleen de laatste output file onder all te zetten, de rest wordt dan automatisch ook gemaakt,
#ook hoef je maar 1 van de outputs van die rule hier neer te zetten.

rule mpileup:
	input:
		bam="../bam/{sample}_indelq_sorted.bam"
	output:
		mpileup=temp("../bam/{sample}.mpileup")
	threads: 1
	resources:
		mem_mb=2000
	params:
		ref=config["all"]["REF"],
	shell:
		"""
		samtools mpileup -f {params.ref} {input.bam} > {output.mpileup}
		"""

rule VarScan_snps:
	input:
		#bam="../bam/{sample}_indelq_sorted.bam"
		mpileup=temp("../bam/{sample}.mpileup")
	output:
		raw_snps=temp("../VarScan/vcf/{sample}_varscan_snps.vcf")
	params:
		ref=config["all"]["REF"],
		#VarScan=config["VarScan_snps"]["VarScan_prog"],
		min_cov=config["VarScan_snps"]["min_cov"],
		min_reads2=config["VarScan_snps"]["min_reads2"],
		min_avg_qual=config["VarScan_snps"]["min_avg_qual"],
		min_var_freq=config["VarScan_snps"]["min_var_freq"],
		p_value=config["VarScan_snps"]["p_value"],
		strand_filter=config["VarScan_snps"]["strand_filter"],
		output_vcf=config["VarScan_snps"]["output_vcf"],
		variants=config["VarScan_snps"]["variants"],
	conda:
		"envs/varscan2.yaml"
	log: "../logs/VarScan/{sample}_varscan.txt"
	shell:
		"""
		varscan mpileup2snp {input.mpileup} \
		--min-coverage {params.min_cov} \
		--min-reads2 {params.min_reads2} \
		--min-avg-qual {params.min_avg_qual} \
		--min-var-freq {params.min_var_freq} \
		--p-value {params.p_value} \
		--strand-filter {params.strand_filter} \
		--output-vcf {params.output_vcf} \
		--variants {params.variants} \
		> {output.raw_snps} 2> {log}
		"""

rule VarScan_indels:
	input:
		#bam="../bam/{sample}_indelq_sorted.bam"
		mpileup=temp("../bam/{sample}.mpileup")
	output:
		raw_indels=temp("../VarScan/vcf/{sample}_varscan_indels.vcf"),
	params:
		ref=config["all"]["REF"],
		min_cov=config["VarScan_indels"]["min_cov"],
		min_reads2=config["VarScan_indels"]["min_reads2"],
		min_avg_qual=config["VarScan_indels"]["min_avg_qual"],
		min_var_freq=config["VarScan_indels"]["min_var_freq"],
		p_value=config["VarScan_indels"]["p_value"],
		strand_filter=config["VarScan_indels"]["strand_filter"],
		output_vcf=config["VarScan_indels"]["output_vcf"],
		variants=config["VarScan_indels"]["variants"],
	conda:
		"envs/varscan2.yaml"
	log: "../logs/VarScan/{sample}_varscan.txt" # TODO: nb this is the same file as in VarScan_snps
	shell:
		"""
		varscan mpileup2indel {input.mpileup} \
		--min-coverage {params.min_cov} \
		--min-reads2 {params.min_reads2} \
		--min-avg-qual {params.min_avg_qual} \
		--min-var-freq {params.min_var_freq} \
		--p-value {params.p_value} \
		--strand-filter {params.strand_filter} \
		--output-vcf {params.output_vcf} \
		--variants {params.variants} \
		> {output.raw_indels} 2> {log}
		"""

rule VarScan_combine:
	input:
		raw_indels="../VarScan/vcf/{sample}_varscan_indels.vcf",
		raw_snps="../VarScan/vcf/{sample}_varscan_snps.vcf"
	output:
		tmp=temp("../VarScan/vcf/{sample}_varscan_tmp.vcf"),
		out="../VarScan/vcf/raw/{sample}_varscan.vcf",
	conda:
		"envs/varscan2.yaml"
	shell:
		"""
		vcfcombine {input.raw_snps} {input.raw_indels} > {output.tmp} &&
		sed 's/\%//' {output.tmp} | sed 's/FREQ\,Number\=1\,Type\=String/FREQ\,Number\=1\,Type\=Float/' > {output.out}
		"""

rule LoFreq:
	input:
		bam="../bam/{sample}_indelq_sorted.bam"
	output:
		raw_snps="../LoFreq/vcf/raw/{sample}_lofreq.vcf"
	params:
		ref=config["all"]["REF"],
		threads=config["all"]["THREADS"],
		min_cov=config["LoFreq"]["min_cov"],
		min_mq=config["LoFreq"]["min_mq"],
		min_bq=config["LoFreq"]["min_bq"],
		min_alt_bq=config["LoFreq"]["min_alt_bq"],
		max_depth=config["LoFreq"]["max_depth"],
		sig=config["LoFreq"]["sig"],
	conda:
		"envs/lofreq.yaml"
	log: "../logs/LoFreq/{sample}_lofreq.txt"
	shell:
		"""
		lofreq call-parallel --pp-threads {params.threads} --call-indels --verbose {input.bam} -f {params.ref} -o {output.raw_snps} \
		--min-cov {params.min_cov} \
		--min-mq {params.min_mq} \
		--min-bq {params.min_bq} \
		--min-alt-bq {params.min_alt_bq} \
		--max-depth {params.max_depth} \
		--sig {params.sig}
		2> {log}
		"""

rule VarScan_readStatFilter:
	input:
		raw_vcf="../VarScan/vcf/raw/{sample}_varscan.vcf",
	output:
		tmp_vcf=temp("../VarScan/vcf/filtered/{sample}_varscan_tmp.vcf"),
		filtered_vcf="../VarScan/vcf/filtered/{sample}_varscan_filt.vcf",
	params:
		SnpSift_filter=config["VarScan_Filter"]["SnpSift_filter"],
		hg19_dict=config["all"]["HG19_DICT"],
	conda:
		"envs/SnpSift.yaml"
	log: "../logs/VarScan/{sample}_varscan_readStatFilter.txt"
	shell:
		"""
		SnpSift filter -f {input.raw_vcf} "(GEN[0].SDP>8) & (GEN[0].AD>3) & (GEN[0].FREQ>4.99) & (GEN[0].PVAL<0.05) & (GEN[0].ADF>0) & (GEN[0].ADR>0)" > {output.tmp_vcf}
		picard SortVcf I={output.tmp_vcf} O={output.filtered_vcf} SD={params.hg19_dict} 2> {log}
		"""


rule LoFreq_readStatFilter:
	input:
		raw_vcf="../LoFreq/vcf/raw/{sample}_lofreq.vcf",
	output:
		tmp_vcf=temp("../LoFreq/vcf/filtered/{sample}_lofreq_tmp.vcf"),
		unsorted_vcf="../LoFreq/vcf/filtered/{sample}_lofreq_unsorted.vcf",
		filtered_vcf="../LoFreq/vcf/filtered/{sample}_lofreq_filt.vcf",
	params:
		af_min=config["LoFreq_Filter"]["af_min"],
		cov_min=config["LoFreq_Filter"]["cov_min"],
		sb_alpha=config["LoFreq_Filter"]["sb_alpha"],
		SnpSift_filter=config["LoFreq_Filter"]["SnpSift_filter"],
		hg19_dict=config["all"]["HG19_DICT"],
	conda:
		"envs/lofreq.yaml"
	log: "../logs/LoFreq/{sample}_lofreq_readStatFilter.txt"
	shell:
		"""
		lofreq filter --verbose --af-min {params.af_min} --cov-min {params.cov_min} --sb-alpha {params.sb_alpha} --sb-incl-indels -i {input.raw_vcf} -o {output.tmp_vcf}
		SnpSift filter -f {output.tmp_vcf} "(DP4[2]>2) & (DP4[3]>2) & ((na HRUN) | (HRUN<8))" > {output.unsorted_vcf}
		picard SortVcf I={output.unsorted_vcf} O={output.filtered_vcf} SD={params.hg19_dict} 2> {log}
		"""

rule VarScan_BlacklistFilter:
	input:
		filt_vcf="../VarScan/vcf/filtered/{sample}_varscan_filt.vcf"
	output:
		blacklisted_vcf="../VarScan/vcf/blacklisted/{sample}_blacklisted.vcf",
		blacklisted_vcf_gz="../VarScan/vcf/blacklisted/{sample}_blacklisted.vcf.gz",
		blacklisted_csv="../VarScan/vcf/blacklisted/{sample}_blacklisted.csv",
		not_blacklisted_vcf="../VarScan/vcf/not_blacklisted/{sample}_not_blacklisted.vcf",
		not_blacklisted_vcf_gz="../VarScan/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz",
		not_blacklisted_csv="../VarScan/vcf/not_blacklisted/{sample}_not_blacklisted.csv",
	params:
		BED_blacklist=config["VarScan_Filter"]["BED_blacklist"],
		Gene_blacklist=config["VarScan_Filter"]["Gene_blacklist"],
		fields='CHROM POS REF ALT DP AF',
	conda:
		"envs/SnpSift.yaml"
	log: "../logs/LoFreq/{sample}_lofreq_readStatFilter.txt"   #TODO: log bestand 1e regel shell; heeft dit zin zo?
	shell:
		"""
		SnpSift intervals -i {input.filt_vcf} {params.BED_blacklist} {params.Gene_blacklist} > {output.blacklisted_vcf} 2> {log}
		SnpSift intervals -i {input.filt_vcf} -x {params.BED_blacklist} {params.Gene_blacklist} > {output.not_blacklisted_vcf}
		SnpSift extractFields -e "." {output.blacklisted_vcf} {params.fields} > {output.blacklisted_csv}
		SnpSift extractFields -e "." {output.not_blacklisted_vcf} {params.fields} > {output.not_blacklisted_csv}
		pbgzip -c {output.blacklisted_vcf} > {output.blacklisted_vcf_gz}
		tabix -s1 -b2 -e2 {output.blacklisted_vcf_gz}
		pbgzip -c {output.not_blacklisted_vcf} > {output.not_blacklisted_vcf_gz}
		tabix -s1 -b2 -e2 {output.not_blacklisted_vcf_gz}
		"""

rule LoFreq_BlacklistFilter:
	input:
		filt_vcf="../LoFreq/vcf/filtered/{sample}_lofreq_filt.vcf"
	output:
		blacklisted_vcf="../LoFreq/vcf/blacklisted/{sample}_blacklisted.vcf",
		blacklisted_vcf_gz="../LoFreq/vcf/blacklisted/{sample}_blacklisted.vcf.gz",
		blacklisted_csv="../LoFreq/vcf/blacklisted/{sample}_blacklisted.csv",
		not_blacklisted_vcf="../LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf",
		not_blacklisted_vcf_gz="../LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz",
		not_blacklisted_csv="../LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.csv",
	params:
		BED_blacklist=config["LoFreq_Filter"]["BED_blacklist"],
		Gene_blacklist=config["LoFreq_Filter"]["Gene_blacklist"],
		fields='CHROM POS REF ALT DP AF',
	conda:
		"envs/SnpSift.yaml"
	log: "../logs/LoFreq/{sample}_lofreq_readStatFilter.txt"  #TODO: log bestand 1e regel shell; heeft dit zin zo?
	shell:
		"""
		SnpSift intervals -i {input.filt_vcf} {params.BED_blacklist} {params.Gene_blacklist} > {output.blacklisted_vcf} 2> {log}
		SnpSift intervals -i {input.filt_vcf} -x {params.BED_blacklist} {params.Gene_blacklist} > {output.not_blacklisted_vcf}
		SnpSift extractFields -e "." {output.blacklisted_vcf} {params.fields} > {output.blacklisted_csv}
		SnpSift extractFields -e "." {output.not_blacklisted_vcf} {params.fields} > {output.not_blacklisted_csv}
		pbgzip -c {output.blacklisted_vcf} > {output.blacklisted_vcf_gz}
		tabix -s1 -b2 -e2 {output.blacklisted_vcf_gz}
		pbgzip -c {output.not_blacklisted_vcf} > {output.not_blacklisted_vcf_gz}
		tabix -s1 -b2 -e2 {output.not_blacklisted_vcf_gz}
		"""

       #TL zijn varscan en lofreq blacklistfilter hetzelfde? zo ja dan kan je hier 1 rule van maken waarbij net als {sample} {program} ook een variable is.
       # ik zou bgzip en tabix daarbij als een losse rule maken, nu zijn het wel veel commandos in 1 rule.

rule Intersect_VariantCalls:
	input:
		lofreq_vcf="../LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz",
		varscan_vcf="../VarScan/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz",
	output:
		intersect_vcf="../vcf/intersect/{sample}_intersect.vcf",
		outersect_vcf="../vcf/outersect/{sample}_outersect",
	conda:
		"envs/SnpSift.yaml"
	log: "../logs/Intersect_variantCalls/{sample}_intersectCalls.txt" 
	shell:
		"""
		bcftools isec -n=2 -w1 {input.lofreq_vcf} {input.varscan_vcf} -o {output.intersect_vcf} -O v
		bcftools isec {input.lofreq_vcf} {input.varscan_vcf} -c all -n +0 -o {output.outersect_vcf} -O v
		"""

rule Annotate_VariantCalls:
	input:
		intersect_vcf="../vcf/intersect/{sample}_intersect.vcf"
	output:
		annotated="../vcf/annotated/{sample}.annotated.vcf"
	params:
		Java_mem=config["all"]["Java_mem"],
		dbsnp=config["all"]["dbsnp"],
		clinvar=config["all"]["clinvar"],
		Cosmic=config["all"]["Cosmic"],
		gnomAD=config["all"]["gnomAD"],
		HMF_PON=config["all"]["HMF_PON"],
	conda:
		"envs/SnpSift.yaml"
	log:
		"../logs/Annotate_variantCalls/{sample}_annotation.txt"
	shell:
		"""
		SnpSift annotate {params.Java_mem} {params.dbsnp} {input.intersect_vcf} 2> {log} |
		SnpSift annotate {params.Java_mem} {params.clinvar} - 2>> {log} |
		SnpSift annotate {params.Java_mem} -v {params.Cosmic} - 2>> {log} |
		SnpSift annotate {params.Java_mem} -v -info 'gnomAD_AF' {params.gnomAD} - 2>> {log} |
		SnpSift annotate {params.Java_mem} -v -info 'PON_COUNT' {params.HMF_PON} - > {output.annotated} 2>> {log}
		"""

rule Effect_Prediction:
	input:
		annotated="../vcf/annotated/{sample}.annotated.vcf"
	output:
		effect="../vcf/annotated/{sample}.annotated_effect.vcf",
		effect_gz="../vcf/annotated/{sample}.annotated_effect.vcf.gz",
		snpEff_stats="../vcf/annotated/stats/{sample}_stats.html",
	conda:
		"envs/SnpSift.yaml"
	params:
		Java_mem=config["all"]["Java_mem"],
		targets=config["all"]["targets"],
		genome_build=config["all"]["genome_build"],
	log:
		"../logs/snpEff/{sample}_snpEff.txt"
	shell:
		"""
		snpEff {params.Java_mem} eff {params.genome_build} -filterInterval {params.targets} -v -hgvs1LetterAa -onlyProtein -strict \
		-stats {output.snpEff_stats} {input.annotated} > {output.effect} 2> {log}
		pbgzip -c {output.effect} > {output.effect_gz}
		tabix -p vcf {output.effect_gz}
		"""

rule Variant_Discrimination:
	input:
		effect="../vcf/annotated/{sample}.annotated_effect.vcf"
	output:
		somatic="../vcf/somatic/all/{sample}_all_somatics.vcf",
		snp="../vcf/snp/{sample}_all_snps.vcf",
		functional_somatic="../vcf/somatic/functional/{sample}_functional_somatics.vcf",
		non_functional_somatic="../vcf/somatic/non_functional/{sample}_non_functional_somatics.vcf",
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift filter -f {input.effect} "((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((na PON_COUNT) | (PON_COUNT<4))" > {output.somatic} &&
		SnpSift filter -f {input.effect} -n "((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((na PON_COUNT) | (PON_COUNT<4))" > {output.snp} &&
		SnpSift filter -f {output.somatic} "(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > {output.functional_somatic} &&
		SnpSift filter -f {output.somatic} -n "(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > {output.non_functional_somatic} 
		"""


rule ExtractVcfFields_csv:
	input:
		somatic="../vcf/somatic/all/{sample}_all_somatics.vcf",
		snp="../vcf/snp/{sample}_all_snps.vcf",
		functional_somatic="../vcf/somatic/functional/{sample}_functional_somatics.vcf",
		non_functional_somatic="../vcf/somatic/non_functional/{sample}_non_functional_somatics.vcf",
	output:
		somatic="../vcf/somatic/all/{sample}_all_somatics.csv",
		snp="../vcf/snp/{sample}_all_snps.csv",
		functional_somatic="../vcf/somatic/functional/{sample}_functional_somatics.csv",
		non_functional_somatic="../vcf/somatic/non_functional/{sample}_non_functional_somatics.csv",
	params:
		fields='CHROM POS REF ALT DP "DP4[2]" "DP4[3]" "SB" "HRUN" AF \
		"ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_P" "ANN[0].HGVS_C" "ANN[0].IMPACT" "ANN[0].EFFECT" \
		"COSM.ID[0]" "GENE[0]" "AA[0]" "CDS[0]" "COMMON" "PON_COUNT" "RS" "CAF" "LOF" "NMD" "MUT" \
		"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" "FATHMM[0]" "MUT.ST[0]"'
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift extractFields -e "." {input.somatic} {params.fields} > {output.somatic} && \
		SnpSift extractFields -e "." {input.snp} {params.fields} > {output.snp} && \
		SnpSift extractFields -e "." {input.functional_somatic} {params.fields} > {output.functional_somatic} && \
		SnpSift extractFields -e "." {input.non_functional_somatic} {params.fields} > {output.non_functional_somatic}
		"""
















































