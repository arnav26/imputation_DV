CONFORM_GT = "/cluster/work/pausch/arnav/imputation/suisag_all_target/conform-gt.24May16.cee.jar"
BEAGLE = "/cluster/work/pausch/arnav/imputation/suisag_all_target/beagle.18Apr22.4de.jar"

rule all:
    input: expand("/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_imputed_merged_ref.vcf.gz.csi",chr = range(1,18))

# Section 1: Prepare REF
rule filter_vcf_ref:
    input:
        vcf = "/cluster/work/pausch/vcf_SSC/2022_07/vcf_beagle/{chr}_beagle.vcf.gz"
    output:
        filtered_vcf = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_filtered.vcf.gz",
        filtered_vcf_index = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_filtered.vcf.gz.csi"
    params:
        bcftools_params = "-q 0.01:minor --threads 4 -O z"
    threads: 4
    shell:'''bcftools view {input} {params.bcftools_params} > {output.filtered_vcf}
    bcftools index {output.filtered_vcf}
    '''
        
rule add_AN_AC_tags_ref:
    input:
        filtered_vcf = rules.filter_vcf_ref.output.filtered_vcf,
        filtered_vcf_index = rules.filter_vcf_ref.output.filtered_vcf_index
    output:
        filtered_vcf_tags = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_filtered_AN_AC_tags.vcf.gz",
        filtered_vcf_tags_index = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_filtered_AN_AC_tags.vcf.gz.csi"
    threads: 4    
    shell: '''bcftools +fill-tags {input.filtered_vcf} --threads 4 -Oz -o {output.filtered_vcf_tags} -- -t AN,AC
    bcftools index {output.filtered_vcf_tags}
    '''

rule shapeit5_ref:
    input:
        filtered_vcf_tags = rules.add_AN_AC_tags_ref.output.filtered_vcf_tags,
        filtered_vcf_tags_index = rules.add_AN_AC_tags_ref.output.filtered_vcf_tags_index,
        pedigree = "/cluster/work/pausch/arnav/imputation/reference_dv_set/pedigree"
    output:
        phased_bcf = temp("/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_phased.bcf")
    params:
         shapeit5_params = "--thread 50 --hmm-ne 100 --region {chr}"
    threads: 50
    resources:
        walltime = "4:00"
    shell:"shapeit5 -I {input.filtered_vcf_tags} --pedigree {input.pedigree} {params.shapeit5_params} -O {output.phased_bcf}"

rule bcf_to_vcf_ref:
    input:
        phased_bcf = rules.shapeit5_ref.output.phased_bcf
    output:
        phased_vcf = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_phased.vcf.gz",
        phased_vcf_index = "/cluster/work/pausch/arnav/imputation/reference_dv_set/{chr}_phased.vcf.gz.csi"
    shell:'''bcftools view {input.phased_bcf} -Oz -o {output.phased_vcf}
    bcftools index {output.phased_vcf}
    '''

# Section 2: Prepare Target
rule conform_target:
    input:
        target_vcf = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first.vcf.gz",
        reference_vcf = rules.bcf_to_vcf_ref.output.phased_vcf
    output:
        conformed_vcf = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed.vcf.gz",
        conformed_vcf_index = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed.vcf.gz.csi"
    resources:
        walltime = "4:00",
        mem_mb = 4000
    params:
        conform_gt = CONFORM_GT,
        conform_gt_params = "chrom={chr} match=POS",
        conform_gt_prefix = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed"
    shell:'''module load openjdk
    java -Xmx4g -jar {params.conform_gt} gt={input.target_vcf} ref={input.reference_vcf} {params.conform_gt_params} out={params.conform_gt_prefix}
    bcftools index {output.conformed_vcf}
    '''

rule add_AN_AC_tags_target:
    input:
        conformed_vcf = rules.conform_target.output.conformed_vcf,
        conformed_vcf_index = rules.conform_target.output.conformed_vcf_index
    output:
        conformed_vcf_tags = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed_AN_AC_tags.vcf.gz",
        conformed_vcf_tags_index = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed_AN_AC_tags.vcf.gz.csi"
    resources:
        walltime = "4:00",
        mem_mb = 4000
    threads: 4
    shell: '''bcftools +fill-tags {input.conformed_vcf} --threads 4 -Oz -o {output.conformed_vcf_tags} -- -t AN,AC
    bcftools index {output.conformed_vcf_tags}
    '''

rule phase_target:
    input:
        conformed_vcf_tags = rules.add_AN_AC_tags_target.output.conformed_vcf_tags,
        conformed_vcf_tags_index = rules.add_AN_AC_tags_target.output.conformed_vcf_tags_index,
        reference_vcf = rules.bcf_to_vcf_ref.output.phased_vcf,
        pedigree = "/cluster/work/pausch/arnav/imputation/reference_dv_set/pedigree"
    output:
        phased_bcf = temp("/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed_phased.bcf")
    threads: 10
    resources:
        walltime = "4:00",
        mem_mb = 4000
    params:
        shapeit5_params = "--thread 10 --hmm-ne 100 --region {chr}"
    shell:"shapeit5 -I {input.conformed_vcf_tags} --pedigree {input.pedigree} --reference {input.reference_vcf} {params.shapeit5_params} -O {output.phased_bcf}"

rule bcf_to_vcf_target:
    input:
        phased_bcf = rules.phase_target.output.phased_bcf
    output:
        phased_vcf = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed_phased.vcf.gz",
        phased_vcf_index = "/cluster/work/pausch/arnav/imputation/suisag_all_target/dup_suppress_first/{chr}_Suisag_all_dup_suppress_first_conformed_phased.vcf.gz.csi"
    shell:'''bcftools view {input.phased_bcf} -Oz -o {output.phased_vcf}
    bcftools index {output.phased_vcf}
    '''

# Section 3: Imputation

rule beagle_imputation:
    input:
        phased_vcf = rules.bcf_to_vcf_target.output.phased_vcf,
        phased_vcf_index = rules.bcf_to_vcf_target.output.phased_vcf_index,
        reference_vcf = rules.bcf_to_vcf_ref.output.phased_vcf
    output: 
        imputed_vcf = "/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_dup_suppress_imputed.vcf.gz",
        imputed_vcf_index = "/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_dup_suppress_imputed.vcf.gz.csi"
    threads: 50
    resources:
        walltime = "4:00",
        mem_mb = 4000
    params:
        beagle = BEAGLE,
        beagle_prefix = "/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_dup_suppress_imputed",
        beagle_params = "chrom={chr} ne=100 nthreads=50"
    shell:'''module load openjdk
    java -Xmx200g -jar {params.beagle} gt={input.phased_vcf} ref={input.reference_vcf} {params.beagle_params} out={params.beagle_prefix}
    bcftools index {output.imputed_vcf}
    '''

rule merge_imputed_ref:
    input:
        imputed_vcf = rules.beagle_imputation.output.imputed_vcf,
        imputed_vcf_index = rules.beagle_imputation.output.imputed_vcf_index,
        reference_vcf = rules.bcf_to_vcf_ref.output.phased_vcf
    output:
        merged_vcf = "/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_imputed_merged_ref.vcf.gz",
        merged_vcf_index = "/cluster/work/pausch/arnav/imputation/suisag_all_imputed/dups_suppressed/{chr}_Suisag_all_imputed_merged_ref.vcf.gz.csi"
    threads: 10
    resources:
        walltime = "4:00",
        mem_mb = 4000
    shell: '''bcftools merge {input.reference_vcf} {input.imputed_vcf} -Oz -o {output.merged_vcf} --force-samples --threads 10
    bcftools index {output.merged_vcf}
    '''
