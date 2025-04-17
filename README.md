# Code for analyzing genome resequencing data 


## Sequence data quality control
 
> Software: FastQC v0.12.1, fastp v0.23.4

Code for filtration

```shell
fastp -i "${RAW}/${per_sample}_raw_1.fq.gz" \
-I "${RAW}/${per_sample}_raw_2.fq.gz" \
-o "${QC_PATH}/${per_sample}_filter_1.fq.gz" \
-O "${QC_PATH}/${per_sample}_filter_2.fq.gz" \
-q 20 \
-l 140 \
--detect_adapter_for_pe \
--html "REPORT_PATH/${per_sample}_report.html" \
--json "REPORT_PATH/${per_sample}_report.json" \
-w 16
```
Code for quality check


```shell
fastqc --outdir ${outdir} \
--noextract \
--threads 10 \
${sample}.fq.gz
```

## Alignment and converting alignment file format


> Software: BWA-MEM2 v2.2.1, Samtools v1.20

Code for alignment and coverting sam to bam

```shell
bwa-mem2 mem -t 68 \
-R "@RG\tID:${sample_name}\tPL:ILLUMINA\tSM:${sample_name}" \ 
-M "${Ref}" \ 
"${QC_PATH}/${sample_name}_filter_1.fq.gz" \ 
"${QC_PATH}/${sample_name}_filter_2.fq.gz | samtools view -S -b - > "${BAM_DIR}/${sample_name}_unsorted.bam"
```

## Sorting, de-duplication and indexing


> Software:Samtools v1.20,picard v3.1.1

Code for sorting the bam file

```shell
samtools sort -@ 12 \
-m 4G \
-o "${BAM_DIR}/${sample_name}_sorted.bam" \ 
"${BAM_DIR}/${sample_name}_unsorted.bam"
```

Code for removing duplication

```shell
picard MarkDuplicates \
--REMOVE_DUPLICATES true \
-I "${input}" \
-M "${BAM_DIR}/${sample_name}_dup.txt" \
-O "${output}"
```

Code for indexing bam file

```shell
samtools index -@ 16 "${output}" "${output}.bai"
```

## Variant calling


> Software: GATK v4.6.1.0

Code for single sample calling

```shell
parallel --jobs 16 \
java -jar $GATK \
HaplotypeCaller \
-R /project_dir/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna \
-I "/project_dir/dup_del_bam/{}_dup_del.bam" \
-ERC GVCF \
--native-pair-hmm-threads 4 \
-O {}.g.vcf \
::: ${sample_name[@]}
```

Code for combining GVCF files

```shell
java -Xmx384G -jar $GATK \
CombineGVCFs \
-R /project_dir/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna \ $(ls *.g.vcf | sed 's/^/-V /') \
-O res/merge_site.g.vcf.gz \
--tmp-dir tmp/
```


Code for joint genotyping

```shell
java -jar -Xmx256G $GATK \
GenotypeGVCFs \
-R /project_dir/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna \
-V res/merge_site.g.vcf.gz \
-O joint_genotype/joint_genotype.vcf.gz
```


## Variants extraction and filtration


> Software: GATK v4.6.1.0, Vcftools v0.1.16 , bcftools v1.20

Code for SNPs selection

```shell
java -jar $GATK SelectVariants \
-V joint_genotype/joint_genotype.vcf.gz \
-O dc_raw_snp.vcf.gz \
--select-type-to-include SNP
```

Code for Indels selection

```shell
java -jar $GATK SelectVariants \
-V joint_genotype/joint_genotype.vcf.gz \
-O dc_raw_snp.vcf.gz \
--select-type-to-include INDEL --max-indel-size 50
```

Code for filtration of SNPs

```shell
java -jar $GATK VariantFiltration \
-R /project_dir/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna \
-V dc_raw_snp.vcf.gz \
-O dc_filtered_snp.vcf.gz \
--filter-name "QD2" --filter-expression "QD < 2.0" \
--filter-name "FS60" --filter-expression "FS > 60.0" \
--filter-name "MQ40" --filter-expression "MQ < 40.0" \
--filter-name "ReadPosRs" --filter-expression "ReadPosRankSum < -8.0" \
--filter-name "MQRankSum" --filter-expression"MQRankSum < -12.5" 
```


Code for filtration of Indels

```shell
java -jar $GATK VariantFiltration \
-R /project_dir/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna  \
-V dc_raw_indel.vcf.gz \
-O dc_filtered_indel.vcf.gz \
--filter-name "QD2" --filter-expression "QD < 2.0" \
--filter-name "FS200" --filter-expression "FS > 200.0" \
--filter-name "MQ40" --filter-expression "MQ < 40.0" \
--filter-name "ReadPosRs" --filter-expression "ReadPosRankSum < -10.0" \
--filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" \
--filter-name "SOR10" --filter-expression "SOR > 10.0" 
```

Code for filtrarion by vcftools

```shell
vcftools \
--gzvcf dc_filtered_snp.vcf.gz \
--max-missing 1 \
--max-alleles 2 \
--min-alleles 2 \
--minDP 15 \
--remove-filtered-all \
--recode \
--stdout | gzip -c > filtered_snp_missing1_alles2_dp15.vcf.gz 

vcftools \
--gzvcf dc_filtered_indel.vcf.gz \
--max-missing 0.9 \
--max-alleles 2 \
--min-alleles 2 \
--minDP 10 \
--remove-filtered-all --recode --stdout | gzip -c > filtered_indel_mis0-9_allele2_dp10.vcf.gz

```

Code for computing af

```shell
bcftools +fill-tags \
filtered_indel_mis0-9_allele2_dp10.vcf.gz \
-Oz -o indels_with_maf.vcf.gz -- -t AF

bcftools +fill-tags \
filtered_snp_missing1_alles2_dp15.vcf.gz \
-Oz -o snp_with_maf.vcf.gz -- -t AF

```

Code for filtration by low maf

this Python script is used to filter variants with a low maf (https://github.com/kunmonster/cattle-genome-analysis.git).



Actually,we can also compute maf using bcftools and filter the variants

```shell

bcftools +fill-tags \
vcf_file \
-Oz -o var_with_maf.vcf.gz -- -t MAF

bcftools view -i 'INFO/MAF>0.05' vcf_file -Oz -o output

```
