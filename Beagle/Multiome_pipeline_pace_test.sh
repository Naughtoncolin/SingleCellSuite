#!/usr/bin/bash

#Description: Pipeline to run GATK on a BAM file from the Multiome analysis (BAM was generated from CellRanger mkfastq)
#GATK Reference Files: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
#BAM files are on the Gibson share: ~/p-ggibson3-shared/Multiome/bamFiles
#Author: Jen Kintzle
#Date: 12/13/2021

sample_list=("PA4")

for item in ${sample_list[@]}; do
	sample_name=$item
	#sample_name="PA4"
	folder_name="/storage/home/hcoda1/0/cnaughton7/scratch/beagle-test"
	#mkdir -p $folder_name
	orig_bam="/storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Emory/p21090-s002_P4/outs/possorted_genome_bam.bam"

	#Depending upon the approach - full file or sample file - comment out the appropriate input_bam line
	input_bam=$folder_name"/possorted_genome_bam_subset.bam" #this version is applicable to the sample file approach
	#input_bam=$orig_bam #this version is applicable to the full file approach


	#This portion of the script is used to grab a subset of the file to test out the pipeline
	#Comment it out if you are planning to run on the full file
	header_line_count=$(samtools view -H $orig_bam | wc -l)
	sample_file_size=$((header_line_count + 10000))
	samtools view -h $orig_bam | head -n $sample_file_size | samtools view -b > $input_bam
	echo "Finished creating sample file"


	#GATK Pipeline - used GATK workflow as https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels/blob/master/rna-germline-variant-calling.wdlv
	ref_fasta="/storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Reference/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta"
	hg38_dbsnp="/storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Reference/Homo_sapiens_assembly38/Homo_sapiens_assembly38.dbsnp138.vcf"
	hg38_known_indels="/storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Reference/Homo_sapiens_assembly38/Homo_sapiens_assembly38.known_indels.vcf"
	hg38_interval_list="/storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Reference/Homo_sapiens_assembly38/Homo_sapiens_assembly38.known_indels.interval_list"

	reverted_bam=$folder_name"/reverted.bam"
	reverted_bam_index=$reverted_bam".bai"
	fastq1=$folder_name"/"$sample_name".1.fastq.gz"
	fastq2=$folder_name"/"$sample_name".2.fastq.gz"
	genome_dir="/storage/home/hcoda1/0/cnaughton7/scratch/beagle-test/STAR_genome_dir"
	outFileNamePrefix=$folder_name"/"$sample_name"_"
	aligned_bam=$outFileNamePrefix"Aligned.sortedByCoord.out.bam"
	merged_bam=$folder_name"/"$sample_name"_merge.bam"
	dupes=$folder_name"/"$sample_name"_dupes.bam"
	dupe_metrics=$folder_name"/"$sample_name".metrics"
	split_bam=$folder_name"/"$sample_name"_split.bam"
	rg_bam=$folder_name"/"$sample_name"_RG.bam"
	recal_data=$folder_name"/"$sample_name".recal_data.csv"
	recalibrated_bam=$folder_name"/"$sample_name".aligned.duplicates_marked.recalibrated.bam"
	scatter=$folder_name"/out_scatter"
	vcf_file=$folder_name"/"$sample_name".vcf.gz"
	filtered_vcf=$folder_name"/"$sample_name".variant_filtered.vcf.gz"

	echo "Starting GATK pipeline for "$sample_name

	#RevertSam
	picard RevertSam --INPUT $input_bam --OUTPUT $reverted_bam --VALIDATION_STRINGENCY SILENT

	#Index bam file
	samtools index -b $reverted_bam $reverted_bam_index

	#SamToFastq
	gatk SamToFastq --INPUT $reverted_bam  --VALIDATION_STRINGENCY SILENT --FASTQ $fastq1 --SECOND_END_FASTQ $fastq2

	#StarGenerateReferences - previously run so commenting out here
	STAR --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $ref_fasta --sjdbGTFfile /storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/Reference/gencode.v38.annotation.gtf --runThreadN 12 --genomeSAsparseD 3 --genomeSAindexNbases 12

	#StarAlign
	STAR --genomeDir $genome_dir --runThreadN 12 --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate --twopassMode Basic --limitBAMsortRAM 45000000000 --outFileNamePrefix $outFileNamePrefix --readFilesIn $fastq1 #second fastq was empty so removing $fastq2

	#MergeBamAlignment 
	gatk MergeBamAlignment --REFERENCE_SEQUENCE $ref_fasta --UNMAPPED_BAM $reverted_bam --ALIGNED_BAM $aligned_bam --OUTPUT $merged_bam --INCLUDE_SECONDARY_ALIGNMENTS false --PAIRED_RUN False --VALIDATION_STRINGENCY SILENT

	#MarkDuplicates
	gatk MarkDuplicates --INPUT $merged_bam --OUTPUT $dupes --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE $dupe_metrics

	#SplitNCigarReads
	gatk SplitNCigarReads -R $ref_fasta -I $dupes -O $split_bam 

	#Add Read Groups (Not part of official pipeline but BaseRecalibrator requires the groups to function; https://gatk.broadinstitute.org/hc/en-us/community/posts/360062229451-A-USER-ERROR-has-occurred-Number-of-read-groups-must-be-1-but-is-0)
	#To validate read groups exist: samtools view -H sample.bam | grep '@RG'
	#If it returns nothing, there are no read groups
	picard AddOrReplaceReadGroups I=$split_bam O=$rg_bam RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=3

	#CallBaseRecalibrator
	gatk --java-options "-Xms4000m" BaseRecalibrator -R $ref_fasta -I $rg_bam --use-original-qualities -O $recal_data -known-sites $hg38_dbsnp -known-sites $hg38_known_indels

	#ApplyBQSR
	gatk --java-options "-Xms3000m" ApplyBQSR --add-output-sam-program-record -R $ref_fasta -I $rg_bam --use-original-qualities -O $recalibrated_bam --bqsr-recal-file $recal_data

	#Create VcfToIntervalList
	picard VcfToIntervalList I=$hg38_known_indels O=$hg38_interval_list

	#ScatterIntervalList
	picard IntervalListTools SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW UNIQUE=true SORT=true INPUT=$hg38_interval_list OUTPUT=$scatter

	#HaplotypeCaller
	gatk --java-options "-Xms6000m" HaplotypeCaller -R $ref_fasta -I $recalibrated_bam -L $hg38_interval_list -O $vcf_file --standard-min-confidence-threshold-for-calling 20 --dbsnp $hg38_dbsnp -ERC GVCF

	#VariantFiltration
	gatk VariantFiltration --R $ref_fasta --V $vcf_file --window 35 --cluster 3 -O $filtered_vcf
done