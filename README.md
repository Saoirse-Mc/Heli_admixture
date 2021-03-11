##Create BWA index##
bwa index -p Hmel2.5bwaidx -a bwtsw Heliconius_melpomene_melpomene_Hmel2.5.fa
##-p (name of index can be whatever you choose) -a(changes depending on long wtsw) or short reads (is))
## bwa index -p whatevername - a bwtsw/is referencegenome.fa


##Convert any files (not all have to be converted) to .fq.gz##
bzcat HM_CS17.read1.bz2 | gzip -c >HM_CS17_R1.fq.gz

##run bwa mem on paired reads - to map to reference## 

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="bwa"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

module load bwa/0.7.15

bwa mem /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/tristero/CAM40013_combined_R1.fq.gz /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/tristero/CAM40013_combined_R2.fq.gz > /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/tristero/CAM40013_R1R2.sam
## bwa mem reference.fa read1.fq.gz read2.fq.gz > r1r2.sam ##


##Convert sam to bam files##

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="SamToBam"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

module load samtools/1.4.1

#loop through directories

for dir in */
do
	#loop through files
        for file in "$dir"/*.sam
        do
		#check if file exists
                if [[ -f $file ]]
                then
			#write message into output, which file is being converted
                        echo "Converting $file to ${file%.*}.bam"
			#convert
                        samtools view -bS $file > ${file%.*}.bam
                fi
        done
done

##Check if bam files are truncated (not necessary but handy for checking)##

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="Checktruncated"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --output=truncatedcheck
#SBATCH --time=1-00:00:00

module load samtools/1.4.1

samtools quickcheck -v *.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'


##Sort bam files for SNP calling## 

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="SortBam"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

module load samtools/1.4.1

#loop through directories
for dir in */
do
	#loop through files
        for file in "$dir"/*.bam
        do
		#check if file exists
                if [[ -f $file ]]
                then
			#write message into output, which file is being converted
                        echo "Converting $file to ${file%.*}.sorted.bam"
			#convert
                samtools sort $file -o ${file%.*}.sorted.bam
                fi
        done
done

##Indexing sorted bam files for calling snps##

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="IndexBam"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

module load samtools/1.4.1

#loop through directories
for dir in */
do
	#loop through files
        for file in "$dir"/*.sorted.bam
        do
		#check if file exists
                if [[ -f $file ]]
                then
			#write message into output, which file is being converted
                        echo "indexing $file"
			#index
                        samtools index $file
                fi
        done
done


##Extract melpomene and numata idividuals from simons data using bcf tools## 
##- tried with gatk super overcomplicated and didn't work##

#!/bin/bash
#SBATCH --partition=usual
#SBATCH --job-name="Bcf_Ext_ind"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --output=Bcf_Ext_ind
#SBATCH --time=1-00:00:00

module load bcftools/1.4.1

bcftools view -s nu_sil.MJ09-4184 -s other_ind_ids Wbar92.DP8MP4BIMAC2HET75.hapFem.minimal.vcf > HNgenomic_combined.vcf

## bcftools view -s samplename variantfile.vcf > combinedoutputfile.vcf##

##Use GATK(Haplotypecaller) to call genotypes##
##First move/copy/download GenomeAnalysisTK.jar into a folder in your directory##

#!/bin/bash
#SBATCH --partition=prevail
#SBATCH --job-name="CallingSnp"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=6-23:59:59


java -jar /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/bellula/24228/dedupl.sorted.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --heterozygosity 0.02 -O /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/bellula/24228/outputNEW.g.vcf

##java -jar path/to/GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fa -I inputindividualfile.sorted.bam  --emitRefConfidence GVCF givesconfidencescores --variant_index_type LINEAR --variant_index_parameter 128000 bothvariantparametersforformatting --heterozygosity levelofdesiredhwterozygosity 0-1 -O outputfile.g.vcf##

##Genotype files to be able to filter by genotype (QUAL etc.). -allSites - keeps invariants also
##used to allow for larger smaple and increase power), Also requires you to have 
##GenomeAnalysisTK.jar in directory (server is not updated)##

#!/bin/bash
#SBATCH -w 'cruncher'
#SBATCH --partition=usual
#SBATCH --job-name="GenotypeAllsites"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

java -jar ~/bellula_tristero_genomic_data/GenomeAnalysisTK.jar -T GenotypeGVCFs \
        -R /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Reference/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa \
        -V ~/bellula_tristero_genomic_data/SNPcalled/CCydno_combined.g.vcf \
        -o CCydnogenotyped_invariants.g.vcf \
        -allSites

##Filter for biallelic, minimunm qualitty of 30, minimum depth of coverage 8, 
##%max missing individuals at any given site,remove indels, recode is specific to vcftools
##to write out, recode info all keeps tag information, stdout -specify output file name##

#!/bin/bash
#SBATCH -w 'cruncher'
#SBATCH --partition=usual
#SBATCH --job-name="vcftools"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

module load  vcftools/0.1.14

vcftools --vcf linaresi_new_genotypedWITHINVARIANT.g.vcf --min-alleles 1 --max-alleles 2 \
--minQ 30 --minDP 8 --max-missing 0.35 --remove-indels --recode \
--recode-INFO-all --stdout > Linvariants_filtered.g.vcf

##remove all */* small invariants because it interferes with admixture analysis###
awk ‘!/*/’ Linvariants_filtered.g.vcf > Linvariants_filteredNEW.g.vcf

##Files need to be bgzippped (bgzip) and tab indexed (tabix) for this step and next)##
bgzip Linvariants_filteredNEW.g.vcf
tabix Linvariants_filteredNEW.g.vcf.gz

##Order all files by chromosome rather than by scaffold (Aligned using reference genome which is 
scaffold organised)## 
#!/bin/bash

#SBATCH -w 'cruncher'
#SBATCH --partition=usual
#SBATCH --job-name="orderbychromes"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

python ~/bellula_tristero_genomic_data/Filtered/vcfChromTransfer.py \
-v /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Filtered/Linvariants_filteredNEW.g.vcf.gz \
-t /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Filtered/coordinates.txt \
>> /data/home/wolfproj/wolfproj-07/bellula_tristero_genomic_data/Filtered/lininvco_NEW.vcf.gz

##Create .geno file for all individuals (combines and converts), Martin (2019) uses this specific 
##format for the rest of the admixture analysis###

#!/bin/bash
#SBATCH -w 'cruncher'
#SBATCH --partition=usual
#SBATCH --job-name="vcftogenoMT"
#SBATCH --mail-user=Saoirse.McMahon@campus.lmu.de
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

python /data/home/wolfproj/wolfproj-07/ABBA_BABA_tut/data/genomics_general-master/VCF_processing/parseVCFs.py \
-i numata.vcf.gz \
-i cydno.vcf.gz \
-i bellinvco_NEW.vcf.gz \
-i lininvco_NEW.vcf.gz \
--threads 25 | bgzip > Belnumlincyd.geno.gz

##Calculate fd values for relevant sites -g .geno -f phased or unphased -o output -w window size (e.g 20kb) -m ?? 
##-s how much window slides by##
## -P1 sister allopatric -P2 sister sympatric -P3 sympatric not sister -O outgroup -T no.of threads 
##--popsFile file needs to be created with individual labels and what population they are e.g CAM40043	mel_bel 
##tab delimited 

#!/bin/bash
#SBATCH -w 'cruncher'
#SBATCH --partition=fat
#SBATCH --job-name="ABBABABA"
#SBATCH --mail-user=113473438@umail.ucc.ie
#SBATCH --mail-type=END
#SBATCH --time=1-00:00:00

python /data/home/wolfproj/wolfproj-07/ABBA_BABA_tut/data/genomics_general-master/ABBABABAwindows.py -g Belnumlincyd.geno.gz -f phased \
-o CLB50kb100m.csv -w 50000 -m 100 -s 5000 \
-P1 cydno -P2 lin -P3 mel_bel -O num -T 2 --popsFile Population_data.pop.txt