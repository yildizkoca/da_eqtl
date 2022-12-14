#!/bin/bash 

#PBS -l qos=biocore
#PBS -l nodes=1:ppn=15
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -l mem=20gb

module load gcc/6.2.0
module load samtools/1.10
module load bcftools/1.16.0
module load tabix/0.2.6
module load qtltools/1.1
module load bedtools
module load htslib/1.10.2


awk '{N=0.0;for(i=7;i<=NF;i++) if($i=="0") N++; if(N/(NF-6)<0.9) print;}' corrected_fctdata_cpm_wgsid.bed > corrected_fctdata_cpm_wgsid_zerofiltered90p.bed

#sort bed file by chrom number and start position
sort -k 1,1 -k2,2n corrected_fctdata_cpm_wgsid_zerofiltered90p.bed > corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed 
bgzip corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed
tabix -p bed corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz

#run pca on phenotypes
QTLtools pca --bed corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --scale --center --out $yk/deji/pca/phenotypes_cpm_wgsid_zerofiltered90p

#run pca on genotypes
QTLtools pca --vcf $yk/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --scale --center --maf 0.05 --distance 50000 --out genotypes


#permutations1 (wo covariates)
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations1;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations1/permutations1_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations1chunk$j.sh;
qsub permutations1chunk$j.sh;
done

#permutations2 (wo covariates, w phenotype correlation structure)
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations2;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --permute 1000 --grp-best --chunk $j 50 --out $yk/deji/stats/permutations2/permutations2_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permumations2chunk$j.sh;
qsub permumations2chunk$j.sh;
done

#permutations3 (w known covariates)
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations3;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/counts/tcovariates.txt.gz --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations3/permutations3_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permumations3chunk$j.sh;
qsub permumations3chunk$j.sh;
done

#permutations4 (w known covariates, w phenotype correlation structure)
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations4;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/counts/tcovariates.txt.gz --permute 1000 --grp-best --chunk $j 50 --out $yk/deji/stats/permutations4/permutations4_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations4chunk$j.sh;
qsub permutations4chunk$j.sh;
done



     

