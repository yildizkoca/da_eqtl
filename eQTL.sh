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
mkdir $deji/jobs/permutations1
mkdir $deji/stats/permutations1
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations1;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations1/permutations1_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations1chunk$j.sh;
qsub permutations1chunk$j.sh;
done

#permutations2 (wo covariates, w phenotype correlation structure)
mkdir $deji/jobs/permutations2
mkdir $deji/stats/permutations2
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations2;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --permute 1000 --grp-best --chunk $j 50 --out $yk/deji/stats/permutations2/permutations2_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations2chunk$j.sh;
qsub permutations2chunk$j.sh;
done

#permutations3 (w known covariates)
mkdir $deji/jobs/permutations3
mkdir $deji/stats/permutations3
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations3;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/counts/covariates_full.txt --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations3/permutations3_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations3chunk$j.sh;
qsub permutations3chunk$j.sh;
done

#permutations4 (w known covariates, w phenotype correlation structure)
mkdir $deji/jobs/permutations4
mkdir $deji/stats/permutations4
for j in $(seq 1 50); do
cd $yk/deji/jobs/permutations4;
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $deji/counts/covariates_full.txt --permute 1000 --grp-best --chunk $j 50 --out $yk/deji/stats/permutations4/permutations4_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh;
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations4chunk$j.sh;
qsub permutations4chunk$j.sh;
done

#permutations for phenotype PC(i-1) (wo known covariates)
for i in $(seq 2 1 21); do
mkdir $yk/deji/jobs/permutations_pc$i
mkdir $yk/deji/stats/permutations_pc$i
cd $yk/deji/jobs/permutations_pc$i
  for j in $(seq 1 50); do
    echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/pca/PCs/PC$i.pca --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations_pc$i/permutations_pc$i"_"$j"_50.txt" --normal --seed 123456" > jobchunk$j.sh;
    (tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations_pc$i"chunk"$j.sh;
    qsub permutations_pc$i"chunk"$j.sh
  done
done

#permutations for phenotype PC(i-1) (w known covariates)
for i in $(seq 2 1 21); do
mkdir $yk/deji/jobs/permutations_cov_and_pc$i
mkdir $yk/deji/stats/permutations_cov_and_pc$i
cd $yk/deji/jobs/permutations_cov_and_pc$i
  for j in $(seq 1 50); do
    echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/pca/cov_and_PCs/cov_and_PC$i.pca --permute 1000 --chunk $j 50 --out $yk/deji/stats/permutations_cov_and_pc$i/permutations_cov_and_pc$i"_"$j"_50.txt" --normal --seed 123456" > jobchunk$j.sh;
    (tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > permutations_cov_and_pc$i"chunk"$j.sh;
    qsub permutations_cov_and_pc$i"chunk"$j.sh
  done
done



#nominals1 (wo covariates)
mkdir $deji/jobs/nominals1
mkdir $deji/stats/nominals1
for j in $(seq 1 50); do
cd $yk/deji/jobs/nominals1
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --nominal 0.05 --chunk $j 50 --out $yk/deji/stats/nominals1/nominals1_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > nominals1chunk$j.sh;
qsub nominals1chunk$j.sh;
done

#nominals3 (w known covariates)
mkdir $deji/jobs/nominals3
mkdir $deji/stats/nominals3
for j in $(seq 1 50); do
cd $yk/deji/jobs/nominals3
echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/counts/covariates_full.txt --nominal 0.05 --chunk $j 50 --out $yk/deji/stats/nominals3/nominals3_$j\_50.txt --normal --seed 123456" > jobchunk$j.sh
(tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > nominals3chunk$j.sh;
qsub nominals3chunk$j.sh;
done

#nominals for phenotype PC(i-1) (wo known covariates)
for i in $(seq 2 1 21); do
mkdir $yk/deji/jobs/nominals_pc$i
mkdir $yk/deji/stats/nominals_pc$i
cd $yk/deji/jobs/nominals_pc$i
  for j in $(seq 1 50); do
    echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/pca/PCs/PC$i.pca --nominal 0.05 --chunk $j 50 --out $yk/deji/stats/nominals_pc$i/nominals_pc$i"_"$j\_50.txt --normal --seed 123456" > jobchunk$j.sh
    (tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > nominals_pc$i"chunk"$j.sh;
    qsub nominals_pc$i"chunk"$j.sh
  done
done


#nominals for phenotype PC(i-1) (w known covariates)
for i in $(seq 2 1 21); do
mkdir $yk/deji/jobs/nominals_cov_and_pc$i
mkdir $yk/deji/stats/nominals_cov_and_pc$i
cd $yk/deji/jobs/nominals_cov_and_pc$i
  for j in $(seq 1 50); do
    echo "QTLtools cis --vcf $yk/deji/fcvcfs/filtered_cleaned_merged_vcfs.vcf.gz --bed $yk/deji/counts/corrected_fctdata_cpm_wgsid_zerofiltered90p.sorted.bed.gz --exclude-samples $yk/deji/counts/file.exc --cov $yk/deji/pca/cov_and_PCs/cov_and_PC$i.pca --nominal 0.05 --chunk $j 50 --out $yk/deji/stats/nominals_cov_and_pc$i/nominals_cov_and_pc$i"_"$j\_50.txt --normal --seed 123456" > jobchunk$j.sh
    (tail -n+1 $yk/deji/jobs/jobheader.sh && cat jobchunk$j.sh) > nominals_cov_and_pc$i"chunk"$j.sh;
    qsub nominals_cov_and_pc$i"chunk"$j.sh
  done
done






     

