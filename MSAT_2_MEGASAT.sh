#!/bin/bash
##REQUIREMENTS:
#Using this script requires the following files in the same directory as the bash script:
#1. Msatcommander primer and microsats output files with "primers.txt" and "sats.txt" as final suffixes on each,  with the information for each locus on the same line for each file  - BEDTOOLS will throw errors if this is not the case

#2. The original genome or sequence file used by Msatcommander to extract microsatellite information in fasta format

##AND
#A version of bedtools downloaded and installed on the Desktop: http://bedtools.readthedocs.io/en/latest/

##TO USE:
#1.Move the script, the genome/sequence and a single pair of primer and microsat files from msatcommander in their own directory.

#2. Ensure the script has permission to run, enter on the command line: chmod +x MSAT_2_MEGASAT.sh

#3. Setup genome variable by entering on the command line: export genome=[genome file name.fasta, fna, fai, etc.]

#4. Run by entering on the command line: ./MSAT_2_MEGASAT.sh

#Set file prefix as mpref
#make bed file for complete primer sequence in bed format  chrom, start, end, and name.
bedtools2=$(find ~/Desktop -name bedtools2)/bin/bedtools
primers=*primers.txt
sats=*sats.txt

#Step 1 Print complete amplified region sequnece
#Get right primer sequence end, complete region end
awk '{ print $14 }' $primers | tr ',' $'\t' > Rightprime
awk '{print $1+$2}' Rightprime > Rightend

#Get left primer sequence beginning
awk '{ print $5 }' $primers | tr ',' $'\t' > Leftprime
awk '{print $1}' Leftprime > Leftstart

#Make a bed input file to output complete sequence
#"Micro" can be replaced with a species name or abbreviation or whatever
awk '{print  "Micro-"$1"-"$3}' $sats > Satnames
awk '{print  $1}' $sats  > Chrom_name
paste -d '\t' Chrom_name Leftstart Rightend Satnames > Compseq.bed

#Get complete sequence
$bedtools2 getfasta -fi $genome -bed Compseq.bed  -fo Comp_seqs -tab
paste -d '\t' Comp_seqs Satnames > Compseq_MSats.txt

#Get 5' and 3' flank sequences
#5'
awk '{print  $5}' $sats > Repstart
awk '{print $1+$2}' Leftprime > Primend
paste -d '\t' Chrom_name Primend Repstart Satnames > 5prime_raw
awk '$2==$3 {$3=$3+1}1' OFS='\t' 5prime_raw > 5prime.bed 
awk '($2==$3)' 5prime_raw | awk '{print $4}' > 5prime_noflank
$bedtools2 getfasta -fi $genome -bed 5prime.bed -fo 5prime -tab

#3'
awk '{print  $6}' $sats > Repend
awk '{print $1}' Rightprime > Rightstart
paste -d '\t' Chrom_name Repend Rightstart Satnames > 3prime_raw
awk '$2==$3 {$3=$3+1}1' OFS='\t' 3prime_raw > 3prime.bed 
awk '($2==$3)' 3prime_raw | awk '{print $4}' > 3prime_noflank
$bedtools2 getfasta -fi $genome -bed 3prime.bed  -fo 3prime -tab


#Get left primer
awk '{ print $6 }' $primers > Leftprimseq
#Reverse complement of right primer 
awk '{ print $15 }' $primers | perl -ne '$seq = $_; $seq =~ tr/atcgATCG/tagcTAGC/;print reverse( $seq )' > Rightprim_rcomp

#Get repeat
#Add repeat length to sequence
awk '{ print $5, length($4)+$5 }' OFS='\t' $sats > Rep_startstop
paste -d '\t' Chrom_name Rep_startstop Satnames > Reps.bed
$bedtools2 getfasta -fi $genome -bed Reps.bed -fo Reps -tab

#Make Megasat Input File
awk '{print $2}' 5prime > 5flank
awk '{print $2}' 3prime > 3flank
awk '{print $2}'  Reps > Rep
paste -d '\t' Satnames Leftprimseq Rightprim_rcomp 5flank 3flank Rep > MegaSat_raw

#Replace missing flanks with X
awk 'NR==FNR{ a[$1]; next }$1 in a{ $4="X"}1' 5prime_noflank MegaSat_raw  > MegaSat_lil_better
awk 'NR==FNR{ a[$1]; next }$1 in a{ $5="X"}1' 3prime_noflank  MegaSat_lil_better > MegaSat_Final.txt

#Clear out the directory
rm Rightprime
rm Rightend
rm Leftprime
rm Leftstart
rm Satnames
rm Chrom_name
rm Compseq.bed
rm Comp_seqs
rm Repstart
rm Primend
rm 5*
rm 3*
rm Repend
rm Rightstart
rm Leftprimseq
rm Rightprim_rcomp
rm Rep_startstop
rm Reps.bed
rm Reps
rm Rep
rm MegaSat_raw
rm MegaSat_lil_better

echo Done!
