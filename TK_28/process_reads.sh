#!/usr/bin/bash

# Load modules
module load bowtie bbtools seqtk blast

#[+] Loading bowtie  2-2.4.5
#[+] Loading samtools 1.15.1  ...
#[+] Loading bbtools  38.96
#[+] Loading seqtk  1.3
#[-] Unloading blast 2.13.0+  ...
#[+] Loading blast 2.13.0+  ...

for i in 1 2 3 4
do        
    SAMPLE=SB-7773_${i}_S${i} 
    echo ********** Processing sample ${SAMPLE} ********** 
    # Merged paired-end reads
    echo -1/4 Merging reads
    bbtools bbmerge in1=./results/01trim/${SAMPLE}/${SAMPLE}_L001_R1_001_val_1.fq.gz  in2=./results/01trim/${SAMPLE}/${SAMPLE}_L001_R2_001_val_2.fq.gz out=${SAMPLE}.merged.fastq  outu1=${SAMPLE}.unmerged.1.fastq  outu2=${SAMPLE}.unmerged.2.fastq adapter=adapters.fasta

    # Adding R1 unpaired and unmerged-1 to fastq file
    gunzip -c ./results/01trim/${SAMPLE}/${SAMPLE}_L001_R1_001_unpaired_1.fq.gz >> ${SAMPLE}.merged.fastq
    cat ${SAMPLE}.unmerged.1.fastq >> ${SAMPLE}.merged.fastq

    # Convert fastq to fasta
    echo -2/4 Converting fastq files to fasta
    seqtk seq -A ${SAMPLE}.merged.fastq > ${SAMPLE}.merged.fasta

    # Prepare blast DB
    echo -3/4 Making blast DB
    makeblastdb -dbtype nucl -in ${SAMPLE}.merged.fasta

    # Run blastn 
    echo -4/4 Running blastn
    CMD=$(blastn -max_target_seqs 10000000 -subject_besthit -query vector_seq.fasta -db ${SAMPLE}.merged.fasta -evalue 1e-10 -dust no -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out ${SAMPLE}.merged.bn)
    echo $CMD
    $CMD

    # Run bowtie2 with giRNAs (it would be possible to use -L 10 -N 1 --no-1mm-upfront to use 10bp seeds with up to 1 mismatch)
    bowtie2-build ${SAMPLE}.merged.fasta ${SAMPLE}.merged
    bowtie2 -x ${SAMPLE}.merged  -f -U crisper_probes.fasta --end-to-end --norc -a |samtools view -F 4 > ${SAMPLE}.sam
done

# Generate bed files with hit coords
echo Creating BED files
for i in 1 2 3 4
do
    grep '^Middle' SB-7773_${i}_S${i}.merged.bn|perl -F'\t' -lane 'print "$F[1]\t$F[8]\t$F[9]\t$F[0]"' > SB-7773_${i}_S${i}.middle.bed
    grep '^3_prime' SB-7773_${i}_S${i}.merged.bn|perl -F'\t' -lane 'print "$F[1]\t$F[8]\t$F[9]\t$F[0]"' > SB-7773_${i}_S${i}.3_prime.bed
    #grep '^giRNA' SB-7773_${i}_S${i}.merged.bn|perl -F'\t' -lane 'print "$F[1]\t$F[8]\t$F[9]\t$F[0]"' > SB-7773_${i}_S${i}.giRNA.bed
done

# Generate fasta files with UMI and giRNA seqs
echo Generating fasta files
for i in 1 2 3 4
do
    seqtk subseq SB-7773_${i}_S${i}.merged.fasta SB-7773_${i}_S${i}.3_prime.bed > SB-7773_${i}_S${i}.3_prime.fasta
    seqtk subseq SB-7773_${i}_S${i}.merged.fasta SB-7773_${i}_S${i}.middle.bed > SB-7773_${i}_S${i}.middle.fasta
    #seqtk subseq SB-7773_${i}_S${i}.merged.fasta SB-7773_${i}_S${i}.giRNA.bed > SB-7773_${i}_S${i}.giRNA.fasta
done

for i in 1 2 3 4
do
    perl ./process_blastn_tab_files.pl -b SB-7773_${i}_S${i}.merged.bn -s SB-7773_${i}_S${i}.sam -f SB-7773_${i}_S${i}.merged.fasta > process_blastn_tab_files.SB-7773_${i}_S${i}.txt
    grep RESULT_SEQ process_blastn_tab_files.SB-7773_${i}_S${i}.txt|perl -ne 'chomp;@x=split /\t/;next if $x[6] != 13 || $x[7] != 20;print "$x[2]:$x[3]:$x[5]\n"'|sort |uniq -c|sort -k1n > process_blastn_tab_files.SB-7773_${i}_S${i}_CLUSTERED.txt
done

module purge

