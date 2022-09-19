#!/home/lorenziha/data//miniconda3/bin/perl
use strict;

my $usage = "$0 -b <blastn.tab> -s <bowtie2 sam file> -f <fasta reads>\n\n";
# input STDIN blastn+ tab file with the following fields:
# query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject length

my %arg = @ARGV;
die $usage unless $arg{-b} && $arg{-s};

my (%prime3, %middle, %readids, %giRNAs );

# BlastN of vector sequences against reads to extract UMI and giRNA coords within reads
open (BLAST, "<$arg{-b}"); 
while(<BLAST>){
    chomp;
    next if m/^#/; # skip comment lines
    my @x = split /\t/;
    my ($qid, $readid, $q_end5, $q_end3, $s_end5,$s_end3) = ($x[0], $x[1], $x[6], $x[7], $x[8], $x[9]);
    if ($qid eq '3_prime'){
        push @{$prime3{$readid}}, $s_end5, $s_end3, $s_end3 + 1 - $s_end5;
    }
    elsif ($qid eq 'Middle'){
        push @{$middle{$readid}}, $s_end5, $s_end3, $s_end3 + 1 - $s_end5;
    }
    $readids{$readid}++;
}
close BLAST;

#Bowtie 2 alignment of giRNAs vs reads to identify giRNAid-ReadID pairs 
open (SAM, "<$arg{-s}");
while(<SAM>){
    chomp;
    my @x = split /\t/;
    my ($giRNAid, $readid, $match) = ($x[0], $x[2], $x[5]);
    push @{$giRNAs{$readid}}, $giRNAid, $match; # I am including up to 1 mismatch  per giRNA
    print "SAM: ($giRNAid, $readid, $match)\n";
}
close SAM;

# Load read seqs from fasta file
my (%seqs, $readid); 
open (FASTA, "<$arg{-f}");
while(<FASTA>){
    chomp;
    if (m/^>(\S+)/){
        $readid = $1;
        print "SEQ:$readid\t"
    } elsif (m/^([ATGCN]{50,})$/i){
        $seqs{$readid} = $1;
        print "$seqs{$readid}\n";
    }
}
close FASTA;

# Estract giRNA and UMI coords within each read, generates UMI and giRNA sequences from their coords in the read and assign giRNA IDs to reads  
for my $readid (keys %readids){
    if ($middle{$readid} && $prime3{$readid}){
        # Both hits exist (likely have giRNA and UMI
        my $umi_len = $middle{$readid}[0] - 1;
        my $middle_len = $middle{$readid}[1] + 1 - $middle{$readid}[0];
        my $girna_len = $prime3{$readid}[0] - $middle{$readid}[1];
        my $girna_id = $giRNAs{$readid}[0] || 'NO_giRNA_HIT';
        my $girnamatch = $giRNAs{$readid}[1] || 'NO_giRNA_HIT';
        print "RESULT: $readid <$umi_len> $middle_len <$girna_len> giRNA=$girna_id ($girnamatch)\n";

        # Get potential UMI and giRNA seqs from read fasta file
        my $umi_seq   = substr($seqs{$readid},0 ,$umi_len);
        my $girna_seq = substr($seqs{$readid},$middle{$readid}[1] + 3, $girna_len - 4);
        print "RESULT_SEQ\t$readid\t$umi_seq\t$girna_seq\t$seqs{$readid}\t$girna_id\t$umi_len\t".($girna_len - 4)."\n";
    }
}
