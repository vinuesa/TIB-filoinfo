#!/usr/bin/env perl
# Author: Pablo Vinuesa
# split a single multi-fasta file into multiple single-sequence FASTA files

my $acc='';

while(<>){
    if(/>/){
       ($acc)=split(/\s+/);
       $acc=~s/[^A-Za-z0-9-_]//g;
       $acc.=".fa";
       close OUT;
       open(OUT,">$acc");
    }
    print OUT $_;
}

