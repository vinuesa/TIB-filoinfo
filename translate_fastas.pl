#!/usr/bin/env perl
#: SCRIPT: translate_fastas.pl
#: Author: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:         http://www.ccg.unam.mx/~vinuesa/
#: Project start: November 2005

#: AIM: batch translate fastas found by file extension in the working directory
#: Assumptions: 
#      1. Bio::SeqIO is installed in the $PATH
#      2. Only DNA fasta sequences with the user-provided extension name 
#         are present in the working directory

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

# Standard modules
use File::Basename;
use Getopt::Std;

# BioPerl modules
use Bio::SeqIO;

# We will work only with fasta-formatted sequences
my $DEF_INPUT_FORMAT    = 'fasta';
my $DEF_OUTPUT_FORMAT   = 'fasta';

my $progname = basename($0); # translate_fastas.pl
my $VERSION = "0.1";

my (%opts,$input_file_extension,$input_format,$output_format, $output_file_extension, $translation_table);
getopts('hi:o:t:e:', \%opts);  # opt p could be a=alphabet [dna|prot]

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "\nusage: $progname version $VERSION [options]\n";
	print   "-h \t this message\n";
	print   "-e \t file extension with nucleotide sequences\n";
	print   "-t \t translation table                                     (optional, default 1)\n";

	print STDOUT <<EOF;
	
 AIM: batch translate fastas found by file extension in the working directory
 Assumptions: 
     1. Bio::SeqIO is installed in the \$PATH
     2. Only DNA fasta sequences with the user-provided extension name 
	are present in the working directory


Usage example:
    $progname -e fa -t 11
	      
EOF
	exit; 
}

# check args
if(defined($opts{'e'})){ $input_file_extension = $opts{'e'}; }
else{ die "# $0 : need a nucleotide file extension\n" }
if(defined($opts{'t'})){ $translation_table = $opts{'t'}; }
else{ $translation_table = 1; }


my @files = < *$input_file_extension >;

foreach my $input_seq_file ( @files ) 
{
   translate_fastaDNA2fastaAA($input_seq_file,$translation_table, $DEF_INPUT_FORMAT, $DEF_OUTPUT_FORMAT);
}

######################################################################################################################

sub translate_fastaDNA2fastaAA 
{
    # takes a dna sequence and translates it to protein using a passed translation table
    # from the list of 16 available in Bioperl 
    # http://www.bioperl.org/Core/Latest/bptutorial.html
    # http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tools/CodonTable.html

    my($dna_fasta, $translation_table, $input_format, $output_format) = @_;
    my ($basename,$ext) = split(/\./, $dna_fasta);
    my $outfile = $basename . '_translated.faa';

    my $sequin;
    my %seqs;
    my $prot;	    
    my $seqin  = Bio::SeqIO->new( -file   => $dna_fasta, 
    				  -format => $input_format);

    my $seqout = Bio::SeqIO->new( -file   => ">$outfile",
    				  -format => $output_format);

    # process each sequence, translating them with the corresponding codon table


    while( my $seq = $seqin->next_seq()) 
    {
	my $pseq = $seq->translate(-codontable_id => $translation_table);
	$seqout->write_seq($pseq);
    }	
    
    # check for translated file
    if(-s $outfile)
    {
       print "> generated file $outfile\n";
    }
    else
    {
       print "> ERROR: could not generate file $outfile\n";
    }
}
