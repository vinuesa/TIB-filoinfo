#!/usr/bin/env perl

# (c) Pablo Vinuesa & Bruno Contreras-Moreira
# CCG-UNAM; http://www.ccg.unam.mx/~vinuesa/
# AAaln2cdnAln.pl Sep 25th, 2009

# reverse translates a protein alignment, given the DNA seqs, to generate the underlying codon alignment!
# use warnings;
use strict;
use File::Basename;

my $progname = basename($0);
my $version  = 0.2; # 2022-11-10; added NOTE in help and checks if cdn_aln was written to disk
                    # v0.1Sep 25th, 2009;

die "\n# $progname V.$version needs two arguments: 
     1) the dna input sequence to align; 
     2) the corresponding protein alignment

     NOTE: requires that the sequences in the protein alignment are in the same order as in the input fasta\n\n" unless @ARGV == 2;

my ($dna_input,$prot_aln) = @ARGV;

my $cdnAln=(convert_AAaln2NTaln($dna_input,$prot_aln));

if (-s $cdnAln){ print "# wrote $cdnAln\n";}
else { print STDERR " ERROR: the expected codon alignment file $cdnAln was not produced !!!\n";}

sub convert_AAaln2NTaln
{
         my ($DNAfile,$AA_aln_file) = @_;

	 my ($basename,$ext) = (split(/\./, $DNAfile))[0,1];

         my $outfile = $basename . "_cdnaln.$ext"; 
	 my ($seq,$sequence,$aa,$n_of_aa,$aligned_nts)=('','','','','',);
	 my ($align_positions,$FASTA,@length,$l);
	 
         # 1) read AA aligned fasta sequences
	 my %AAalign = read_FASTA_sequence($AA_aln_file);

         # 2) read NT fasta sequence
	 my %NTsequences = read_FASTA_sequence($DNAfile);
	 
         # 3) replicate alignment	 
	 foreach $seq (sort {$a<=>$b} (keys(%NTsequences)))
	 {
		  $FASTA .= $NTsequences{$seq}{'NAME'};
		  
		  # 3.1) dissect codons
		  my @codons;
		  $sequence = $NTsequences{$seq}{'SEQ'};
		  while($sequence)
		  {
		   	push(@codons,substr($sequence,0,3));
			$sequence = substr($sequence,3);
		  }
		  
		  # 3.2) loop through aminoacids in $AAalign{$seq}
		  $sequence = $AAalign{$seq}{'SEQ'};
		  $n_of_aa = length($sequence);
		  $align_positions=0;
	 	  $aligned_nts = '';
		  for($aa=0;$aa<=$n_of_aa;$aa++)
		  {
		   	if(substr($sequence,$aa,1) ne '-')
			{
				 $aligned_nts .= $codons[$align_positions];
				 $align_positions++;
			}
			else
			{
				 $aligned_nts .= '---';
			}
		  } 	
	 
	 	  $l = length($aligned_nts);
	 	  if(!grep(/$l/,@length)){ push(@length,$l); }
		  
		  $FASTA .= $aligned_nts."\n";
	 }
	 
	 if(scalar(@length) > 1)
	 {
	 	  print "# convert_AAaln2NTaln : warning, sequence lengths differ\n";
	 	  
	 }
	 
         # write the output FASTA file
	 open(FASTA,">$outfile") || die "# convert_AAaln2NTaln : cannot write to $outfile\n";
	 print FASTA $FASTA;
	 close FASTA;
    
	 return $outfile;
}

sub read_FASTA_sequence
{
    my ($infile) = @_;
	 
    my (%FASTA,$name,$seq,$n_of_sequences);

    $n_of_sequences = 0;
    open(FASTA,$infile) || die "# read_FASTA_sequence: cannot read $infile\n";
    while(<FASTA>)
    {
    	     if(/^\>/)
	     {
	   	   $name = $_; 
	   	   $n_of_sequences++;
	   	   $FASTA{$n_of_sequences}{'NAME'} = $name;
	     }  		   
	     else
	     {
	   	   $FASTA{$n_of_sequences}{'SEQ'} .= $_;
	   	   $FASTA{$n_of_sequences}{'SEQ'} =~ s/[\s|\n]//g;
	     }
    }
    close(FASTA);

    return %FASTA;
}
