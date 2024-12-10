#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#: AUTHOR: Pablo Vinuesa, CCG-UNAM, Mexico; https://www.ccg.unam.mx/~vinuesa/; @pvinmex
#: AIM: extract a user-defined number of N- or C-terminal residues from a multi-FASTA file,
#:    and compute basic stats and counts for a certain residue. 
#:    This is useful, for example, to verify and analyze the distribution of cysteines
#:      among the five C-terminal residues in Rab GTPases

my $progname = basename($0); # extract_N-or_C-terminal_regions_and_compute_counts.pl
my $version = '0.3.1_2024-12-09'; # v0.3.1_2024-12-09 added missing # before comment following version; 
                                # v0.3_2024-11-05; corrected version number and improved description of the script's aim and output
                                # v0.2_2024-11-02; added single letter options
                                # v0.1_2024-11-01 first commit

# Variables for command-line options
my ($fasta_file, $output_prefix, $residue_count, $region, $count_residue, $help);

# Getopt::Long for command-line options
GetOptions(
    'fasta|f=s'        => \$fasta_file,
    'output|o=s'       => \$output_prefix,
    'residues|n=i'     => \$residue_count,
    'region|r=s'       => \$region,
    'count_residue|a=s' => \$count_residue,
    'help'           => \$help,
);

# Display help menu
if ($help || !$fasta_file || !$output_prefix || !$residue_count || !$region || !$count_residue) {
    print <<"HELP";
Usage for $progname v$version:

$progname --fasta <input.fasta> --output <output_prefix> --residues <num_residues> --region <N|C> --count_residue <residue>

Options:
    --fasta|-f          Input multi-FASTA file of protein sequences
    --output|-o         Output file prefix (for generated FASTA and tabular files)
    --residues|-n       Number of N-terminal or C-terminal residues to extract
    --region|-r         Region to extract: 'N' for N-terminal or 'C' for C-terminal
    --count_residue|-a  Residue to count in the extracted region (e.g., M for methionine, C for cysteine)
    --help              Show this help message

AIM:
   Extracts a specified number of N- or C-terminal residues from each protein sequence in a multi-FASTA file,
     computes basic overal residue frequencies, and counts the occurrences of a certain focal amino acid. 
   This is useful, for example, to analyze the distribution of cysteines found among the five C-terminal 
    residues of Rab GTPases, which require at least one cysteine as a prennylation site.

OUTPUT:
    - multi-FASTA file with the extracted N- or C-terminal regions
    - A three-col table with the seq_id, sequence, counts for user-defined residue in the extracted sequence
    - A Matrix with overall proportions of each amino acid at each site of the extracted sequences 

HELP
    exit;
}

# Validate the region option
die "Invalid region option. Use 'N' for N-terminal or 'C' for C-terminal.\n" unless $region =~ /^[NC]$/i;

# Define proteinogenic amino acids
my @amino_acids = qw(A R N D C Q E G H I L K M F P S T W Y V);

# Open input and output files
open my $in_fh, '<', $fasta_file or die "Could not open '$fasta_file': $!";
open my $fasta_out, '>', "${output_prefix}_${region}term.fasta" or die "Could not write to FASTA output file: $!";
open my $tab_out, '>', "${output_prefix}_${region}term.tsv" or die "Could not write to TSV output file: $!";
open my $freq_out, '>', "${output_prefix}_freq_table_${region}term.tsv" or die "Could not write to frequency table output file: $!";

my %frequency_table;
my %residue_counts;
my %accessions;
my $sequence_count = 0;

# Process each sequence in the FASTA file
{
    local $/ = "\n>"; # FASTA record separator
    while (my $record = <$in_fh>) {
        $record =~ s/^>//;
        my ($header, @seq_lines) = split /\n/, $record;
        my $sequence = join '', @seq_lines;
        $sequence =~ s/\s//g;  # Remove any whitespace
        $sequence =~ s/>//g; #Clean any potential '>' characters from the end of the sequence

        # Extract accession from the header
        my ($accession) = split /\s/, $header;
        $accessions{$sequence_count} = $accession;

        # Extract specified region (N- or C-terminal)
        my $term_seq;
        if ($region =~ /N/i) {
            $term_seq = substr($sequence, 0, $residue_count);
        } else {
            # Ensure exactly the requested number of C-terminal residues are extracted
            $term_seq = substr($sequence, -$residue_count) if length($sequence) >= $residue_count;
        }

        # Count the specified residue in the extracted region
	#print "DEBUG: count_residue: $count_residue for term_seq: $term_seq\n";
        #my $residue_count_specific = ($term_seq =~ tr/$count_residue/$count_residue/);
	#my $residue_count_specific = $term_seq =~ tr/$count_residue//;
	my $residue_count_specific = () = $term_seq =~ /$count_residue/g;
	#print "DEBUG: residue_count_specific: $residue_count_specific\n";
        $residue_counts{$sequence_count} = $residue_count_specific;

        # Write to FASTA output
        print $fasta_out ">$accession\n$term_seq\n";

        # Write to tabular output
        print $tab_out "$accession\t$term_seq\t$residue_count_specific\n";

        # Update frequency table for proteinogenic amino acids
        for my $i (0 .. length($term_seq) - 1) {
            my $residue = substr($term_seq, $i, 1);
            $frequency_table{$i}{$residue}++ if grep { $_ eq $residue } @amino_acids;
        }
        $sequence_count++;
    }
}

# Print frequency table for each amino acid at each site of the extracted sequences
print $freq_out "Position\t", join("\t", @amino_acids), "\n";
for my $seq_idx (0 .. 0) {
    for my $pos (0 .. $residue_count - 1) {
        print $freq_out "\t$pos";
        for my $residue (@amino_acids) {
            my $count = $frequency_table{$pos}{$residue} // 0;
	    my $freq = sprintf("%.2f", $count/$sequence_count); 
            print $freq_out "\t$freq";
        }
        print $freq_out "\n";
    }
}

# Close files
close $in_fh;
close $fasta_out;
close $tab_out;
close $freq_out;

print "Processing complete. Results saved with prefix '$output_prefix'.\n";
