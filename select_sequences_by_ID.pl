#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#: AUTHOR: Pablo Vinuesa, CCG-UNAM, Mexico; https://www.ccg.unam.mx/~vinuesa/; @pvinmex
#: AIM: select [include|exclude] sequences listed in --ids IDs.list from the
#:      user-provided multi-FASTA file, and write the filtered FASTA file to disk.


my $progname = basename($0); # select_sequences_by_ID.pl
my $version = '0.2_2024-11-01'; # V0.2_2024-11-05, improved description of the script's aim
                                # v0.1_2024-11-02; first commit

# Variables for command-line options
my ($fasta_file, $id_file, $include, $exclude, $output_file, $help);

# Getopt::Long interface
GetOptions(
    'fasta|f=s'   => \$fasta_file,
    'ids|l=s'     => \$id_file,
    'include|i'   => \$include,
    'exclude|e'   => \$exclude,
    'output|o=s'  => \$output_file,
    'help|-h'     => \$help,
) or die "Incorrect usage!\n";

# Display help if requested or if necessary parameters are missing
if ($help || !$fasta_file || !$id_file || !$output_file || (!$include && !$exclude)) {
print << "END_HELP";
Usage for $progname v$version:

$progname --fasta <input_fasta> --ids <id_file> --output <output_fasta> <--include | --exclude> [--help]

OR

$progname -f <input_fasta> -l <id_file> -o <output_fasta> <-i | -e> [--help|-h]


Options:
    --fasta|-f     (Path to the) multi-FASTA file containing sequences
    --ids|-l       (Path to the) file containing sequence identifiers (one per line)
    --output|-o    (Path to the) output file for the filtered sequences
    --include|-i   Include only the sequences listed in the ID file
    --exclude|-e   Exclude the sequences listed in the ID file
    --help|-h      Show this help message

Note:
    Either --include or --exclude must be specified, but not both.
    
AIM: select [include|exclude] sequences listed in --ids IDs.list from the
     user-provided multi-FASTA file, and write the filtered FASTA file to disk.
         
END_HELP
    exit;
}

# Read sequence IDs from the ID file into a hash for fast lookup
open my $id_fh, '<', $id_file or die "Could not open ID file '$id_file': $!\n";
my %ids;
while (<$id_fh>) {
    chomp;
    $ids{$_} = 1;
}
close $id_fh;

# Open the FASTA file and output file
open my $fasta_fh, '<', $fasta_file or die "Could not open FASTA file '$fasta_file': $!\n";
open my $out_fh, '>', $output_file or die "Could not open output file '$output_file': $!\n";

# Variables to store sequence data
my $write_seq = 0;
my $header = '';
my $sequence = '';

# Process the FASTA file
while (<$fasta_fh>) {
    chomp;
    if (/^>(\S+)/) {
        # Print the previous sequence if applicable
        if ($header && $write_seq) {
            print $out_fh "$header\n$sequence\n";
        }
        
        # Reset for the new sequence
        $header = $_;
        $sequence = '';
        
        # Determine if we should include or exclude this sequence
        my $seq_id = $1;
        if ($include) {
            $write_seq = exists $ids{$seq_id};
        } elsif ($exclude) {
            $write_seq = !exists $ids{$seq_id};
        }
    } else {
        # Append sequence lines
        $sequence .= $_;
    }
}

# Print the last sequence if applicable
if ($header && $write_seq) {
    print $out_fh "$header\n$sequence\n";
}

# Close file handles
close $fasta_fh;
close $out_fh;

print "Filtering complete. Results saved to '$output_file'.\n";
