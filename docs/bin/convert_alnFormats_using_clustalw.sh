#!/usr/bin/env bash

#: AUTHOR: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:         http://www.ccg.unam.mx/~vinuesa/
#: AIM: a very simple bash script to teach basic shell scripting to my students at the
#:      Bachelor's Program in Genomic Scinces, UNAM-Mexico (http://www.lcg.unam.mx). 
#:      It calls clustalw to perform sequence format conversions.

# Capture the program name
progname=$(basename $0)
VERSION=0.1

# save command line arguments to named variables
infile=$1
outformat=$2

# define the help function
function help
{
    # This is a "here document"; very convenient for this matter: 
    # just print text interpolating variables
    cat << EOF

    usage for $progname v.$VERSION:
    
    $progname <input_sequence_file2convert> <outfile_format>

    example: $progname my_alignment.fasta phylip
             (converts fasta file to phylip file)

    NOTE: Clustalw can only read: clustalw and fasta formats
                    It can write: clustalw fasta nexus phylip
EOF

}

# Check that the user provides the expected number of arguments to the script or die
if [ $# -ne 2 ]
then
      # print help and exit
      help
      exit
fi

# Prepare the output filename using some basic variable string manipulations 
if [ "$outformat" = clustalw ]
then
   short_format=aln
else   
   short_format=${outformat:0:3}
fi

outfile=${infile%.*}.$short_format

# This is the clustalw command that does the actual format conversion; 
# work silently: send STOUT and STDERR to bit bucket
clustalw -infile=$infile -convert -output=$outformat -outfile=$outfile &> /dev/null

if [ -s $outfile ]
then
    echo
    echo "# Converted file $infile to $outfile in $outformat format!"
    echo

else
    echo
    echo "# ERROR: the expected output file $outfile was not produced!"
    echo
fi
