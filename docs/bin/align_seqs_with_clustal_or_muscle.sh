#!/usr/bin/env bash

# Author: Pablo Vinuesa; 
# January 14th, 2011

# program: align_seqs_with_clustal_or_muscle.sh
# sample shell script written for TLEM11 
# to teach basic bash scripting constructs

version=1.2 # Sept 17th, 2012

#>>> functions <<<
# 1)function to check that the required binaries are in $PATH
function check_binaries
{
   # the script assumes that clustalw and muscle are in found in $PATH; lets check it
  for program in clustalw muscle
  do
     bin=$(type -P $program )
     
     if [ -z $bin ]
     then
        echo "# $0 ERROR: can't find $program in \$PAH:"
        echo $PATH
     else 
        echo "# looks ok, found $program here: $bin"
     fi
  done
}

# 2) check that two arguments are passed to the script from the commmand line
if [ $# -ne 2 ]; then
   echo
   echo "# $0 vers.$version needs two arguments: "
   echo -e "#\t1) the input fasta file extension_name <[fas|fasta|fna|faa]>"
   echo -e "#\t2) the alignment program to use <[muscle|clustalw]>"
   echo "# usage example: $0 fna muscle"
   echo
   echo "# NOTE: the script assumes that clustalw and muscle are in found in \$PATH"
   echo "# will check now for the presence of both binaries in \$PATH" 
   check_binaries
   echo
   exit 1
fi

#------------END FUNCTIONS --------------------#

# 3) assign positional parameters to named variables
ext="$1"
program="$2"

# 4) align all *.$ext fasta files in the current directory, 
for file in $(ls *.$ext); do
   if [ -s $file ]; then
        if [ $program = clustalw ]; then
            echo "# running clustalw -infile=$file -align -output=fasta -outfile=${file%.$ext}_cluw_aln.$ext"
	    clustalw -infile=$file -align -output=fasta -outfile=${file%.$ext}_cluw_aln.$ext
      elif [ $program = muscle ]; then
            echo "# running muscle -in $file -out ${file%.$ext}_mus_aln.$ext"
	    muscle -in $file -out ${file%.$ext}_mus_aln.$ext
      else
          echo "# sorry, don't know alignment program $algorithm"
	  echo "# $0 only knows clustal and muscle"
	  exit 1
      fi
   else
      echo "# $0 Error: could not find a valid *.$ext file in dir $(pwd)"
      exit 2
   fi
done
