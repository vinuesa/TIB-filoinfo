#!/usr/bin/bash -

#: progname: compute_RBH_clusters.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#
#: AIM: Edit FASTA header of proteomes not fetched from standard databases (sp|tr ...)
#:      This is a companion script for compute_RBH_clusters.sh
#: - for example, for ACAC3.faa, >FUN_000001-T1 FUN_000001 is transformed to: 
#: 			      >lcl|FUN_000001-T1_ACAC3 FUN_000001
#----------------------------------------------------------------------------------------
#: LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
#----------------------------------------------------------------------------------------
#: DISCLAIMER 
#: THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
#: APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
#: HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
#: OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
#: THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#: PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
#: IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
#: ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#----------------------------------------------------------------------------------------
#: GitHub repo: you can fetch the latest version of the script from:
#   https://github.com/vinuesa/TIB-filoinfo/blob/master/compute_blastp_RBH_orthologous_clusters.sh
# wget -c https://raw.githubusercontent.com/vinuesa/TIB-filoinfo/master/fix_FASTA_headers4blastdb.sh
#----------------------------------------------------------------------------------------


# set Bash's unofficial strict mode
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

progname=${0##*/}
vers='0.1_2023-10-28'  # first commit 


function print_help {

   cat << HELP

   Usage: $progname <fasta2edit>
   
   * AIM: Edit FASTA header of proteomes not fetched from standard databases (sp|tr ...)
          This is a companion script for compute_RBH_clusters.sh
     - for example, for ACAC3.faa, >FUN_000001-T1 FUN_000001 is transformed to: 
                                   >lcl|FUN_000001-T1_ACAC3 FUN_000001
   
   # To process specific files, list them in a for loop like shown below:
   for f in ACAC3.faa Neff2.faa; do $progname \$f; done

   # To process all *.faa files, list them in a for loop like shown below:
   for f *.faa; do $progname \$f; done

HELP

   exit 1

}

[ $# -ne 1 ] && print_help

perl -pe 'if(/^>(\S+)\s+(.*)?$/){ $id1=$1; $id2=$2; $f=(split(/\./, $ARGV))[0]; s/>.*$/>lcl|$id1\_$f $id2/}' $1 > ${1%.*}ed.faa
