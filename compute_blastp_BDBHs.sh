#!/usr/bin/env bash

#: progname: compute_pw_blastp_BDBHs.sh

progname=${0##*/}
vers=0.1_2022-12-10

set -eo pipefail

# DECLARE & INITIALIZE GLOBALS
#declare -a args
DEBUG=0
qcov=80
num_aln=10
n=''

RED='\033[0;31m'
NC='\033[0m' # No Color => end color

#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#

function print_start_time()
{
   printf '%(%T)T %s' '-1'
}
#----------------------------------------------------------------------------------------- 

function filter_best_hits()
{
    local infile n s 

    declare -a names 
    declare -A max best

    infile=$1

    while read -r n f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 s; 
    do
        if [[ -z "${max[$n]}" ]]
        then
             names+=( "$n" )
        fi

        # float comparisons using bc
        if [[ ! "${max[$n]}" ]] || [[ 1 -eq "$(echo "$s > ${max[$n]}" | bc)" ]]
        then
            max["$n"]="$s"
            unset 'best[$n]'
        fi

        # float comparisons using bc
        if [[ 1 -eq "$(echo "$s == ${max[$n]}" | bc)" ]]
        then
             best[$n]="${n}\t$f2\t$f3\t$f4\t$f5\t$f6\t$f7\t$f8\t$f9\t${f10}\t${f11}\t$s"
        fi
    done < "$infile"


    for n in "${names[@]}"; do
        [[ -z "${best[$n]}" ]] && continue
        echo -e "${best[$n]}"
    done
}
#----------------------------------------------------------------------------------------- 

function check_dependencies()
{
    for programname in fas2tab.pl tab2fas.pl blastp makeblastdb 
    do
       #if which $programname >/dev/null; then <== avoid which
       # see: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script

       bin=$(type -P $programname)
       if [ -z "$bin" ]; then
          echo
          echo -e "${RED}# ERROR: $programname not in place!${NC}\n"
          echo "# ... you will need to install \"$programname\" first or include it in \$PATH"
          echo "# ... exiting"
          exit 1
       fi
    done

    echo
    echo '# Run check_dependencies() ... looks good: all required binaries and perl scripts are in place.'
    echo
}
#----------------------------------------------------------------------------------------- 

function print_version()
{
   cat <<EOF

$progname v.$vers

EOF
  exit 1
}
#----------------------------------------------------------------------------------------- 

function print_help()
{
   cat <<EOF
   $progname v.$vers usage:
   
   REQUIRED:
    -a <string> input proteome 1 (FASTA)
    -b <string> input proteome 2 (FASTA)
    
   OPTIONAL
   -D <flag> print debugging info
   -h <flag> print this help
   -n <int> number of blast alignments [def:$num_aln]
   -q <int> minimum query coverage percentage cutoff [def:$qcov]
   -v <flag> print version
   
   AIM: Computes blastp-based BDBHs between a pair of proteome 
        (protein FASTA) files
   
   SOURCE: the latest version can be fetched from 
           https://github.com/vinuesa/TIB-filoinfo
	   
   LICENSE: GPL v3.0. 
      See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
        
EOF

   check_dependencies
   
   exit 1  
}
#----------------------------------------------------------------------------------------- 


#------------------------------------#
#----------- GET OPTIONS ------------#
#------------------------------------#

a_file=''
b_file=''
runmode=''

args=("$@")

while getopts ':a:b:n:q:hDv?:' OPTIONS
do
   case $OPTIONS in

   a)   a_file=$OPTARG
        ;;
   b)   b_file=$OPTARG
        ;;
   n)   num_aln=$OPTARG
        ;;
   q)   qcov=$OPTARG
        ;;
   h)   print_help
        ;;
   v)   print_version
        ;;
   D)   DEBUG=1
        ;;
   :)   printf "argument missing from -%s option\n" "$OPTARG"
   	 print_help
     	 exit 2 
     	 ;;
   ?)   echo "need the following args: "
   	 print_help
         exit 3
	 ;;
   *)   echo "An  unexpected parsing error occurred"
         echo
         print_help
	 exit 4
	 ;;	 
   esac >&2   # print the ERROR MESSAGES to STDERR
done

shift $((OPTIND - 1))

if [ -z "$a_file" ]
then
       echo "# ERROR: no input fasta file a!"
       print_help
       exit 1    
fi

if [ -z "$b_file" ]
then
       echo "# ERROR: no input fasta file b!"
       print_help
       exit 1    
fi

###>>> Exported variables !!!
#declare -x skip_seqs_gt=$skip_seqs_gt perl # export only2perl!!!  $ENV{skip_seqs_gt}

wkdir=$(pwd)

start=$(print_start_time)
hostn=$(hostname)

# print run parameters
echo "
===================================================================================================== 
$progname vers. $vers 
-----------------------------------------------------------------------------------------------------
 run on $hostn at ${start/ /} with the following parameters: 
 - wkdir=$wkdir | num_aln=$num_aln | qcov=$qcov | DEBUG=$DEBUG 
 - invocation: $progname ${args[*]}
===================================================================================================== 
 " >&2

#
# >>>> MAIN CODE <<<<
#

# 1. edit headers for makeblastdb (blast+)
print_start_time && echo '# formatting input FASTA files for indexing'
perl -pe 'if(/^>/){$c++; s/>/>lcl\|A_DB_$c /}' "$a_file" > "${a_file}ed"
perl -pe 'if(/^>/){$c++; s/>/>lcl\|B_DB_$c /}' "$b_file" > "${b_file}ed"

# 2. run makeblastdb on both input proteomes
ext="${a_file##*.}"

print_start_time && echo "# Generating indexed blastp databases"
for f in *"${ext}"ed
do 
    #Note: -parse_seqids results in query and subject IDs with different structure => lcl|A_DB_203	B_DB_166
    makeblastdb -in "$f" -dbtype prot -parse_seqids &> /dev/null
    
    fas2tab.pl "$f" | sed '/^$/d' > "${f}"tab
done

# 3. Run pairwise blastp  
print_start_time && echo "# running: blastp -query ${a_file}ed -db ${b_file}ed -qcov_hsp_perc $qcov -outfmt 6 -num_alignments $num_aln -num_threads 10 > AvsB" >&2
blastp -query "${a_file}"ed -db "${b_file}"ed -qcov_hsp_perc "$qcov" -outfmt 6 -num_alignments "$num_aln" -num_threads 10 > AvsB

print_start_time && echo "# runnign: blastp -query ${b_file}ed -db ${a_file}ed -qcov_hsp_perc $qcov -outfmt 6 -num_alignments $num_aln -num_threads 10 > BvsA" >&2
blastp -query "${b_file}"ed -db "${a_file}"ed -qcov_hsp_perc "$qcov" -outfmt 6 -num_alignments "$num_aln" -num_threads 10 > BvsA 


# 4. Retrieve the highest-scoring hit out of the -num_alignments $num_aln hits
# Note: makeblastdb with -parse_seqids produces sobject cols without the lcl| prfix: => lcl|A_DB_203	B_DB_166
#        and therefore we remove them to have subject and query IDs with the same structure,
#        required for hash key comparisons later in the code
print_start_time && echo "# filter_best_hits AvsB > AvsB.best" >&2
filter_best_hits AvsB | sed 's#lcl|##' > AvsB.best

print_start_time && echo "# filter_best_hits BvsA > BvsA.best" >&2
filter_best_hits BvsA | sed 's#lcl|##' > BvsA.best

# 5. Compute A vs B reciprocal best hits
# the following hashes will hold query=>subject and subject=>query results
#    that will be used to to identify the reicprocal best hits (RBHs)
declare -A AB_hits BA_hits 
print_start_time && echo "# Computing A vs B reciprocal best hits @ qcov=$qcov" >&2
while read -r l
do
     q=$(echo "$l" | awk '{print $1}')
     subj=$(echo "$l" | awk '{print $2}')
     AB_hits["$q"]="$subj"
done < AvsB.best

while read -r l
do
     q=$(echo "$l" | awk '{print $1}')
     subj=$(echo "$l" | awk '{print $2}')
     BA_hits["$subj"]="$q"
done < BvsA.best

# 5.1 find the RBHs by comparing the two hashes and sort output by A_DB_# key
RBH_outfile="${a_file%.*}"_vs_"${b_file%.*}"_qcov"${qcov}"_RBHs_2col.tsv
for kA in "${!AB_hits[@]}"
do
    if [[ "${BA_hits[$kA]}" == "${AB_hits[$kA]}" ]]
    then
         printf "%s\t%s\n" "$kA" "${AB_hits[$kA]}"
    fi
done | sort -k1.6g > "$RBH_outfile"

# 5.2 make a more informative RBH_outfile name, including number of RBHs found
num_RBHs=$(wc -l "$RBH_outfile" | awk '{print $1}')
RBH_out="${a_file%.*}"_vs_"${b_file%.*}"_qcov"${qcov}"_${num_RBHs}RBHs_2col.tsv
mv "$RBH_outfile" "$RBH_out"
RBH_outfile="$RBH_out"
print_start_time && echo "# Found $num_RBHs reciprocal best hits between $a_file and $b_file at qcov=${qcov}%" >&2

# 6. Write cluster fastas
print_start_time && echo "# Generating cluster fastas"
c=0
while read -r A B; do 
    c=$((c + 1)) # note that simply c=$((c++)) or ((c++)) does not work; it won't even enter the loop
    grep -w "$A" "${a_file}"edtab > cluster_"${c}".fastab
    grep -w "$B" "${b_file}"edtab >> cluster_"${c}".fastab
done < "$RBH_outfile"

# 6.1 reconstitute fasta
for f in cluster_*.fastab
do
    tab2fas.pl "$f" > "${f%tab}"
done

clusters_dir="${a_file%.*}"_vs_"${b_file%.*}"_qcov"${qcov}"_"${num_RBHs}"_clusters
rm ./*fastab

[[ -d "$clusters_dir" ]] && rm -rf "$clusters_dir"
[[ ! -d "$clusters_dir" ]] && mkdir "$clusters_dir" || { echo "ERROR: could not generate dir $clusters_dir" && exit 1 ; }

print_start_time && echo "# moving cluster fastas to $clusters_dir"
mv cluster*.fas "$clusters_dir"

# 7. final cleanup
print_start_time && echo "# final cleanup"
((DEBUG == 0)) && rm AvsB BvsA AvsB.best BvsA.best "${a_file}ed" "${b_file}ed" ./*ed.p* ./*edtab

exit 0
