#!/usr/bin/env bash

#: progname: compute_RBH_clusters.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#
#: AIM: Wrapper script around NCBI's blastp, to compute reciprocal best hits (RBHs) between a
#:      reference genome (manually or automatically selected) and a set of additional
#:      proteomes (protein fasta files). It also computes the core and non_core RBH sets, 
#:      writing the corresponding clusters (FASTA files) to disk
#
#: Design: The code is partially parallelized:
#:         - the user can specify a number of threads to parallelize blastp
#:         - the final blastdbcmd is called from xargs to parallelize the 
#:           retrieval and writing to disk of RBH clusters 
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
# wget -c https://raw.githubusercontent.com/vinuesa/TIB-filoinfo/master/compute_blastp_RBH_orthologous_clusters.sh
#----------------------------------------------------------------------------------------

progname=${0##*/}
vers='1.1.6_2023-10-28' # compute_RBH_clusters.sh v1.1.6_2023-10-28 
		       #  - minor fixes: quoted 'EOF' in print_notes; changed invocation: $progname ${args[*]} to ${args[@]}
		       #  - additional cleanup of old snippets in the main code section; some added to notes()

min_bash_vers=4.4 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

# set Bash's unofficial strict mode
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# fixed custom blastp header
cols='6 qseqid sseqid pident gaps length qlen slen qcovs evalue score'

# Initialize variables with default values; no undefined/unbound variables allowed in strict mode 
DEBUG=0
qcov=60
num_aln=1
threads=4  #$(nproc) # all cores available
ext=faa
ref=''

task='blastp-fast'
mat=BLOSUM62
Eval=0.00001
seg=yes      # yes|no
mask=true    # true|false
best_hit_overhang=0.1    # recommended value in blastp -help
best_hit_score_edge=0.1  # recommended value in blastp -help
#fix_header=0 # do not fix FASTA header for makeblastdb; disabled in vers. 1.1.3_2023-10-26

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
#BLUE='\033[0;34m'
LBLUE='\033[1;34m'
#CYAN='\033[0;36m'
NC='\033[0m' # No Color => end color

#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#

# Function to print error messages
function error {
    echo -e "${RED}Error: $1 ${NC}" 1>&2
    exit 1
}
#-----------------------------------------------------------------------------------------

# Function to check that files were generated
function check_file {
    if [ -s "$1" ]
    then
        echo -e "${GREEN} wrote $1 ${NC}"
    elif [ ! -s "$1" ] && [ -n "$2" ] # pass a second arg, like warn
    then
         echo -e "${YELLOW}WARNING: could not write $1 ${NC}" 1>&2
    else
        echo -e "${RED}Error: could not write $1 ${NC}" 1>&2
        exit 1
    fi
}
#-----------------------------------------------------------------------------------------

# Function to check Bash version
function check_bash_version {
   # Checks if the Bash version meets the minimum required version.
   # the more modern bash array syntax and functions require bash version >= 4.4
   local bash_vers min_bash_vers
   min_bash_vers=$1
   bash_vers=$(bash --version | awk 'NR==1{print $4}' | sed 's/(.*//' | cut -d. -f1,2)
   
   if (( $(bc <<< "$bash_vers < $min_bash_vers") )); then
       error "You need Bash version ${min_bash_vers} or higher to run this script."
   fi
   
   echo "${bash_vers}"
}
#-----------------------------------------------------------------------------------------

function print_start_time {
   # Prints the current time in the format "HH:MM:SS".
   printf '%(%T)T %s' '-1'
}
#----------------------------------------------------------------------------------------- 

function print_start_date {
   # Prints the current date in the format "YYYY-MM-DD".
   printf '%(%F)T %s' '-1'
}
#----------------------------------------------------------------------------------------- 

# Function to check dependencies
function check_dependencies {
    local dependencies=("blastp" "makeblastdb" "blastdbcmd" "blastdb_aliastool")
    
    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            error "$dep not found in PATH. Install it or add it to PATH."
        fi
    done

    echo -e "${GREEN} All required dependencies are in place.${NC}"
}

#----------------------------------------------------------------------------------------- 

function select_ref { 
     # Automatically selects the smallest genome as the reference.
     local ext
     ext=$1
     
     # find the genome with the smallest number of gene|protein sequences
     for f in *."${ext}"
     do 
          echo -ne "$f\t"
	  grep -c '>' "$f"
     done | sort -k2g | awk 'NR == 1{print $1}'
}
#----------------------------------------------------------------------------------------- 

function print_n_processors {
    # Prints the number of processors/cores on the system.
    # could also use nproc, which is part of GNU core utils ;)
    awk '/processor/{p++}END{print p}' /proc/cpuinfo
}
#----------------------------------------------------------------------------------------- 

function print_end_message {
   # Prints a message to acknowledge the use of the script.
   cat <<EOF
  
  Done!
  
  ========================================================================================
  If you use $progname v.$vers for your research, then please:
  
  
  1. Cite the code in your work as:   
  Pablo Vinuesa. $progname v.$vers 
       https://github.com/vinuesa/TIB-filoinfo/blob/master/$progname
  
  2. Give it a like on the https://github.com/vinuesa/TIB-filoinfo/ repo
  
  Thanks!

EOF

exit 0

}
#----------------------------------------------------------------------------------------- 

function print_version {
   #  Prints the version information.
   cat <<EOF

$progname v.$vers

EOF
  exit 1
}
#----------------------------------------------------------------------------------------- 

function run_blastp {
   local task q db outfile
   task=$1 # default: blastp-fast
   q=$2
   db=$3
   outfile=$4
   
   # parameterizations in part following suggestions by Julie E. Hernandez-Salmeronóand Gabriel Moreno-Hagelsieb 
   # in BMC Genomics volume 21, Article number: 741 (2020); https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07132-6
   blastp -task "$task" -query "$q" -db "$db" -matrix "$mat" -outfmt "$cols" -num_alignments $num_aln -num_threads $threads -qcov_hsp_perc "$qcov" \
    -seg "$seg" -soft_masking "$mask" -best_hit_overhang "$best_hit_overhang" -best_hit_score_edge "$best_hit_score_edge" -use_sw_tback \
    -evalue "$Eval" -out "$outfile" 2> /dev/null
    
   check_file "$outfile"
}
#----------------------------------------------------------------------------------------- 

function print_notes {

   cat << 'EOF'
   
    1. The auxiliary scritp fix_FASTA_headers4blastdb.sh might be helpful
        to generate properly-formatted FASTA headers of local/user proteomes.
	It is available from https://github.com/vinuesa/TIB-filoinfo

    2. A detailed tutorial on using blast efficiently on the Linux command line can found here:
         https://vinuesa.github.io/TIB-filoinfo/sesion3_BLAST/
    
    3. For a versatile and highly customizable pan-genome analysis software package, 
         consider using GET_HOMOLOGUES
         https://doi.org/10.1128%2FAEM.02411-13
	 https://github.com/eead-csic-compbio/get_homologues
	 https://hub.docker.com/u/eeadcsiccompbio
	 
    4. To perform core- and/or pan-genome phylogenomic analyses,
       consider using the GET_PHYLOMARKERS package
       https://doi.org/10.3389/fmicb.2018.00771 
       https://github.com/vinuesa/get_phylomarkers
       https://hub.docker.com/r/vinuesa/get_phylomarkers  
 
 
  ====== DEV NOTES =====
#-----------------------------
# 6. Write cluster FASTA files
#-----------------------------
# Below is the code for three strategies to write RBH cluster FASTA files, two of them commented out
#  The final (uncommented) one is the fastest of the benchmarked strategies. 

## >>> Take 1 6.1 read all source FASTA file into memory (as the hash seqs) for later filtering
# print_start_time && echo "# reading all source FASTA files into memory \(as the hash seqs\) ..."
#
#declare -A seqs; 
#
#for f in "${infiles[@]}"
#do
#   while read -r l
#   do 
#      if [[ "$l" =~ ^\> ]] # line contains the FASTA header
#      then 
#           key=$l
#      else               # concatenate sequence lines
#           seqs["$key"]="${seqs[$key]}""${l}"
#      fi
#   done < "$f"
#done


## 6.2 filter the hash seqs with cluster keys saved in each line of RBHs_matrix.tsv
# for k in "${!seqs[@]}"; do if [[ $k =~ \>Q8MYF2 ]]; then echo -e "$k\n${seqs[$k]}"; fi;  done < test.faa
#print_start_time && echo "# Filtering the seqs hash and writing RBH cluster FASTA files ..."

#initialize cluster counter
#c=0

## read each line of the RBHs_matrix.tsv
##  into the ids array using the read -r -a idiom,
##  and filter the seqs hash with the corresponding IDs as keys
##  to write out the cluster_x.fas files
## NOTE: the comment block below implementing
##  a hash traversing strategy below is very slow, 
##  even with the continue 2 statement; 
##    should benchmark the following options: 
##          - filtering fastab with grep
##          - use blastdbcmd and 
#
#while read -r -a ids
#do 
#    ((c++)) 
#    # iterate over all indexes of the idx array
#    #  and append them to the growing cluster_"${c}".fastab file
#    for (( idx=0; idx <= ((${#ids[@]} -1)); idx++))
#    do
#	# iterate of all source sequence FASTA headers (k)
#	#  and print out only sequences matching the ids 
#	#  in cluster c, concatenating them '>>' to cluster_"${c}".fastab
#	for k in "${!seqs[@]}"
#        do 
#            if [[ "$k" =~ "${ids[$idx]}" ]]
#	    then 
#	         ((DEBUG > 0)) && echo "DEBUG (write clusters): k:$k; idx:$idx; ids[$idx]}:${ids[$idx]}; c=$c" >&2
#		 echo -e "$k\n${seqs[$k]}"
#		 continue 2
#	    fi
#        done >> cluster_"${c}".fas
#    done 
#done < RBHs_matrix.tsv

           
EOF

   exit 1  

}
#----------------------------------------------------------------------------------------- 

function print_help {
   # Prints the help message explaining how to use the script.
   cat <<EOF
   
   Usage: $progname -d <dir> [-e <ext>] [-r <reference proteome>] [-E <Eval>] [-m <matrix>] [-n <num_aln>] [-q <qcov>]  
             [-t <threads>] [-T <blastp|blastp-fast>] [-S <yes|no>] -M [<false|true>] [-D] [-h] [-v]
   
   REQUIRED:
    -d <string> path to directory containing the input proteomes (protein FASTA)
    
   OPTIONAL
   -D <flag> print debugging info
   -e <string> fasta file extension name [def:$ext]
   -E <integer or float -ge 0> E-value [def:$Eval]
   -h <flag> print this help
   -n <int> number of blast alignments [def:$num_aln]
   -m <string> matrix name <BLOSUM45|BLOSUM62|BLOSUM80> [def:$mat]
   -M <false|true> soft_masking [def:$mask]
   -q <int> minimum query coverage percentage cutoff [def:$qcov]
   -r <string> name of user-selected reference proteome
   -S <yes|no> SEG filtering [def:$seg]
   -t <int> number of threads for blastp runs [def:$threads]
   -T <string> blastp task <blastp|blastp-fast> [def:$task]
   -v <flag> print version
   
   EXAMPLES:
     $progname -d . -m BLOSUM80 -T blastp -q 85 -r my_ref.fa -e fa
     $progname -d proteome_files -n 5 -t \$(nproc)
   
   AIM: 
      Wrapper script around NCBI\'s blastp, to compute reciprocal best hits (RBHs) between a
      reference genome (manually or automatically selected) and a set of additional proteomes
      (protein fasta files). It also computes the core and non_core RBH sets, writing the 
      corresponding clusters (FASTA files) to disk in the core_clusters and nonCore_clusters dirs.
   
   SOURCE: the latest version can be fetched from https://github.com/vinuesa/TIB-filoinfo
	   
   LICENSE: GPL v3.0. 
      See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
      
   NOTES: 
    1. Assumes that the input FASTA sequences haver properly-formatted headers 
         for indexing with makeblastdb; locally/user generated proteomes should have
	 the following FASTA header structure: '>lcl|uniqueID_ORGNemonic'
	 like in '>lcl|FUN_005793_ACAC3' or '>lcl|000762_Sm18'. 
    2. run $progname -N for additional notes	 
	 
EOF

   check_dependencies
   
   exit 1  
}
#----------------------------------------------------------------------------------------- 


#------------------------------------#
#----------- GET OPTIONS ------------#
#------------------------------------#
# This section uses getopts to process command-line arguments 
#    and set the corresponding variables accordingly. 

[ $# -eq 0 ] && print_help

args=("$@")


while getopts ':d:e:E:m:M:n:q:r:S:t:T:hDvN?:' OPTIONS
do
   case $OPTIONS in

   d)   proteome_dir=$OPTARG
        ;;
   e)   ext=$OPTARG
        ;;
   E)   Eval=$OPTARG
        ;;
   h)   print_help
        ;;
   m)   mat=$OPTARG
        ;;
   M)   mask=$OPTARG
        ;;
   n)   num_aln=$OPTARG
        ;;
   N)   print_notes
        ;;
   q)   qcov=$OPTARG
        ;;
   r)   ref=$OPTARG
        ;;
   S)   seg=$OPTARG
        ;;
   t)   threads=$OPTARG
        ;;
   T)   task=$OPTARG
        ;;
   v)   print_version
        ;;
   D)   DEBUG=1
        ;;
   :)   printf "argument missing from -%s option\n" "$OPTARG"
   	 print_help
     	 ;;
   ?)   echo "need the following args: "
   	 print_help
	 ;;
   *)   echo "An  unexpected parsing error occurred"
         echo
         print_help
	 ;;	 
   esac >&2   # print the ERROR MESSAGES to STDERR
done

shift $((OPTIND - 1))

# Check required params were provided
if [[ -z "$proteome_dir" ]]; then
       error "Missing required option: -d <dir>"
       print_help
fi

if [[ -z "$ref" ]]
then
     ref_selection="auto"
else
     ref_selection="$ref"
fi

today=$(print_start_date)
hostn=$(hostname)

bash_vers=$(check_bash_version "$min_bash_vers")
nprocs=$(print_n_processors)
os=$(uname -o)
blast_vers=$(blastp -version | awk 'NR == 1{print $2}')

# OK, ready to start the analysis ...
start_time=$SECONDS

# print run parameters
echo "
===================================================================================================== 
$progname vers. $vers 
-----------------------------------------------------------------------------------------------------
 run on $hostn using $os with $nprocs processors and bash v.${bash_vers} on ${today/ /} 
   with the following parameters: 
 - proteome_dir=$proteome_dir | fasta_extension=$ext
 - BLASTP params: blastp v.${blast_vers} | task=$task | num_aln=$num_aln | qcov=$qcov | 
                 Eval=$Eval | mat=$mat | seg=$seg | mask=$mask | threads=$threads 
 - reference=$ref_selection | DEBUG=$DEBUG 
 - invocation: $progname ${args[@]}
===================================================================================================== 
 " >&2

# -------------------
# >>>> MAIN CODE <<<<
# -------------------
# Main script logic starts here

# ---------------------------------------
# 0. setup pipeline and select reference
# ---------------------------------------

# Move to proteome directory
print_start_time && echo "# working in $proteome_dir"
cd "$proteome_dir" || error "Could not cd into $proteome_dir"

wkdir=$(pwd)

print_start_time && echo '# Selecting the smallest genome as the reference'
# automatically select the smallest reference, if not provided as ARG
if [[ "$ref_selection" == "auto" ]]; then
    ref=$(select_ref "$ext")
    echo "  - Selected $ref as the reference genome"
fi

printf '%s\n' '-----------------------------------------------------------------------------------------------------'

# ----------------------------------------------------------
# 1. Read proteome files into the non_ref and infiles arrays
# ----------------------------------------------------------

# keep the ref as first element in fasta arrays
declare -a non_ref infiles 
non_ref=( $(ls *"${ext}" | grep -v "$ref") )
infiles=("$ref" "${non_ref[@]}")


# -------------------------------------------------------------------
# 2. run makeblastdb on all input proteomes with edited FASTA headers
# -------------------------------------------------------------------
# Each proteome is used to create a blastp database using makeblastdb.
print_start_time && echo "# Generating indexed blastp databases"

for f in "${infiles[@]}"
do 
    makeblastdb -in "$f" -dbtype prot -parse_seqids &> /dev/null
done


print_start_time && echo "# Generating the aliased blastp database allDBs ..."
blastdb_aliastool -dblist_file <(printf '%s\n' "${infiles[@]}") -dbtype prot -out allDBs -title allDBs

printf '%s\n' '-----------------------------------------------------------------------------------------------------'


#-----------------------------------
# 3. Run and process pairwise blastp 
#-----------------------------------
# For each non-reference proteome, blastp is run against the reference proteome 
#   (REFvsGENO) and vice versa (GENOvsREF).

for f in "${non_ref[@]}"
do
    ref_vs_geno_blastout=${ref%.*}vs${f%.*}_best_hits.tmp
    geno_vs_ref_blastout=${f%.*}vs${ref%.*}_best_hits.tmp
    
    print_start_time && echo "# Running: run_blastp $task ${ref}ed ${f}ed $ref_vs_geno_blastout" 
    run_blastp "$task" "$ref" "$f" "$ref_vs_geno_blastout"

   # Retrieve the best nonREF proteome database hits using blastdbcmd, onlfy if qcov > \$qcov
   print_start_time && echo "# Retrieving the best hits from $ref_vs_geno_blastout with blastdbcmd ... "
   blastdbcmd -entry_batch <(awk -F"\t" -v qcov="$qcov" '$8 > qcov{print $2}' "$ref_vs_geno_blastout" | sort -u) -db "$f" > "${ref%.*}vs${f%.*}"_besthits.faa
   check_file "${ref%.*}vs${f%.*}"_besthits.faa
   
   num_hits=$(grep -c '^>' "${ref%.*}vs${f%.*}"_besthits.faa)
   
   if ((num_hits == 0))
   then
        echo "WARNING: no hits in ${ref%.*}vs${f%.*}_besthits.faa"
	rm "${ref%.*}vs${f%.*}"_besthits.faa
	continue
   fi
   
    print_start_time && printf '%s\n' "# Running: run_blastp $task ${ref%.*}vs${f%.*}_besthits.faa ${ref}ed $geno_vs_ref_blastout ..."    
    run_blastp "$task" "${ref%.*}vs${f%.*}"_besthits.faa "$ref" "$geno_vs_ref_blastout"

    # Sort the blastp output table from the preceding search by increasing E-values (in column 9) and decreasing scores (col 10)
    #    & filter out unique REF vs nonREF RBHs using AWK hashes from the sorted blast output table with qcov > $qcov
    print_start_time && echo "# Filtering out unique REF vs nonREF RBHs from the sorted blast output table with qcov > $qcov"
    for GENOid in $(cut -f1 "$geno_vs_ref_blastout" | sort -u)
    do 
       grep "$GENOid" "$ref_vs_geno_blastout"
    done | sort -gk9,9 -gk10,10 | \
           awk -v qcov="$qcov" 'BEGIN{FS=OFS="\t"}{REFid[$1]++; GENOid[$2]++; if(REFid[$1] == 1 && GENOid[$2] == 1 && $8 > qvov) print }' > \
             "${ref_vs_geno_blastout%.*}"_RBHs_qcov_gt"${qcov}".tsv
    
    check_file "${ref_vs_geno_blastout%.*}"_RBHs_qcov_gt"${qcov}".tsv
done

printf '%s\n' '-----------------------------------------------------------------------------------------------------'


#-----------------------------------------------------------------------
# 4. Identify REF proteins shared by all tsv files holding pairwise RBHs
#-----------------------------------------------------------------------
# Find the intersections of REFs in all tsv files
print_start_time && echo "# Computing the intersections of REF proteins in all tsv files holding pairwise RBHs ... "
awk '{r[$1]++; if(r[$1] == ARGC-1) print $1}' ./*.tsv > REF_RBH_IDs.list
[[ ! -s REF_RBH_IDs.list ]] && error "could not write REF_RBH_IDs.list"

intersection_size=$(wc -l REF_RBH_IDs.list | awk '{print $1}')
((intersection_size > 0)) && echo -e "${LBLUE}  Found $intersection_size core RBHs shared by ${#non_ref[*]} nonREF proteomes with the $ref reference proteome${NC}"
((intersection_size == 0)) && error "# ERROR: found $intersection_size core orhtologous genes among ${#infiles[*]} input proteomes ..."

printf '%s\n' '-----------------------------------------------------------------------------------------------------'


#-------------------------------------------------------------------------------------------------------------------------
# 5. Loop over tsv files and generate RBHs_matrix.tsv core_genome_clusters.tsv and nonCore_genome_clusters.tsv tables
#-------------------------------------------------------------------------------------------------------------------------
# Cluster Computation

print_start_time && echo "# Computing clusters of homologous sequences ..."

# The core_ref hash counts the instances of the REFERNCE_IDs in the RBH tables
declare -A core_ref
core_ref=()

# The all_clusters hash is indexed by REFERNCE_IDs 
#   and as its value holds the RBHs as a tab-separated
#   string of IDs from nonREF proteomes    
declare -A all_clusters
all_clusters=()

# 5.1 Construct the all_clusters hash, indexed by reference proteome, 
#  containing as value a string of tab-separated nonREF proteome RBH IDs.
print_start_time && echo "# Populating the all_clusters hash ..."
for t in *RBHs_*.tsv; do
    [ ! -s "$t" ] && error "file: $t does not exist or is empty"
    while read -r REF QUERY rest
    do
        # count the instances of REFERNCE_IDs in each RBHs_*.tsv tables
	(( core_ref["$REF"]++ ))
	if (( ${core_ref["$REF"]} == 1 ))
	then
	    all_clusters["$REF"]="$QUERY" 
	else
	    all_clusters["$REF"]="${all_clusters[$REF]}\t$QUERY"
	fi
    done < "$t" || { echo "Failed to process file: $t"; exit 1; } # required test for set -e compliance
done


# 5.2 print the RBHs_matrix.tsv
print_start_time && echo "# Printing the RBHs_matrix, core_genome_clusters, and nonCore_genome_clusters files ..."
for key in "${!all_clusters[@]}"
do
      echo -e "${key}\t${all_clusters[$key]}"
done > RBHs_matrix.tsv
check_file RBHs_matrix.tsv

# 5.3 print the core_genome_clusters.tsv
awk -v ninfiles="${#infiles[@]}" 'NF == ninfiles' RBHs_matrix.tsv > core_genome_clusters.tsv
check_file core_genome_clusters.tsv

# 5.4 print the nonCore_genome_clusters.tsv
awk -v ninfiles="${#infiles[@]}" 'NF != ninfiles' RBHs_matrix.tsv > nonCore_genome_clusters.tsv
check_file nonCore_genome_clusters.tsv warn

printf '%s\n' '-----------------------------------------------------------------------------------------------------'


#-----------------------------
# 6. Write cluster FASTA files
#-----------------------------

print_start_time && printf '%s\n' '# Extracting RBH cluster FASTA files ...'

# 6.1 (take 3; see notes) blastdbcmd is called from parallel or xargs 
#  - write the each line of the RBHs_matrix.tsv IDs to a tmpfile
#    to pass the list of tmpfiles to a parallel call of blastdbcmd
print_start_time && echo "  - Writing each line of the RBHs_matrix.tsv IDs to an idstmp file por parallelization ..."

#initialize cluster counter
c=0
while read -r -a ids
do 
    ((c++)) 
    # write each line of the RBHs_matrix.tsv to a temporal file
    printf '%s\n' "${ids[@]}" > cluster_"${c}".idstmp
done < RBHs_matrix.tsv || { echo "Failed to process file: RBHs_matrix.tsv"; exit 1; } # required test for set -e compliance

# 6.2 Pass the list of tmpfiles to a parallel call of blastdbcmd, or, if not available, call it from xargs
if command -v parallel &> /dev/null 
then 
    print_start_time && echo "  - Running blastdbcmd through parallel with -j \$(nproc) ..."
    find . -name '*.idstmp' | parallel --gnu -j $(nproc) 'blastdbcmd -db allDBs -dbtype prot -entry_batch {} -out {.}.fas'  
else
    # use the more portable find | xargs idiom to parallelize the blastdbcmd call of parallel is not found on host
    print_start_time && echo "  - Running blastdbcmd in parallel with xargs using all available cores \$(nproc) ..."
    find . -name '*.idstmp' -print0 | xargs -0 -P $(nproc) -I % blastdbcmd -db allDBs -dbtype prot -entry_batch % -out %.fas

    # rename *.idstmp.fas cluster files with rename, if available
    if command -v rename &> /dev/null; 
    then 
        rename 's/\.idstmp//' *.fas
    fi
fi

# 6.3 filter out core and nonCore clusters
print_start_time && echo "# Moving core and nonCore clusters to their respective directories ..."

mkdir core_clusters || error "could not mkdir core_clusters"
mkdir nonCore_clusters || error "could not mkdir nonCore_clusters"
for f in cluster_*.fas
do 
    if [[ $(grep -c '^>' "$f") -eq "${#infiles[@]}" ]]
    then 
         mv "$f" core_clusters/core_"${f}" || error "could not mv $f to core_clusters/core_${f}"
    else
         mv "$f" nonCore_clusters/nonCore_"${f}" || error "could not mv $f to nonCore_clusters/nonCore_${f}"
    fi
done

printf '%s\n' '-----------------------------------------------------------------------------------------------------'


#-----------------
# 7. final cleanup
#-----------------
#Unnecessary files generated during the process are removed.
print_start_time && echo "# Tidying up $wkdir ..."

printf '%s\n' '-----------------------------------------------------------------------------------------------------'

((DEBUG == 0)) && rm *."${ext}".* ./*best_hits.tmp ./REF_RBH_IDs.list ./*besthits.faa ./*.idstmp allDBs.pal
gzip *RBHs_qcov_gt*tsv

elapsed=$(( SECONDS - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days, %H hr, %M min, %S sec')"

print_end_message

exit 0



