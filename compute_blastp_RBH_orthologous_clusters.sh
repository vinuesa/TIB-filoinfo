#!/usr/bin/env bash

#: progname: compute_RBH_clusters.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#
#: AIM: Simple script around NCBI's blastp, to compute reciprocal best hits between a
#:      reference genome (manually or automatically selected) and a set of additional
#:      proteomes (protein fasta files). It also computes the core set, writing
#:      the corresponding clusters (FASTA files) to disk
#
#: Design: all blastp and other computations are run sequentially, 
#:         although the user can specify a number of threads to parallelize blastp 
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
#:   NOTES: 
#:    1. For a versatile and highly customizable pan-genome analysis 
#:       software package, consider using GET_HOMOLOGUES
#:         https://doi.org/10.1128%2FAEM.02411-13
#:	 https://github.com/eead-csic-compbio/get_homologues
#:	 https://hub.docker.com/u/eeadcsiccompbio
#:	 
#:    2. To perform core- and/or pan-genome phylogenomic analyses
#:       consider using the GET_PHYLOMARKERS package
#:       https://doi.org/10.3389/fmicb.2018.00771 
#:       https://github.com/vinuesa/get_phylomarkers
#:       https://hub.docker.com/r/vinuesa/get_phylomarkers          
#----------------------------------------------------------------------------------------
#: GitHub repo: you can fetch the latest version of the script from:
#   https://github.com/vinuesa/TIB-filoinfo/blob/master/compute_blastp_RBH_orthologous_clusters.sh
# wget -c https://raw.githubusercontent.com/vinuesa/TIB-filoinfo/master/compute_blastp_RBH_orthologous_clusters.sh
#----------------------------------------------------------------------------------------

progname=${0##*/}
vers=1.1_2023-10-25 # v1.1_2023-10-25 significant speedup by parallelizing cluster writing with xargs call of blastdbcmd
    # - Major rewrite and simplification of the RBH filtering code, now using AWK hashes. 
    # - Major rewrite and simplification of the code to write out clusters
    # - Customized and more informative BLAST results table fields
    # - higher customizability of BLAST runs
		         

min_bash_vers=4.4 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

set -o pipefail


# fixed custom blastp header
cols='6 qseqid sseqid pident gaps length qlen slen qcovs evalue score'

# Initialize variables with default values
DEBUG=0
qcov=70
num_aln=5
threads=$(nproc) # all cores available
ext=faa

mat=BLOSUM62
Eval=0.001
seg=yes      # yes|no
mask=true    # true|false
fix_header=0 # do not fix FASTA header for makeblastdb 

# Color codes for output
RED='\033[0;31m'
NC='\033[0m' # No Color => end color

#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#

# Function to print error messages
function error() {
    echo -e "${RED}Error: ${NC} $1" >&2
    exit 1
}
#-----------------------------------------------------------------------------------------

# Function to check that files were generated
function check_file() {
    if [ -s "$1" ]
    then
        echo " wrote $1"
    elif [ ! -s "$1" ] && [ -n "$2" ] # pass a second arg, like warn
    then
         echo -e "${RED}WARNING: could not write ${NC} $1" >&2
    else
        echo -e "${RED}Error: could not write ${NC} $1" >&2
        exit 1
    fi
}
#-----------------------------------------------------------------------------------------

# Function to check Bash version
function check_bash_version()
{
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

function print_start_time()
{
   # Prints the current time in the format "HH:MM:SS".
   printf '%(%T)T %s' '-1'
}
#----------------------------------------------------------------------------------------- 

function print_start_date()
{
   # Prints the current date in the format "YYYY-MM-DD".
   printf '%(%F)T %s' '-1'
}
#----------------------------------------------------------------------------------------- 

# Function to check dependencies
function check_dependencies() {
    local dependencies=("fas2tab.pl" "tab2fas.pl" "blastp" "makeblastdb" "blastdbcmd")
    
    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            error "$dep not found in PATH. Install it or add it to PATH."
        fi
    done

    echo "All required binaries and scripts are in place."
}

#----------------------------------------------------------------------------------------- 

function select_ref()
{ 
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

function print_n_processors()
{
    # Prints the number of processors/cores on the system.
    awk '/processor/{p++}END{print p}' /proc/cpuinfo
}
#----------------------------------------------------------------------------------------- 

function print_end_message()
{
   # Prints a message to acknowledge the use of the script.
   cat <<EOF
  
  Done!
  
  ========================================================================================
  If you use $progname v.$vers for your research,
  I would appreciate that you:
  
  1. Cite the code in your work as:   
  Pablo Vinuesa. $progname v.$vers 
       https://github.com/vinuesa/TIB-filoinfo/blob/master/$progname
  
  2. Give it a like on the https://github.com/vinuesa/TIB-filoinfo/ repo
  
  Thanks!

EOF

exit 0

}
#----------------------------------------------------------------------------------------- 

function print_version()
{
   #  Prints the version information.
   cat <<EOF

$progname v.$vers

EOF
  exit 1
}
#----------------------------------------------------------------------------------------- 

function run_blastp()
{
   local q db outfile
   q=$1
   db=$2
   outfile=$3
   
   blastp -query "$q" -db "$db" -matrix "$mat" -outfmt "$cols" -num_alignments $num_aln -num_threads $threads -qcov_hsp_perc "$qcov" \
    -seg "$seg" -soft_masking "$mask" -evalue "$Eval" > "$outfile"
}
#----------------------------------------------------------------------------------------- 


function print_help()
{
   # Prints the help message explaining how to use the script.
   cat <<EOF
   
   Usage: $progname -d <dir> [-e <ext>] [-E <Eval>] [-m <matrix>] [-n <num_aln>] [-q <qcov>] [-t <threads>] 
                              [-S <yes|no> SEG filtering] -M [<false|true> soft_masking] [-D] [-v]
   
   REQUIRED:
    -d <string> path to directory containing the input proteomes (protein FASTA)
    
   OPTIONAL
   -D <flag> print debugging info
   -e <string> fasta file extension name [def:$ext]
   -E <integer or float -ge 0> E-value [def:$Eval]
   -F <0|1> fix FASTA header for makeblastdb [def:$fix_header]
   -h <flag> print this help
   -n <int> number of blast alignments [def:$num_aln]
   -m <string> matrix name <BLOSUM45|BLOSUM62|BLOSUM80> [def:$mat]
   -M <false|true> soft_masking [def:$mask]
   -q <int> minimum query coverage percentage cutoff [def:$qcov]
   -S <yes|no> SEG filtering [def:$seg]
   -t <int> number of threads for blastp runs [def:$threads]
   -v <flag> print version
   
   EXAMPLES:
     $progname -d .
     $progname -d proteome_files -n 5 -q 85 -t 8
   
   AIMS: 
     1. Compute blastp reciprocal best hits (RBHs) 
        between a set of proteomes (protein FASTA files) 
	and a single reference 	proteome, which is 
	automatically selected as the smallest one. 
     2. Clusters (FASTA FILES) of core genome loci
	are computed
   
   SOURCE: the latest version can be fetched from 
           https://github.com/vinuesa/TIB-filoinfo
	   
   LICENSE: GPL v3.0. 
      See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
      
   NOTES: 
    1. For a versatile and highly customizable pan-genome analysis 
       software package, consider using GET_HOMOLOGUES
         https://doi.org/10.1128%2FAEM.02411-13
	 https://github.com/eead-csic-compbio/get_homologues
	 https://hub.docker.com/u/eeadcsiccompbio
	 
    2. To perform core- and/or pan-genome phylogenomic analyses
       consider using the GET_PHYLOMARKERS package
       https://doi.org/10.3389/fmicb.2018.00771 
       https://github.com/vinuesa/get_phylomarkers
       https://hub.docker.com/r/vinuesa/get_phylomarkers          
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

args=("$@")

while getopts ':d:e:E:F:m:M:n:q:r:S:t:hDv?:' OPTIONS
do
   case $OPTIONS in

   d)   proteome_dir=$OPTARG
        ;;
   e)   ext=$OPTARG
        ;;
   E)   Eval=$OPTARG
        ;;
   F)   fix_header=$OPTARG
        ;;
   h)   print_help
        ;;
   m)   mat=$OPTARG
        ;;
   M)   mask=$OPTARG
        ;;
   n)   num_aln=$OPTARG
        ;;
   q)   qcov=$OPTARG
        ;;
   r)   ref=$OPTARG
        ;;
   S)   seg=$OPTARG
        ;;
   t)   threads=$OPTARG
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

###>>> Exported variables !!!
declare -x g=$genome_ID perl # export only2perl!!!; call as $ENV{genome_ID}

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
 run on $hostn running $os with $nprocs processors and bash v.${bash_vers} at ${today/ /} 
   using the following parameters: 
 - proteome_dir=$proteome_dir | fasta_extension=$ext
 - BLASTP params: blastp v.${blast_vers} | num_aln=$num_aln | qcov=$qcov | threads=$threads |
                      Eval=$Eval | mat=$mat | seg=$seg | mask=$mask
 - reference=$ref_selection | fix_header=$fix_header | DEBUG=$DEBUG 
 - invocation: $progname ${args[*]}
===================================================================================================== 
 " >&2

# -------------------
# >>>> MAIN CODE <<<<
# -------------------
# Main script logic starts here

# ---------------------------------------
# 0. Validate user input & setup pipeline
# ---------------------------------------
# Check dependencies and Bash version
check_dependencies
v=$(check_bash_version "$min_bash_vers")
echo "# running under bash v.$v"

# Move to proteome directory
print_start_time && echo "# working in $proteome_dir"
cd "$proteome_dir" || error "Could not cd into $proteome_dir"

wkdir=$(pwd)

print_start_time && echo '# Selecting the smallest genome as the reference'
# automatically select the smallest reference, if not provided as ARG
if [[ -z "$ref" ]]; then
    ref=$(select_ref "$ext")
    echo "  - Selected $ref as the reference genome"
fi

echo '-----------------------------------------------------------------------------------------------------'

# keep the ref as first element in fasta arrays
declare -a non_ref infiles 
non_ref=( $(ls *"${ext}" | grep -v "$ref") )
infiles=("$ref" "${non_ref[@]}")

# ---------------------------------------------------------------------------
# 1. edit headers for makeblastdb (blast+): Formatting and Indexing Proteomes
# ---------------------------------------------------------------------------
# The FASTA headers of all proteomes are edited using Perl to add a unique identifier.
if ((fix_header == 1))
then
    print_start_time && echo "# Formatting ${#infiles[@]} input FASTA files for indexing"
    perl -pe 'if(/^>/){$c++; s/>/>lcl\|REF_$c /}' "$ref" > "${ref}ed"

    genome_ID=0
    g=0
    for f in "${non_ref[@]}"; do
        genome_ID=$(( genome_ID + 1 ))
        g="$genome_ID"
        perl -pe 'if(/^>/){$c++; s/>/>lcl\|GENO$ENV{g}\_$c /}' "$f" > "${f}ed"
    done
else
    for f in *"${ext}"
    do
        ln -s "$f" "${f}ed"   
    done
fi

declare -a faaed_files
faaed_files=($(ls *.faaed))

# -------------------------------------------------------------------
# 2. run makeblastdb on all input proteomes with edited FASTA headers
# -------------------------------------------------------------------
# Each proteome is used to create a blastp database using makeblastdb.
print_start_time && echo "# Generating indexed blastp databases"

for f in *"${ext}"ed
do 
    #Note: -parse_seqids results in query and subject IDs with different structure => lcl|A_DB_203  B_DB_166
    makeblastdb -in "$f" -dbtype prot -parse_seqids &> /dev/null
done


print_start_time && echo "# Generating the aliased blastp database allDBs ..."
blastdb_aliastool -dblist_file <(printf '%s\n' "${faaed_files[@]}") -dbtype prot -out allDBs -title allDBs

echo '-----------------------------------------------------------------------------------------------------'

#-----------------------------------
# 3. Run and process pairwise blastp 
#-----------------------------------
# For each non-reference proteome, blastp is run against the reference proteome 
#   (REFvsGENO) and vice versa (GENOvsREF).

genome_ID=0

for f in "${non_ref[@]}"; do
    genome_ID=$(( genome_ID + 1 ))
    
    ref_vs_geno_blastout=${ref%.*}vs${f%.*}_best_hits.tmp
    geno_vs_ref_blastout=${f%.*}vs${ref%.*}_best_hits.tmp
    
    print_start_time && echo "# Running: blastp -seg yes -soft_masking true -query ${ref}ed -db ${f}ed -qcov_hsp_perc $qcov -outfmt $cols -num_alignments $num_aln -num_threads $threads -evalue $Eval -matrix $mat > $ref_vs_geno_blastout"
    
    blastp -seg yes -soft_masking true -query "${ref}"ed -db "${f}"ed -qcov_hsp_perc "$qcov" -outfmt "$cols" -num_alignments "$num_aln" \
      -num_threads "$threads" -evalue "$Eval" -matrix "$mat" > "$ref_vs_geno_blastout"
    
    check_file "$ref_vs_geno_blastout"

   # Retrieve the best nonREF proteome database hits using blastdbcmd, onlfy if qcov > \$qcov
   print_start_time && echo "# Retrieving the best hits from $ref_vs_geno_blastout with blastdbcmd ..."
   blastdbcmd -entry_batch <(awk -F"\t" -v qcov=$qcov '$8 > qcov{print $2}' "$ref_vs_geno_blastout" | sort -u) -db "${f}"ed > "${ref%.*}vs${f%.*}"_besthits.faa
   
   check_file "${ref%.*}vs${f%.*}"_besthits.faa
   
   num_hits=$(grep -c '^>' "${ref%.*}vs${f%.*}"_besthits.faa)
   
   if ((num_hits == 0))
   then
        echo "WARNING: no hits in ${ref%.*}vs${f%.*}"_besthits.faa"
	rm ${ref%.*}vs${f%.*}"_besthits.faa
	continue
   fi

    print_start_time && echo "# Running: blastp -seg yes -soft_masking true -query ${ref%.*}vs${f%.*}_besthits.faa -db ${ref}ed -qcov_hsp_perc $qcov -outfmt $cols -num_alignments $num_aln -num_threads $threads -evalue $Eval -matrix $mat > $geno_vs_ref_blastout"
    
    blastp -seg yes -soft_masking true -query "${ref%.*}vs${f%.*}"_besthits.faa -db "${ref}"ed -qcov_hsp_perc "$qcov" -outfmt "$cols" -num_alignments "$num_aln" \
      -num_threads "$threads" -evalue "$Eval" -matrix "$mat"> "$geno_vs_ref_blastout"
      
    check_file "$geno_vs_ref_blastout"

    # Sort the blastp output table from the preceding search by increasing E-values (in column 9) and decreasing scores (col 10)
    #    & filter out unique REF vs nonREF RBHs using AWK hashes from the sorted blast output table with qcov > $qcov
    print_start_time && echo "# Filtering out unique REF vs nonREF RBHs from the sorted blast output table with qcov > $qcov"
    for GENOid in $(cut -f1 "$geno_vs_ref_blastout" | sort -u)
    do 
       grep "$GENOid" "$ref_vs_geno_blastout"
    done | sort -gk9,9 -gk10,10 | \
           awk -v qcov=$qcov 'BEGIN{FS=OFS="\t"}{REFid[$1]++; GENOid[$2]++; if(REFid[$1] == 1 && GENOid[$2] == 1 && $8 > qvov) print }' > \
             "${ref_vs_geno_blastout%.*}"_RBHs_qcov_gt"${qcov}".tsv
    
    check_file "${ref_vs_geno_blastout%.*}"_RBHs_qcov_gt"${qcov}".tsv ]]
done

echo '-----------------------------------------------------------------------------------------------------'


#-----------------------------------------------------------------------
# 4. Identify REF proteins shared by all tsv files holding pairwise RBHs
#-----------------------------------------------------------------------
# Find the intersections of REFs in all tsv files
print_start_time && echo "# Computing the intersections of REF proteins in all tsv files holding pairwise RBHs"
awk '{r[$1]++; if(r[$1] == ARGC-1) print $1}' ./*.tsv > REF_RBH_IDs.list
[[ ! -s REF_RBH_IDs.list ]] && error "could not write REF_RBH_IDs.list"

intersection_size=$(wc -l REF_RBH_IDs.list | awk '{print $1}')
((intersection_size > 0)) && echo "# Found $intersection_size core RBHs shared by ${#non_ref[*]} nonREF proteomes with the $ref reference proteome"
((intersection_size == 0)) && error "# ERROR: found $intersection_size core orhtologous genes among ${#infiles[*]} input proteomes ..."

echo '-----------------------------------------------------------------------------------------------------'

#-------------------------------------------------------------------------------------------------------------------------
# 5. Loop over tsv files and generate RBHs_matrix.tsv core_genome_clusters.tsv and nonCore_genome_clusters.tsv tables
#-------------------------------------------------------------------------------------------------------------------------
# Cluster Computation

print_start_time && echo "# Computing clusters of homologous sequences"

# The core_ref hash counts the instances of the REFERNCE_IDs in the RBH tables
declare -A core_ref

# The pangenome_clusters hash is indexed by REFERNCE_IDs 
#   and as its value holds the RBHs as a tab-separated
#   string of IDs from nonREF proteomes    
declare -A pangenome_clusters

# Construct the pangenome_clusters hash, indexed by reference proteome, 
#  containing a value a string of tab-separated nonREF proteome RBH IDs.
print_start_time && echo "# Populating the pangenome_clusters hash ..."
for t in *RBHs_*.tsv; do
    while read -r REF QUERY rest
    do
        # count the instances of REFERNCE_IDs in each RBHs_*.tsv tables
	(( core_ref["$REF"]++ ))
	if [[ ${core_ref["$REF"]} -eq 1 ]]
	then
	    pangenome_clusters["$REF"]="$QUERY" 
	else
	    pangenome_clusters["$REF"]="${pangenome_clusters[$REF]}\t$QUERY"
	fi
    done < "$t"
done


# 5.1 print the RBHs_matrix.tsv
print_start_time && echo "# Printing the pangenome_matrix, core_genome_clusters and nonCore_genome_clusters files"
for key in "${!pangenome_clusters[@]}"
do
      echo -e "${key}\t${pangenome_clusters[$key]}"
done > RBHs_matrix.tsv
check_file RBHs_matrix.tsv

# 5.2 print the core_genome_clusters.tsv
awk -v ninfiles="${#infiles[@]}" 'NF == ninfiles' RBHs_matrix.tsv > core_genome_clusters.tsv
check_file core_genome_clusters.tsv

# 5.3 print the nonCore_genome_clusters.tsv
awk -v ninfiles="${#infiles[@]}" 'NF != ninfiles' RBHs_matrix.tsv > nonCore_genome_clusters.tsv
check_file nonCore_genome_clusters.tsv warn

echo '-----------------------------------------------------------------------------------------------------'

#-----------------------------
# 6. Write cluster FASTA files
#-----------------------------
# Below is the code for three strategies, two of them commented out
#  The final (uncommented) one is the fastest of the benchmarked strategies. 
print_start_time && echo "# Extracting RBH cluster FASTA files ..."

## >>> Take 1 6.1 read all source FASTA file into memory (as the hash seqs) for later filtering
# print_start_time && echo "# reading all source FASTA files into memory (as the hash seqs) ..."
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


print_start_time && echo "# Reading RBHs_matrix.tsv and extracting cluster proteins with blastdbcmd ..."
## >>> Take 2: 6.1 This is the blastdb-based approach, which is a bit faster than the 
##   hash traversing and filtering approach shown above, but not much more
## read each line of the RBHs_matrix.tsv
##  into the ids array using the read -r -a idiom,
##  and extract the ids from the aliased allDBs blastdb

#declare -a entries
#while read -r -a ids
#do 
#    ((c++)) 
#    # iterate over all indexes of the idx array
#    #  and append them to the entries array
#    for (( idx=0; idx <= ((${#ids[@]} -1)); idx++))
#    do
#        entries+=("${ids[$idx]}")
#    done 
#    # extract the entries from the aliased allDBs input proteomes
#    #  note the use of stdbuf to unbuffer the blasdbcmd output, in an attempt to speed it up,
#    #  but seems to have no effect
#    #  https://stackoverflow.com/questions/3465619/how-to-make-output-of-any-shell-command-unbuffered
#    stdbuf -i0 -o0 -e0 blastdbcmd -db allDBs -dbtype prot -entry_batch <(printf '%s\n' "${entries[@]}") -out cluster_"${c}".fas
#    unset -v entries
#done < RBHs_matrix.tsv

# >>> Take 3, 6.1 uses blastdbcmd but calling it in parallel with xargs 
# write the each line of the RBHs_matrix.tsv IDs to a tmpfile
#  to pass the list of tmpfiles to a parallel call of blastdbcmd
print_start_time && echo "# Writing each line of the RBHs_matrix.tsv IDs to an idstmp file por parallelization ..."

#initialize cluster counter
c=0
while read -r -a ids
do 
    ((c++)) 
    # write each line of the RBHs_matrix.tsv to a temporal file
    printf '%s\n' "${ids[@]}" > cluster_${c}.idstmp
done < RBHs_matrix.tsv

## 6.2 pass the list of tmpfiles to a parallel call of blastdbcmd
##  This works nicely, but must ensure that parallel is available on host
#ls *.idstmp | parallel --gnu -j "$threads" 'blastdbcmd -db allDBs -dbtype prot -entry_batch {} -out {.}.fas'  

# 6.2 use the more portable find | xargs idiom to parallelize the blastdbcmd call
print_start_time && echo "# Running blastdbcmd in parallel with xargs using -P $threads ..."
find . -name '*.idstmp' -print0 | xargs -0 -P "$threads" -I % blastdbcmd -db allDBs -dbtype prot -entry_batch % -out %.fas


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
         mv "$f" nonCore_clusters/nonCore_"${f}" || error "could mv $f to nonCore_clusters/nonCore_${f}"
    fi
done

echo '-----------------------------------------------------------------------------------------------------'


#----------------
# 7. final cleanup
#----------------
#Unnecessary files generated during the process are removed.
print_start_time && echo "# Tidying up $wkdir ..."

echo '-----------------------------------------------------------------------------------------------------'

((DEBUG == 0)) && rm ./*faaed ./*faaed.* ./*best_hits.tmp ./REF_RBH_IDs.list ./*besthits.faa ./*.idstmp

elapsed=$(( SECONDS - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days, %H hr, %M min, %S sec')"

print_end_message

exit 0



