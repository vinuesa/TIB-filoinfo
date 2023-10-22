#!/usr/bin/env bash

#: progname: compute_blastp_RBH_orthologous_clusters.sh
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
vers=0.6_2023-08-26 # improved code modularization and comments.       

min_bash_vers=4.4 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

set -o pipefail

# Initialize variables with default values
DEBUG=0
qcov=80
num_aln=10
threads=4
ext=faa

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

function filter_best_hits()
{
    # Processes blastp output to retain best hits for each query.
    # sorts the hits for each query, to identify the one with the 
    # hightest bitscore (column 12 of standard tabular blastp -m 6),
    # which is passed to filter_best_hits as its only parameter
    
    local infile n s 
    
    # Declare arrays to store the names and best hits
    declare -a names 
    declare -A max best

    infile=$1

    # Read the input file and parse the fields
    while read -r n f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 s; 
    do
        if [[ -z "${max[$n]}" ]]
        then
             names+=( "$n" )
        fi

        # Compare the score for the current name to the current max value
        # float comparisons using bc
        if [[ ! "${max[$n]}" ]] || [[ 1 -eq "$(echo "$s > ${max[$n]}" | bc)" ]]
        then
            max["$n"]="$s"
            unset 'best[$n]'
        fi

        # If the current score is equal to the current max value, store the line
        # float comparisons using bc
        if [[ 1 -eq "$(echo "$s == ${max[$n]}" | bc)" ]]
        then
             best[$n]="${n}\t$f2\t$f3\t$f4\t$f5\t$f6\t$f7\t$f8\t$f9\t${f10}\t${f11}\t$s"
        fi
    done < "$infile"

    # Print the best hits for each name
    for n in "${names[@]}"; do
        [[ -z "${best[$n]}" ]] && continue
        echo -e "${best[$n]}"
    done
}
#----------------------------------------------------------------------------------------- 

# Function to check dependencies
function check_dependencies() {
    local dependencies=("fas2tab.pl" "tab2fas.pl" "blastp" "makeblastdb")
    
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

function print_help()
{
   # Prints the help message explaining how to use the script.
   cat <<EOF
   
   Usage: $progname -d <dir> [-e <ext>] [-n <num_aln>] [-q <qcov>] [-t <threads>] [-D] [-v]
   
   REQUIRED:
    -d <string> path to directory containing the input proteomes (protein FASTA)
    
   OPTIONAL
   -D <flag> print debugging info
   -e <string> fasta file extension name [def:$ext]
   -h <flag> print this help
   -n <int> number of blast alignments [def:$num_aln]
   -q <int> minimum query coverage percentage cutoff [def:$qcov]
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

while getopts ':d:e:n:q:r:t:hDv?:' OPTIONS
do
   case $OPTIONS in

   d)   proteome_dir=$OPTARG
        ;;
   e)   ext=$OPTARG
        ;;
   h)   print_help
        ;;
   n)   num_aln=$OPTARG
        ;;
   q)   qcov=$OPTARG
        ;;
   r)   ref=$OPTARG
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
 - BLASTP params: blastp v.${blast_vers} | num_aln=$num_aln | qcov=$qcov | threads=$threads 
 - reference=$ref | DEBUG=$DEBUG 
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
check_bash_version "$min_bash_vers"

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
print_start_time && echo "# formatting ${#infiles[@]} input FASTA files for indexing"
perl -pe 'if(/^>/){$c++; s/>/>lcl\|REF_$c /}' "$ref" > "${ref}ed"

genome_ID=0
g=0
for f in "${non_ref[@]}"; do
    genome_ID=$(( genome_ID + 1 ))
    g="$genome_ID"
    perl -pe 'if(/^>/){$c++; s/>/>lcl\|GENO$ENV{g}\_$c /}' "$f" > "${f}ed"
done

# -------------------------------------------------------------------
# 2. run makeblastdb on all input proteomes with edited FASTA headers
# -------------------------------------------------------------------
# Each proteome is used to create a blastp database using makeblastdb.
print_start_time && echo "# Generating indexed blastp databases"
for f in *"${ext}"ed
do 
    #Note: -parse_seqids results in query and subject IDs with different structure => lcl|A_DB_203	B_DB_166
    makeblastdb -in "$f" -dbtype prot -parse_seqids &> /dev/null
    
    # The tabular output of the fas2tab.pl script is saved for each proteome.
    fas2tab.pl "$f" | sed '/^$/d' > "${f}"tab
done
echo '-----------------------------------------------------------------------------------------------------'

# declar hashes holding [query]-target and [target]-query hits
declare -A AB_hits BA_hits 

#-----------------------------------
# 3. Run and process pairwise blastp 
#-----------------------------------
genome_ID=0

# declare hashes to  hold non-redundant query and target hits
declare -A seen seen2

# For each non-reference proteome, blastp is run against the reference proteome 
#   (REFvsGENO) and vice versa (GENOvsREF).
for f in "${non_ref[@]}"; do
    genome_ID=$(( genome_ID + 1 ))
    print_start_time && echo "# running: blastp -seg yes -soft_masking true -query ${ref}ed -db ${f}ed -qcov_hsp_perc $qcov -outfmt 6 -num_alignments $num_aln -num_threads $threads > REFvsGENO${genome_ID}"
    blastp -seg yes -soft_masking true -query "${ref}"ed -db "${f}"ed -qcov_hsp_perc "$qcov" -outfmt 6 -num_alignments "$num_aln" -num_threads "$threads" > REFvsGENO"${genome_ID}"

    print_start_time && echo "# running: blastp -seg yes -soft_masking true -query ${f}ed -db ${ref}ed -qcov_hsp_perc $qcov -outfmt 6 -num_alignments $num_aln -num_threads $threads > GENO${genome_ID}vsREF"
    blastp -seg yes -soft_masking true -query "${f}"ed -db "${ref}"ed -qcov_hsp_perc "$qcov" -outfmt 6 -num_alignments "$num_aln" -num_threads "$threads" > GENO"${genome_ID}"vsREF


    # 3.1 The best hits are extracted using the filter_best_hits function.
    #      Retrieve the highest-scoring hit out of the -num_alignments $num_aln hits
    # Note: makeblastdb with -parse_seqids produces subject cols without the lcl| prfix: => lcl|A_DB_203	B_DB_166
    #        and therefore we remove them to have subject and query IDs with the same structure,
    #        required for hash key comparisons later in the code
    print_start_time && echo "# filter_best_hits REFvsGENO${genome_ID} > REFvsGENO${genome_ID}.best"
    filter_best_hits REFvsGENO"${genome_ID}" | sed 's#lcl|##' > REFvsGENO"${genome_ID}".best

    print_start_time && echo "# filter_best_hits GENO${genome_ID}vsREF > GENO${genome_ID}vsREF.best"
    filter_best_hits GENO"${genome_ID}"vsREF | sed 's#lcl|##' > GENO"${genome_ID}"vsREF.best


    # 3.2. Reciprocal best hits (RBHs) are determined by comparing the hits between the two blastp runs.
    #      Compute A vs B reciprocal best hits: 
    #      the following hashes will hold query=>subject and subject=>query results
    #      that will be used to identify the reicprocal best hits (RBHs)
    print_start_time && echo "# Computing REF vs GENO reciprocal best hits @ qcov=$qcov"
    
    #REF_1	GENO1_92	100.000	234	0	0	5	238	1	234	1.18e-180	488
    #REF_2	GENO1_141	100.000	278	0	0	1	278	1	278	0.0	569
    while read -r q subj rest
    do
	AB_hits["$q"]="$subj"
    done < REFvsGENO"${genome_ID}".best

    
    #GENO1_1	REF_52	99.441	179	1	0	10	188	1	179	1.65e-131	359
    #GENO1_2	REF_53	100.000	95	0	0	1	95	1	95	6.70e-68	191
    while read -r q subj rest
    do
        (( seen[$q]++ ))
        if (( ${seen[$q]} == 1 )); then
	    BA_hits["$subj"]="$q"
	fi
    done < GENO"${genome_ID}"vsREF.best
    
    # 3.3 find the RBHs by comparing the two hashes and sort output by A_DB_# key
    RBH_outfile="${ref%.*}"_vs_"${f%.*}"_qcov"${qcov}"_RBHs_2col.tsv
    for kA in "${!AB_hits[@]}"
    do
        (( seen2[${AB_hits[$kA]}]++ ))
	if [[ "${BA_hits[$kA]}" == "${AB_hits[$kA]}" ]] && (( "${seen2[${AB_hits[$kA]}]}" == 1 ))
        then
             printf "%s\t%s\n" "$kA" "${AB_hits[$kA]}"
        else
	     continue
	fi
    done | sort -k1.6g > "$RBH_outfile"


    # 3.4 make a more informative RBH_outfile name, including number of RBHs found
    num_RBHs=$(wc -l "$RBH_outfile" | awk '{print $1}')
    RBH_out="${ref%.*}"_vs_"${f%.*}"_qcov"${qcov}"_${num_RBHs}RBHs_2col.tsv
    mv "$RBH_outfile" "$RBH_out"
    RBH_outfile="$RBH_out"
    print_start_time && echo "# Found $num_RBHs reciprocal best hits between $ref and $f at qcov=${qcov}%"
    echo '==='
done

echo '-----------------------------------------------------------------------------------------------'

#--------------------------------------------------------------------
# 4. Loop over tsv files and count the number of hits for each REF_ID
#--------------------------------------------------------------------
# Cluster Computation and Output

print_start_time "# Computing clusters of homologous sequences"

# A hash (core_ref) counts the occurrences of each reference ID in the RBH tables.
declare -A core_ref # hash counting the instances of the REFERNCE_IDs in the RBH tables


# For each RBH table, the RBH clusters are extracted and stored in cluster_X.fastab files.
for f in *RBHs_2col.tsv; do
    while read -r REF QUERY
    do
        # count the instances of REFERNCE_IDs in each RBHs_2col.tsv tables
	(( core_ref["$REF"]++ ))
    done < "$f"
done

((DEBUG > 0 )) && echo "# DEBUG: Loop over tsv files and count the number of hits for each REF_ID; \${#non_ref[@]}: ${#non_ref[@]}" >&2 \
               &&  for k in "${!core_ref[@]}"; do (("${core_ref[$k]}" == "${#non_ref[@]}")) && echo -e "$k\t${core_ref[$k]}"; done >&2

# A list of reference IDs shared with all genomes is created (ref_orth_IDs).
declare -a ref_orth_IDs # array holding the REFERNCE_IDs found in all RBH tables

# ref_orth_IDs contains the reference_IDs of orthologous loci
#   shared by the REFERENCE with all genomes; 
ref_orth_IDs=( $(for k in "${!core_ref[@]}"; do (("${core_ref[$k]}" == "${#non_ref[@]}")) && printf '%s\n' "$k"; done) )

# sort the ref_orth_IDs array; note sort -k1.5g for IDs like REF_58
ref_orth_IDs=( $(printf '%s\n' "${ref_orth_IDs[@]}" | sort -k1.5g) ) 
((DEBUG > 0 )) && echo "# DEBUG: the sorted  \${ref_orth_IDs[@]} array:" >&2 \
               && echo "${ref_orth_IDs[@]}" >&2


echo "# found ${#ref_orth_IDs[@]} RBH clusters for ${#non_ref[@]} genomes against $ref"
printf '%s\n' "${ref_orth_IDs[@]}" > "${ref%.*}"_core_IDs.list

for f in *RBHs_2col.tsv; do
    base="${f%RBHs_2col.tsv}"
    for id in "${ref_orth_IDs[@]}"
    do
          grep -w "$id" "$f" | cut -f2
    done > "${base}"_core_IDs.list
done

print_start_time && echo "# Generating a table of ortholog IDs shared with $ref" 
paste ./*core_IDs.list > ORTHOLOG_IDs_SHARED_WITH_REF.tsv

# ORTHOLOG_IDs_SHARED_WITH_REF.tsv
#REF_1	GENO1_496	GENO2_843
#REF_3	GENO1_498	GENO2_845
#REF_4	GENO1_500	GENO2_847

#--------------------------
# 5. Write cluster fastas
#--------------------------
# The clusters are split into "core" and "non-core" based on the number of genomes sharing the locus.
print_start_time && echo '# Generating FASTA files for orthologous clusters'

declare -a non_ref_faaedtab_files faaedtab_files
ref_faaedtab=$(ls "${ref}"edtab)
non_ref_faaedtab_files=( $(ls *."${ext}"edtab | grep -v "${ref}"edtab) )
# put the ref_faaedtab in the first position idx[0] of the faaedtab_files array
faaedtab_files=( "$ref_faaedtab" "${non_ref_faaedtab_files[@]}" )


# 5.1 reconstitute cluster fastas
# NOTE: calling blasdbcmd multiple times in a loop to retrieve 
#       the cluster sequences is quite slow. 
#       Therefore they're grepped out of the *edtab files
declare -A seen3
c=0
while read -r -a ids; do 
    c=$((c + 1))
    # grep out the reference sequence, stored at idx[0]
    grep -w "${ids[0]}" "${faaedtab_files[0]}" > cluster_"${c}".fastab
    
    # next process the non ref sequences stored at idxs >= 1
    #  and append them to the growing cluster_"${c}".fastab file
    for (( idx=1; idx <= ((${#ids[@]} -1)); idx++)); do
	(( seen3[${ids[$idx]}]++ ))
	((DEBUG > 0)) && echo "DEBUG: idx:$idx; c=$c; grep -w ${ids[$idx]} ${faaedtab_files[$idx]}" >&2
        if (( ${seen3[${ids[$idx]}]} == 1 )); then
	    grep -w "${ids[$idx]}" "${faaedtab_files[$idx]}" >> cluster_"${c}".fastab
	fi    
    done
done < ORTHOLOG_IDs_SHARED_WITH_REF.tsv

((DEBUG == 0)) && [[ -s ORTHOLOG_IDs_SHARED_WITH_REF.tsv ]] && rm ./*core_IDs.list
[[ ! -s ORTHOLOG_IDs_SHARED_WITH_REF.tsv ]] && error "could not write ORTHOLOG_IDs_SHARED_WITH_REF.tsv"

# Convert fastab files back to FASTA
for f in cluster_*.fastab
do
     tab2fas.pl "$f" > "${f%tab}"
done


# Cluster (RBH) FASTA files are moved to appropriate directories (core and non_core).
clusters_dir=${ref%.*}_vs_${#non_ref_faaedtab_files[@]}genomes_qcov"${qcov}"_"${#ref_orth_IDs[@]}"_clusters
[[ -d "$clusters_dir" ]] && rm -rf "$clusters_dir"
[[ ! -d "$clusters_dir" ]] && mkdir "$clusters_dir" || error "could not generate dir $clusters_dir"

print_start_time && echo "# moving cluster FASTA files and tables to $clusters_dir"
mv cluster*.fas GENO*.best REF*.best ORTHOLOG_IDs_SHARED_WITH_REF.tsv "$clusters_dir"

cd "$clusters_dir" || error "could not cd into $clusters_dir"
mkdir core || error "could not mkdir core"
mkdir non_core || error "could not mkdir non_core"

core_loci=0
non_core_loci=0
for f in *fas; do 
   if grep -c '^>' "$f" | grep "${#infiles[@]}" &> /dev/null; then 
      core_loci=$((core_loci + 1))
      mv "$f" core
   else
      non_core_loci=$((non_core_loci + 1))
      mv "$f" non_core
   fi
done
print_start_time && echo "# Moved $core_loci core FASTA files to $clusters_dir/core"
print_start_time && echo "# Moved $non_core_loci non-core FASTA files to $clusters_dir/non_core"


# ----------------
# 6. final cleanup
# ----------------
# Unnecessary files generated during the process are removed.
cd "$wkdir" || error "could not cd into $wkdir"
print_start_time && echo "# final cleanup in $wkdir"
((DEBUG == 0)) && rm ./*fastab GENO*vsREF REFvsGENO[[:digit:]]* ./*"${ext}"ed* ./*RBHs_2col.tsv

echo

elapsed=$(( SECONDS - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days, %H hr, %M min, %S sec')"

print_end_message