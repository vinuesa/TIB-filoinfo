#!/usr/bin/env bash

#: phyml_protModelFinder.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex
#: AIM: simple wraper script around phyml, to select a good model for protein alignments
#:      compute AIC, BIC, delta_BIC and BICw, and estimate a ML phylogeny using the best-fitting model
#: LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
 
#: Desgin:  phyml_protModelFinder.sh evaluates a set of the named empirical substitution matices 
#      currently implemented in phml v3.*, combining them or not with +G and/or +f
# - Amino-acid based models : LG (default) | WAG | JTT | MtREV | Dayhoff | DCMut | RtREV | CpREV | VT | AB
#		              Blosum62 | MtMam | MtArt | HIVw |  HIVb | custom

# set bash strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
#set -euo pipefail # fails in some calls on tepeu with old phyml and bash; on buluc or tepeu use set -uo pipefail
set -euo pipefail

host=$(hostname)

progname=${0##*/}
version=0.3 # 2022-11-25;
min_bash_vers=4.3 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

declare -a models       # array holding the base models (empirical substitution matrices to be evaluated)
declare -A model_cmds   # hash holding the commands for each model 
declare -A model_scores # hash holding the model lnL scores and AICi values 
declare -A model_options

# array of models to evaluate
nuclear_genome_models=(AB BLOSUM62 DAYHOFF DCMut JTT LG VT WAG)
organelle_genemome_models=(CpREV MTMAM MtREV MtArt)
nuclear_and_organellar=(AB BLOSUM62 DAYHOFF DCMut JTT LG VT WAG CpREV MTMAM MtREV MtArt)
viral_genome_models=(HIVw HIVb RtREV)
all_models=(AB BLOSUM62 DAYHOFF DCMut JTT LG VT WAG CpREV MTMAM MtREV MtArt HIVw HIVb RtREV)

# hash mapping option => model_set
model_options['1']='nuclear_genome_models'
model_options['2']='organelle_genemome_models'
model_options['3']='nuclear_and_organellar_models'
model_options['4']='viral_genome_models'
model_options['5']='all_models'

#==============================#
# >>> FUNCTION DEFINITIONS <<< #
#------------------------------#
function check_bash_version()
{
   local bash_vers min_bash_vers
   min_bash_vers=$1
   bash_vers=$(bash --version | head -1 | awk '{print $4}' | sed 's/(.*//' | cut -d. -f1,2)
   awk -v bv="$bash_vers" -v mb="$min_bash_vers" \
     'BEGIN { if (bv < mb){print "FATAL: you are running acient bash v"bv, "and version >=", mb, "is required"; exit 1}else{print "# Bash version", bv, "OK"} }'
}
#-----------------------------------------------------------------------------------------

function check_is_phylip()
{
   local phylip
   phylip=$1
   
   if ! awk 'NR==1 && NF==2' "$phylip" &> /dev/null; then 
       echo "FATAL ERROR: input file $phylip does not seem to by a canonical phylip alingment"
       print_help
   fi
}
#-----------------------------------------------------------------------------------------

function print_start_time()
{
   #echo -n "[$(date +%T)] "
   printf '%(%T )T' '-1' # requires Bash >= 4.3
}
#-----------------------------------------------------------------------------------------

function compute_AICi()
{
   local score n_branches extra_params total_params
   
   score=$1
   n_branches=$2
   extra_params=$3
   
   total_params=$((n_branches + extra_params))
 
   # AICi=-2*lnLi + 2*Ni
   echo "(-2 * $score) + (2 * $total_params)" | bc -l
}
#-----------------------------------------------------------------------------------------

function compute_AICc()
{
   local AIC score n_branches extra_params total_params
   
   score=$1
   n_branches=$2
   extra_params=$3
   n_sites=$4
   AIC=$5
   
   total_params=$((n_branches + extra_params))
 
   # AICi=-2*lnLi + 2*Ni
   #AIC=$( echo "(-2 * $score) + (2 * $total_params)" | bc -l )   
   #echo $AIC + (2 * $total_params($total_params + 1)/($n_sites - $total_params -1)) | bc
   
   echo "$AIC + ( 2 * ($total_params * ($total_params + 1))/($n_sites - $total_params -1) )" | bc -l
}
#-----------------------------------------------------------------------------------------

function compute_BIC()
{
   local score n_branches extra_params total_params n_sites
   
   score=$1
   n_branches=$2
   extra_params=$3
   n_sites=$4
 
   total_params=$((n_branches + extra_params))

   # BICi= k*ln(n) -2*lnLi
   awk -v lnL="$score" -v k="$total_params" -v n="$n_sites" 'BEGIN{ BIC= (-2 * lnL) + (k * log(n)); printf "%0.5f", BIC }'
}
#-----------------------------------------------------------------------------------------

function print_help(){

   cat <<-EoH

   $progname v${version} requires two arguments provided on the command line:
   
   $progname <string [input phylip file (aligned PROTEIN sequences)> <int [1-4; default=1]>
    
      # model sets to choose from: 
      1 -> nuclear genes (AB BLOSUM62 DAYHOFF DCMut JTT LG VT WAG)
      2 -> organellar genes (CpREV MTMAM MtREV MtArt)
      3 -> nuclear and organellar (1 + 2)
      4 -> retroviral genes (HIVw HIVb RtREV)
      5 -> all (1+2+3+4)
   
   AIM: $progname v${version} will evaluate the fit of the the seleced 
   	 model set, combined or not with +G and/or +f, computing AIC, BIC
	 deltaBIC, BICw and inferring the ML tree under the BIC-selected
	 model	
    	 
   LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE     

EoH
   
   exit 0

}
#-----------------------------------------------------------------------------------------
#============================= END FUNCTION DEFINITIONS ==================================
#=========================================================================================

# ============ #
# >>> MAIN <<<
# ------------ #

# 0. Check that the input file was provided, and that the host runs bash >= v4.3
(( $# != 2 )) && print_help

infile="$1"
model_set="$2"

wkd=$(pwd)

# verify input & bash vesion
[[ ! -s "$infile" ]] && echo "FATAL ERROR: could not find $infile in $wkd" && exit 1
(( model_set < 1 )) || ((model_set > 5 )) && print_help
check_is_phylip "$infile"

# OK, ready to start the analysis ...
start_time=$SECONDS
echo "================================================================================"
check_bash_version "$min_bash_vers"
echo -n "# $progname v$version run on $host started on: "; printf '%(%F at %T)T\n' '-1'
echo "================================================================================"
echo

# 1. get sequence stats
no_seq=$(awk 'NR == 1{print $1}' "$infile") 
no_sites=$(awk 'NR == 1{print $2}' "$infile") 
no_branches=$((2 * no_seq - 3))

# 2. set the selected model set, making a copy of the set array into the models array
case "$model_set" in
   1) models=( "${nuclear_genome_models[@]}" ) ;;
   2) models=( "${organelle_genemome_models[@]}" ) ;;
   3) models=( "${nuclear_and_organellar[@]}" );;
   4) models=( "${viral_genome_models[@]}" );;
   5) models=( "${all_models[@]}" );;
   *) echo "unknown model set!" && print_help ;;
esac
   
# 3. Compute a fast NJ tree estimating distances with the LG matrix 
print_start_time 
echo "1. Computing NJ-LG tree for input file $infile with $no_seq sequences"
echo '--------------------------------------------------------------------------------'
phyml -i "$infile" -d aa -m LG -c 1 -b 0 -o n &> /dev/null

# 2. rename the outfile for future use as usertree
if [[ -s "${infile}"_phyml_tree.txt ]]; then
   mv "${infile}"_phyml_tree.txt "${infile}"_LG-NJ.nwk
else
    echo "FATAL ERROR: could not compute ${infile}_phyml_tree.txt" && exit 1
fi

# 4. run a for loop to combine all base models with (or not) +G and or +f
#     and fill the model_scores and model_cmds hashes
echo "2. running in a for loop to combine all base models in model_set ${model_set}=>${model_options[$model_set]}, with (or not) +G and or +f"
for mat in "${models[@]}"; do
     print_start_time && echo "# running: phyml -i $infile -d aa -m $mat -u ${infile}_LG-NJ.nwk -c 1 -v 0 -o lr &> /dev/null"
     phyml -i "$infile" -d aa -m "$mat" -u "${infile}"_LG-NJ.nwk -c 1 -o lr &> /dev/null 
     extra_params=0 
     total_params=$((no_branches + extra_params))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$extra_params")
     AICc=$(compute_AICc "$score" "$no_branches" "$extra_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$no_branches" "$extra_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat"

     print_start_time && echo "# running: phyml -i $infile -d aa -m $mat -c 4 -a e -u ${infile}_LG-NJ.nwk -o lr &> /dev/null"
     phyml -i "$infile" -d aa -m "${mat}" -c 4 -a e -u "${infile}"_LG-NJ.nwk -o lr &> /dev/null     
     extra_params=1 
     total_params=$((no_branches + extra_params))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$extra_params")
     AICc=$(compute_AICc "$score" "$no_branches" "$extra_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$no_branches" "$extra_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -c 4 -a e"

     print_start_time && echo "# running: phyml -i $infile -d aa -m $mat -f e -c 1 -u ${infile}_LG-NJ.nwk -o lr &> /dev/null"
     phyml -i "$infile" -d aa -m "$mat" -f e -c 1 -u "${infile}"_LG-NJ.nwk -o lr &> /dev/null     
     extra_params=19 #19 from AA frequencies
     total_params=$((no_branches + extra_params))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$extra_params")
     AICc=$(compute_AICc "$score" "$no_branches" "$extra_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$no_branches" "$extra_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+f"]="$model_string"
     model_cmds["${mat}+f"]="$mat -f e"

     print_start_time && echo "# running: phyml -i $infile -d aa -m $mat -u ${infile}_LG-NJ.nwk -f e -a e -o lr &> /dev/null"
     phyml -i "$infile" -d aa -m "$mat" -u "${infile}"_LG-NJ.nwk -f e -a e -c 4 -o lr &> /dev/null     
     extra_params=20 #19 from AA frequencies + 1 gamma 
     total_params=$((no_branches + extra_params))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$extra_params")
     AICc=$(compute_AICc "$score" "$no_branches" "$extra_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$no_branches" "$extra_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+f+G"]="$model_string"
     model_cmds["${mat}+f+G"]="$mat -f e -c 4 -a e"
done
echo '--------------------------------------------------------------------------------------------------'
echo

# 5. print a sorted summary table of model fits from the model_scores hash
print_start_time

echo "# writing ${infile}_sorted_model_set_${model_set}_fits.tsv, sorted by BIC"
for m in "${!model_scores[@]}"; do
    echo -e "$m\t${model_scores[$m]}"
done | sort -nk7 > "${infile}"_sorted_model_set_"${model_set}"_fits.tsv


# 6 compute delta_BIC and BICw, based on "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
declare -a BIC_a
declare -a BIC_deltas_a
declare -a BICw_a
BIC_a=( $(awk '{print $7}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv) )
min_BIC="${BIC_a[0]}"

# 6.1 fill BIC_deltas_a array
BIC_deltas_a=()
for i in "${BIC_a[@]}"
do
     BIC_deltas_a+=( $( echo "$i" - "$min_BIC" | bc -l) )
done

# 6.2 paste the BIC_deltas_a values as a new column to "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
paste "${infile}"_sorted_model_set_"${model_set}"_fits.tsv <(for i in "${BIC_deltas_a[@]}"; do echo "$i"; done) > t
[[ -s t ]] && mv t "${infile}"_sorted_model_set_"${model_set}"_fits.tsv

# 6.3 Compute the BICw_sums (denominator) of BICw
BICw_sums=0
for i in "${BIC_deltas_a[@]}"; do 
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }')  
   #echo "num:$BICw_numerator"
   BICw_sums=$(bc <<< $BICw_sums'+'$BICw_numerator)
done
#echo BICw_sums:$BICw_sums

# 6.4 fill the BICw_a array
BICw_a=()
for i in "${BIC_deltas_a[@]}"; do
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }' 2> /dev/null)   
   BICw=$(echo "$BICw_numerator / $BICw_sums" | bc -l)
   BICw_a+=( $(printf "%.2f" "$BICw") )
done

# 6.5 paste the BICw_a values as a new column to "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
paste "${infile}"_sorted_model_set_"${model_set}"_fits.tsv <(for i in "${BICw_a[@]}"; do echo "$i"; done) > t
[[ -s t ]] && mv t "${infile}"_sorted_model_set_"${model_set}"_fits.tsv


# 6.6 Display  the final "${infile}"_sorted_model_set_"${model_set}"_fits.tsv and extract the best model name
if [[ -s "${infile}"_sorted_model_set_"${model_set}"_fits.tsv ]]; then
    # display models sorted by BIC
    best_model=$(awk 'NR == 1{print $1}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv)
    [[ -z "$best_model" ]] && echo "FATAL ERROR: unbound \$best_model at $LINENO" && exit 1

    # print table with header to STDOUT and save to file
    awk 'BEGIN{print "model\tK\tsites/K\tlnL\tAIC\tAICc\tBIC\tdeltaBIC\tBICw"}{print}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv | column -t
    awk 'BEGIN{print "model\tK\tsites/K\tlnL\tAIC\tAICc\tBIC\tdeltaBIC\tBICw"}{print}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv > t
    mv t "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
else
    echo "ERROR: could not write ${infile}_sorted_model_set_${model_set}_fits.tsv"
fi

# cleanup: remove phyml output files from the last pass through the loop
[[ -s "${infile}"_phyml_stats.txt ]] && rm "${infile}"_phyml_stats.txt
[[ -s "${infile}"_phyml_tree.txt ]] && rm "${infile}"_phyml_tree.txt


# 7. compute ML tree under best-fitting model
echo '--------------------------------------------------------------------------------------------------'
echo "* NOTE: when sites/K < 40, the AICc is recommended over AIC."
echo
echo "... will estimate the ML tree under best-fitting model $best_model selected by BIC"

print_start_time

echo "# running: phyml -i $infile -d aa -m ${model_cmds[$best_model]} -o tlr -s BEST  &> /dev/null"

# note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the comand fail ??!!
phyml -i "$infile" -d aa -m "${model_cmds[$best_model]}" -o tlr -s BEST  &> /dev/null

# 7.1 Check and rename final phyml output files
if [[ -s "${infile}"_phyml_stats.txt ]]; then
     mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_phyml_stats.txt
     echo "# wrote ${infile}_${best_model}_phyml_stats.txt"
else
     echo "FATAL ERROR: ${infile}_phyml_stats.txt was not generated!"
fi

if [[ -s "${infile}"_phyml_tree.txt ]]; then
     mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_phyml_tree.txt
     echo "# wrote ${infile}_${best_model}_phyml_tree.txt"
else
     echo "FATAL ERROR: ${infile}_phyml_tree.txt was not generated!"
fi
echo '--------------------------------------------------------------------------------------------------'

echo

elapsed=$(( SECONDS - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days, %H hr, %M min, %S sec')"

echo 'Done!'
echo
exit 0
