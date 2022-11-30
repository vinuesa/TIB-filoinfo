#!/usr/bin/env bash

#: phyml_DNAmodelFinder.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#: AIM: simple wraper script around phyml, to select a reasonable substitution model for DNA alignments
#:      compute AIC, BIC, delta_BIC and BICw, and estimate a ML phylogeny using the best-fitting model
#: LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE
 
#: Design: phyml_DNAmodelFinder.sh evaluates the named parametric substitution models
#     currently implemented in phml v3.*, combining them or not with +G and/or +I
#      - Nucleotide-based models : HKY85 (default) | JC69 | K80 | F81 | F84 | TN93 | GTR 
#     Under runmode 2, the set is significantly expanded, by adding equal|unequal frquency 
#         model sets, which are automatically selected based on delta_BIC between JC69 and F81
#     The models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates, in order
#        to calulate each model's AIC, BIC, delta_BIC and BICw. 
#     The best model is selected by BIC


#: GitHub repo: you can fetch the latest version of the script from:
#   https://github.com/vinuesa/TIB-filoinfo/blob/master/phyml_DNAmodelFinder.sh
# wget -c https://raw.githubusercontent.com/vinuesa/TIB-filoinfo/master/phyml_DNAmodelFinder.sh

# set Bash's unofficial strict mode
set -euo pipefail
host=$(hostname)

progname=${0##*/}
version='0.3_2022-11-28' 
min_bash_vers=4.3 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

# declare array and hash variables
declare -a models             # array holding the base models (empirical substitution matrices to be evaluated)
declare -A model_cmds         # hash holding the commands for each model 
declare -A model_scores       # hash holding the model lnL scores and AICi values 
declare -A model_options      # hash mapping option => model_set
declare -A model_free_params  # hash mapping base models with their free parameters

# standard_models
model_free_params['JC69']=0
model_free_params['K80']=1
model_free_params['F81']=3
model_free_params['F84']=4
model_free_params['HKY85']=4
model_free_params['TN93']=5
model_free_params['GTR']=8

# extended_models_ef
model_free_params['010020ef']=2  # TNef
model_free_params['012230ef']=3  # TIMef
model_free_params['012314ef']=4  # TVMef
model_free_params['012310ef']=3  # TVM1ef 
model_free_params['010213ef']=3  # TVM2ef
model_free_params['012213ef']=3  # TVM3ef
model_free_params['012345ef']=5  # SYMef
model_free_params['012013ef']=3  # TVM4ef
model_free_params['010012ef']=2  # TVM5ef 
model_free_params['012012ef']=2  # TVM6ef 
model_free_params['010212ef']=2  # TVM7ef 
model_free_params['012313ef']=3  # TVM8ef 
model_free_params['010023ef']=3  # TIM1ef
model_free_params['012234ef']=4  # TIM2ef
model_free_params['012232ef']=2  # TIM3ef

# extended_models_uf
model_free_params['012210']=5  # K81uf 
model_free_params['012230']=6  # TIM
model_free_params['012314']=7  # TVM
model_free_params['012310']=6  # TVM1
model_free_params['010213']=6  # TVM2
model_free_params['012213']=6  # TVM3
model_free_params['012013']=6  # TVM4
model_free_params['010012']=5  # TVM5 
model_free_params['012012']=5  # TVM6 
model_free_params['010212']=5  # TVM7 
model_free_params['012313']=6  # TVM8 
model_free_params['010023']=6  # TIM1
model_free_params['012234']=7  # TIM2
model_free_params['012232']=6  # TIM3

# array of models to evaluate
standard_models=(JC69 K80 F81 HKY85 TN93 GTR)

#                    TNef   TIMef  TVMef TVM1ef TVM2ef TVM3ef SYMef  TVM4ef TVM5ef TVM6ef TVM7ef TVM8ef TIM1ef TIM2ef TIM3ef
extended_models_ef=(010020 012230 012314 012310 010213 012213 012345 012013 010012 012012 010212 012313 010023 012234 012232)

#                   K81uf   TIM    TVM   TVM1   TVM2   TVM3   TVM4   TVM5   TVM6   TVM7   TVM8   TIM1   TIM2   TIM3
extended_models_uf=(012210 012230 012314 012310 010213 012213 012013 010012 012012 010212 012313 010023 012234 012232)   
                                                       
test_models=(HKY85 GTR)

# hash mapping option => model_set
model_options['1']='standard_models'
model_options['2']='extended_models'
model_options['3']='test_models'

#==============================#
# >>> FUNCTION DEFINITIONS <<< #
#------------------------------#

function check_dependencies()
{
    declare -a progs required_binaries optional_binaries
    local p programname
    
    required_binaries=(awk bc sed perl phyml)
    optional_binaries=(mpirun phyml-mpi)
    
    for p in "${optional_binaries[@]}"
    do
          if type -P "$p" >/dev/null
	  then
	      progs=("${optional_binaries[@]}")
	  else
	      progs=()
	  fi
    done
    
    progs+=("${required_binaries[@]}")
    
    for programname in "${progs[@]}"
    do
       if ! type -P "$programname"; then  # NOTE: will print paths of binaries to STDOUT (no >/dev/null)
          echo
          echo "$# ERROR: $programname not in place!"
          echo "  ... you will need to install \"$programname\" first, or include it in \$PATH"
          echo "  ... exiting"
          exit 1
       else
          continue
       fi
    done
    
    echo
    echo '# Run check_dependencies() ... looks good: all required binaries are in place.'
    echo
}
#-----------------------------------------------------------------------------------------

function check_bash_version()
{
   local bash_vers min_bash_vers
   min_bash_vers=$1
   bash_vers=$(bash --version | head -1 | awk '{print $4}' | sed 's/(.*//' | cut -d. -f1,2)
   
   echo "$bash_vers"
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

function compute_nt_freq_in_phylip()
{
  local phylip
  phylip=$1
 
  awk '
  BEGIN{print "idx\tNT\tobs_freq"}
  {
    # ignore first row and column
    if( NR > 1 && NF > 1){
       # remove empty spaces
       gsub(/[ ]+/," ")
       l=length($0)

       for(i=1; i<=l; i++){
          c = substr($0, i, 1)
         
	  # count only standard amino acids
	  if (c ~ /[ACGT]/){
              ccounts[c]++
              letters++
          }
       }
    }
  }
  # print relative frequency of each residue
  END {
     for (c in ccounts){ 
        aa++ 
        printf "%i\t%s\t%.4f\n", aa, c, (ccounts[c] / letters )
     }	
  }' "$phylip"
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
   # AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
   local score n_branches total_params
   
   score=$1
   n_branches=$2
   total_params=$3
    
   # AICi=-2*lnLi + 2*Ni
   echo "(-2 * $score) + (2 * $total_params)" | bc -l
}
#-----------------------------------------------------------------------------------------

function compute_AICc()
{
   # AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
   local score extra_params total_params n_sites AIC
   
   score=$1
   total_params=$2
   n_sites=$3
   AIC=$4
 
   # AICi=-2*lnLi + 2*Ni
   #AIC=$( echo "(-2 * $score) + (2 * $total_params)" | bc -l )   
   #echo $AIC + (2 * $total_params($total_params + 1)/($n_sites - $total_params -1)) | bc
   
   echo "$AIC + ( 2 * ($total_params * ($total_params + 1))/($n_sites - $total_params -1) )" | bc -l
}
#-----------------------------------------------------------------------------------------

function compute_BIC()
{
   # BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
   local score extra_params total_params n_sites
   
   score=$1
   total_params=$2
   n_sites=$3
 
   # BICi= k*ln(n) -2*lnLi
   awk -v lnL="$score" -v k="$total_params" -v n="$n_sites" 'BEGIN{ BIC= (-2 * lnL) + (k * log(n)); printf "%0.5f", BIC }'
}
#-----------------------------------------------------------------------------------------

function check_compositional_heterogeneity()
{
    # returns 1|0 if there is or not significant compositional_heterogeneity based on delta_AIC > 2
    local JC_BICi F81_BICi uf_models_flag
    
    JC_BICi=$1
    F81_BICi=$2
    
    if [[ $(echo "$JC_BICi - $F81_BICi" | bc -l | cut -d. -f1) -gt 2 ]]; then 
        uf_models_flag=1
    else
        uf_models_flag=0
    fi

    echo "$uf_models_flag"
}
#-----------------------------------------------------------------------------------------

function print_help(){

   bash_vers=$(check_bash_version "$min_bash_vers")
   bash_ge_5=$(awk -v bv="$bash_vers" 'BEGIN { if (bv >= 5.0){print 1}else{print 0} }')
   
   if ((bash_ge_5 > 0)); then 
   cat <<-EoH

   $progname v${version} requires two arguments provided on the command line:
   
   $progname <string [input phylip file (aligned DNA sequences)> <int [model sets:1-3]>
    
      # model sets to choose from: 
      1 -> (JC69 K80 F81 HKY85 TN93 GTR)
     
                                    TNef     TIMef    TVMef    TVM1ef   TVM2ef   TVM3ef   SYMef    TVM4ef   TVM5ef   TVM6ef
      2 -> 1 + extended_ef_models: (010020ef 012230ef 012314ef 123421ef 121324ef 123324ef 012345ef 012013ef 010012ef 012012ef 
      				   TVM7ef   TVM8ef   TIM1ef   TIM2ef   TIM3ef
                                   010212ef 012313ef 010023ef 012234ef 012232ef)
             OR  
	                                    K81uf   TIM    TVM   TVM1   TVM2   TVM3   TVM4   TVM5   TVM6   TVM7   TVM8   TIM1   TIM2   TIM3
           1 + extended_uf_models: (012210 012230 012314 012310 010213 012213 012013 010012 012012 010212 012313 010023 012234 012232)   

           $progname automatically chooses the proper extended set (ef|uf) to evaluate, 
                 based on delta_BIC evaluation of compositional bias (JC69 vs F81)
       
      3 -> minimal test set (K80 F81 HKY85 GTR)
   
   AIM: * $progname v${version} will evaluate the fit of the seleced model set,
            combined or not with +G and/or +I, computing AIC, BICm, delta_BIC, BICw 
	    and inferring the ML tree under the BIC-selected model	
        * The models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates, 
	    in order to calculate each model\'s AIC, BIC, delta_BIC and BICw. 
        * The best model is selected by BIC
    	 
   SOURCE: the latest version of the program is available on GitHub:
            https://github.com/vinuesa/TIB-filoinfo
   
   LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 
   
EoH

   else
   cat <<-EoH

   $progname v${version} requires two arguments provided on the command line:
   
   $progname <string [input phylip file (aligned DNA sequences)> <int [model sets:1|3]>
    
      # model sets to choose from: 
      1 -> (JC69 K80 F81 HKY85 TN93 GTR)
      2 -> WILL NOT RUN properly on Bash < v5.0, sorry (see NOTE below) 
      3 -> minimal test set (K80 F81 HKY85 GTR)
   
   AIM: * $progname v${version} will evaluate the fit of the seleced model set,
            combined or not with +G and/or +I, computing AIC, BICm, delta_BIC, BICw 
	    and inferring the ML tree under the BIC-selected model	
        * The models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates, 
	    in order to calculate each model\'s AIC, BIC, delta_BIC and BICw. 
        * The best model is selected by BIC
    	 
   SOURCE: the latest version of the program is available on GitHub:
            https://github.com/vinuesa/TIB-filoinfo
   
   LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 
   
   NOTE: you are running the old Bash version $bash_vers. 
         Update to v 5.0 or greater, to profit from the full set of models 
	         and features implemented in $progname
   
EoH
   
   fi
   
   exit 0

}
#-----------------------------------------------------------------------------------------
#============================= END FUNCTION DEFINITIONS ==================================
#=========================================================================================

# ============ #
# >>> MAIN <<<
# ------------ #

## Check environment
# 0. Check that the input file was provided, and that the host runs bash >= v4.3
(( $# < 2 )) || (( $# > 3 )) && print_help

infile="$1"
model_set="$2"
#num_threads=${3:-6}

wkd=$(pwd)

check_is_phylip "$infile"

# OK, ready to start the analysis ...
start_time=$SECONDS
echo "========================================================================================="
bash_vers=$(check_bash_version "$min_bash_vers")
awk -v bv="$bash_vers" -v mb="$min_bash_vers" \
  'BEGIN { if (bv < mb){print "FATAL: you are running acient bash v"bv, "and version >=", mb, "is required"; exit 1}else{print "# Bash version", bv, "OK"} }'

echo -n "# $progname v$version running on $host. Run started on: "; printf '%(%F at %T)T\n' '-1'
echo "# workding directory: $wkd"
check_dependencies
echo "# infile:$infile model_set:$model_set => ${model_options[$model_set]}" 
echo "========================================================================================="
echo

# 1. get sequence stats
print_start_time
echo " # 1. Computing sequence stats for ${infile}:"

no_seq=$(awk 'NR == 1{print $1}' "$infile") 
echo "- number of sequences: $no_seq"

no_sites=$(awk 'NR == 1{print $2}' "$infile") 
echo "- number of sites: $no_sites"

no_branches=$((2 * no_seq - 3))
echo "- number of branches: $no_branches"

echo "- observed nucleotide frequencies:"
compute_nt_freq_in_phylip "$infile"
echo '--------------------------------------------------------------------------------'
echo 


# 2. set the selected model set, making a copy of the set array into the models array
case "$model_set" in
   1) models=( "${standard_models[@]}" ) ;;
   2) models=( "${standard_models[@]}" ) ;; # plus automatic selection of extended_models to evaluate
   3) models=( "${test_models[@]}" ) ;;
   *) echo "unknown model set!" && print_help ;;
esac


# 3. Compute a fast NJ tree estimating distances with the JC model
print_start_time 
echo "1. Computing NJ-JC tree for input file $infile with $no_seq sequences"
echo '--------------------------------------------------------------------------------'
phyml -i "$infile" -d nt -m JC69 -c 1 -b 0 -o n &> /dev/null

# 4. rename the outfile for future use as usertree
if [[ -s "${infile}"_phyml_tree.txt ]]; then
   mv "${infile}"_phyml_tree.txt "${infile}"_JC-NJ.nwk
else
    echo "FATAL ERROR: could not compute ${infile}_phyml_tree.txt" && exit 1
fi

# 5.1 run a for loop to combine all base models with (or not) +G and or +I
#     and fill the model_scores and model_cmds hashes
echo "2.1. running in a for loop to combine all base models in model_set ${model_set}=>${model_options[$model_set]},
      with (or not) +G and or +I, and compute the model lnL, after optimizing branch lengths and rates"

# globals for compositional_heterogenety check
JC_BICi=0
F81_BICi=0
compositional_heterogenety=''
freq_cmd=''
    
for mat in "${models[@]}"; do
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     phyml -i "$infile" -d nt -m "$mat" -u "${infile}"_JC-NJ.nwk -c 1 -o lr &> /dev/null 
     extra_params=0 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     
     # save the JC_BICi to test for significant compositional heterogeneity
     if ((model_set == 2)); then
         if [[ "$mat" == 'JC69' ]]; then
	     JC_BICi="$AICi"
	     echo "JC_BICi: $JC_BICi"
	 elif [[ "$mat" == 'F81' ]]; then
	     F81_BICi="$AICi"
	     echo "F81_BICi: $F81_BICi"
	 fi
     fi
     
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat"

     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "${mat}" -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr &> /dev/null
     extra_params=1 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -c 4 -a e"

     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "$mat" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr &> /dev/null
     extra_params=1 # 1 pInv
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I"]="$model_string"
     model_cmds["${mat}+I"]="$mat -v e"

     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
     phyml -i "$infile" -d nt -m "$mat" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr &> /dev/null
     extra_params=2 #19 from AA frequencies + 1 gamma 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I+G"]="$model_string"
     model_cmds["${mat}+I+G"]="$mat -v e -c 4 -a e"
done

# cleanup: remove phyml output files from the last pass through the loop
[[ -s "${infile}"_phyml_stats.txt ]] && rm "${infile}"_phyml_stats.txt
[[ -s "${infile}"_phyml_tree.txt ]] && rm "${infile}"_phyml_tree.txt


##### 5.2 if extended set, then compute and set the compositional_heterogenety flag 
#           and loop over corresponding set of extended models

      
# check_compositional_heterogeneity and set compositional_heterogenety flag (1|0), accordingly
if ((model_set == 2)); then
   if [[ -n "$JC_BICi" ]] && [[ -n "$F81_BICi" ]]; then
       compositional_heterogenety=$(check_compositional_heterogeneity "$JC_BICi" "$F81_BICi") 
       echo '--------------------------------------------------------------------------------'
       print_start_time
       echo '# Starting evaluation and automatic selection of extended model set'
       echo "# setting compositional_heterogenety flag to: $compositional_heterogenety"
   fi

   # fill the models array with the proper ones, 
   #   based on the compositional_heterogenety flag
   models=()
   if ((compositional_heterogenety == 1)); then
   	echo '# will evaluate models with unequal frequencies'
	models=( "${extended_models_uf[@]}" )
        freq_cmd="-f e"
   else
   	echo '# will evaluate models with equal frequencies'
   	models=( "${extended_models_ef[@]}" )	
   	freq_cmd="-f 0.25,0.25,0.25,0.25"
   fi

   echo '--------------------------------------------------------------------------------'
   print_start_time

   echo "2.2. running a loop to combine the extended models in model_set ${model_set}=>${model_options[$model_set]},
      with (or not) +G and or +I, and compute the model lnL, after optimizing branch lengths and rates"

   # 5.2 loop over the set of extended models, passing the proper freq_cmd, based on the compositional_heterogenety flag
   for mat in "${models[@]}"; do
     ((compositional_heterogenety == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -c 1 -o lr &> /dev/null 
     extra_params=0 
     ((compositional_heterogenety == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")     
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat"

     ((compositional_heterogenety == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "${mat}" "$freq_cmd" -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr &> /dev/null
     extra_params=1 # 1 gamma 
     ((compositional_heterogenety == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -c 4 -a e"

     ((compositional_heterogenety == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr &> /dev/null
     extra_params=1 # 1 pInv
     ((compositional_heterogenety == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I"]="$model_string"
     model_cmds["${mat}+I"]="$mat -v e"

     ((compositional_heterogenety == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
     phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr &> /dev/null
     extra_params=2 # 1 pInv + 1 gamma 
     ((compositional_heterogenety == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-/{print $3}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I+G"]="$model_string"
     model_cmds["${mat}+I+G"]="$mat -v e -c 4 -a e"
   done
fi # extended set

echo '--------------------------------------------------------------------------------'

echo
##### 

# 6. print a sorted summary table of model fits from the model_scores hash
print_start_time

echo "# writing ${infile}_sorted_model_set_${model_set}_fits.tsv, sorted by BIC"
echo '--------------------------------------------------------------------------------'
for m in "${!model_scores[@]}"; do
    echo -e "$m\t${model_scores[$m]}"
done | sort -nk7 > "${infile}"_sorted_model_set_"${model_set}"_fits.tsv

# 7. compute delta_BIC and BICw, based on "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
declare -a BIC_a
declare -a BIC_deltas_a
declare -a BICw_a
BIC_a=( $(awk '{print $7}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv) )
min_BIC="${BIC_a[0]}"

# 7.1 fill BIC_deltas_a array
BIC_deltas_a=()
for i in "${BIC_a[@]}"
do
     BIC_deltas_a+=( $( echo "$i" - "$min_BIC" | bc -l) )
done

# 7.2 Compute the BICw_sums (denominator) of BICw
BICw_sums=0
for i in "${BIC_deltas_a[@]}"; do 
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }')  
   #echo "num:$BICw_numerator"
   BICw_sums=$(bc <<< "$BICw_sums"'+'"$BICw_numerator")
done
#echo BICw_sums:$BICw_sums

# 7.3 fill the BICw_a array
BICw_a=()
for i in "${BIC_deltas_a[@]}"; do
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }' 2> /dev/null)   
   BICw=$(echo "$BICw_numerator / $BICw_sums" | bc -l)
   BICw_a+=( $(printf "%.2f" "$BICw") )
done

# 7.4 paste the BIC_deltas_a & BICw_a values as a new column to "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
paste "${infile}"_sorted_model_set_"${model_set}"_fits.tsv <(for i in "${BIC_deltas_a[@]}"; do echo "$i"; done) \
                                                           <(for i in "${BICw_a[@]}"; do echo "$i"; done) > t
[[ -s t ]] && mv t "${infile}"_sorted_model_set_"${model_set}"_fits.tsv

# 7.5 Display  the final "${infile}"_sorted_model_set_"${model_set}"_fits.tsv and extract the best model name
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


#--------------------------------------------
# 8. compute ML tree under best-fitting model
#--------------------------------------------
echo '--------------------------------------------------------------------------------------------------'
echo "* NOTE 1: when sites/K < 40, the AICc is recommended over AIC."
echo "* NOTE 2: Best model selected by BIC, because AIC is biased in favour of parameter-rich models."
echo
echo '=================================================================================================='
echo '#  Will estimate the ML tree under best-fitting model $best_model selected by BIC'
echo '--------------------------------------------------------------------------------------------------'

print_start_time

echo "# running: phyml -i $infile -d nt -m ${model_cmds[$best_model]} -o tlr -s BEST"

# note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the comand fail
phyml -i "$infile" -d nt -m ${model_cmds[$best_model]} -o tlr -s BEST &> /dev/null

# 8.1 Check and rename final phyml output files
if [[ -s "${infile}"_phyml_stats.txt ]]; then
     mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_phyml_stats.txt
     echo "# Your results:"
     echo "  - ${infile}_${best_model}_phyml_stats.txt"
else
     echo "FATAL ERROR: ${infile}_phyml_stats.txt was not generated!"
fi

if [[ -s "${infile}"_phyml_tree.txt ]]; then
     mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_phyml_tree.txt
     echo "  - ${infile}_${best_model}_phyml_tree.txt"
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
