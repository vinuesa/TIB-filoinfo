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
#   https://github.com/vinuesa/TIB-filoinfo/blob/master/phyml_DNAmodelFinder.sh
# wget -c https://raw.githubusercontent.com/vinuesa/TIB-filoinfo/master/phyml_DNAmodelFinder.sh
#----------------------------------------------------------------------------------------

# set Bash's unofficial strict mode
set -euo pipefail
host=$(hostname)

progname=${0##*/}
version='0.9.1_2023-11-12' # v0.9.1_2023-11-12; fixed unbound variable vers (should be version) in print_end_message
min_bash_vers=4.4 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

n_starts=1         # seed trees
delta_BIC_cutoff=2 # to set compositional_heterogeneity, transitional_heterogeneity and pInv flags 

# declare array and hash variables
declare -a models             # array holding the base models (empirical substitution matrices to be evaluated)
declare -A model_cmds         # hash holding the commands for each model 
declare -A model_scores       # hash holding the model lnL scores and AICi values 
declare -A model_options      # hash mapping option => model_set
declare -A model_free_params  # hash mapping base models with their free parameters

declare -a TIMef TVMef TIMuf TVMuf # array holding the model codes

# To be implemented
#named_model_codes[JC]=000000
#named_model_codes[F81]=000000
#named_model_codes[K80]=010010
#named_model_codes[HKY]=010010
#named_model_codes[TN]=010020
#named_model_codes[TNe]=010020
#named_model_codes[K81]=012210  # TPM1
#named_model_codes[K81u]=012210 # TPM1
#named_model_codes[TPM2]=010212
#named_model_codes[TPM2u]=010212
#named_model_codes[TPM3]=012012
#named_model_codes[TPM3u]=012012
#named_model_codes[TIM]=012230
#named_model_codes[TIMe]=012230
#named_model_codes[TIM2]=010232
#named_model_codes[TIM2e]=010232
#named_model_codes[TIM3]=012032
#named_model_codes[TIM3e]=012032
#named_model_codes[TVM]=012314
#named_model_codes[TVMe]=012314
#named_model_codes[SYM]=012345
#named_model_codes[GTR]=012345

# standard_models # 6
model_free_params['JC69']=0
model_free_params['K80']=1
model_free_params['F81']=3
model_free_params['HKY85']=4
model_free_params['TN93']=5
model_free_params['GTR']=8

### extended_models_ef - TVM* # 15
model_free_params['012210ef']=2  # K81
model_free_params['012314ef']=4  # TVMef
model_free_params['012310ef']=3  # TVM1ef 
model_free_params['010213ef']=3  # TVM2ef
model_free_params['012213ef']=3  # TVM3ef
model_free_params['012013ef']=3  # TVM4ef
model_free_params['010012ef']=2  # TVM5ef 
model_free_params['012012ef']=2  # TVM6ef 
model_free_params['010212ef']=2  # TVM7ef 
model_free_params['010210ef']=2  # TVM8ef 
model_free_params['012313ef']=3  # TVM9ef 
model_free_params['010011ef']=1  # TVM10ef
model_free_params['012212ef']=2  # TVM11ef
model_free_params['011010ef']=1  # TVM12ef
model_free_params['001101ef']=1  # TVM13ef

TVMef[0]=012210
TVMef[1]=012314
TVMef[2]=012310
TVMef[3]=010213
TVMef[4]=012213
TVMef[5]=012013
TVMef[6]=010012
TVMef[7]=012012
TVMef[8]=010212
TVMef[9]=010210
TVMef[10]=012313
TVMef[11]=010011
TVMef[12]=012212
TVMef[13]=011010
TVMef[14]=001101

### extended_models_ef - TNef|TIM*|SYM # 26
model_free_params['010020ef']=2  # TNef
model_free_params['012230ef']=3  # TIMef
model_free_params['012345ef']=5  # SYMef
model_free_params['010023ef']=3  # TIM1ef
model_free_params['010232ef']=4  # TIM2ef
model_free_params['012232ef']=3  # TIM3ef
model_free_params['012332ef']=3  # TIM4ef
model_free_params['012342ef']=4  # TIM5ef
model_free_params['012343ef']=4  # TIM6ef
model_free_params['012340ef']=4  # TIM7ef
model_free_params['010021ef']=2  # TIM8ef
model_free_params['010022ef']=2  # TIM9ef
model_free_params['011123ef']=3  # TIM10ef
model_free_params['012223ef']=3  # TIM11ef
model_free_params['010120ef']=2  # TIM12ef
model_free_params['000120ef']=2  # TIM13ef
model_free_params['000121ef']=2  # TIM14ef
model_free_params['001021ef']=2  # TIM15ef
model_free_params['012234ef']=4  # TIM16ef
model_free_params['010231ef']=3  # TIM17ef
model_free_params['011230ef']=3  # TIM18ef
model_free_params['011020ef']=2  # TIM19ef
model_free_params['012130ef']=3  # TIM20ef
model_free_params['010121ef']=2  # TIM21ef
model_free_params['010122ef']=2  # TIM22ef
model_free_params['010123ef']=3  # TIM23ef

TIMef[0]=010020
TIMef[1]=012230
TIMef[2]=012345
TIMef[3]=010023
TIMef[4]=010232
TIMef[5]=012232
TIMef[6]=012332
TIMef[7]=012342
TIMef[8]=012343
TIMef[9]=012340
TIMef[10]=010021
TIMef[11]=010022
TIMef[12]=011123
TIMef[13]=012223
TIMef[14]=010120
TIMef[15]=000120
TIMef[16]=000121
TIMef[17]=001021
TIMef[18]=012234
TIMef[19]=010231
TIMef[20]=011230
TIMef[21]=011020
TIMef[22]=012130
TIMef[23]=010121
TIMef[24]=010122
TIMef[25]=010123

## extended_models_uf TVM* # 16
model_free_params['012210']=5  # K81uf 
model_free_params['012314']=7  # TVM
model_free_params['012310']=6  # TVM1
model_free_params['010213']=6  # TVM2
model_free_params['012213']=6  # TVM3
model_free_params['012013']=6  # TVM4
model_free_params['010012']=5  # TVM5 
model_free_params['012012']=5  # TVM6 
model_free_params['010212']=5  # TVM7 
model_free_params['012313']=6  # TVM8 
model_free_params['010011']=4  # TVM9
model_free_params['012212']=5  # TVM10 
model_free_params['010210']=5  # TVM11 
model_free_params['011010']=4  # TVM12
model_free_params['001101']=4  # TVM13 

TVMuf[0]=012210
TVMuf[1]=012314
TVMuf[2]=012310
TVMuf[3]=010213
TVMuf[4]=012213
TVMuf[5]=012013
TVMuf[6]=010012
TVMuf[7]=012012
TVMuf[8]=010212
TVMuf[9]=012313
TVMuf[10]=010011
TVMuf[11]=012212
TVMuf[12]=010210
TVMuf[13]=011010
TVMuf[14]=001101

## extended_models_uf TIM* # 25
model_free_params['012230']=6  # TIM
model_free_params['010023']=6  # TIM1
model_free_params['010232']=7  # TIM2
model_free_params['012232']=6  # TIM3
model_free_params['012332']=6  # TIM4
model_free_params['012342']=7  # TIM5
model_free_params['012343']=7  # TIM6
model_free_params['012340']=7  # TIM7
model_free_params['010021']=5  # TIM8
model_free_params['010022']=5  # TIM9
model_free_params['010120']=5  # TIM10
model_free_params['011123']=6  # TIM11
model_free_params['012223']=6  # TIM12
model_free_params['012222']=5  # TIM13 
model_free_params['000120']=5  # TIM14
model_free_params['000121']=5  # TIM15
model_free_params['001021']=5  # TIM16
model_free_params['012234']=7  # TIM17
model_free_params['010231']=6  # TIM18 
model_free_params['011020']=5  # TIM19
model_free_params['012130']=6  # TIM20
model_free_params['010121']=5  # TIM21
model_free_params['010122']=5  # TIM22
model_free_params['010123']=6  # TIM23
model_free_params['012345']=8  # GTR

TIMuf[0]=012230
TIMuf[1]=010023
TIMuf[2]=010232
TIMuf[3]=012232
TIMuf[4]=012332
TIMuf[5]=012342
TIMuf[6]=012343
TIMuf[7]=012340
TIMuf[8]=010021
TIMuf[9]=010022
TIMuf[10]=010120
TIMuf[11]=011123
TIMuf[12]=012223
TIMuf[13]=012222
TIMuf[14]=000120
TIMuf[15]=000121
TIMuf[16]=001021
TIMuf[17]=012234
TIMuf[18]=010231
TIMuf[19]=011020
TIMuf[20]=012130
TIMuf[21]=010121
TIMuf[22]=010122
TIMuf[23]=010123
TIMuf[24]=012345


# array of models to evaluate
standard_models=(JC69 K80 F81 HKY85 TN93 GTR)

extended_models_ef=( "${TIMef[@]}" "${TVMef[@]}" ) # 26

extended_models_uf=( "${TIMuf[@]}" "${TVMuf[@]}" ) # 26
                                                       
test_models=(JC69 K80 F81 HKY85 TN93)

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
       if ! type -P "$programname"; then  # print paths of binaries to STDOUT
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
   bash_vers=$(bash --version | awk 'NR==1{print $4}' | sed 's/(.*//' | cut -d. -f1,2)

   # float comparisons using bc
   if [[ 1 -eq "$(echo "$bash_vers < $min_bash_vers" | bc)" ]]
   then
      echo "# FATAL: you are using the old bash v.${bash_vers}
               but $progname requires bash >= v${min_bash_vers}
	       to use hashes and other goodies"
      exit 1	       
   fi
   
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
    
    #if [[ $(echo "$JC_BICi - $F81_BICi" | bc -l | cut -d. -f1) -gt "$delta_BIC_cutoff" ]]; then 
    if [[ "$(echo "if (${JC_BICi} - ${F81_BICi} >= ${delta_BIC_cutoff}) 1" | bc)" -eq 1 ]]; then 
        uf_models_flag=1
    else
        uf_models_flag=0
    fi

    echo "$uf_models_flag"
}
#-----------------------------------------------------------------------------------------

function check_transitional_heterogeneity()
{
    # returns 1|0 if there is or not significant transitional_heterogeneity based on delta_AIC > 2
    local HKY85_BICi TN93_BICi ti_models_flag
    
    HKY85_BICi=$1
    TN93_BICi=$2
    
    ti_models_flag=''
        
    #if [[ $(echo "if (${HKY85_BICi - $TN93_BICi" | bc | cut -d. -f1) -gt "$delta_BIC_cutoff" ]]; then 
    if [[ "$(echo "if (${HKY85_BICi} - ${TN93_BICi} >= ${delta_BIC_cutoff}) 1" | bc)" -eq 1 ]]; then 
        ti_models_flag=1
    else
        ti_models_flag=0
    fi

    echo "$ti_models_flag"
}
#-----------------------------------------------------------------------------------------

function check_pInv()
{
    # returns 1|0 if there are or not a significant proportion of invariant sites
    #  based on delta_AIC > 2
    local HKY85I_BICi HKY85IG_BICi pInv_flag
    
    HKY85I_BICi=$1
    HKY85IG_BICi=$2
    
    #if [[ $(echo "$HKY85I_BICi - $HKY85IG_BICi" | bc -l | cut -d. -f1) -gt "$delta_BIC_cutoff" ]]; then 
    if [[ "$(echo "if (${HKY85I_BICi} - ${HKY85IG_BICi} >= ${delta_BIC_cutoff}) 1" | bc)" -eq 1 ]]; then 
        pInv_flag=1
    else
        pInv_flag=0
    fi

    echo "$pInv_flag"
}
#-----------------------------------------------------------------------------------------

function print_model_details()
{
   cat <<EoH

1 -> standard models (JC69 K80 F81 HKY85 TN93 GTR)

2 -> 1 + extended_ef_models and extended_uf_models
     the corresponding model subset being automatically selected
EoH

  echo "# number of models in extended_models_ef = ${#extended_models_ef[@]}"
  for m in "${!extended_models_ef[@]}"; do
      echo "$m => ${extended_models_ef[$m]}"
  done
  echo '--------------------------'
  echo "# number of models in extended_models_uf = ${#extended_models_uf[@]}"
  for m in "${!extended_models_uf[@]}"; do
      echo "$m => ${extended_models_uf[$m]}"
  done
  exit 0

}
#-----------------------------------------------------------------------------------------

function print_end_message()
{
   cat <<EOF
  ========================================================================================
  If you use $progname v.$version for your research,
  I would appreciate that you:
  
  1. Cite the code in your work as:   
  Pablo Vinuesa. $progname v.$version 
       https://github.com/vinuesa/TIB-filoinfo/blob/master/$progname
  
  2. Give it a like on the https://github.com/vinuesa/TIB-filoinfo/ repo
  
  Thanks!

EOF
}
#----------------------------------------------------------------------------------------- 

function print_help(){

   bash_vers=$(check_bash_version "$min_bash_vers")
   bash_ge_5=$(awk -v bv="$bash_vers" 'BEGIN { if (bv >= 5.0){print 1}else{print 0} }')
   
if ((bash_ge_5 > 0)); then 
   cat <<EoH

$progname v${version} requires two arguments provided on the command line, the third one being optional:

$progname <string [input phylip of aligned DNA sequences)> <int [model sets:1-3]> <int [no_rdm_starts; default:$n_starts]>
 
# model sets to choose from: 
1 -> standard models (JC69 K80 F81 HKY85 TN93 GTR)

2 -> standard + ${#extended_models_ef[@]} extended_ef_models, OR
     standard + ${#extended_models_uf[@]} extended_uf_models 
     			     
NOTE: $progname automatically chooses the proper extended set (ef|uf) to evaluate, 
        based on delta_BIC evaluation of compositional bias (JC69 vs F81)
 
3 -> minimal test set (JC69 F81 HKY85 TN93)

EXAMPLE: $progname primates.phy 2
 
AIM:  $progname v${version} will evaluate the fit of the the seleced model set,
	combined or not with +G and/or +f, computing AICi, BICi, deltaBIC, BICw 
     	  and inferring the ML tree under the BIC-selected model  

PROCEDURE:
 - Models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates 
      to calculate their AICi, BICi, delta_BIC and BICw
 - Only relevant matrices among the extended set are evaluated, based on delta_BIC
      comparisons between JC69-F81, to decide if ef|uf models should be evaluated
      and comparisons between KHY85-TN93, to determine if models with two Ti rates should 
      be evaluated
 - pInv is automatically excluded in the extended model set, 
	if the delta_BICi_HKY+G is =< $delta_BIC_cutoff when compared with the BIC_HKY+G+I
 - The best model is selected by BIC
 - SPR searches can be launched starting from multiple random trees
 - Default single seed tree searches use a BioNJ with BEST moves
     
SOURCE: the latest version of the program is available on GitHub:
	 https://github.com/vinuesa/TIB-filoinfo

LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 
   
EoH
      
else
   cat <<EoH

$progname v${version} requires two arguments provided on the command line, the third one being optional:

$progname <string [input phylip file (aligned DNA sequences)> <int [model sets:1|3]> <int [no_rdm_starts; default:$n_starts]>
 
   # model sets to choose from: 
   1 -> (JC69 K80 F81 HKY85 TN93 GTR)
   2 -> WILL NOT RUN properly on Bash < v5.0, sorry (see NOTE below) 
   3 -> minimal test set (JC K80 F81 HKY85 TN93)

AIM:  $progname v${version} will evaluate the fit of the the seleced model set,
	combined or not with +G and/or +f, computing AICi, BICi, deltaBIC, BICw 
     	  and inferring the ML tree under the BIC-selected model  

PROCEDURE
  - Models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates 
       to calculate their AICi, BICi, delta_BIC and BICw
  - The best model is selected by BIC
  - SPR searches can be launched starting from multiple random trees
  - Default single seed tree searches use a BioNJ with BEST moves     
      
SOURCE: the latest version of the program is available on GitHub:
	 https://github.com/vinuesa/TIB-filoinfo

LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 

NOTE: you are running the old Bash version $bash_vers. 
      Update to version >=5.0 to profit from the full set of models 
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
n_starts="${3:-1}"

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
echo "# infile:$infile; model_set:$model_set => ${model_options[$model_set]}; seed trees: $n_starts; delta_BIC_cutoff=$delta_BIC_cutoff" 
echo "========================================================================================="
echo ''

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
echo ''


# 2. set the selected model set, making a copy of the set array into the models array
case "$model_set" in
   1) models=( "${standard_models[@]}" ) ;;
   2) models=( "${test_models[@]}" ) ;; # plus automatic selection of extended_models to evaluate
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

# globals for compositional_heterogeneity check
declare -A seen
seen["HKY85+G"]=0
JC_BICi=0
F81_BICi=0
TN93_BICi=0
HKY85_BICi=0
HKY85I_BICi=0
compositional_heterogeneity=''
transitional_heterogeneity=''
use_pInv=''
freq_cmd=''
    
for mat in "${models[@]}"; do
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     phyml -i "$infile" -d nt -m "$mat" -u "${infile}"_JC-NJ.nwk -c 1 -o lr --leave_duplicates --no_memory_check &> /dev/null 
     extra_params=0 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat"

     # save the JC_BICi to test for significant compositional heterogeneity (JC vs F81)
     #  and transitional bias (HKY vs TN93
     if ((model_set == 2)); then
         if [[ "$mat" == 'JC69' ]]; then
	     JC_BICi="$BICi"
	     echo "JC_BICi: $JC_BICi"
	 elif [[ "$mat" == 'F81' ]]; then
	     F81_BICi="$BICi"
	     echo "F81_BICi: $F81_BICi"
	 elif [[ "$mat" == 'HKY85' ]]; then
	     HKY85_BICi="$BICi"
	     echo "HKY85_BICi: $HKY85_BICi"
	 elif [[ "$mat" == 'TN93' ]]; then
	     TN93_BICi="$BICi"
	     echo "TN93_BICi: $TN93_BICi"
	 fi
     fi
     
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "${mat}" -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     extra_params=1 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -c 4 -a e"
     
     if ((model_set == 2)); then
	 if [[ $mat == 'HKY85' ]] && [[ -n ${model_cmds["${mat}+G"]} ]] && [[ ${seen["${mat}+G"]} -eq 0 ]]; then
	     HKY85G_BICi="$AICi"
	     echo "HKY85+G_BICi: $HKY85G_BICi"	
             seen['HKY85+G']=1
	 fi
     fi
     
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "$mat" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     extra_params=1 # 1 pInv
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I"]="$model_string"
     model_cmds["${mat}+I"]="$mat -v e"
          
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
     phyml -i "$infile" -d nt -m "$mat" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --leave_duplicates --no_memory_check &> /dev/null
     extra_params=2 #19 from AA frequencies + 1 gamma 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I+G"]="$model_string"
     model_cmds["${mat}+I+G"]="$mat -v e -c 4 -a e"

     if ((model_set == 2)); then
	 
	 if [[ $mat == 'HKY85' ]] && [[ -n ${model_cmds["${mat}+I+G"]} ]] && [[ ${seen["${mat}+G"]} -eq 1 ]]; then
	     HKY85IG_BICi="$AICi"
	     echo "HKY85+I+G_BICi: $HKY85IG_BICi"
	 fi
     fi
done

# cleanup: remove phyml output files from the last pass through the loop
[[ -s "${infile}"_phyml_stats.txt ]] && rm "${infile}"_phyml_stats.txt
[[ -s "${infile}"_phyml_tree.txt ]] && rm "${infile}"_phyml_tree.txt


##### 5.2 if extended set, then compute and set the compositional_heterogeneity flag 
#           and loop over corresponding set of extended models

# check_compositional_heterogeneity and set compositional_heterogeneity flag (1|0), accordingly
if ((model_set == 2)); then
   if [[ -n "$JC_BICi" ]] && [[ -n "$F81_BICi" ]]; then
       compositional_heterogeneity=$(check_compositional_heterogeneity "$JC_BICi" "$F81_BICi") 
       echo '--------------------------------------------------------------------------------'
       print_start_time
       echo '# Starting evaluation and automatic selection of extended model set'
       echo "# setting compositional_heterogeneity flag to: $compositional_heterogeneity"
   fi

   if [[ -n "$HKY85_BICi" ]] && [[ -n "$TN93_BICi" ]]; then
       transitional_heterogeneity=$(check_transitional_heterogeneity "$HKY85_BICi" "$TN93_BICi") 
       echo '--------------------------------------------------------------------------------'
       print_start_time
       echo '# ... evaluation and automatic selection of extended model set'
       echo "# setting transitional_heterogeneity flag to: $transitional_heterogeneity"
   fi

   if [[ -n "$HKY85G_BICi"  ]] && [[ -n "$HKY85IG_BICi" ]]; then
       use_pInv=$(check_pInv "$HKY85G_BICi" "$HKY85IG_BICi") 
       echo '--------------------------------------------------------------------------------'
       print_start_time
       echo '# ... evaluation and automatic selection of extended model set'
       echo "# setting use_pInv flag to: $use_pInv"
   fi
  
   # fill the models array with the proper ones, 
   #   based on the compositional_heterogeneity flag
   models=()
   if ((compositional_heterogeneity == 1)); then
   	echo '# will evaluate models with unequal frequencies'
	models=( "${extended_models_uf[@]}" )
        freq_cmd=' -f m '
   else
   	echo '# will evaluate models with equal frequencies'
   	models=( "${extended_models_ef[@]}" )	
   	freq_cmd=' -f 0.25,0.25,0.25,0.25 '
   fi

   echo '--------------------------------------------------------------------------------'
   print_start_time

   echo "2.2. running a loop to combine the extended models in model_set ${model_set}=>${model_options[$model_set]},
      with (or not) +G and or +I, and compute the model lnL, after optimizing branch lengths and rates"

   # 5.2 loop over the set of extended models, passing the proper freq_cmd, based on the compositional_heterogeneity flag
   for mat in "${models[@]}"; do
     # skip models with two transition rates if transitional_heterogeneity == 0
     ((transitional_heterogeneity == 0)) \
       && [[ "$mat" =~ (01[0-6][0-6]2[0-6]|01[0-6][0-6]3[0-6]|01[0-6][0-6]4[0-6]|00[0-1][0-6][1-6][1-6]|000[0-6][1-6][0-6]) ]] \
       && echo "skipping TN|TIM|SYM matrix $mat" && continue
     ((transitional_heterogeneity == 1)) \
       && [[ "$mat" =~ (01[0-6][0-6]1[0-6]|00[0-6][0-6]0[0-6] ) ]] && echo "skipping TVM matrix $mat" && continue
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -c 1 -v 0 -o lr --leave_duplicates --no_memory_check &> /dev/null 
     extra_params=0 
     ((compositional_heterogeneity == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")     
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat -c 1 -v 0"

     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     phyml -i "$infile" -d nt -m "${mat}" "$freq_cmd" -v 0 -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     extra_params=1 # 1 gamma 
     ((compositional_heterogeneity == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -v 0 -c 4 -a e"

     if ((use_pInv > 0)); then
         ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
         phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
         extra_params=1 # 1 pInv
         ((compositional_heterogeneity == 0)) && mat="${mat}ef"
         total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
         sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
         score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
         AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
         AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
         BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
         printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
	 ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         model_scores["${mat}+I"]="$model_string"
         model_cmds["${mat}+I"]="$mat -v e -c 1"

         ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
         phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --leave_duplicates --no_memory_check &> /dev/null
         extra_params=2 # 1 pInv + 1 gamma 
         ((compositional_heterogeneity == 0)) && mat="${mat}ef"
         total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
         sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
         score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
         AICi=$(compute_AICi "$score" "$no_branches" "$total_params")
         AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
         BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
         printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
	 ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         model_scores["${mat}+I+G"]="$model_string"
         model_cmds["${mat}+I+G"]="$mat -v e -c 4 -a e"
     else
         echo "skipping  ${mat}+I and ${mat}+I+G" && continue
     fi
   done
fi # extended set

echo '--------------------------------------------------------------------------------'

echo ''
##### 

# 6. print a sorted summary table of model fits from the model_scores hash
print_start_time

((compositional_heterogeneity == 1)) && echo "# writing ${infile}_sorted_model_set_${model_set}_fits.tsv, sorted by BIC; extmodels+F true"
((compositional_heterogeneity == 0)) && echo "# writing ${infile}_sorted_model_set_${model_set}_fits.tsv, sorted by BIC; extmodels+F false"
echo '--------------------------------------------------------------------------------'

# Add the +F decoration to model codes when compositional_heterogeneity == 1
#  Note: need to remove it later on, when calling phyml for the final phylogeny
#   -m ${model_cmds[${best_model%+F}]}
if ((compositional_heterogeneity == 1)) ; then
     for m in "${!model_scores[@]}"; do
          if [[ "$m" =~ (JC69|K80|F81|HKY85|TN93) ]]; then
	        echo -e "$m\t${model_scores[$m]}"
	  else
	        echo -e "$m+F\t${model_scores[$m]}"
	  fi	  
    done | sort -gk7 > "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
else
    for m in "${!model_scores[@]}"; do
        #m="${m%ef}"
        echo -e "$m\t${model_scores[$m]}"
    done | sort -gk7 > "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
fi


# 7. compute delta_BIC and BICw, based on "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
declare -a BIC_a
declare -a BIC_deltas_a
declare -a BICw_a
declare -a BICcumW_a
#BIC_a=( $(awk '{print $7}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv) )
mapfile -t BIC_a < <(awk '{print $7}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv)
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
   #echo "DEBUG: delta=$i 'BEGIN{printf '%.10f', exp(-1/2 * delta) }')"
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }')  
   #echo "DEBUG num:$BICw_numerator"
   BICw_sums=$(bc <<< "$BICw_sums"'+'"$BICw_numerator")
done
#echo BICw_sums:$BICw_sums

# 7.3 fill the BICw_a and BICcumW_a arrays
BICw_a=()
BICcumW_a=()
BICcumW=0
for i in "${BIC_deltas_a[@]}"; do
   BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }' 2> /dev/null)   
   BICw=$(echo "$BICw_numerator / $BICw_sums" | bc -l)
   BICw_a+=( $(printf "%.2f" "$BICw") )
   BICcumW=$(echo "$BICcumW + $BICw" | bc)
   BICcumW_a+=( $(printf "%.2f" "$BICcumW") )
done

# 7.4 paste the BIC_deltas_a & BICw_a values as a new column to "${infile}"_sorted_model_set_"${model_set}"_fits.tsv
paste "${infile}"_sorted_model_set_"${model_set}"_fits.tsv <(for i in "${BIC_deltas_a[@]}"; do echo "$i"; done) \
                                                           <(for i in "${BICw_a[@]}"; do echo "$i"; done) \
							   <(for i in "${BICcumW_a[@]}"; do echo "$i"; done) > t
							   
[[ -s t ]] && mv t "${infile}"_sorted_model_set_"${model_set}"_fits.tsv


# 7.5 Display  the final "${infile}"_sorted_model_set_"${model_set}"_fits.tsv and extract the best model name
if [[ -s "${infile}"_sorted_model_set_"${model_set}"_fits.tsv ]]; then
    # display models sorted by BIC
    best_model=$(awk 'NR == 1{print $1}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv)
    [[ -z "$best_model" ]] && echo "FATAL ERROR: unbound \$best_model at $LINENO" && exit 1

    # print table with header to STDOUT and save to file
    awk 'BEGIN{print "model\tK\tsites/K\tlnL\tAIC\tAICc\tBIC\tdeltaBIC\tBICw\tBICcumW"}{print}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv | column -t
    awk 'BEGIN{print "model\tK\tsites/K\tlnL\tAIC\tAICc\tBIC\tdeltaBIC\tBICw\tBICcumW"}{print}' "${infile}"_sorted_model_set_"${model_set}"_fits.tsv > t
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
echo "* NOTE 1: when sites/K < 40, the AICc is recommended over AIC"
echo "* NOTE 2: The best model is selected by BIC, because AIC is biased, favoring parameter-rich models"
echo ''
echo '=================================================================================================='
echo '#  Will estimate the ML tree under best-fitting model $best_model selected by BIC'
echo '--------------------------------------------------------------------------------------------------'

print_start_time

# Note: need to remove the +F or ef decoration from model codes, if present, to index the hash
((compositional_heterogeneity == 1)) && best_model=${best_model%+F}
((compositional_heterogeneity == 0)) && best_model=${best_model%ef}

if ((n_starts == 1)); then
    echo "# running: phyml -i $infile -d nt -m ${model_cmds[${best_model}]} -o tlr -s BEST"

    # note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the comand fail (Bash v4.4)
    phyml -i "$infile" -d nt -m ${model_cmds[${best_model}]} -o tlr -s BEST &> /dev/null
else
    echo "# running: phyml -i $infile -d nt -m ${model_cmds[${best_model}]} -o tlr -s SPR --rand_start --n_rand_starts $n_starts"

    # note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the comand fail (Bash v4.4)
    phyml -i "$infile" -d nt -m ${model_cmds[${best_model}]} -o tlr -s SPR --rand_start --n_rand_starts "$n_starts" &> /dev/null
fi


# 8.1 Check and rename final phyml output files
if [[ -s "${infile}"_phyml_stats.txt ]]; then
     
     if ((n_starts == 1)); then
         mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_BESTmoves_phyml_stats.txt
         echo "# Your results:"
         echo "  - ${infile}_${best_model}_BESTmoves_phyml_stats.txt"
     else
         mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_"${n_starts}"rdmStarts_SPRmoves_phyml_stats.txt
         echo "# Your results:"
         echo "  - ${infile}_${best_model}_${n_starts}rdmStarts_SPRmoves_phyml_stats.txt"
     fi
else
     echo "FATAL ERROR: ${infile}_phyml_stats.txt was not generated!"
fi

if [[ -s "${infile}"_phyml_tree.txt ]]; then
     if ((n_starts == 1)); then
         mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_BESTmoves_phyml_tree.txt
         echo "  - ${infile}_${best_model}_BESTmoves_phyml_tree.txt"
     else
         mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_"${n_starts}"rdmStarts_SPRmoves_phyml_tree.txt
         echo "  - ${infile}_${best_model}_${n_starts}rdmStarts_SPRmoves_phyml_tree.txt"
     fi
else
     echo "FATAL ERROR: ${infile}_phyml_tree.txt was not generated!"
fi

if ((n_starts > 1)) && [[ -s "${infile}"_phyml_rand_trees.txt ]]; then
    mv "${infile}"_phyml_rand_trees.txt "${infile}"_phyml_"${n_starts}"rand_trees.txt
    echo "  - ${infile}_phyml_${n_starts}rand_trees.txt"
fi 

echo '--------------------------------------------------------------------------------------------------'

echo ''

elapsed=$(( SECONDS - start_time ))

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days, %H hr, %M min, %S sec')"

echo 'Done!'

print_end_message

exit 0
