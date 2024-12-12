#!/usr/bin/env bash

#: phyml_DNAmodelFinder.sh
#: Author: Pablo Vinuesa, CCG-UNAM, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#: AIM: wraper script around phyml, to select a reasonable substitution model for DNA alignments.
#:      Computes AIC, BIC, delta_BIC and BICw, and estimate a ML phylogeny using the best-fitting model
#:         selected under the BIC.
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
version='2.1.0_2024-12-12_GUADALUPE' # v2.1.0_2024-12-12_GUADALUPE: fixed awk exp() out of range error message, by replacing with 2.7183^()
   # v2.0_2024-12-12_GUADALUPE; major upgrade
   # - added getopts interface with key options -m -s -b to control model set, search and branch-support methods.
   # - improved checking of user-provided input data and parameters.
   # v1.2.2_2024-12-3' strict check of phylip format
   # phyml_DNAmodelFinder.sh v1.2.0_2023-11-18; 
   # - implements extended (base) model set: 
   #	* 6 standard named models
   #	* 64 extended equal-frequency (ef) models (TIM and TVM sets)
   #	* 62 extended unequal-frequency (uf) models (TIM and TVM sets)
   # - minor code cleanup (removed unused variable)

min_bash_vers=4.4 # required to write modern bash idioms:
                  # 1.  printf '%(%F)T' '-1' in print_start_time; and 
                  # 2. passing an array or hash by name reference to a bash function (since version 4.3+), 
		  #    by setting the -n attribute
		  #    see https://stackoverflow.com/questions/16461656/how-to-pass-array-as-an-argument-to-a-function-in-bash

min_phyml_version=3.3 # this corresponds to 3.3.year; check_phyml_version also extracts the year, suggesting to update to v2022
phymlv=''
phymlyr=''
model_set=''
bash_ge_5=''

DEBUG=0

n_starts=1         # seed trees
delta_BIC_cutoff=2 # to set compositional_heterogeneity, transitional_heterogeneity and pInv flags 
search_method=BEST # BEST of NNI and SPR at each move
boot=-5


# declare array and hash variables
declare -a models             # array holding the base models (empirical substitution matrices to be evaluated)
declare -A model_cmds         # hash holding the commands for each model 
declare -A model_scores       # hash holding the model lnL scores and AICi values 
declare -A model_options      # hash mapping option => model_set
declare -A model_free_params  # hash mapping base models with their free parameters

declare -a TIMef TVMef TIMuf TVMuf # array holding the model codes

# To be implemented
# http://www.iqtree.org/doc/Substitution-Models#dna-models
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

### extended_models_ef - TVM* # 21

model_free_params['012210ef']=2  # K81
model_free_params['012314ef']=4  # TVMef
model_free_params['012310ef']=3  # TVM1ef 
model_free_params['010213ef']=3  # TVM2ef
model_free_params['012213ef']=3  # TVM3ef
model_free_params['012013ef']=3  # TVM4ef
model_free_params['010012ef']=2  # TVM5ef 
model_free_params['012012ef']=2  # TVM6ef   # TPM3
model_free_params['010212ef']=2  # TVM7ef   # TPM2
model_free_params['010210ef']=2  # TVM8ef 
model_free_params['012313ef']=3  # TVM9ef 
model_free_params['010011ef']=1  # TVM10ef
model_free_params['012212ef']=2  # TVM11ef
model_free_params['011010ef']=1  # TVM12ef
model_free_params['001101ef']=1  # TVM13ef
#==
model_free_params['000100ef']=1  # TVM14ef
model_free_params['000101ef']=1  # TVM15ef
model_free_params['010110ef']=1  # TVM16ef
model_free_params['000102ef']=2  # TVM17ef
model_free_params['010112ef']=2  # TVM18ef
model_free_params['010211ef']=2  # TVM19ef

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
#== 
TVMef[15]=000100
TVMef[16]=000101
TVMef[17]=010110
TVMef[18]=000102
TVMef[19]=010112
TVMef[20]=010211


### extended_models_ef - TNef|TIM*|SYM # 43
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
model_free_params['012222ef']=2  # TIM11ef
model_free_params['012223ef']=3  # TIM12ef
model_free_params['010120ef']=2  # TIM13ef
model_free_params['000120ef']=2  # TIM14ef
model_free_params['000121ef']=2  # TIM15ef
model_free_params['001021ef']=2  # TIM16ef
model_free_params['012234ef']=4  # TIM17ef
model_free_params['010231ef']=3  # TIM18ef
model_free_params['011230ef']=3  # TIM19ef
model_free_params['011020ef']=2  # TIM20ef
model_free_params['012130ef']=3  # TIM21ef
model_free_params['010121ef']=2  # TIM22ef
model_free_params['010122ef']=2  # TIM23ef
model_free_params['010123ef']=3  # TIM24ef
#==
model_free_params['001020ef']=2  # TIM25ef
model_free_params['000123ef']=3  # TIM26ef
model_free_params['010203ef']=3  # TIM27ef
model_free_params['010223ef']=3  # TIM28ef
model_free_params['010230ef']=3  # TIM29ef
model_free_params['012032ef']=3  # TIM30ef # < TIM3
model_free_params['010233ef']=3  # TIM31ef 
model_free_params['001234ef']=4  # TIM32ef
model_free_params['010234ef']=4  # TIM33ef
model_free_params['011234ef']=4  # TIM34ef
model_free_params['012234ef']=4  # TIM35ef
model_free_params['012134ef']=4  # TIM36ef
model_free_params['012304ef']=4  # TIM37ef
model_free_params['012324ef']=4  # TIM38ef
model_free_params['012334ef']=4  # TIM39ef
model_free_params['012341ef']=4  # TIM40ef
model_free_params['012344ef']=4  # TIM41ef

TIMef[0]=010020   # TNef
TIMef[1]=012230   # TIMef
TIMef[2]=012345
TIMef[3]=010023
TIMef[4]=010232   # TIM2
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
#==
TIMef[26]=001020
TIMef[27]=000123
TIMef[28]=010203
TIMef[29]=010223
TIMef[30]=010230
TIMef[31]=012032
TIMef[32]=010233
TIMef[33]=001234
TIMef[34]=010234
TIMef[35]=011234
TIMef[36]=012134
TIMef[37]=012222
TIMef[38]=012304
TIMef[39]=012324
TIMef[40]=012334
TIMef[41]=012341
TIMef[42]=012344


## extended_models_uf TVM* # 21
model_free_params['012210']=5  # K81uf # K81=TPM1
model_free_params['012314']=7  # TVM   # TVM
model_free_params['012310']=6  # TVM1
model_free_params['010213']=6  # TVM2
model_free_params['012213']=6  # TVM3
model_free_params['012013']=6  # TVM4
model_free_params['010012']=5  # TVM5 
model_free_params['012012']=5  # TVM6   # TPM3
model_free_params['010212']=5  # TVM7   # TPM2
model_free_params['012313']=6  # TVM8 
model_free_params['010011']=4  # TVM9
model_free_params['012212']=5  # TVM10 
model_free_params['010210']=5  # TVM11 
model_free_params['011010']=4  # TVM12
model_free_params['001101']=4  # TVM13 
#==
model_free_params['000100']=4  # TVM14
model_free_params['000101']=4  # TVM15
model_free_params['010110']=4  # TVM16
model_free_params['000102']=5  # TVM17
model_free_params['010112']=5  # TVM18
model_free_params['010211']=5  # TVM19

TVMuf[0]=012210
TVMuf[1]=012314   # TVM
TVMuf[2]=012310
TVMuf[3]=010213
TVMuf[4]=012213
TVMuf[5]=012013
TVMuf[6]=010012
TVMuf[7]=012012   # TPM3
TVMuf[8]=010212   # TPM2
TVMuf[9]=012313
TVMuf[10]=010011
TVMuf[11]=012212
TVMuf[12]=010210
TVMuf[13]=011010
TVMuf[14]=001101
#== 
TVMuf[15]=000100
TVMuf[16]=000101
TVMuf[17]=010110
TVMuf[18]=000102
TVMuf[19]=010112
TVMuf[20]=010211


## extended_models_uf TIM* # 41
model_free_params['012230']=6  # TIM   # < TIM1
model_free_params['010023']=6  # TIM1
model_free_params['010232']=7  # TIM2  # < TIM2
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
model_free_params['011230']=6  #
model_free_params['012130']=6  # TIM20
model_free_params['010121']=5  # TIM21
model_free_params['010122']=5  # TIM22
model_free_params['010123']=6  # TIM23
#==
model_free_params['001020']=5  # TIM24
model_free_params['000123']=6  # TIM25
model_free_params['010203']=6  # TIM26
model_free_params['010223']=6  # TIM27
model_free_params['010230']=6  # TIM28
model_free_params['012032']=6  # TIM29 # < TIM3
model_free_params['010233']=6  # TIM30
model_free_params['001234']=7  # TIM31
model_free_params['010234']=7  # TIM32
model_free_params['011234']=7  # TIM33
model_free_params['012234']=7  # TIM34
model_free_params['012134']=7  # TIM35
model_free_params['012304']=7  # TIM36
model_free_params['012324']=7  # TIM37
model_free_params['012334']=7  # TIM38
model_free_params['012341']=7  # TIM39
model_free_params['012344']=7  # TIM40

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
TIMuf[20]=011230
TIMuf[21]=012130       
TIMuf[22]=010121       
TIMuf[23]=010122       
TIMuf[24]=010123       
#==		      
TIMuf[25]=001020
TIMuf[26]=000123
TIMuf[27]=010203
TIMuf[28]=010223
TIMuf[29]=010230
TIMuf[30]=012032
TIMuf[31]=010233
TIMuf[32]=001234
TIMuf[33]=010234
TIMuf[34]=011234
TIMuf[35]=012134
TIMuf[36]=012304
TIMuf[37]=012324
TIMuf[38]=012334
TIMuf[39]=012341
TIMuf[40]=012344

# array of models to evaluate
standard_models=(JC69 K80 F81 HKY85 TN93 GTR)      # 6

extended_models_ef=( "${TIMef[@]}" "${TVMef[@]}" ) # 64

extended_models_uf=( "${TIMuf[@]}" "${TVMuf[@]}" ) # 62
                                                       
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

function print_version()
{
   echo "$progname v$version"
   exit 0
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

function check_phyml_version {
   local phyml_version min_phyml_version phyml_version_year
   
   min_phyml_version=$1
   phyml_version_year=''
   phyml_version=''
   
   # This version extracts 3.3 from '3.3.20220408.'; but the series goes back to 2017 ...
   #   v3.3.20220408 is required to have the --leave_duplicates option
   phyml_version=$(phyml --version | awk '/This is PhyML version/{print substr($NF, 0, 4)}' | sed 's/\.$//') 
   #phyml_OK=$(awk -v v="$phyml_version" -v m="$min_phyml_version" 'BEGIN{if(v < m || v > 2000) { print 0}else{print 1}  }')
   
   if [[ 1 -eq "$(echo "$phyml_version == $min_phyml_version" | bc)" ]]
   then
        phyml_version_year=$(phyml --version | awk '/This is PhyML version/{print substr($NF, 0, 8)}' | sed 's/.*\.//g')
   else
        phyml_version_year="$phyml_version"
   fi
      
   printf '%s\n' "${phyml_version}_${phyml_version_year}"
}
#-----------------------------------------------------------------------------------------

function check_is_phylip(){
    local phylip_file="$1"

    if [[ ! -f "$phylip_file" ]]; then
        echo "Error: File not found: $phylip_file"
        return 1
    fi

    awk '
    BEGIN {
        valid = 1
    }
    NR == 1 {
        # First line: check for two integers separated by whitespace
        if (!($1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && NF == 2)) {
            print "Error: First line must contain two integers separated by whitespace."
            valid = 0
            exit 1
        }
        num_sequences = $1
        sequence_length = $2
    }
    NR > 1 && /^[^[:space:]]+/{ 
        # Other lines: verify alignment and proper separation
	id = substr($0, 1, 10)  # PHYLIP identifiers are typically up to 10 chars
	idcounts++
        seq = $2
        #gsub(/^ +| +$/, "", seq)  # Trim spaces around the sequence part
        if (!match(seq, /^[-A-Za-z.]+$/)) {
            print "Error: Invalid sequence format on line " NR "."
            valid = 0
            exit 1
        }
        # Check alignment (ensure space between ID and sequence)
        #if (substr($0, 10, 1) != " ") {
	if (! /^[^[:space:]]+[[:space:]+]/){
            print "Error: No space separating ID and sequence on line " NR "."
            valid = 0
            exit 1
        }
    }
    NR > 1 && /^[[:space:]+][-A-Za-z.]+/{
        seq = $0
        gsub(/^ +| +$| /, "", seq)  # Trim spaces around the sequence part
        if (!match(seq, /^[-A-Za-z.]+$/)) {
            print "Error: Invalid sequence format on line " NR "."
            valid = 0
            exit 1
        }
    }
    END {
        if (idcounts != num_sequences) {
            print "Error: Number of sequences " idcounts " does not match specified count (" num_sequences ")."
            valid = 0
        }
        exit valid ? 0 : 1
    }
    ' "$phylip_file"
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
   # AICi=$(compute_AICi "$score" "$total_params")
   local score total_params
   
   score=$1
   total_params=$2
    
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

function print_help {
   # Prints the help message explaining how to use the script.

   bash_vers=$(check_bash_version "$min_bash_vers")
   bash_ge_5=$(awk -v bv="$bash_vers" 'BEGIN { if (bv >= 5.0){print 1}else{print 0} }')   
   
if (( bash_ge_5 > 0)); then 

 cat <<EOF
 USAGE for $progname v$version:
 $progname -i <infile> -m <model set to evaluate> [-b <(-)int>] [ -h <print help>] [-r <number of random seed trees>] [-s <SEARCH METHOD>] [ -v <print version>]
 
 REQUIRED:
  -i <string> input alignment in PHYLIP format
  -m <int> model set to evaluate by BIC
       1 -> standard models (JC69 K80 F81 HKY85 TN93 GTR)
       2 -> WILL NOT RUN properly on Bash < v5.0, sorry (see NOTE below) 
       3 -> minimal test set (JC69 F81 HKY85 TN93)

 OPTIONAL
  -h <flag> print help (this message)
  -b <int>
      int > 0: int is the number of bootstrap replicates.
       0: neither approximate likelihood ratio test nor bootstrap values are computed.
      -1: approximate likelihood ratio test returning aLRT statistics.
      -2: approximate likelihood ratio test returning Chi2-based parametric branch supports.
      -4: SH-like branch supports alone.
      -5: (default) approximate Bayes branch supports.  
  -s <NNI|SPR|BEST> search (branch-swapping) method; default:$search_method
  -r <integer> number of random start trees to use for sequential searches; default rand_starts:$n_starts
  -v <flag> print version
 
 EXAMPLE:
   $progname -i my_PHYLIP_alignment.phy -m 2 -b -4 -r 10 -s NNI
 
 NOTES: 1. Assumes DNA input sequences are ALIGNED and in (relaxed) PHYLIP format 
	 
 AIM:  $progname v${version} will evaluate the fit of the selected model set,
	combined or not with +G and/or +f, computing AICi, BICi, deltaBIC and BICw, 
     	  inferring the ML tree under the BIC-selected model  

 PROCEDURE
  - Models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates 
       to calculate their AICi, BICi, delta_BIC and BICw
  - The best model is selected by BIC
  - SPR searches can be launched starting from multiple random trees
  - Default single seed tree searches use a BioNJ with BEST moves     
      
 SOURCE: the latest version of the program is available on GitHub:
	 https://github.com/vinuesa/TIB-filoinfo

 LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 

EOF

else

 cat <<EOF
 USAGE for $progname v$version: 
 $progname -i <infile> -m <model set to evaluate> [-b <(-)int>] [ -h <print help>] [-r <number of random seed trees>] [-s <SEARCH METHOD>] [ -v <print version>]
 
 REQUIRED:
  -i <string> input alignment in PHYLIP format
  -m <int> model set to evaluate by BIC
       1 -> standard models (JC69 K80 F81 HKY85 TN93 GTR)
       2 -> WILL NOT RUN properly on Bash < v5.0, sorry (see NOTE below) 
       3 -> minimal test set (JC69 F81 HKY85 TN93)

 OPTIONAL
  -h <flag> print help (this message)
  -b <int>
      int > 0: int is the number of bootstrap replicates.
       0: neither approximate likelihood ratio test nor bootstrap values are computed.
      -1: approximate likelihood ratio test returning aLRT statistics.
      -2: approximate likelihood ratio test returning Chi2-based parametric branch supports.
      -4: SH-like branch supports alone.
      -5: (default) approximate Bayes branch supports.  
  -s <NNI|SPR|BEST> search (branch-swapping) method; default:$search_method
  -r <integer> number of random start trees to use for sequential searches; default rand_starts:$n_starts
  -v <flag> print version
 
 EXAMPLE:
   $progname -i my_PHYLIP_alignment.phy -m 1 -b -4 -s NNI -r 10

 NOTES: 1. you are running the old Bash version $bash_vers. 
           Update to version >=5.0 to profit from the full set of models 
             and features implemented in $progname
        2. Assumes DNA input sequences are ALIGNED and in (relaxed) PHYLIP format 


 AIM:  $progname v${version} will evaluate the fit of the selected model set,
	combined or not with +G and/or +f, computing AICi, BICi, deltaBIC and BICw, 
     	  inferring the ML tree under the BIC-selected model  

 PROCEDURE
  - Models are fitted using a fixed NJ-JC tree, optimizing branch lenghts and rates 
       to calculate their AICi, BICi, delta_BIC and BICw
  - The best model is selected by BIC
  - SPR searches can be launched starting from multiple random trees
  - Default single seed tree searches use a BioNJ with BEST moves     
      
 SOURCE: the latest version of the program is available on GitHub:
	 https://github.com/vinuesa/TIB-filoinfo

 LICENSE: GPL v3.0. See https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE 

EOF

fi
   
   exit 1  
}

#-----------------------------------------------------------------------------------------
#============================= END FUNCTION DEFINITIONS ==================================
#=========================================================================================
#------------------------------------#
#----------- GET OPTIONS ------------#
#------------------------------------#
# This section uses getopts to process command-line arguments 
#    and set the corresponding variables accordingly. 

[ $# -eq 0 ] && print_help

args=("$@")

while getopts ':b:i:m:r:s:hv?:' OPTIONS
do
   case $OPTIONS in

   b)   boot=$OPTARG
        ;;
   i)   infile=$OPTARG
        ;;
   m)   model_set=$OPTARG
        ;;
   r)   n_starts=$OPTARG
        ;;
   s)   search_method=$OPTARG
        ;;	
   h)   print_help
        ;;
   v)   print_version
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
if [[ -z "$infile" ]]; then
       echo "Missing required input PHYLIP file name: -i <input>"
       print_help
fi

if [[ -z "$model_set" ]]; then
       echo "Missing required model_set option: -m <1|2|3>" warn
       print_help
fi

# Validate user-provided data and params
[[ ! -s "$infile" ]] && echo "FATAL ERROR: could not find $infile in $wkd" && exit 1

# Check input file is a relaxed or canonical PHYLIP file
if ! check_is_phylip "$infile" &> /dev/null; then 
   echo "FATAL ERROR: input file ${infile} does not seem to be a canonical phylip file. Will exit now!"
   exit 1
fi 

if (( model_set < 1 )) || (( model_set > 3 )); then
   print_help
fi

if (( model_set == 2 )) && (( bash_ge_5 < 0 )); then
   echo "You are using bash version $bash_vers, which is incopatible with model_set == 2; upgrade your Bash, or use model_set == 1"
   print_help
fi


if (( boot < -5 )); then
   echo "# you are setting branch support value mode to $boot; values < -5 are not allowed!"
   print_help
fi

if (( boot > 10000 )); then
   echo "# you are setting bootstrapping to $boot; that will take a very long time to compute. Use a number of pseudoreplicates <= 10000"
   print_help
fi


# Check that the provided search_method matches the allowed set
regex='NNI|SPR|BEST'

if ! [[ "$search_method" =~ $regex ]]; then
   echo "ERROR: the provided search_method '-S $search_method' does not match the allowed set $regex"
   exit 1
fi

# ============ #
# >>> MAIN <<<
# ------------ #
# Main script logic starts here

## Check environment
# 0. Check that the input file was provided, and that the host runs bash >= v4.3
echo "# invocation: $progname " "${args[@]}"

wkd=$(pwd)


## Check Bash & PhyMLs vesions
# check PhyML version
echo "# checking PhyML version"
phymlv=$(check_phyml_version "$min_phyml_version")
phymlyr=$(echo "$phymlv" | cut -d_ -f2)
   
   
# OK, ready to start the analysis ...
start_time=$SECONDS
echo "========================================================================================="
bash_vers=$(check_bash_version "$min_bash_vers")
awk -v bv="$bash_vers" -v mb="$min_bash_vers" \
  'BEGIN { if (bv < mb){print "FATAL: you are running acient bash v"bv, "and version >=", mb, "is required"; exit 1}else{print "# Bash version", bv, "OK"} }'

echo "# running with phyml v.${phymlv}"
((phymlyr < 2022)) && printf '%s\n%s\n' "# Warning: running old PhyML version from $phymlyr!" "   Update to the latest one, using the phyml's GitHub repo: https://github.com/stephaneguindon/phyml/releases" 

check_dependencies

echo -n "# $progname v$version running on $host. Run started on: "; printf '%(%F at %T)T\n' '-1'
echo "# working directory: $wkd"

# Print run parameters
echo "# infile:$infile; model_set:$model_set => ${model_options[$model_set]}; seed trees: $n_starts; delta_BIC_cutoff=$delta_BIC_cutoff; branch_support_type=$boot"
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
(($? > 0)) && { echo "FATAL ERROR: input file ${infile} does not seem to be a canonical phylip file. Will exit now!"; exit 1 ; }
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
    echo "FATAL ERROR: could not compute ${infile}_phyml_tree.txt (${infile}_JC-NJ.nw)" && exit 1
fi

# 5.1 run a for loop to combine all base models with (or not) +G, +I +G+I
#     and fill the model_scores and model_cmds hashes
echo "2.1. running in a for loop to combine all base models in model_set ${model_set}=>${model_options[$model_set]},
     with (or without) +G and or +I, and compute the model lnL, after optimizing branch lengths and rates"

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
    
# loop to test standard models, with or without +G, +I, +G+I; 
# Note the use of -f m to estimate base frequencies under ML. PhyML takes care of using or not -f m, given the standard model name
for mat in "${models[@]}"; do
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -f m -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -u "${infile}"_JC-NJ.nwk -c 1 -o lr --leave_duplicates --no_memory_check &> /dev/null 
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -u "${infile}"_JC-NJ.nwk -c 1 -o lr --no_memory_check &> /dev/null 
     extra_params=0 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat"

     # save the JC_BICi to test for significant compositional heterogeneity (JC vs F81)
     #  and transitional bias between purines and pyrimidines (HKY vs TN93)
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
     
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -f m -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "${mat}" -f m -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "${mat}" -f m -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --no_memory_check &> /dev/null
     extra_params=1 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")
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
     
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -f m -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --no_memory_check &> /dev/null
     extra_params=1 # 1 pInv
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     model_scores["${mat}+I"]="$model_string"
     model_cmds["${mat}+I"]="$mat -v e"
          
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat -f m -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --leave_duplicates --no_memory_check &> /dev/null
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" -f m -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --no_memory_check &> /dev/null
     extra_params=2 #19 from AA frequencies + 1 gamma 
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")
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
       && [[ "$mat" =~ (00[0-5][0-5]2[0-5]|01[0-5][0-5]2[0-6]|01[0-5][0-5]0[0-6]|01[0-5][0-5]3[0-5]|01[0-5][0-5]4[0-5]|01[0-5][0-5]5[0-5]|00[0-1][0-5][1-5][1-5]|000[0-5][1-5][0-5]) ]] \
       && echo "skipping TN|TIM|SYM matrix $mat" && continue
     ((transitional_heterogeneity == 1)) \
       && [[ "$mat" =~ (01[0-5][0-5]1[0-5]|00[0-5][0-5]0[0-5]|000[0-5]0[0-5] ) ]] && echo "skipping TVM matrix $mat" && continue
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -c 1 -v 0 -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -c 1 -v 0 -o lr --leave_duplicates --no_memory_check &> /dev/null 
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -c 1 -v 0 -o lr --no_memory_check &> /dev/null 
     extra_params=0 
     ((compositional_heterogeneity == 0)) && mat="${mat}ef"
     ((DEBUG)) && echo "### DEBUG:  mat=$mat"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")     
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     model_scores["${mat}"]="$model_string"
     model_cmds["${mat}"]="$mat -c 1 -v 0"

     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -c 4 -a e -u ${infile}_JC-NJ.nwk -o lr"
     ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "${mat}" "$freq_cmd" -v 0 -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
     ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "${mat}" "$freq_cmd" -v 0 -c 4 -a e -u "${infile}"_JC-NJ.nwk -o lr --no_memory_check &> /dev/null
     extra_params=1 # 1 gamma 
     ((compositional_heterogeneity == 0)) && mat="${mat}ef"
     total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
     sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
     score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
     AICi=$(compute_AICi "$score" "$total_params")
     AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
     BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
     printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
     ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
     model_scores["${mat}+G"]="$model_string"
     model_cmds["${mat}+G"]="$mat -v 0 -c 4 -a e"

     if ((use_pInv > 0)); then
         ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -v e -c 1 -u ${infile}_JC-NJ.nwk -o lr"
         ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --leave_duplicates --no_memory_check &> /dev/null
         ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -v e -c 1 -u "${infile}"_JC-NJ.nwk -o lr --no_memory_check &> /dev/null
         extra_params=1 # 1 pInv
         ((compositional_heterogeneity == 0)) && mat="${mat}ef"
         total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
         sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
         score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
         AICi=$(compute_AICi "$score" "$total_params")
         AICc=$(compute_AICc "$score" "$total_params" "$no_sites" "$AICi")
         BICi=$(compute_BIC "$score" "$total_params" "$no_sites")
         printf -v model_string "%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f" "$total_params" "$sites_by_K" "$score" "$AICi" "$AICc" "$BICi"
	 ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         model_scores["${mat}+I"]="$model_string"
         model_cmds["${mat}+I"]="$mat -v e -c 1"

         ((compositional_heterogeneity == 0)) && mat="${mat%ef}"
         print_start_time && echo "# running: phyml -i $infile -d nt -m $mat $freq_cmd -u ${infile}_JC-NJ.nwk -v e -a e -o lr"
         ((phymlyr >= 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --leave_duplicates --no_memory_check &> /dev/null
	 ((phymlyr < 2022)) && phyml -i "$infile" -d nt -m "$mat" "$freq_cmd" -u "${infile}"_JC-NJ.nwk -v e -a e -c 4 -o lr --no_memory_check &> /dev/null
         extra_params=2 # 1 pInv + 1 gamma 
         ((compositional_heterogeneity == 0)) && mat="${mat}ef"
         total_params=$((no_branches + extra_params + ${model_free_params[$mat]}))
         sites_by_K=$(echo 'scale=2;'"$no_sites/$total_params" | bc -l)
         score=$(awk '/Log-l/{print $NF}' "${infile}"_phyml_stats.txt)
         AICi=$(compute_AICi "$score" "$total_params")
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
# e (Euler's number; constant | e = 2.718281828459) base of the natural logarithm and exponential function
#  use e^(-1/2 * delta) instead of exp(-1/2 * delta) to avoid exponential out of range messages like:
#  awk: cmd. line:1: warning: exp: argument -1233.51 is out of range # on chaac
# https://stackoverflow.com/questions/69285812/awk-error-argument-is-out-of-the-range-how-to-solve-this-problem
# https://www.rapidtables.com/math/number/e_constant.html
e=2.718281828459
BICw_sums=0
for i in "${BIC_deltas_a[@]}"; do 
   #echo "DEBUG: delta=$i 'BEGIN{printf '%.10f', exp(-1/2 * delta) }')"
   #BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }' 2> /dev/null)  
   BICw_numerator=$(awk -v delta="$i" -v e="$e" 'BEGIN{printf "%.10f", e^(-1/2 * delta) }' 2> /dev/null)
   #echo "DEBUG num:$BICw_numerator"
   BICw_sums=$(bc <<< "$BICw_sums"'+'"$BICw_numerator")
done
#echo BICw_sums:$BICw_sums

# 7.3 fill the BICw_a and BICcumW_a arrays
BICw_a=()
BICcumW_a=()
BICcumW=0
for i in "${BIC_deltas_a[@]}"; do
   #BICw_numerator=$(awk -v delta="$i" 'BEGIN{printf "%.10f", exp(-1/2 * delta) }' 2> /dev/null)   
   BICw_numerator=$(awk -v delta="$i" -e="$e" 'BEGIN{printf "%.10f", e^(-1/2 * delta) }' 2> /dev/null) 
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
echo "#  Will estimate the ML tree under best-fitting model $best_model selected by BIC"
echo '--------------------------------------------------------------------------------------------------'

print_start_time

# Note: need to remove the +F or ef decoration from model codes, if present, to index the hash
((compositional_heterogeneity == 1)) && best_model=${best_model%+F}
((compositional_heterogeneity == 0)) && best_model=${best_model%ef}

if ((n_starts == 1)); then
    echo "# running: phyml -i $infile -d nt -m ${model_cmds[${best_model}]} -o tlr -s $search_method -b $boot"

    # note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the command fail (Bash v4.4)
    phyml -i "$infile" -d nt -m ${model_cmds[${best_model}]} -o tlr -s $search_method -b "$boot" &> /dev/null
else
    echo "# running: phyml -i $infile -d nt -m ${model_cmds[${best_model}]} -o tlr -s $search_method --rand_start --n_rand_starts $n_starts -b $boot"

    # note that on tepeu, the quotes around "${model_cmds[$best_model]}" make the comand fail (Bash v4.4)
    phyml -i "$infile" -d nt -m ${model_cmds[${best_model}]} -o tlr -s "$search_method" --rand_start --n_rand_starts "$n_starts" -b "$boot" &> /dev/null
fi


# 8.1 Check and rename final phyml output files
if [[ -s "${infile}"_phyml_stats.txt ]]; then
     
     if ((n_starts == 1)); then
         mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_"${search_method}"moves_phyml_stats.txt
         echo "# Your results:"
         echo "  - ${infile}_${best_model}_${search_method}moves_phyml_stats.txt"
     else
         mv "${infile}"_phyml_stats.txt "${infile}"_"${best_model}"_"${n_starts}"rdmStarts_"${search_method}"moves_phyml_stats.txt
         echo "# Your results:"
         echo "  - ${infile}_${best_model}_${n_starts}rdmStarts_${search_method}moves_phyml_stats.txt"
     fi
else
     echo "FATAL ERROR: ${infile}_phyml_stats.txt was not generated!"
fi

if [[ -s "${infile}"_phyml_tree.txt ]]; then
     if ((n_starts == 1)); then
         mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_"${search_method}"moves_phyml_tree.txt
         echo "  - ${infile}_${best_model}_${search_method}moves_phyml_tree.txt"
     else
         mv "${infile}"_phyml_tree.txt "${infile}"_"${best_model}"_"${n_starts}"rdmStarts_"${search_method}"moves_phyml_tree.txt
         echo "  - ${infile}_${best_model}_${n_starts}rdmStarts_${search_method}moves_phyml_tree.txt"
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
