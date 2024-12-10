#!/usr/bin/awk -f

# fasta_toolkit.awk VERSION:0.1 released Dec 21, 2020
# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
# Source: https://github.com/vinuesa/intro2linux
# AIM: munge fasta file containing one or more CDSs using one of 5 runmodes:
#       -R 1 [filter sequences matching -m string]
#       -R 2 [reverse complement DNA sequence] 
#       -R 3 [extract sequence by -s start -e end coordinates]
#       -R 4 [translate CDSs, using universal genetic code]
#       -R 5 [print basic DNA sequence stats]
#       -R 6 [print sequence lengths (also for proteins!]
# USAGE: call the script without arguments to print the help menu
# NOTES:  
#   1. uses Arnold Robbin's getopt() function from the gawk distro to deal with options and arguments
#   2. can read FASTA files from file or STDIN
#   3. pass only single FASTA files to fasta_toolkit.awk
#   4. prints results to STDOUT

#---------------------------------------------------------------------------------------------------------#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
#---------------------------------------------------------------------------------------------------------#

function read_fasta(      i, h, s)
{ 
   # fill hash seqs with the DNA/CDS sequences
   s=""
   h=$1
   for (i=2; i<=NF; i++) s=s$i
   seqs[h]=s
}
#---------------------------------------------------------------------------------------------------------

function rev_compl(header, dna_seq,        i, k, x) # reverse complement
{  # receives two arguments: the fasta header and the CDS sequence to be translated
    k=x=""
    
    dna_seq = toupper(dna_seq) # to match the keys of compl_nucl
    for( i = length(dna_seq); i !=0; i-- ) { # note that the loop reads the sequence from end to beginnig
         k = substr(dna_seq, i, 1)
         x = x compl_nucl[k]
    }
    printf ">%s\n%s\n", header, x
}
#---------------------------------------------------------------------------------------------------------

function extract_sequence_by_coordinates(header, dna_seq, start, end) {
       if(start == 1) print ">"header"\n"substr(dna_seq, start, end)
       if(start > 1)  print ">"header"\n"substr(dna_seq, start, end-start+1)
}
#---------------------------------------------------------------------------------------------------------

function translate_dna(header, dna_seq,      s, i, p)
{  # receives two arguments: the fasta header and the CDS sequence to be translated
   
   # Initialize variables: 
   #  do-while loop control variable i (nt counter) 
   #   and p, which will hold the translation product
   {i=1; p=""; triplet_counter=0}

   # Here we run a do-while loop; the do loop is a variation of the while looping statement. 
   #  The do loop always executes the body once and then repeats the body as long as the condition is true
   # We use the do-while loop, to get a first triplet string saved in s; 
   #  then the while loop keeps going until substr() got the last triplet, resulting in an empty s="".
   do {
   	  # First check that the script got some input
   	  #   if not, exit with an error message
   	  if(length(dna_seq) == 0) {
   	      print "ERROR: need a DNA sequence string to translate (valid DNA sequence, divisible by 3) "
   	      exit 1
       
   	  # Check that the DNA sequence string is divisible by 3 using the modulo operator
   	  #   if not, exit with an error message
   	  } else if(length(dna_seq)%3) { 
   	      printf "# WARNING: input DNA sequence for %s not divisible by 3. Will skip it!\n", header
	      break
   	  }

   	  # use substr() to split the input sequence (dna_seq) into triplets saved in s 	
   	  s=substr(dna_seq, i, 3)
       
   	  # keep track of processed triplets (codons)
   	  triplet_counter++
       
   	  # check that the input corresponds to valid nucleotides
   	  if ( s ~ /[^acgtACGT]+/ ) { 
   	      print "ERROR: input triplet", triplet_counter, "=", s, 
   		    "contains a non-valid nucleotide symbol ..."
   	      exit 3
   	  }

   	  # make sure that input nt symbols are uppercase to match the hash keys
   	  s=toupper(s)
   		
   	  # use the codon hash c as lookup table to translate the s triplet
   	  #   appending codons[s] to the growing peptide p
   	  { 
   	      # break out of loop if we get no more triplets 
   	      #   out of the input DNA string with substr()
   	      if (codons[s]=="") { 
   		 break
   	      }
   	      else if (s in codons == 0) { 
   		 # if the triplet is not contained in c, append "X" to p
   		 p=p unknown
   	      } else { 
   		 # append aminoacid codons[s] to growing peptide
   		 p = p codons[s]
   	     }
   	 }
   	 i=i+3 # increment the counter of processed dna nucleotides by 3 
    }
    # run while loop until substring cannot retrieve any more triplets
    while (s!="")
    prots[header]=p
}
#---------------------------------------------------------------------------------------------------------
function print_sequence_stats(header, dna_seq,         i,l) 
{  # receives two arguments: the fasta header and the CDS sequence to be translated
   sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=""
   l=length(dna_seq)
   
   for (i=1;i<=l;i++) {
             if (substr(dna_seq,i,1)=="T") sumT+=1
        else if (substr(dna_seq,i,1)=="A") sumA+=1
        else if (substr(dna_seq,i,1)=="G") sumG+=1
        else if (substr(dna_seq,i,1)=="C") sumC+=1
        else if (substr(dna_seq,i,1)=="N") sumN+=1
   }
   # print stats
   printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%3.2f\n", header, sumA, sumC, sumG, sumT, sumN, l, (sumC+sumG)/l*100
}
#---------------------------------------------------------------------------------------------------------

function print_sequence_lengths(header, dna_seq,         l) 
{  # receives two arguments: the fasta header and the CDS sequence to be translated
   l=length(dna_seq)
   
   #length_summary_stats(l)   
   # print sequence length to STDOUT
   printf "%s\t%d\n", header, l 
}
#---------------------------------------------------------------------------------------------------------
#function length_summary_stats (l,   sum,mean_len,min,max,counter)
#{
#   sum += l
#   counter++
#   
#   # print sequence length to STDOUT
#   mean_len =sum/counter
#         
#         if (l < min){min = l} else if (l > max) {max = l}
#
#           printf "%d\t%d\t%3.2f\n", min, max, mean_len 
#   
#}
#---------------------------------------------------------------------------------------------------------

function print_help(prog, vers) # (program, version)
{   
   print "OPTIONS for " prog " v"vers > "/dev/stderr" 
   print " # Required options" > "/dev/stderr" 
   print "    -R <int> [runmode]" > "/dev/stderr" 
   print "         1 [filter DNA|PROT sequences matching -m string]" > "/dev/stderr" 
   print "         2 [reverse complement DNA sequence]" > "/dev/stderr" 
   print "         3 [extract DNA|PROT sequence by -s start -e end coordinates]" > "/dev/stderr" 
   print "         4 [translate CDSs]" > "/dev/stderr" 
   print "         5 [print basic DNA sequence stats]" > "/dev/stderr" 
   print "         6 [print DNA|PROT sequence lenghts]" > "/dev/stderr" 

   print "\n # Runmode-specific options" > "/dev/stderr" 
   print "     -m <string> [match 'string'] for -R 1" > "/dev/stderr" 
   print "     -M <string> [reverse match 'string'] for -R 1" > "/dev/stderr" 
   print "     -s <int> [start coord] for -R 3" > "/dev/stderr" 
   print "     -e <int> [end coord] for -R 3" > "/dev/stderr" 

   print "\n # Other options" > "/dev/stderr" 
   print "     -d [FLAG; sets DEBUG=1 to print extra debugging info]" > "/dev/stderr" 
  
   print "\n # Usage examples:"
   print "    ./" prog " -R 1 -m match_string input_DNA_or_PROT.fasta" > "/dev/stderr" 
   print "    ./" prog " -R 1 -M RevMatch_string input_DNA.fasta" > "/dev/stderr" 
   print "    ./" prog " -R 2 input_DNA.fasta" > "/dev/stderr" 
   print "    ./" prog " -R 3 -s 2 -e 5 input_DNA_or_PROT.fasta" > "/dev/stderr" 
   print "    ./" prog " -R 4 input_DNA.fasta" > "/dev/stderr" 
   print "    ./" prog " -R 5 input_DNA.fasta" > "/dev/stderr" 
   print "    cat input.fasta | ./"prog " -R 5" > "/dev/stderr" 
   print "    ./" prog " -R 6 input_DNA_or_PROT.fasta | sort -t$'\\t' -nk2,2" > "/dev/stderr" 
   print "    ./" prog " -R 6 input_DNA_or_PROT.fasta | cut -f2 | col_sumStats.sh" > "/dev/stderr" 

   print "\n # Notes:" > "/dev/stderr" 
   print "    1. Pass only single FASTA files to " prog > "/dev/stderr" 
   print "    2. prints results to STDOUT" > "/dev/stderr"

   exit 1
}
#---------------------------------------------------------------------------------------------------------

# available in /usr/share/awk/
# getopt.awk --- Do C library getopt(3) function in awk
#
# Arnold Robbins, arnold@skeeve.com, Public Domain
#
# Initial version: March, 1991
# Revised: May, 1993

# External variables:
#    Optind -- index in ARGV of first nonoption argument
#    Optarg -- string value of argument to current option
#    Opterr -- if nonzero, print our own diagnostic
#    Optopt -- current option letter

# Returns:
#    -1     at end of options
#    "?"    for unrecognized option
#    <c>    a character representing the current option

# Private Data:
#    _opti  -- index in multiflag option, e.g., -abc
function getopt(argc, argv, options,    thisopt, i)
{
    if (length(options) == 0)    # no options given
        return -1

    if (argv[Optind] == "--") {  # all done
        Optind++
        _opti = 0
        return -1
    } else if (argv[Optind] !~ /^-[^:[:space:]]/) {
        _opti = 0
        return -1
    }
    if (_opti == 0)
        _opti = 2
    thisopt = substr(argv[Optind], _opti, 1)
    Optopt = thisopt
    i = index(options, thisopt)
    if (i == 0) {
        if (Opterr)
            printf("%c -- invalid option\n", thisopt) > "/dev/stderr"
        if (_opti >= length(argv[Optind])) {
            Optind++
            _opti = 0
        } else
            _opti++
        return "?"
    }
    if (substr(options, i + 1, 1) == ":") {
        # get option argument
        if (length(substr(argv[Optind], _opti + 1)) > 0)
            Optarg = substr(argv[Optind], _opti + 1)
        else
            Optarg = argv[++Optind]
        _opti = 0
    } else
        Optarg = ""
    if (_opti == 0 || _opti >= length(argv[Optind])) {
        Optind++
        _opti = 0
    } else
        _opti++
    return thisopt
}
#---------------------------------------------------------------------------------------------------------

BEGIN {
    # Initializations
    Debug  = 0
    Opterr = 1    # default is to diagnose
    Optind = 1    # skip ARGV[0]
    
    progname = "fasta_toolkit.awk"
    version  = 0.7  # v0.7 Dec 9, 2024, added improved documentation to parse sequence-length stats with col_sumStats.sh
                    # v0.6 Nov 2, 2024, added -R 6, to compute only sequence lengths for DNA or Protein sequences; improved documentation
                    # v0.5 dec 30, 2020, added -M string for -R 1, to select sequences not matching string (revMatch)
                    # v0.4 dec 25, 2020. added PROCINFO["sorted_in"] = "@ind_num_asc" to print sorted results
                    # v0.3 dec 23, 2020. Prints warning and does not exit, if dna_seq not divisible by 3
                    # v0.2 dec 22, 2020, improved layout; fixed typos
                    # v0.1 dec 21, 2020, first commit

    # print the FASTA sequences stored in hashes in ascending order by locus_tag number 
    PROCINFO["sorted_in"] = "@ind_num_asc"

    # check that the script receives input either from file or pipe
    if ( ARGC < 2 ) print_help(progname, version)

    while ((c = getopt(ARGC, ARGV, "dm:M:e:R:s:")) != -1) {
        if (c == "R") {
	    runmode = Optarg
	} 
	else if (c == "m") {
	    string = Optarg
	} 
	else if (c == "M") {
	    RevMatch = Optarg
	} 
	else if (c == "s") {
	    start = Optarg
	} 
	else if (c == "e") {
	    end = Optarg
	} 
	else if (c == "d") {
	    Debug = 1
	}
	else {
            print "ERROR: option -" Optopt, "is not defined"
	    print_help(progname, version)
	}    
    }

    #input_fasta = ARGV[ARGC-1]

    #Debug=1 --> check ARGS
    if(Debug) # activate with -d option
    {
        print "ARGC="ARGC
        print "ARGV[ARGC-1]="ARGV[ARGC-1]

        for(i=1; i<=ARGC; i++)
        {
	   print "ARGV["i"]=" ARGV[i]
        }
    }
    
    # clear out options
    for (i = 1; i < Optind; i++)
       ARGV[i] = ""

  #------------------------- END GETOPT --------------------------#

    # Model FASTA file
    RS=">"
    FS=OFS="\n"

  #-----------------------  END MODEL FASTA  ---------------------#
   
    # check that the user provided the required options and arguments, depending on runmode
    if (runmode == 1 && ( ! string && ! RevMatch)) {
        print "ERROR: -R 1 requires -m match_string to filter the input FASTA" > "/dev/stderr"
        print_help(progname, version)
    }   

    if (runmode == 3 && length(start) == 0) {
        print "ERROR: -R 3 requires -s <int> to provide the START coordinate to extract the sequence string from the input FASTA" > "/dev/stderr"
        print_help(progname, version)
    }   
    
    if (runmode == 3 && length(end) == 0) {
        print "ERROR: -R 3 requires -e <int> to provide the END coordinate to extract the sequence string from the input FASTA" > "/dev/stderr"
        print_help(progname, version)
    }   
 
    # set runmode-specific hashes
    if (runmode == 2) {
        # complement sequences
        compl_nucl["T"]="A"
        compl_nucl["A"]="T"
        compl_nucl["C"]="G"
        compl_nucl["G"]="C"
        compl_nucl["N"]="N"
    }

    if(runmode == 4) {
        # initialize a hash named "codons" holding the codon-aminoacid pairs, 
        #   based on the universal genetic code
        codons["ATA"]="I"; codons["ATC"]="I"; codons["ATT"]="I"; codons["ATG"]="M";
        codons["ACA"]="T"; codons["ACC"]="T"; codons["ACG"]="T"; codons["ACT"]="T";
        codons["AAC"]="N"; codons["AAT"]="N"; codons["AAA"]="K"; codons["AAG"]="K";
        codons["AGC"]="S"; codons["AGT"]="S"; codons["AGA"]="R"; codons["AGG"]="R";
        codons["CTA"]="L"; codons["CTC"]="L"; codons["CTG"]="L"; codons["CTT"]="L";
        codons["CCA"]="P"; codons["CCC"]="P"; codons["CCG"]="P"; codons["CCT"]="P";
        codons["CAC"]="H"; codons["CAT"]="H"; codons["CAA"]="Q"; codons["CAG"]="Q";
        codons["CGA"]="R"; codons["CGC"]="R"; codons["CGG"]="R"; codons["CGT"]="R";
        codons["GTA"]="V"; codons["GTC"]="V"; codons["GTG"]="V"; codons["GTT"]="V";
        codons["GCA"]="A"; codons["GCC"]="A"; codons["GCG"]="A"; codons["GCT"]="A";
        codons["GAC"]="D"; codons["GAT"]="D"; codons["GAA"]="E"; codons["GAG"]="E";
        codons["GGA"]="G"; codons["GGC"]="G"; codons["GGG"]="G"; codons["GGT"]="G";
        codons["TCA"]="S"; codons["TCC"]="S"; codons["TCG"]="S"; codons["TCT"]="S";
        codons["TTC"]="F"; codons["TTT"]="F"; codons["TTA"]="L"; codons["TTG"]="L";
        codons["TAC"]="Y"; codons["TAT"]="Y"; codons["TAA"]="*"; codons["TAG"]="*";
        codons["TGC"]="C"; codons["TGT"]="C"; codons["TGA"]="*"; codons["TGG"]="W";
    }

    if (runmode == 5) {
         # print table header for sequence stats (-R 5)
         print "seq_name\tA\tC\tG\tT\tN\tlength\tGC%"
    }
}
# -------------------- # 
# >>> MAIN PROGRAM <<< #  
# -------------------- # 

NR > 1 { 
   # read the CDS (DNA) sequence into the global seqs hash


   if (runmode > 6) {
      printf "ERROR: runmomde %d is not defined\n", runmode
      print_help(progname, version)
   }
   else {
      read_fasta()
   }
}

END { 
    for (h in seqs) { 
	if( length(seqs[h] >= 2) && runmode == 1 && string) 
	{ 
	    if(h ~ string) print ">"h, seqs[h] 
	}

	if( length(seqs[h] >= 2) && runmode == 1 && RevMatch) 
	{ 
	    if(h !~ RevMatch) print ">"h, seqs[h] 
	}

	if( length(seqs[h] >= 2) && runmode == 2 ) { rev_compl(h, seqs[h]) }
        if( length(seqs[h] >= 3) && runmode == 3 ) { extract_sequence_by_coordinates(h, seqs[h], start, end) }
        if( length(seqs[h] >= 3) && runmode == 4 ) 
	{
           # 1. translate the CDS
	   translate_dna(h, seqs[h]) 
           
	   # 2. print tranlsated fasta
           for (h in prots) {
	       # make sure we print only non-redundant sequences
	       hcount[h]++
	       if ($1 != "" && hcount[h] == 1) printf ">%s\n%s\n", h, prots[h] 
	   }    
        }
	if( length(seqs[h] >= 2) && runmode == 5 ) { print_sequence_stats(h, seqs[h]) }
	if( length(seqs[h] >= 2) && runmode == 6 ) { print_sequence_lengths(h, seqs[h]) }
    }
}
