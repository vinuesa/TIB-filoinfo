#!/usr/bin/env bash

progname=${0##*/} # simplify_domtblout.sh
vers=2024-11-06 # v2024-11-06 USA post-election trauma; first commit;  

function print_help(){
    cat << HELP
       USAGE for $progname v$vers
       
       $progname <tblout_file> to parse
       
   AIM: prints: target\taccession\ttlen\tquery_name\tqlen\tdE-value\tdScore\tDbias\ti-Evalue\tscore\tbias\thmmStart\thmmEnd\tenvStart\tenvEnd\tacc 
         table with columnt -t formatting, for convenient analysis of hmmer3 domtblout tables to identify proper ga, tc and nc cutoffs  
       
HELP

   exit 1

}

# 1. Check and capture required user input
(($# < 1)) && print_help

tblout_file=$1
                           #  1       2         3       4         6       7         8      9    13        14    15    16        17       20        21     22
awk 'BEGIN{OFS="\t"; print "target\taccession\ttlen\tquery_name\tqlen\tdE-value\tdScore\tDbias\ti-Evalue\tscore\tbias\thmmStart\thmmEnd\tenvStart\tenvEnd\tacc"} 
      !/^#/ {print $1,$2,$3,$4,$6,$7,$8,$9,$13,$14,$15,$16,$17,$20,$21,$22}' "${tblout_file}" | column -t -s $'\t'


# Some fields to simplify tblout
#awk 'BEGIN{OFS="\t"; print "target\tsE-value\tsScore\tDescr" }$0 !~/^#/{print $1,$5,$6,$8,$9,$19,$20,$21,$22,$23,$24,$25}' "${tblout_file}" | column -t -s $'\t'
