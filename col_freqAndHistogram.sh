#!/usr/bin/env bash

# AUTHOR: Pablo Vinuesa; https://www.ccg.unam.mx/~vinuesa/ @pvinmex
# little awk script to compute basic frequency analysis and raw dotted histogram
#  receives the output of cut for a numeric column in a table 

progname=${0##*/}
min_count=1
scale=100

function print_help()
{
   cat << _EOF_
     $progname usage: expects numeric data from a table column extracted for example with cut,
         excluding the header!
     
   EXAMPLE:  
     cut -fx table | sed '1d' | $progname [optional: min_count, default:<$min_count>, to be considered] [scaling factor, def:<$scale>]
     
   AIM:
      little awk script to compute basic frequency analysis and raw dotted histogram
         receives the output of cut for a numeric column in a table   
_EOF_

  exit 1
}

# test if stdin is a tty, and not a pipe
[ -t 0 ] && print_help

# receive optional user input
min_count=${1:-1}
scale=${2:-100}

awk -v minCount="$min_count" -v scale="$scale" '
{l=$1; a[l]++; if (a[l]>max) max=a[l]} 

END {printf("Length\tFrequency"); 
  for (i in a) 
    if(a[i] >= minCount){printf("\n%s\t%s\t",i,a[i]); 
      for (j=0;j<(int(a[i]*(scale/max)));j++)
        printf("*")} print ""} 
'

