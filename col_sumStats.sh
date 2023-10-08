#!/usr/bin/env bash

# little awk script to compute basic summary statistics
#  receives the output of cut for a numeric column in a table 


progname=${0##*/}


function print_help()
{
   cat << _EOF_
   
     $progname usage: expects numeric data from a table column extracted for example with cut
     
     cut -fx table | $progname
_EOF_

  exit 1
}

# test if stdin is a tty, and not a pipe
[ -t 0 ] && print_help

sort -n | awk '
BEGIN{ max = min = Mode = 'NaN' }
{     
     sum+=$1
     a[x++]=$1
     b[$1]++
     if(b[$1]>hf){hf=b[$1]}
}
NR == 1 { min=$1; max=$1 }
$1 < min { min=  $1 }
$1 > max { max=  $1 }

END{ n = asort(a);idx=int((x+1)/2)
    print "Min: " min
    print "Mean: " sum/x 
    print "Median: " ((idx==(x+1)/2) ? a[idx] : (a[idx]+a[idx+1])/2) 
    for (i in b){if(b[i]==hf){(k=="") ? (k=i):(k=k FS i)}{FS=","}}
    print "Mode: " k
    print "Max: " max
}
'
