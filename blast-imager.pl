#!/usr/bin/env perl
# blast-imager.pl
use strict;
use warnings;
use File::Basename;
use GD;

# Adapted by Pablo Vinuesa from O'Reilly's BLAST book
# By Ian Korf, Mark Yandell, Joseph Bedell
# Publisher: O'Reilly Media
# Final Release Date: July 2003 

my $progname = basename($0);
my $VERSION  = '0.3_26Jul19'; # v0.3_26Jul19: added blastx (blast+) invocation example
  # 8Sep15; added the -t in if ( -t and scalar @ARGV < 1 ){ print_help($progname,$VERSION)  }
  #   to allow the script to work both with pre-existing blast_output.tab files and directly 
  #   reading blast output from STDIN!

# see the -t test for STDIN at http://www.perlmonks.org/?node_id=202564
# http://stackoverflow.com/questions/518992/how-can-i-check-peek-stdin-for-piped-data-in-perl-without-using-select
if ( -t and scalar @ARGV < 1 ){ print_help($progname,$VERSION)  }

# Parse tabular -m 8 blast data (or standard -outfmt 6 of blast+ programs)
my ($Q, %S, @S);
my ($MAX, $MIN, $HSPs) = (0, 1e20, 0);
while (<>) {
	my ($q, $id, $p, $l, $m, $g, $qb, $qe, $sb, $se, $e, $b) = split;
	$Q =$q;
	if ($qb > $qe) {($qb, $qe, $sb, $se) = ($qe, $qb, $se, $sb)}
	$MAX = $qe if $qe > $MAX;
	$MIN = $qb if $qb < $MIN;
	push @S, $id if not defined $S{$id};
	push @{$S{$id}}, [$qb, $qe, $sb, $se, $p];
	$HSPs++;
}

# Setup graph params
#my ($L, $B, $R, $H, $F) = (150, 600, 50, 40, 20); # graph regions
#my ($W, $Hsep, $Ssep)   = (3, 14, 18);            # line width and spacing
my ($L, $B, $R, $H, $F) = (150,550, 50, 40, 20);   # graph regions
my ($W, $Hsep, $Ssep)   = (3, 14, 18);             # line width and spacing


my $vsize = $H + $F + $Hsep * $HSPs + $Ssep * (keys %S);
$vsize = 100 if $vsize < 100;
$vsize = 4000 if $vsize > 4000;
my $hsize = $L + $B + $R;
my $SCALE = $B / ($MAX - $MIN + 1);  
my $image = new GD::Image($hsize, $vsize);

# Colors
my @Color;
my @data = ([0,0,0], [196,0,255], [0,0,255], [0,255,255], [0,255,0], 
	[255,255,0], [255,196,0], [255,0,0], [128,128,128]);
for (my $i = 0; $i < @data; $i++) {
	$Color[$i] = $image->colorAllocate(@{$data[$i]});
}
my $White = $image->colorAllocate(255,255,255);
my $Black = $Color[0];

# Header
$image->filledRectangle(0, 0, $hsize, $vsize, $White);
$image->string(gdMediumBoldFont, 5, $H-8, substr($Q,0,18), $Black);
$image->line($L, $H, $L+$B, $H, $Black);
$image->string(gdSmallFont, $L, $H-20, $MIN, $Black);
$image->string(gdSmallFont, $L+$B, $H-20, $MAX, $Black);

# Percent identity key
$image->string(gdSmallFont, 670, 5, "% Identity", $Black);
#for (my $i = 20; $i <= 100; $i += 10) {
for (my $i = 92; $i <= 100; $i += 1) {   #change scale from 92%-100% ID
	my $x = $L+$B/2 + $i*2;
	#$image->filledRectangle($x, 5, $x+10, 15, colormap($i));
	$image->filledRectangle($x, 1, $x+1, 10, colormap($i));
	#$image->string(gdTinyFont, $x, 17, $i, $Black);
	$image->string(gdTinyFont, $x, 20, $i, $Black); # $x, px from top, $i, $Black
}

# Alignments
my @Depth;
my $v = 0;
foreach my $id (@S) {
	$v += $Ssep;
	$image->string(gdSmallFont, 10, $H+$v+9, substr($id,0,18), $Black);
	foreach my $hsp (@{$S{$id}}) {
		$v += $Hsep;
		my ($qb, $qe, $sb, $se, $pct) = @$hsp;
		my $strand = $sb < $se ? '+' : '-';
		my ($x1, $x2, $y) = (scale($qb), scale($qe), $H+$v);
		foreach my $x ($x1..$x2) {$Depth[$x]++}
		my $c = colormap($pct);
		$image->filledRectangle($x1, $y, $x2, $y+$W, $c);		
		$image->string(gdTinyFont, $x1 -(5*length($qb)), $y-5, $qb, $Black);
		$image->string(gdTinyFont, $x2+2, $y-5, $qe, $Black);
		$image->string(gdTinyFont, $x1 -(5*length($sb)), $y+2, $sb, $Black);
		$image->string(gdTinyFont, $x2+2, $y+2, "$se $strand", $Black);
	}
}

# Alignment depth
my $MaxDepth = 0;
foreach my $d (@Depth) {$MaxDepth = $d if defined $d and $d > $MaxDepth}
my $Dscale = int($MaxDepth/10) +1;
$image->string(gdTinyFont, $L+$B+2, $H+2, "$Dscale/line", $Black);
for (my $i = 0; $i < @Depth; $i++) {
	next unless defined $Depth[$i];
	my $level = $Depth[$i]/$Dscale +1;
	for (my $j = 0; $j < $level; $j++) {
		$image->line($i, $H+$j*2, $i, $H+$j*2, $Black);
	}
}

# Output (edit this for your installation/taste)
print $image->png;
#print $image->jpg;
# print $image->gif;

sub colormap {
	my ($value) = @_;
	#my $n = ($value >= 100) ? 0: int((109 - $value) / 10);
	my $n = ($value >= 100) ? 0: int((100 - $value) );   # change scale to have values from 92-100
	return defined $Color[$n] ? $Color[$n] : $Color[@Color-1];
}

sub scale {
	my ($x) = @_;
	my $scale = ($x - $MIN) * $SCALE + $L;
	return $scale;
}

sub print_help
{
   my ($prog, $vers) = @_;
   
   print <<HELP;

    $prog v.$vers USAGE:
    
    $prog requires a blast tabular output file provied as single argument or piped in directly from blast
    
    The output should be redirected to a file with png extension.
    
    USAGE EXAMPLES:
    
      1) Calling $prog on a single blast output table
        blast-imager.pl blast_m8_output.tab > blast_graph.png
     
      2) in a shell for loop:
        for file in *fa; do blastall -p blastx -i \$file -d integron_cassettes4blastdb.faaed -m 8 -b 40 | \
	   blast-imager.pl > \${file%.fa}_blast_graph.png; done
      
      3) using blastx (from the blast+ suite; note -outfmt 6 usage, wich is equivalent to blastall -p blastx -m 8)
          for file in *fa; do blastx -query \$file -db integron_cassettes4blastdb.faaED -evalue 1e-20 -outfmt 6 -num_alignments 100 | \
          blast-imager.pl > \${file%.*}.png; done	   
      

HELP

    exit(1);
}
