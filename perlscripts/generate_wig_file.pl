#!/usr/bin/perl
use strict;
use warnings;

## generate_wig_file.pl input.srf cell_type.xls
## default bs = 1000
my $bs = 1;

open IN, $ARGV[1] or die $!;
my %cell_types;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$cell_types{$tmp[0]} = $tmp[1];
}
close IN;


my %hash;
open IN, "zcat $ARGV[0]|" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $cell_id = $tmp[0];
	next if not exists $cell_types{$cell_id};
	my $type = $cell_types{$cell_id};
	my $chr = $tmp[1];
	my $epos = $tmp[2];
	my $pos = int($tmp[2] / $bs);
	my $value = $tmp[5];
	if(not exists $hash{$type}{$chr}{$pos}{"total"}){
		$hash{$type}{$chr}{$pos}{"total"} = 0;
		$hash{$type}{$chr}{$pos}{"methy"} = 0;
	}
	$hash{$type}{$chr}{$pos}{"total"}++;
	$hash{$type}{$chr}{$pos}{"methy"} += $value;
}
close IN;


my $prefix = substr($ARGV[0], 0, length($ARGV[0]) - 4);

foreach my $type (keys %hash){
	open OUT, "|gzip - >$prefix\_$bs\_$type.wig.gz" or die $!;
	print OUT "track type wiggle_0\n";
	foreach my $chr (sort keys %{$hash{$type}}){
		print OUT "variableStep chrom=$chr span=$bs\n";
		foreach my $pos (sort {$a<=>$b} keys %{$hash{$type}{$chr}}){
			my $spos = $pos*$bs+1;
			my $m = $hash{$type}{$chr}{$pos}{"methy"};
			my $t = $hash{$type}{$chr}{$pos}{"total"};
			my $value = int($m / $t * 100);
			print OUT "$spos\t$value\n";
		}
	}
	close OUT;
}



