#!/usr/bin/perl
use strict;
use warnings;

# ss61.readInfor2Mtx.pl input.readinfo.txt.gz binsize
# mtx output
# value = percent methylation [0,100]

if($#ARGV<1){
	print "perl perlscripts/04.srf2mtx.pl input.rsf binsize\n";
	exit();
}

my $bs = $ARGV[1];
my $input = $ARGV[0];
my $prefix = substr($ARGV[0], 0, length($ARGV[0])-4);

my %blk;

system("mkdir $prefix\_$bs\_blk");
my %data;
my %cid;
my %bid;

open IN, "zcat $input|" or die $!;
my $n_values = 0;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $cell_id = $tmp[0];
	my $chr = $tmp[1];
	my $pos = int($tmp[2] / $bs);
	my $value = $tmp[5];
	my $bin = "$chr\_$pos";
	$cid{$cell_id} = 1;
	$bid{$chr}{$pos} = 1;
	if(not exists $data{$cell_id}{$chr}{$pos}{"m"}){
		$data{$cell_id}{$chr}{$pos}{"m"} = 0;
		$data{$cell_id}{$chr}{$pos}{"u"} = 0;
		$n_values++;
	}
	$data{$cell_id}{$chr}{$pos}{"m"} ++ if $value == 1;
	$data{$cell_id}{$chr}{$pos}{"u"} ++ if $value == 0;
}
close IN;

my $index = 0;
open OUT, ">$prefix\_$bs\_blk\/barcodes.tsv" or die $!;
foreach my $cell_id (sort keys %cid){
	$index++;
	$cid{$cell_id} = $index;
	print OUT $cell_id."\n";
}
close OUT;
my $n_cells = $index;

$index = 0;
open OUT1, ">$prefix\_$bs\_blk\/genes.tsv" or die $!;
open OUT2, ">$prefix\_$bs\_blk\/peaks.bed" or die $!;
foreach my $chr (sort keys %bid){
	foreach my $pos (sort keys %{$bid{$chr}}){
		$index++;
		$bid{$chr}{$pos} = $index;
		my $pos_s = $pos*$bs;
		my $pos_e = $pos_s + $bs;
		print OUT1 "$chr\:$pos_s\-$pos_e\t$chr\:$pos_s\-$pos_e\n";
		print OUT2 "$chr\t$pos_s\t$pos_e\t$chr\:$pos_s\-$pos_e\n";
	}
}
close OUT1;
close OUT2;
my $n_bins = $index;

open OUT, ">$prefix\_$bs\_blk\/matrix.mtx" or die $!;
print OUT "\%\%MatrixMarket matrix coordinate real general\n\%\n$n_bins $n_cells $n_values\n";
foreach my $cell_id (sort keys %data){
	my $index_cell = $cid{$cell_id};
	foreach my $chr (sort keys %{$data{$cell_id}}){
		foreach my $pos (sort keys %{$data{$cell_id}{$chr}}){
			my $index_bin = $bid{$chr}{$pos};
			my $mC = $data{$cell_id}{$chr}{$pos}{"m"};
			my $uC = $data{$cell_id}{$chr}{$pos}{"u"};
			my $cur_value = 0;
			$cur_value = int(100* $mC / ($mC+$uC) ) if ($mC + $uC)>0;
			my $output = "$index_bin $index_cell $cur_value\n";
			print OUT $output;
		}
	}
}
close OUT;







