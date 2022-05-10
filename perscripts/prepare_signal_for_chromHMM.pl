#!/usr/bin/perl
use strict;
use warnings;

my $bs = 5000;

#### Cell \t chr1
#### Mark1 \t Mark2 \t Mark3

#### mCPG mCHG mCHH hCPG hCHG hCHH
my %cells1;
my %cells2;
open IN, "/projects/ren-transposon/home/chz272/05.mESC/02.PBMC/01.clustering/2022_01_20_Seurat_annotation.xls" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $type = substr($tmp[0],0,2);
	my $cid = substr($tmp[0],3,8);
	my $anno = $tmp[1];
	if($type == "01"){
		$cells1{$cid} = $anno;
	}
	if($type == "02"){
		$cells2{$cid} = $anno;
	}
}
close IN;


my @cells = qw/T_naive T_reg B NK Monocytes/;
my @marks = qw/5mC_CPG 5mC_CHG 5mC_CHH 5hmC_CPG 5hmC_CHG 5hmC_CHH/;

my $chrom_size = "/projects/ps-renlab/chz272/genome_ref/hg38/hg38.sim.sizes.xls";
my %chroms;
open IN, $chrom_size or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[0];
	foreach my $pos (1 .. int($tmp[1]/$bs)+1){
		$chroms{$chr}{$pos} = 1;
	}
}
close IN;

my %data;
foreach my $mark (@marks){
	open IN, "zcat /projects/ren-transposon/home/chz272/transposon/02.others/04.pbmc/01.source/PBMC_1_hg38_sorted.bam_$mark\.srf|" or die $!;
	print STDERR "Processing sub 1 $mark...\n";
	my $line = 0;
	my $first = 0;
	while(<IN>){
		$line++;
		if($line > 1000000){
			#last;
			$first++;
			print STDERR  "$first M lines processed\n";
			$line=0;
		}
		chomp;
		my @tmp = split/\s+/, $_;
		next if not exists $cells1{$tmp[0]};
		my $anno = $cells1{$tmp[0]};
		next if not exists $chroms{$tmp[1]};
		my $pos = int($tmp[2]/$bs)+1;
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"} = 0 if not exists $data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"};
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"}++;
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"} = 0 if not exists $data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"};
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"}+=$tmp[5];
	}
	close IN;
	open IN, "zcat /projects/ren-transposon/home/chz272/transposon/02.others/04.pbmc/01.source/PBMC_2_hg38_sorted.bam_$mark\.srf|" or die $!;
	print STDERR "Processing sub 2 $mark...\n";
	$line = 0;
	$first = 0;
	while(<IN>){
		$line++;
		if($line > 1000000){
			#last;
			$first++;
			print STDERR  "$first M lines processed\n";
			$line=0;
		}
		chomp;
		my @tmp = split/\s+/, $_;
		next if not exists $cells2{$tmp[0]};
		my $anno = $cells2{$tmp[0]};
		next if not exists $chroms{$tmp[1]};
		my $pos = int($tmp[2]/$bs)+1;
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"} = 0 if not exists $data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"};
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Total"}++;
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"} = 0 if not exists $data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"};
		$data{$anno}{$mark}{$tmp[1]}{$pos}{"Value"}+=$tmp[5];
	}
	close IN;
}

system("mkdir 2022_01_22_methylStates");
print STDERR "Outputing states...\n";
foreach my $anno (@cells){
	open OUT, ">2022_01_22_methylStates/$anno\_signal" or die $!;
	foreach my $chr (sort keys %chroms){
		print OUT "$anno\t$chr\n";
		print OUT "5mC_CPG\t5mC_CHG\t5mC_CHH\t5hmC_CPG\t5hmC_CHG\t5hmC_CHH\n";
		foreach my $pos (sort {$a<=>$b} keys %{$chroms{$chr}}){
			my $output = "";
			foreach my $mark (@marks){
				my $value = 0;
				if(exists $data{$anno}{$mark}{$chr}{$pos}{"Value"}){
					my $t0 = $data{$anno}{$mark}{$chr}{$pos}{"Total"};
					my $t1 = $data{$anno}{$mark}{$chr}{$pos}{"Value"};
					$value = int($t1/$t0 * 1000);
				}
				$output .= "$value\t";
			}
			chomp($output);
			$output .= "\n";
			print OUT $output;
		}
	}
	close OUT;
}