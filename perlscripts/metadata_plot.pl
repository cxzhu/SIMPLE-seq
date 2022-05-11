#!/usr/bin/perl
use strict;
use warnings;
my $type1 = "naive";
#my $c = "5mC";

open IN, "RNA_seq_RPKM.xls" or die $!;
print STDERR "Reading gene expression...\n";
my %genes;
my %hi;
my %me;
my %lo;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $gene_id = $tmp[1];
	#next if substr($gene_id, 0, 2) ne "EN";
	my $fpkm = $tmp[14];
	$fpkm = $tmp[15] if $type1 eq "primed";
	$hi{$gene_id} = 1 if $fpkm>=50;
	$me{$gene_id} = 1 if $fpkm > 1 && $fpkm < 50;
	$lo{$gene_id} = 1 if $fpkm <=1; #&& $fpkm > 0.1;
}
close IN;


# 1-25	26-75		76-100
# TSS		Genbody	TES
my $bs = 50;
my %anno_hi;
my %anno_me;
my %anno_lo;
open IN, "zcat /projects/ren-transposon/home/chz272/05.mESC/01.mESC/03.TAPS_WGBS/gencode.vM21.basic.annotation.gtf.gz|" or die $!;
print STDERR "Reading gene annotation...\n";
while(<IN>){
	next if substr($_, 0, 3) ne "chr";

	chomp;
	my @tmp = split/\t+/, $_;
	next if $tmp[2] ne "transcript";
	#$tmp[8] =~ m/gene_id "(.*)"; transcript_id/;
	$tmp[8] =~ m/gene_name "(.*)"; transcript_type/;
	next if not defined $1;
	#print "$1\n";
	#exit();
	my $gid=$1;
	if(exists $hi{$gid}){
		if($tmp[6] eq "+"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[3] - (26-$cur_step)*$step_size;
				my $pos_end = $tmp[3] - (25-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 75){
				my $pos_start = $tmp[3]+($cur_step - 26)*$step_size;
				my $pos_end = $tmp[3]+($cur_step - 25)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[4]+($cur_step - 76)*$step_size;
				my $pos_end = $tmp[4]+($cur_step - 75)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}
		}
		if($tmp[6] eq "-"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[4] + (25-$cur_step)*$step_size;
				my $pos_end = $tmp[4] + (26-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 76){
				my $pos_start = $tmp[4]-($cur_step - 25)*$step_size;
				my $pos_end = $tmp[4]-($cur_step - 26)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[3]-($cur_step - 100)*$step_size;
				my $pos_end = $tmp[3]-($cur_step - 101)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}	
		}
	}
	if(exists $me{$gid}){
		if($tmp[6] eq "+"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[3] - (26-$cur_step)*$step_size;
				my $pos_end = $tmp[3] - (25-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_me{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 75){
				my $pos_start = $tmp[3]+($cur_step - 26)*$step_size;
				my $pos_end = $tmp[3]+($cur_step - 25)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_hi{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[4]+($cur_step - 76)*$step_size;
				my $pos_end = $tmp[4]+($cur_step - 75)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_me{$tmp[0]}{$spos} = $cur_step;
				}
			}
		}
		if($tmp[6] eq "-"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[4] + (25-$cur_step)*$step_size;
				my $pos_end = $tmp[4] + (26-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_me{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 76){
				my $pos_start = $tmp[4]-($cur_step - 25)*$step_size;
				my $pos_end = $tmp[4]-($cur_step - 26)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_me{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[3]-($cur_step - 100)*$step_size;
				my $pos_end = $tmp[3]-($cur_step - 101)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_me{$tmp[0]}{$spos} = $cur_step;
				}
			}	
		}
	}
	if(exists $lo{$gid}){
		if($tmp[6] eq "+"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[3] - (26-$cur_step)*$step_size;
				my $pos_end = $tmp[3] - (25-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 75){
				my $pos_start = $tmp[3]+($cur_step - 26)*$step_size;
				my $pos_end = $tmp[3]+($cur_step - 25)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[4]+($cur_step - 76)*$step_size;
				my $pos_end = $tmp[4]+($cur_step - 75)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}
		}
		if($tmp[6] eq "-"){
			my $step_size = 200;
			foreach my $cur_step(1 .. 25){
				my $pos_start = $tmp[4] + (25-$cur_step)*$step_size;
				my $pos_end = $tmp[4] + (26-$cur_step)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}
			$step_size = ($tmp[4] - $tmp[3]) / 50;
			foreach my $cur_step(26 .. 76){
				my $pos_start = $tmp[4]-($cur_step - 25)*$step_size;
				my $pos_end = $tmp[4]-($cur_step - 26)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}			
			$step_size = 200;
			foreach my $cur_step(76 .. 100){
				my $pos_start = $tmp[3]-($cur_step - 100)*$step_size;
				my $pos_end = $tmp[3]-($cur_step - 101)*$step_size;
				foreach my $spos (int($pos_start / $bs) .. int($pos_end / $bs)){
					$anno_lo{$tmp[0]}{$spos} = $cur_step;
				}
			}	
		}
	}
}

my %data_hi;
my %data_me;
my %data_lo;

print STDERR "Reading gene bed...\n";
foreach my $type (qw/5mC 5hmC/){
	open IN, "zcat 2022_01_22_SIMPLE_mESC_$type1\_$type\_CPG.bed.gz|" or die $!;
	while(<IN>){
		chomp;
		my @tmp = split/\s+/, $_;
		next if(($tmp[4]+$tmp[5])<5);
		my $spos = int($tmp[1]/$bs);
		if(exists $anno_hi{$tmp[0]}{$spos}){
			my $step = $anno_hi{$tmp[0]}{$spos};
			if(not exists $data_hi{$type}{$step}{"total"}){
				$data_hi{$type}{$step}{"total"}=0;
				$data_hi{$type}{$step}{"value"}=0;
			}
			$data_hi{$type}{$step}{"total"}++;
			$data_hi{$type}{$step}{"value"}+=$tmp[3];
		}
		if(exists $anno_me{$tmp[0]}{$spos}){
			my $step = $anno_me{$tmp[0]}{$spos};
			if(not exists $data_me{$type}{$step}{"total"}){
				$data_me{$type}{$step}{"total"}=0;
				$data_me{$type}{$step}{"value"}=0;
			}
			$data_me{$type}{$step}{"total"}++;
			$data_me{$type}{$step}{"value"}+=$tmp[3];
		}
		if(exists $anno_lo{$tmp[0]}{$spos}){		
			my $step = $anno_lo{$tmp[0]}{$spos};
			if(not exists $data_lo{$type}{$step}{"total"}){
				$data_lo{$type}{$step}{"total"}=0;
				$data_lo{$type}{$step}{"value"}=0;
			}
			$data_lo{$type}{$step}{"total"}++;
			$data_lo{$type}{$step}{"value"}+=$tmp[3];
		}
	}
}
close IN;


open OUT, ">2022_01_23_gene_metaplot\_$type1\.xls" or die $!;
foreach my $type (qw/5mC 5hmC/){
	print OUT "$type\_High";
	foreach my $step (1..100){
		my $value = 0;
		if(exists $data_hi{$type}{$step}{"total"}){
			$value = $data_hi{$type}{$step}{"value"}/$data_hi{$type}{$step}{"total"};
		}
		print OUT "\t$value";
	}
	print OUT "\n";
	print OUT "$type\_Medium";
	foreach my $step (1..100){
		my $value = 0;
		if(exists $data_me{$type}{$step}{"total"}){
			$value = $data_me{$type}{$step}{"value"}/$data_me{$type}{$step}{"total"};
		}
		print OUT "\t$value";
	}
	print OUT "\n";
	print OUT "$type\_Low";
	foreach my $step (1..100){
		my $value = 0;
		if(exists $data_lo{$type}{$step}{"total"}){
			$value = $data_lo{$type}{$step}{"value"}/$data_lo{$type}{$step}{"total"};
		}
		print OUT "\t$value";
	}
	print OUT "\n";
}
close OUT;




