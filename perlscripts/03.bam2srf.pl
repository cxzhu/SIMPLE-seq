#!/usr/bin/perl
use strict;
use warnings;

### read refs
my $genome_ref = "path to your reference fasta file"; ### modify the path here.
if($genome_ref eq "path to your reference fasta file"){
	print "Please modify the genome reference in the script file\n";
	exit(0);
}
open IN, $genome_ref or die $!;
my $chr = "";
my %fa;
print STDERR "Read reference file...\n";
while(<IN>){
	chomp;
	if(m/\>/){
		$chr = substr($_, 1, length($_)-1);
		$fa{$chr} = "";
		next;
	}
	$fa{$chr} .= $_;
}
close IN;
print STDERR "Read reference file finished...\n";

open IN, "samtools view $ARGV[0]|" or die $!;
my $prefix = substr($ARGV[0], 0, length($ARGV[0])-4);



open OUT_CPG, "|gzip - >$prefix\_CPG.srf" or die $!;
open OUT_CHG, "|gzip - >$prefix\_CHG.srf" or die $!;
open OUT_CHH, "|gzip - >$prefix\_CHH.srf" or die $!;

my $cur_chr = "";
my %cpg;
my %chg;
my %chh;
my $meta = 0;
my $chr_length = 0;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp [2];
	my $p = $tmp[3];
	my $cigar = $tmp[5];
	my $seq = $tmp[9];
	my @sp = split/[A-Z]/, $cigar;
	my $total_cp = 0;
	my $total_len = 0;

	### correct seq
	foreach my $num (@sp){
		my $mark = substr($total_cp + length($num), 1);
		if($mark eq "S" && $total_cp == 0){
			$seq = substr($seq, $num, length($seq) - $num);
		}
		if($mark eq "S" && $total_cp !=0){
			$seq = substr($seq, 0, length($seq)- $num);
		}
		if($mark eq "M"){
			$total_len += $num;
		}
		if($mark eq "I"){
			my $tr_seq = substr($seq, 0, $total_len).substr($seq, $total_len + $num, length($seq)-$total_len - $num);
			$seq = $tr_seq;
		}
		if($mark eq "D"){
			my $tr_seq = substr($seq, 0, $total_len);
			my $i = 0;
			while($i<$num){
				$tr_seq.="N";
				$i++;
			}
			$tr_seq.=substr($seq, $total_len, length($seq)-$total_len);
			$total_len += $num;
			$seq = $tr_seq;
		}
		$total_cp+=length($num);
		$total_cp++;
	}

	$meta = $p + 1000000 if $meta == 0;
	if( ($cur_chr ne $chr or $p > $meta)){
		if($cur_chr ne "" or $p > $meta){
			## output hash
			# 1 output cpg;
			foreach my $pos (keys %cpg){
				foreach my $bc (keys %{$cpg{$pos}}){
					my $c2t = $cpg{$pos}{$bc}{"c2t"};
					my $tot = $cpg{$pos}{$bc}{"tot"};
					my $methylated = 0;
					$methylated = 1 if $c2t/$tot > 0.5;
					my $output = "$bc\t$chr\t$pos\t$tot\t$c2t\t$methylated\n";
					print OUT_CPG $output;
				}
			}
			# 2  output chg;
			foreach my $pos (keys %chg){
				foreach my $bc (keys %{$chg{$pos}}){
					my $c2t = $chg{$pos}{$bc}{"c2t"};
					my $tot = $chg{$pos}{$bc}{"tot"};
					my $methylated = 0;
					$methylated = 1 if $c2t/$tot > 0.5;
					my $output = "$bc\t$chr\t$pos\t$tot\t$c2t\t$methylated\n";
					print OUT_CHG $output;
				}
			}			
			# 3  output chh;
			foreach my $pos (keys %chh){
				foreach my $bc (keys %{$chh{$pos}}){
					my $c2t = $chh{$pos}{$bc}{"c2t"};
					my $tot = $chh{$pos}{$bc}{"tot"};
					my $methylated = 0;
					$methylated = 1 if $c2t/$tot > 0.5;
					my $output = "$bc\t$chr\t$pos\t$tot\t$c2t\t$methylated\n";
					print OUT_CHH $output;
				}
			}			
			## clear hash
			%cpg = ();
			%chg = ();
			%chh = ();
		}
		print STDERR "Processing $chr to $meta 000000 bp... \n";
		$cur_chr = $chr;
		$meta = $p + 1000000;
		$chr_length = length($fa{$chr})-100;
	}

	#my $seq = $tmp[9];
	my $len = length($seq) - 1;
	next if $p < 100;
	next if $len > $chr_length;
	my $cid = substr($tmp[0], -8, 8);
	foreach my $cur_p (0 .. $len){
		my $cur_pos = $cur_p + $p;
		my $base = substr($seq, $cur_p, 1);
		my $cur_base = uc(substr($fa{$chr}, $cur_pos-1, 1));
		if($cur_base eq 'C'){
			my $p1 = uc(substr($fa{$chr}, $cur_pos, 1));
			if($p1 eq 'G'){
				$cpg{$cur_pos}{$cid}{"tot"} = 0 if not exists $cpg{$cur_pos}{$cid}{"tot"};
				$cpg{$cur_pos}{$cid}{"tot"}++;
				$cpg{$cur_pos}{$cid}{"c2t"} = 0 if not exists $cpg{$cur_pos}{$cid}{"c2t"};
				$cpg{$cur_pos}{$cid}{"c2t"}++ if $base eq 'T';
				next;
			}
			else{
				my $p2 = uc(substr($fa{$chr}, $cur_pos+1, 1));
				if($p2 eq 'G'){
					$chg{$cur_pos}{$cid}{"tot"} = 0 if not exists $chg{$cur_pos}{$cid}{"tot"};
					$chg{$cur_pos}{$cid}{"tot"}++;
					$chg{$cur_pos}{$cid}{"c2t"} = 0 if not exists $chg{$cur_pos}{$cid}{"c2t"};
					$chg{$cur_pos}{$cid}{"c2t"}++ if $base eq 'T';
					next;
				}
				else{
					$chh{$cur_pos}{$cid}{"tot"} = 0 if not exists $chh{$cur_pos}{$cid}{"tot"};
					$chh{$cur_pos}{$cid}{"tot"}++;
					$chh{$cur_pos}{$cid}{"c2t"} = 0 if not exists $chh{$cur_pos}{$cid}{"c2t"};
					$chh{$cur_pos}{$cid}{"c2t"}++ if $base eq 'T';
					next;					
				}
			}
		}
		elsif($cur_base eq 'G'){
			my $p1 = uc(substr($fa{$chr}, $cur_pos-2, 1));
			if($p1 eq 'C'){
				$cpg{$cur_pos}{$cid}{"tot"} = 0 if not exists $cpg{$cur_pos}{$cid}{"tot"};
				$cpg{$cur_pos}{$cid}{"tot"}++;
				$cpg{$cur_pos}{$cid}{"c2t"} = 0 if not exists $cpg{$cur_pos}{$cid}{"c2t"};
				$cpg{$cur_pos}{$cid}{"c2t"}++ if $base eq 'A';
				next;
			}
			else{
				my $p2 = uc(substr($fa{$chr}, $cur_pos-3, 1));
				if($p2 eq 'C'){
					$chg{$cur_pos}{$cid}{"tot"} = 0 if not exists $chg{$cur_pos}{$cid}{"tot"};
					$chg{$cur_pos}{$cid}{"tot"}++;
					$chg{$cur_pos}{$cid}{"c2t"} = 0 if not exists $chg{$cur_pos}{$cid}{"c2t"};
					$chg{$cur_pos}{$cid}{"c2t"}++ if $base eq 'A';
					next;
				}
				else{
					$chh{$cur_pos}{$cid}{"tot"} = 0 if not exists $chh{$cur_pos}{$cid}{"tot"};
					$chh{$cur_pos}{$cid}{"tot"}++;
					$chh{$cur_pos}{$cid}{"c2t"} = 0 if not exists $chh{$cur_pos}{$cid}{"c2t"};
					$chh{$cur_pos}{$cid}{"c2t"}++ if $base eq 'A';
					next;					
				}
			}		
		}
	}
}

foreach my $pos (keys %cpg){
	foreach my $bc (keys %{$cpg{$pos}}){
		my $c2t = $cpg{$pos}{$bc}{"c2t"};
		my $tot = $cpg{$pos}{$bc}{"tot"};
		my $methylated = 0;
		my $pos1 = $pos;
		$methylated = 1 if $c2t/$tot > 0.5;
		my $output = "$bc\t$chr\t$pos1\t$tot\t$c2t\t$methylated\n";
		print OUT_CPG $output;
	}
}
foreach my $pos (keys %chg){
	foreach my $bc (keys %{$chg{$pos}}){
		my $c2t = $chg{$pos}{$bc}{"c2t"};
		my $tot = $chg{$pos}{$bc}{"tot"};
		my $methylated = 0;
		my $pos1 = $pos;
		$methylated = 1 if $c2t/$tot > 0.5;
		my $output = "$bc\t$chr\t$pos1\t$tot\t$c2t\t$methylated\n";
		print OUT_CHG $output;
	}
}			
foreach my $pos (keys %chh){
	foreach my $bc (keys %{$chh{$pos}}){
		my $c2t = $chh{$pos}{$bc}{"c2t"};
		my $tot = $chh{$pos}{$bc}{"tot"};
		my $methylated = 0;
		my $pos1 = $pos;
		$methylated = 1 if $c2t/$tot > 0.5;
		my $output = "$bc\t$chr\t$pos1\t$tot\t$c2t\t$methylated\n";
		print OUT_CHH $output;
	}
}	
close OUT_CPG;
close OUT_CHG;
close OUT_CHH;
close IN;
