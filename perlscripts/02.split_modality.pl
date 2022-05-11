#!/usr/bin/perl
use strict;
use warnings;

my $mc = 0;
my $h1 = 0;
my $h2 = 0;
my $h3 = 0;

my $i = 0;

# TTGACAGGCG:TTGATG:56:25:09


open IN, "samtools view -h $ARGV[0]|" or die $!;


open OUT1, "|samtools view - -b > $ARGV[0]\_5mC.bam" or die $!;
open OUT2, "|samtools view - -b > $ARGV[0]\_5hmC.bam" or die $!;
open OUT3, "|samtools view - -b > $ARGV[0]\_other.bam" or die $!;

while(<IN>){
	if(m/^\@/){
		print OUT1 $_;
		print OUT2 $_;
		print OUT3 $_;
	}
	else{
		my @tmp = split/\s+/, $_;
		my $ta = substr($tmp[0], -15, 1);
		my $tb = substr($tmp[0], -11, 1);
		if($ta eq "C" && $tb eq "C"){
			print OUT1 $_;
		}
		elsif($ta eq "T" && $tb eq "C"){
			print OUT3 $_;
		}
		elsif($ta eq "C" && $tb eq "T"){
			print OUT3 $_;
		}
		elsif($ta eq "T" && $tb eq "T"){
			print OUT2 $_;
		}
	}
}
close IN;
close OUT1;
close OUT2;
close OUT3;

