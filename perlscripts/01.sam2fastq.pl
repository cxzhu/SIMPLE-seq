#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0]."_BC.sam";
open IN, "cat $file|" or die $!;
my $prefix = substr($file, 0, length($file)-7);
open OUT, "|gzip - > $prefix\_cov.fq.gz" or die $!;
while(<IN>){
	next if $_ =~ m/^\@/;
	chomp;
	my @tmp = split/\s+/, $_;
	my $cid = $tmp[2];
	next if $cid eq '*';
	my @sp = split/\:/, $tmp[0];
	my $read1 = substr($tmp[0], 42, 150); ### for BGI-seq data only; if illumina, modify the length based on read name format
	my $qual1 = substr($tmp[0], 193, 150); ### for BGI-seq data only; if illumina, modify the length based on read name format
	my $pre = substr($read1, 0, 6);
	my $read1_s = substr($read1, 6, length($read1)-6);
	my $qual1_s = substr($qual1, 6, length($qual1)-6);
	my $rname = "\@$sp[0]\:$sp[1]\:$pre\:$cid";
	print OUT "$rname\n$read1_s\n+\n$qual1_s\n";
}
close IN;
close OUT;
