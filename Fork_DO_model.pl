#!/usr/bin/perl -w

use strict;
use warnings;
use Parallel::ForkManager;

my $filename = $ARGV[0];

open FILE, "<$filename" or die $!;
my @seqs;
while (<FILE>){
	if ($_ =~ /^>(\S+)/){
		push @seqs,$1;
	}
}

my $manager = new Parallel::ForkManager( 20 );
my $t = scalar(@seqs);
my $o = 1;
foreach my $seq (@seqs){
	print "$o from $t\n";
	$o++;
		$manager->start and next;
		my $command = "perl DO_model.pl $seq up > ../logs/log_$seq"."_$o.txt";
      	system( $command );
      	$manager->finish;
}