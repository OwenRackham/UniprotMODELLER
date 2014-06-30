#!/usr/bin/perl -w

use strict;
use warnings;
use Parallel::ForkManager;

my $filename = $ARGV[0];
my $species = $ARGV[1];
my $LOGFILE = "/tmp/";
open FILE, "<$filename" or die $!;
my @seqs;
while (<FILE>){
	if ($_ =~ /^>(\S+)/){
		push @seqs,$1;
	}
}

my $manager = new Parallel::ForkManager( 18 );
my $t = scalar(@seqs);
my $o = 1;
foreach my $seq (@seqs){
	print "$o from $t\n";
	$o++;
		#unless (-e  "$LOGFILE"."up_r2_$seq"."_$o.out") {
		$manager->start and next;
		my $command = "perl DO_model.pl $seq up $species > $LOGFILE"."up_r2_$seq"."_$o.out";
      	system( $command );
      	$manager->finish;
		#}
}
