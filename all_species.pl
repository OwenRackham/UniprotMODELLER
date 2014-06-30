#!/usr/bin/perl -w

use strict;
use warnings;

print "bakers:1/10 \n";
my $cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/bakersyeast.fa bakersyeast";
system($cmd);
print "fission:2/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/fissionyeast.fa fissionyeast";
system($cmd);
print "human:3/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/human.fa human";
system($cmd);
print "plas:4/10 \n";
my $cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/plasmodium.fa plasmodium";
system($cmd);
print "worm:5/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/worm.fa worm";
system($cmd);
print "arab:6/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/arabidopsis.fa arabidopsis";
system($cmd);
print "ecoli:7/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/ecoli.fa ecoli";
system($cmd);
print "fly:8/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/fly.fa fly";
system($cmd);
print "mouse:9/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/mouse.fa mouse";
system($cmd);
print "staph:10/10 \n";
$cmd = "perl Fork_DO_model.pl /home/rackham/workspace/data/r2genomes/Genome3D/staphylococcus.fa staphylococcus";
system($cmd);

