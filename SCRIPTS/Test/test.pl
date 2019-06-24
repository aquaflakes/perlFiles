#!/usr/bin/perl -w
#20150814

use strict;
use warnings;

use Mylib1::seqFile;

my @allreads;
#&seqFile::getSeqbyLine(\@allreads,$ARGV[0],$ARGV[1]);
my $liglength=&seqFile::liglength($ARGV[0]);
use Data::Dumper;
#print Dumper(@allreads);
#print "\n".($#allreads+1)."\n";
print "$liglength\n";