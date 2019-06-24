#!/usr/bin/perl -w
#2016-05-03
#used seed a guide file automatically using yimeng's guide "Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper"

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;


		
######### global variables to adjust #############

######### time and usage ###
my $this_programme = basename($0);
my $usage = "\n\$$this_programme -g xx.txt [-s seedfile --m] \n
	seed a guide file automatically using yimeng's guide 'Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper'
	
	***options***
	-g     target guide file
	-s     file with stdSeed ('/var/www/kaz/guide/Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper.txt' as default)
	--m    using methylated seeds (default using non methylated seeds)
";


########## Info from cmd line #################
	##### options
	GetOptions ('g=s' => \my $guideFile, 's=s' => \my $seedFile, 'm' => \my $methylated);
	unless ($guideFile) {print $usage; exit;}
	$seedFile||="/var/www/kaz/guide/Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper.txt";

# backup the guide first
my $bk_guide_dir="/wrk/data/fangjie/DefaultFilesforProc/guide_bk_seedStd/";
use Mylib1::fileIO; my $name_only=fileIO::getfile($guideFile);
`cp -a $guideFile $bk_guide_dir$name_only`;	
my @stdSeeds= `cat $seedFile`; chomp @stdSeeds;

	# grep only methylated or unmethylated
unless (defined $methylated)
{	@stdSeeds= map {my @currLine=split /\t/, $_; $currLine[0] =~ m/GROUP/ || $currLine[7] =~ m/methylated/i? ():\@currLine; } @stdSeeds; }#returning null list means delete
else
{	@stdSeeds= map {my @currLine=split /\t/, $_; $currLine[0] =~ m/GROUP/ || $currLine[7] !~ m/methylated/i? ():\@currLine; } @stdSeeds;}

	#convert into hash
my %stdSeeds; 
	#substitute if the defined entry is dimer motif or if the seed is longer than current
map { defined($stdSeeds{uc $$_[0]}) && $stdSeeds{uc $$_[0]}->{isdimer}!~m/dimer/? ():($stdSeeds{uc $$_[0]} = {FL=>$$_[1],barcode=>$$_[2],batch=>$$_[3],seed=>$$_[4],multinomial=>$$_[5],cyc=>$$_[6],family=>$$_[7],isdimer=>$$_[8]}); } @stdSeeds;

open(my $in, "<", $guideFile) or die; my $currline=1; my ($seedcol,$TFcol);
open(my $out, ">", $guideFile."_test.txt") or die;
while (<$in>)
{
	chomp; $_=~tr/\r//d;
	if ($currline==1) {
		use List::MoreUtils; my @title=split /\t/,$_; $seedcol= List::MoreUtils::first_index {$_ eq "Inplantedseed"} @title; $TFcol= List::MoreUtils::first_index {$_ eq "TF"} @title;
		print $out $_."\n";
	}
	else{
		my @currLine=split /\t/, $_; $currLine[$seedcol]=$stdSeeds{uc $currLine[$TFcol]}->{seed}||"AA";
		print $out join "\t", @currLine; print $out "\n";		
	}	
	$currline++;
}
close $in;
close $out;
#`cp -a $guideFile* /wrk/data/fangjie`;
system ("mv $guideFile"."_test.txt $guideFile");
`chmod 666 $guideFile`;

exit;