#!/usr/bin/perl -w
##perl programme #2
# rename and collect all read files in subforlders called JustRead as links
#2015-10-26

use FindBin;               # locate this script
use lib "$FindBin::Bin/../.."; 
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use autodie;
#use perllib::Me::Nu_cnt::checkfile;
#use POSIX qw(ceil strftime);
use File::Spec qw(abs2rel);
#use Parallel::ForkManager;
#use Cwd; # get current dir
# use File::Temp qw(tempdir);
		
######### global variables to adjust #############
my $wkPATH= "/wrk/data/fangjie/"; # set the rootpath of working

######### time and usage
my $tiimaalku =localtime(); my $start_time=time();
my $this_programme = basename($0);
my $usage = "\n\$$this_programme [folder] (--jr -b s)
             $this_programme  . (-b s)

   modify batch name and to all txt file under subfolders named JustRead
 
   
   --jr for folders containing read files but without 'Justread' in name
   -b	specify batch name, preserve original if no batchName provided.

   ";

########## Info from cmd line #################
	if (!$ARGV[0]) { print "usage: $usage\n\n";exit;}
	##### options
	GetOptions ('jr' => \my $isJustReadDir,'pc' => \my $cycPresrv, 'pb' => \my $batchPresrv,'c=i' => \my $forcedCyc, 'b=s' => \my $batchName); # i:int   f:float   s:string

	##### remaining param
	
	my $specified_path_absolute = File::Spec->rel2abs($ARGV[0]);
	##### extract the name list of reads files to be processed
	my @JustReads_files;
		# if the current path is "JustRead then take the files included"
		if(($specified_path_absolute =~ m/justread/i)||$isJustReadDir) {@JustReads_files = <$ARGV[0]/*.txt>;}  
		# if current path is not JustRead then search for all JustRead subfolder
		else{ 
				my @JustReads_file_folders =`find $specified_path_absolute -iname "*Justread*" -type d 2>/dev/null | grep -v "copied"`;
				foreach (@JustReads_file_folders)
				{
					chomp; # eliminate the \n
					my @filestmp=glob("$_/*.txt");
					push(@JustReads_files,@filestmp);
				}
			}
	##### construct output path (extract only the path info from $ARGV[0])
	my $output_pathAdd = "/allReadsFiles_ln/";
		#folder name of the output folder under the current
	my $output_path;    if ( -f $specified_path_absolute ) {$specified_path_absolute =~ m/^(.*\/)/;  my $matched = $1 || "."; $output_path = $matched.$output_pathAdd;} else { $output_path = $specified_path_absolute.$output_pathAdd;}
		# if it is a file, get its path or use ".", if input is directory than use it for concatanation directly
	if (-d "$output_path") {system ("rm -r $output_path");}
		# if the output folder already exists delete it

######### Loop through all JustReads files gen cut files
my $temp_cnt=1;


foreach my $JustReadsFile (@JustReads_files)
{   

	##### get abs path for each JR file
	my $JustReadsFile = File::Spec->rel2abs($JustReadsFile); 
	##### construct the fileName to be renamed
	my $JRfile_newName;
		
	$JustReadsFile=~ m/\/(?<batch>[^\/]*)_(?<addinfo>[^\/]*?)(?<barcode>[A-Z]*[\d]*N[A-Z]*)_(?<cyc>.?)\.txt$/i;
	$batchName || ($batchName = $+{batch});
	$JRfile_newName=$batchName.$+{addinfo}."_".$+{barcode}."_".$+{cyc}.".txt";
	#$JustReadsFile =~s/(\/)[^\/]*(_[^\/]*_)\d*(\.txt$)/$1$batchName$2$cyc$3/ir; # -r: not modifying but just return the modified value
	#system ("mv $JustReadsFile $JRfile_newName");  # rename the original file under JustRead folder
	
	# create links to all justread files under outputDir
	if (!(-d "$output_path")) {system ("mkdir -p $output_path");}
	system ("ln -fs $JustReadsFile $output_path/$JRfile_newName");


	
	##### extract info form the name of JustReadfile			

    #------------- fixed message for ending a JustReadfile proc
    #print $JustReadsFile." processing finished";
	#print "\n";
	#print "############ file $temp_cnt #############\n\n";
	#$pm->finish;#PPPPPP
    $temp_cnt++;
}
#$pm->wait_all_children;	#PPPPPP

#------------- fixed ending of program 
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
exit;







