#!/usr/bin/perl -w
##perl programme #2
#used to cut long reads into 40bp shorter ones for autoseed analysis

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
my $usage = "\n\$$this_programme [folder] (--jr -c i -b s)
             Autoseed_rename.pl  . -pc -pb

   add Batch&Cyc to all txt file under subfolders named JustRead
   forlder name of each cyc should include cyc0-cyc4

   
   --jr for folders of JustRead files but without 'Justread' in name
   -c	specify cyc number for files without cyc No.
   -b	specify batch name for files without cyc No.
   -pc	preserve original cyc number of file
   -pb  preserve original batch name
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

#my $pm = Parallel::ForkManager->new($cores||6); #PPPPPP
#Process:#PPPPPP
foreach my $JustReadsFile (@JustReads_files)
{   
   	#print "\n########### file $temp_cnt ##############\n";
	#$pm->start and next Process;#PPPPPPPP
	
	##### get abs path for each JR file
	my $JustReadsFile = File::Spec->rel2abs($JustReadsFile); 
	##### construct the fileName to be renamed
	my $JRfile_newName;
	
	$JustReadsFile =~ m/cyc(\d)/i; my $cyc = $1; #extract cyc number
	if (defined($forcedCyc)) {$cyc= $forcedCyc;}
	
	$JustReadsFile=~ m/\/(?<batch>[^\/]*)_[^\/]*_(?<cyc>.?)\.txt$/;
	if ($batchPresrv && $+{batch})
	{$batchName=$+{batch};}
	if ($cycPresrv && $+{cyc})
	{$cyc=$+{cyc};}	
	$JRfile_newName= $JustReadsFile =~s/(\/)[^\/]*(_[^\/]*_)\d*(\.txt$)/$1$batchName$2$cyc$3/ir; # -r: not modifying but just return the modified value
	
	

	system ("mv $JustReadsFile $JRfile_newName");
	
	# creat links to all justread files under outputDir
	my ($volume,$pathOfDir,$JRfileName_noPath) = File::Spec->splitpath($JRfile_newName);
	my $JRLinkName = $output_path.$JRfileName_noPath;
	if (!(-d "$output_path")) {system ("mkdir -p $output_path");}
	system ("ln -fs $JRfile_newName $JRLinkName");


	
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







