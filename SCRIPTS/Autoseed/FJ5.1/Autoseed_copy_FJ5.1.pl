#!/usr/bin/perl -w
##perl programme #1
#used to cut long reads into 40bp shorter ones for autoseed analysis
#v2 for FJ5.1KSEL with yimeng ligand, c1-3

use FindBin;               # locate this script
use lib "$FindBin::Bin/../../.."; 
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use autodie;
use perllib::Me::Nu_cnt::checkfile;
use POSIX qw(ceil strftime);
use File::Spec qw(abs2rel);
use Parallel::ForkManager;
#use Cwd; # get current dir
# use File::Temp qw(tempdir);
		
######### global variables to adjust #############
my $wkPATH= "/wrk/data/fangjie/"; # set the rootpath of working

######### time and usage
my $tiimaalku =localtime(); my $start_time=time();
my $this_programme = basename($0);
my $usage = "\n\$$this_programme [JustReads file | folder]\n
	Processes all files undr JustRead subfolders
	Example: $this_programme -g guidefileName
	         Autoseed_copy.pl -g FJ5.1KSELt_r1_Fangj201510.txt

	***options***
	-g  specify the guide file to copy and run for Autoseed
	";


########## Info from cmd line #################
	##### options
	GetOptions ('l=i' => \my $lengthDef, 'guidegen' => \my $guideGen,'test' => \my $isTest, 'copy' => \my $toCopy, 'c=i' => \my $cores, 'g=s' => \my $guideName,
				'jr' => \my $isJustReadDir, 'b=s' => \my $bk_cyc,);#, 'c=i' => \$cores, 'n=s' => \$appended_foldername, 'p=f' => \$p_value ); # i:int   f:float   s:string
	

########## copy for autoseed ##########
my $output_path= File::Spec->rel2abs(".");
$guideName=~ m{([^/]*).txt};
my  $guideFileNameOnly= $1;

	system("cp -a $output_path/*.seq /var/www/kaz/data/sequence");
	if ($guideName)
	{
		system("cp -a $guideName /var/www/kaz/guide");
		chdir("/var/www/kaz");
		system("nohup php asks_gen.php $guideFileNameOnly 8 thread=10 order=0123 >$output_path/nohup_ask_gen.txt&"); # "&" to run at background
		# cyc like 4b0 will not run unless delete "all"
		
		#if (defined($bk_cyc))
		#{	# cyc like 4b0 will not run unless delete "all"
		#	system("nohup php asks_gen.php $guideFileNameOnly kmerspace 8 thread=10 >$output_path/nohup_ask_gen_bk.txt&"); 
		#}
		print("nohup log to $output_path/nohup_ask_gen.txt\n");
	}



#------------- fixed ending of program 
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
exit;







