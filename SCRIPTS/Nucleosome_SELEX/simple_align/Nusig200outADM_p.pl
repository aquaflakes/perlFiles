#!/usr/bin/perl -w
#BEGIN
#{
#	push (@INC, "/wrk/data/fangjie/lib");
#	push (@INC, "/Volumes/nutcase.biosci.ki.se -aea No num/lib/"); #for local debug
#}
use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../..";  # use this directory
use warnings;
use strict;
use perllib::Me::Nu_cnt::ADM;
use File::Basename;
use perllib::Me::parallel;
use perllib::Me::Nu_cnt::Nusig200outADM;
#use Data::Dumper;


 #initialize stuff
##################
my $title = "<ligand nt frequency calc of 200bp full>";
my $this_programme = basename($0);
my $usage = "\$$this_programme <seq reads file | all> <row number | all> <file containing std adm | d> <core numbers> ";
print "$title\n";
if(!defined($ARGV[0]) or !defined($ARGV[1]) or !defined($ARGV[2])){ print "$usage\n\n"; exit;}

######### variables specified from cmd #############
my $file_param = $ARGV[0];
my $row_param = $ARGV[1];
my $cores=$ARGV[3]||10;
my $filterfile = $ARGV[2]; if (($filterfile eq "d") or ($filterfile eq "default") ) {$filterfile="/wrk/data/fangjie/DefaultFilesforProc/adm_FAC_14704GGTC94N_4_nt_counts.txt";}

######### variables to adjust #############
my $tmpFile_path="./analysis/$this_programme"."_tmp";                 my $tmpFile_path_ori=$tmpFile_path;
my $final_output_path="./analysis/$this_programme"; #."_".$cores."c";

my $reads_anal_methodRef = \&Nusig200outADM::Nusig200outADM; # process reads and store result in hash tmp, write after all reads finished
my $tmpHash_output_methodRef = \&ADM::wtADM_200outADM; # write hash in the process of each core into a tmp file
my $Summed_hash_output_methodRef = \&ADM::WriteSummedADM;

my $filterHashRef;
my $filter_file_process_methodRef = \&ADM::rdADMv1;
my $tmpreads_output_methodRef;  # process reads and write to tmp file immediately
my $Summed_seq_output_methodRef;
my $isFastq =0; # 1 if input file is fastq format, 0 if is JustReads
################################################


my $tiimaalku =localtime(); my $start_time=time();
print "Program started $tiimaalku\n";


my @files;
if ($file_param eq "all")
	{@files = <*200*>; } #all files in the folder into the array
	else
	{push (@files, $file_param);}

my $loopcnt=0; #number of files processed

FileLoop:
foreach my $filename ( @files )
{   #loop through all files
		my $temp_cnt=$loopcnt+1;
		my $liglength = &checkfile::formatCheck_and_liglength($filename,$isFastq); #($filename, 0 || 1(check if only fastq file required))
		unless ($liglength) {print "!! skip \n"; next FileLoop;}
		$filterHashRef=$filter_file_process_methodRef->($filterfile,$liglength);
		
		print "\n########### file $temp_cnt ##############\n";
		print "length of the reads set to $liglength bp for $filename \n";
	
		my $filename_no_ext=$filename; $filename_no_ext =~ tr/\.txt//d; 
		my $tmpFile_path =$tmpFile_path."/$filename_no_ext"; #modify temp file path for each file
		
		########################################
		my $filename_final_anal= $filename_no_ext."_out_ADM.txt"; #		<<<<<<<<<<< set final file name >>>>>>>>>>>>
		my $filename_final_seq= $filename_no_ext."_seq.txt"; #		<<<<<<<<<<< set final file name >>>>>>>>>>>>
		########################################
														
	my ($tmpFile_path_r, $tmp_hash_files_arrRef, $temp_seq_files_arrRef) =
    &parallel::calc($filename, $liglength, $row_param, $cores, $tmpFile_path, $isFastq,      $reads_anal_methodRef, $tmpHash_output_methodRef,      $tmpreads_output_methodRef, $filterHashRef ); #process for the files <<<<<<<<<<<<<<<
	
	############# after added the number of N ###########
	$liglength+=53;
	#####################################################
		
	if (defined($tmpHash_output_methodRef)) {$Summed_hash_output_methodRef->($tmpFile_path_r, $tmp_hash_files_arrRef, $liglength, $final_output_path, $filename_final_anal);}
	if (defined($tmpreads_output_methodRef)) {$Summed_seq_output_methodRef->($tmpFile_path_r, $tmp_hash_files_arrRef, $liglength, $final_output_path, $filename_final_anal);}
	
		print $filename." processing finished";
		print "\n";
		print "############ file $temp_cnt #############\n\n";
		$loopcnt++;	
}

print "$loopcnt files processed in total\n";
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
if (-d "$tmpFile_path_ori") {system ("rm -r $tmpFile_path_ori");} # if the temp folder exist delete it


########## plot ##########
my $ADM_output_path=$final_output_path."/ADM";
#print "$ADM_output_path \n";
chdir $ADM_output_path;
system ('R --no-save </wrk/data/fangjie/lib/Rlib/Nusel_align147.R');

exit;







	


