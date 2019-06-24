#!/usr/bin/perl -w
#20150814 edited from Nusig147_freq
use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../..";  # use this directory
use warnings;
use strict;
use perllib::Me::Nu_cnt::ADM;
use File::Basename;
use perllib::Me::parallel;
use Data::Dumper;
use List::Util qw(sum);
use perllib::Me::FT::fft_win_phase_test;
use Parallel::ForkManager;

 #initialize stuff
##################
my $title = "<ligand nt frequency calc>";
my $this_programme = basename($0);
my $usage = "\$$this_programme <seq reads file | all> <row number | all> <core numbers(opt)>";
print "$title\n";
if(!defined($ARGV[0]) or !defined($ARGV[1])){ print "$usage\n\n"; exit;}

######### variables specified from cmd #############
my $file_param = $ARGV[0];
my $row_param = $ARGV[1];
my $cores=$ARGV[2]||10;
my $filterfile;

######### variables to adjust #############
my $tmpFile_path="./analysis/$this_programme"."_tmp";                 my $tmpFile_path_ori=$tmpFile_path;
my $final_output_path="./analysis/$this_programme";

my $reads_anal_methodRef = \&ADM::addSR_simple_align; # process reads and store result in hash tmp, return a hash, write after all reads finished
my $tmpHash_output_methodRef = \&ADM::wtADM; # write hash in the process of each core into a tmp file
my $Summed_hash_output_methodRef = \&ADM::WriteSummedADM;

my $filterHashRef;
my $filter_file_process_methodRef;
my $tmpreads_output_methodRef;  # process reads and write to tmp file immediately
my $Summed_seq_output_methodRef;
my $isFastq =0; # 1 if input file is fastq format, 0 if is JustReads
################################################

my $tiimaalku =localtime(); my $start_time=time();
print "Program started $tiimaalku\n";

if (-d "$final_output_path") {system ("rm -r $final_output_path");}

my @files;
if ($file_param eq "all")
	{@files = <*.txt>; } #all files in the folder into the array
	else
	{push (@files, $file_param);}

my $loopcnt=0; #number of files processed

FileLoop:
foreach my $filename ( @files )
{   #loop through all files
		my $temp_cnt=$loopcnt+1;
		my $liglength = &checkfile::formatCheck_and_liglength($filename,$isFastq); #($filename, 0 || 1(check if only fastq file required))
		#$liglength =51;
		unless ($liglength) {print "!! skip \n"; next FileLoop;}
		
		print "\n########### file $temp_cnt ##############\n";
		print "length of the reads set to $liglength bp for $filename \n";
	
		my $filename_no_ext=$filename; $filename_no_ext =~ tr/\.txt//d; 
		my $tmpFile_path =$tmpFile_path."/$filename_no_ext"; #modify temp file path for each file
		
		########################################
		my $filename_final_anal= $filename_no_ext."_ADM.txt"; #		<<<<<<<<<<< set final file name >>>>>>>>>>>>
		my $filename_final_seq= $filename_no_ext."_seq.txt"; #		<<<<<<<<<<< set final file name >>>>>>>>>>>>
		########################################
														
	my ($tmpFile_path_r, $tmp_hash_files_arrRef, $temp_seq_files_arrRef) =
    &parallel::calc($filename, $liglength, $row_param, $cores, $tmpFile_path, $isFastq,      $reads_anal_methodRef, $tmpHash_output_methodRef,      $tmpreads_output_methodRef, $filterHashRef ); #process for the files <<<<<<<<<<<<<<<

	if (defined($tmpHash_output_methodRef)) {$Summed_hash_output_methodRef->($tmpFile_path_r, $tmp_hash_files_arrRef, $liglength, $final_output_path, $filename_final_anal);}
	if (defined($tmpreads_output_methodRef)) {$Summed_seq_output_methodRef->($tmpFile_path_r, $tmp_hash_files_arrRef, $liglength, $final_output_path, $filename_final_anal);}
	
	#-----bk correction---
		#---read in summed ADM just generated--- 
	my $Summed_ADM_path= $final_output_path."/ADM/".$filename_final_anal;
	my $Summed_ADM_bf_bkCorr= &ADM::rdADMv1_notlg($Summed_ADM_path); system ("rm $Summed_ADM_path");
	my $Summed_ADM_bf_bkCorr_lineCnt = $Summed_ADM_bf_bkCorr->{"A"}->[4]+$Summed_ADM_bf_bkCorr->{"T"}->[4]+
						$Summed_ADM_bf_bkCorr->{"C"}->[4]+$Summed_ADM_bf_bkCorr->{"G"}->[4]; #use the sum of ATCG of the 5th position as the line count 
	
		#---read in background ADM----
		"#########"=~ m/(.*)/; #reset $1
	my @cyc0_bkfiles = </wrk/data/fangjie/seqFiles/cyc0/analysis/Nusig147_p.pl/ADM/*>;
	my $target_bkfile;
	#my $target_bkfile="/wrk/data/fangjie/seqFiles/cyc0/analysis/Nusig147_p.pl/ADM/_14704GGTC94N_0_ADM.txt"; 
	$filename_no_ext =~ m/_(147[0-9]*[A-Z]*)[0-9A-Z]*_/;  # m/_(147[0-9A-Z]*)_/;
	my $match_pattern= $1;
	foreach my $bkADM (@cyc0_bkfiles)
		{
				if ($bkADM =~ m/\Q$match_pattern/) {$target_bkfile=$bkADM;last;}
								
		}
		
	#-------substract------
	my $Summed_ADM_af_bkCorr;
	if ($target_bkfile) #when bk exists correct for bk
	{
		print ">>>bk file<<<   ".$target_bkfile."\n";
		my $bck_ADM_c0= &ADM::rdADM_asIs($target_bkfile);
		$bck_ADM_c0= &ADM::trim_multint($bck_ADM_c0);	
		$Summed_ADM_bf_bkCorr= &ADM::trim_multint($Summed_ADM_bf_bkCorr);
		
		my $bck_ADM_c0_lineCnt = $bck_ADM_c0->{"A"}->[4]+$bck_ADM_c0->{"T"}->[4]+
				$bck_ADM_c0->{"C"}->[4]+$bck_ADM_c0->{"G"}->[4]; #use the sum of ATCG of the 5th position as the line count 
				foreach my $keynt (keys(%{$Summed_ADM_bf_bkCorr}))
			{	
				my $meanOfKey= &mean(@{$bck_ADM_c0->{$keynt}}); # substract mean in correction, retaining information of counts in the resulting ADM
				for (my $i = 0; $i < scalar(@{$Summed_ADM_bf_bkCorr->{$keynt}}); ++$i)
				{
					#no warnings; 
					$Summed_ADM_af_bkCorr->{$keynt}->[$i] = $Summed_ADM_bf_bkCorr->{$keynt}->[$i] - ($bck_ADM_c0->{$keynt}->[$i]-$meanOfKey) * ($Summed_ADM_bf_bkCorr_lineCnt/$bck_ADM_c0_lineCnt);
				}
			}
	
	}else{print "no matching bk file found \n";$Summed_ADM_af_bkCorr=$Summed_ADM_bf_bkCorr;}
		
	my $Summed_ADM_freq = &ADM::cntToFreq($Summed_ADM_af_bkCorr); # convert counts to frequency
	my ($fftPower,$peakPercent,$phaseAngle,$peakArea) = &fft_win_phase_test::ftADM($Summed_ADM_freq);
	&ADM::wtADM_asIs_capLine($fftPower,$filename_final_anal,$final_output_path."/FT","freq");#final param is key for cap row
	&ADM::wtADM_asIs_capLine($peakPercent,$filename_final_anal,$final_output_path."/FTS","keys");
	&ADM::wtADM_asIs_capLine($phaseAngle,$filename_final_anal,$final_output_path."/phaseAng","keys");
	&ADM::wtADM_asIs_capLine($peakArea,$filename_final_anal,$final_output_path."/Area","keys");
	
	my $Summed_ADM_freq_nomalized=&normalize($Summed_ADM_freq); # devided by average and minus average
	$Summed_ADM_freq->{pos}= \my @emptyArr; push $Summed_ADM_freq->{pos},(1..scalar(@{$Summed_ADM_freq->{A}}));
	&ADM::wtADM_asIs_capLine($Summed_ADM_freq,$filename_final_anal,$final_output_path."/ADM","pos","--prefill","--fillLong"); #fill diNt with 0 in the end bit 


	#-----------bk corr end--------	
	
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
my $ADM_output_path=$final_output_path;#."/ADM";
#print "$ADM_output_path \n";
chdir $ADM_output_path;
system ('R --no-save </wrk/data/fangjie/lib/SCRIPTS/Nucleosome_SELEX/FT_related/R_NuselFT_align_stdphase.R');
#print "!!!!!!!!!!!!! a few bp trimmed in each end for FT analysis !!!!!!!!!!!!\n";
exit;


#------------------- func used---------
sub mean
{
    return sum(@_)/@_;
}

sub normalize #(hash without cap line)
{
	my $ntFreqHash=shift;
	foreach my $keynt (keys($ntFreqHash))
    {
		my $meanOfKey= &mean(@{$ntFreqHash->{$keynt}});
        for (my $i = 0; $i < scalar(@{$ntFreqHash->{$keynt}}); ++$i)
		{
			$ntFreqHash->{$keynt}->[$i] = ($ntFreqHash->{$keynt}->[$i]-$meanOfKey)/$meanOfKey;
		}
	}
	return $ntFreqHash;
}




	


