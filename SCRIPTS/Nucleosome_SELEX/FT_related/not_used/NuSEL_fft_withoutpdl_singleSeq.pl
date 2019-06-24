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
use perllib::Me::FT::fft_win;
use perllib::FileIO::FileIO;

 #initialize stuff
##################

my $this_programme = basename($0);
my $usage = "\$$this_programme <seq reads file | all> <row number | all> --bkcrr [core numbers(opt)] --NameAdd";
if(!defined($ARGV[0]) or !defined($ARGV[1])){ print "$usage\n\n"; exit;}

######### variables specified from cmd #############
my $file_param = $ARGV[0];
my $row_param = $ARGV[1];
my $cores=$ARGV[3]||10;
	my $pm = Parallel::ForkManager->new($cores);
my $bkcrr=$ARGV[2];
my $finalNameAdd=$ARGV[4]||"";

######### variables to adjust #############
my $tmpFile_path="./analysis/$this_programme"."_tmp";                 my $tmpFile_path_ori=$tmpFile_path;
my $final_output_path="./analysis/$this_programme._$finalNameAdd/";

my $reads_anal_methodRef = \&ADM::addSR_simple_align; # process reads and store result in hash tmp, return a hash, write after all reads finished
my $tmpHash_output_methodRef = \&ADM::wtADM; # write hash in the process of each core into a tmp file
my $Summed_hash_output_methodRef = \&ADM::WriteSummedADM;

################################################

my $tiimaalku =localtime(); my $start_time=time();
print "Program started $tiimaalku\n";

if (-d "$final_output_path") {system ("rm -r $final_output_path");}

my @files;
if ($file_param eq "all")
	{@files = <*.txt>; } #all files in the folder into the array
	else
	{push (@files, $file_param);}

#my $loopcnt=0; #number of files processed

FileLoop:
foreach my $filename ( @files )
{   #loop through all files
		#my $temp_cnt=$loopcnt+1;
		my $liglength = &checkfile::formatCheck_and_liglength($filename); #($filename, 0 || 1(check if only fastq file required))
		unless ($liglength) {print "!! skip \n"; next FileLoop;}
		
		$pm->start and next FileLoop;
		
		#print "\n########### file $temp_cnt ##############\n";
		print "length of the reads set to $liglength bp for $filename \n";
		my $filename_no_ext=$filename; $filename_no_ext =~ tr/\.txt//d;
		
		open(my $seq,"<",$filename) or die "cannot open sequence file $filename";
		my %ntCount;
		my %peakHeight; # peak height of 10.2bp
		my @keys= qw (A T C G AA AT AC AG CA CT CC CG GA GT GC GG TA TT TC TG); #(AA GG TA);
		
		$row_param = $ARGV[1];
		while (my $seq = <$seq>)
		{
			chomp $seq;			
			my $Nufreq= \my %Nufreq; for (@keys) {$Nufreq->{$_}=[split (",",("0," x $liglength))]};
			
			for (@keys) {  while ($seq=~ m/(?=$_)/g) {$Nufreq{$_}->[$-[0]]++;}  }
			# since by default search begins from the end of last pos, using ?= matches the 0 length position before homodinucleotides
			# calculate the NT occurance in each position, $-[0] is startPos of the match
						
			my ($peakHeightHR,$ntCountHR) = &fft_win::ftSRfreq($Nufreq);
			for (@keys) {$peakHeight{$_}+=$peakHeightHR->{$_}; $ntCount{$_}+=$ntCountHR->{$_};}
			$row_param--; $row_param||last;
			
		}
		close $seq;
		
		my %normalizedPeak= map {$_, $peakHeight{$_}/$ntCount{$_}."\t".$peakHeight{$_}."\t"."$ntCount{$_}"} @keys;
		my $addedCap;
		if ($bkcrr)
		{
				my @cyc0_bkfiles = glob("/wrk/data/fangjie/seqFiles/cyc0/analysis/Nusig147_fft_singleSeq.pl/*.txt");
				my $target_bkfile; $filename_no_ext =~ m/_([0-9]*[A-Z]*)[0-9A-Z]*_/;  # m/_(147[0-9A-Z]*)_/;
				my $match_pattern= $1;
				foreach my $bkADM (@cyc0_bkfiles)
				{ if ($bkADM =~ m/\Q$match_pattern/) {$target_bkfile=$bkADM;last;} }
				print ">>>bk file<<<   ".$target_bkfile."\n";
				if ($target_bkfile)
				{		
				my $nt_freq_bk= &FileIO::readTable_1_lvHash($target_bkfile);
				%normalizedPeak= map{$_,($peakHeight{$_}/$ntCount{$_}-$nt_freq_bk->{$_}->[0])."\t".$normalizedPeak{$_}} @keys;
				$addedCap=["key","Peak10.2crr","normalizedPeak10,2H","peakHeight","ntAbsCounts"];
				}
				#print Dumper ($nt_freq_bk);die;

		}
		
		
		@keys=sort @keys;
		&FileIO::writeTable_1_LvHR("$final_output_path/$filename_no_ext"."_SRft_.txt",\%normalizedPeak,\@keys,$addedCap||["key","normalizedPeak10,2H","peakHeight","ntAbsCounts" ]);#outFilepath data_HR, subkeys_to_print(ordered)_AR, capline_AR
		#for (@keys) {print $_."\t".$peakHeight{$_}/$ntCount{$_}."\n"} ;

	#-------------------	
	
		print $filename." processing finished";
		print "\n";
		#print "############ file $temp_cnt #############\n\n";
		
		$pm->finish;
		#$loopcnt++;	
}
$pm->wait_all_children;


print scalar(@files)."files processed in total\n";
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
if (-d "$tmpFile_path_ori") {system ("rm -r $tmpFile_path_ori");} # if the temp folder exist delete it

########## plot ##########

chdir $final_output_path;
system ('R --no-save </wrk/data/fangjie/lib/SCRIPTS/Nucleosome_SELEX/FT_related/Nusel_singleread_ft.R');
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




	


