#!/usr/bin/perl -w
#20150824
use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../..";  # use this directory
use warnings;
use strict;
use perllib::Me::Nu_cnt::ADM();
use perllib::Me::Nu_cnt::checkfile();
use File::Basename;
use Parallel::ForkManager;
#use Data::Dumper;
use List::Util qw(sum);
#use perllib::Me::FT::fft_win;
use perllib::FileIO::FileIO();
use PDL;
use PDL::FFT qw/fft/;
use POSIX qw(dup2);
use FileHandle ();



 #initialize stuff
##################
my $this_programme = basename($0);
my $usage = "\$$this_programme <seq reads file | all> <row number | all> --bkcrr(0) [core numbers(opt)] --NameAdd
no bk correction available yet";
if(!defined($ARGV[0]) or !defined($ARGV[1])){ print "$usage\n\n"; exit;}

######### variables specified from cmd #############
my $file_param = $ARGV[0];
my $row_param = $ARGV[1];
my $cores=$ARGV[3]||10;
	my $pm = Parallel::ForkManager->new($cores);
my $bkcrr=$ARGV[2]||0;
my $finalNameAdd=$ARGV[4]||"";

######### variables to adjust #############
my $tmpFile_path="./analysis/$this_programme"."_tmp";                 my $tmpFile_path_ori=$tmpFile_path;
my $final_output_path="./analysis/$this_programme"."$finalNameAdd/";

my $reads_anal_methodRef = \&ADM::addSR_simple_align; # process reads and store result in hash tmp, return a hash, write after all reads finished
my $tmpHash_output_methodRef = \&ADM::wtADM; # write hash in the process of each core into a tmp file
my $Summed_hash_output_methodRef = \&ADM::WriteSummedADM;

################################################

my $tiimaalku =localtime(); my $start_time=time();
print "Program started $tiimaalku\n";

if (-d "$final_output_path") {system ("rm -r $final_output_path");}
system ("mkdir -p $final_output_path");

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
		# my $lineCnt; #used for normalization of peak height
		my @keys= sort (qw (A T C G AA AT AC AG GA GT GC GG CA CT CC CG TA TT TC TG)); 
		my @FFTsum;
		$liglength--;
		$row_param = &checkfile::SetRowNum ($filename,$row_param);
		while (my $seq = <$seq>)
		{
			chomp $seq;			
			for (0..(@keys-1))
			{
				my @FFTtemp_r=(0) x ($liglength); my @FFTtemp_i=(0) x ($liglength); #fill the blank position with "0"s
				while ($seq=~ m/$keys[$_]/g) {$FFTtemp_r[$-[0]]++;} #not counting neiboring dinucleotides to avoid baseline distortion,otherwise use (?=$keys[$_])
				splice @FFTtemp_r, $liglength, 1; # delete the final element when using single nucleotied as key, because output wcol requires all $FFT[$_] having the same length
				my $FFTtemp_r= pdl @FFTtemp_r; my $FFTtemp_i= pdl @FFTtemp_i;
				fft($FFTtemp_r,$FFTtemp_i);
				$FFTsum[$_]+= sqrt($FFTtemp_r**2+$FFTtemp_i**2);
				#$FFTsum[$_]/=PDL::average($FFTsum[$_]); # normalize against average (performed in R now)
			}
			$row_param--; $row_param||last;		
		}
		close $seq;
		
		# generate the according axis points
		my @FFTaxis_gen=(0) x ($liglength); my $FFTaxis_gen = pdl @FFTaxis_gen; my $axisPoints;
		my $N=$liglength; # points of data for FFT
		my $D=1; # step distance between data points
		@FFTaxis_gen & 1 ? ($axisPoints = $FFTaxis_gen->xlinvals(-($N/2-0.5)/$N/$D,($N/2-0.5)/$N/$D)->rotate(-($N-1)/2)) : ($axisPoints = $FFTaxis_gen->xlinvals(-($N/2-1)/$N/$D,1/2/$D)->rotate(-($N/2 -1))); 
		
		PDL::wcols $axisPoints, @FFTsum, "$final_output_path/$filename_no_ext"."pdlfft.txt",{ Header => "keys\t".join ("\t",@keys), Colsep => "\t" };
		
	
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
system ('R --no-save </wrk/data/fangjie/lib/SCRIPTS/Nucleosome_SELEX/FT_related/R_NuselFFT_SingleRead.R');
print "\n\n 93 bit of the ligand used in FFT for both dinucleotides and single nucleotides!!!!!!!!!!!!!!\n\n";
exit;


#------------------- func used---------
sub mean
{
    return main::sum(@_)/@_;
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




	


