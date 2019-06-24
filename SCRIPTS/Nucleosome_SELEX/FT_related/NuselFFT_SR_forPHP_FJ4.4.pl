#!/usr/bin/perl -w
#2016-02-04  support of .gz file
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
use Getopt::Long;



#*** INIT ***
    my $this_programme = basename($0);
    my $cmdInput=$this_programme." ".join(" ", @ARGV);
    my $usage = "\n\$$this_programme -f '../all_ln/FJ6.1_OriMe1/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-*IIIc5*' -t 30 -l 100000 -o ./SR_FFT/MeOricyc5/ -k 2 -c 5 -g /var/www/kaz/guide/FJ6.1_OriMe1_Fangj201605.txt -phpdir /var/www/kaz/zhu_data/FJ6.1_OriMe1_Fangj201605/SRFFT
or \$$this_programme -f '../all_ln/FJ6.1_OriMe1/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-*IIIc5*' -t 30 -l 100000 -o ./SR_FFT/MeOricyc5/ -k 2


  ## general ##
   -f 'str' read files, META CHAR OK
   -o       output folder
   -h       help
   -t int   total thread No. (1 /file)
   
  ## specific ##

   -l int   proc only first -l lines /file
   -k int	specify kmer length 
   
   -g str	specify guide file to edit
   -c int	specify the cycle
   -phpdir str  specify the folder to link output images to (for php display)
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob, 'o=s' => \my $outputPath, 'g=s' => \my $guide, 'phpdir=s' => \my $phpdir,'c=s' => \my $cyc,  'h' => \my $helpInfo, 'l=i' => \my $lines, 'k=i' => \my $k,'t=i' => \our $threads, 'st=i' => \my $sub_threads ); # i:int   f:float   s:string
    if ($helpInfo||(!$fileGlob)) { print "\nusage: $usage\n";exit;}
    $outputPath||="./analysis/$this_programme/";
    $lines||=100000000000000000;
    $threads||=2; $sub_threads||=10; 
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath);
    qx(echo ">> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd.txt); # memo input cmd
    my @files = glob($fileGlob); chomp @files;

# -----------

# ----interface----
	my $pm = Parallel::ForkManager->new($threads);
	my $final_output_path=$outputPath;


################################################

my $tiimaalku =localtime(); my $start_time=time();
print "Program started $tiimaalku\n";

FileLoop:
foreach my $filename ( @files )
{   		
		$pm->start and next FileLoop;
		
		#print "\n########### file $temp_cnt ##############\n";
		my $filename_no_ext=$filename; $filename_no_ext =~ tr/\.txt//d;
		
		my @allreads;
		use Mylib1::seqFile qw(); &seqFile::getallSeq(\@allreads,$filename);
		@allreads=splice(@allreads,0,$lines) if @allreads>$lines; # take specified amount of lines
		my $liglength =length($allreads[0]);
		print "\nligand length for $filename is $liglength\n";
				
		use Mylib1::seqFile; my $Aref = &seqFile::kmerArrGen($k);
		my @keys= ($k==2 ? sort  (qw (A T C G AA AT AC AG GA GT GC GG CA CT CC CG TA TT TC TG)) : @$Aref ); 
		my @FFTsum;
		$liglength--;
		foreach my $seq (@allreads)
		{		
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
		}
		
		# generate the according axis points
		my @FFTaxis_gen=(0) x ($liglength); my $FFTaxis_gen = pdl @FFTaxis_gen; my $axisPoints;
		my $N=$liglength; # points of data for FFT
		my $D=1; # step distance between data points
		
		# determine freq of each point according to odd or even, if odd=1
		@FFTaxis_gen & 1 ? ($axisPoints = $FFTaxis_gen->xlinvals(-($N/2-0.5)/$N/$D,($N/2-0.5)/$N/$D)->rotate(-($N-1)/2)) : ($axisPoints = $FFTaxis_gen->xlinvals(-($N/2-1)/$N/$D,1/2/$D)->rotate(-($N/2 -1))); 
		
		use Mylib1::fileIO; my $filename_out=&fileIO::getfile($filename);
		system ("mkdir -p $final_output_path/FFT_SR/");
		PDL::wcols $axisPoints, @FFTsum, "$final_output_path/FFT_SR/$filename_out"."fft.txt",{ Header => "keys\t".join ("\t",@keys), Colsep => "\t" };
		
	
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

########## plot ##########

chdir $final_output_path;
#plot individual well, link for display
system ('Rscript /wrk/data/fangjie/lib/Rlib/FJ4.4/plotFFT_forPHP.R \'./FFT_SR/*.txt\' ./forPHP/');
qx(mkdir -p $phpdir) if (!(-d $phpdir));
use File::Spec qw(rel2abs); my $link_source= File::Spec->rel2abs("./forPHP");
qx(ln -sf $link_source $phpdir/SRFFTc$cyc);


#calc peak intensity of individual well and write guide
system ("Rscript /wrk/data/fangjie/lib/Rlib/FJ4.4/FFTSR_peakH_writeGuide.R 'FFT_SR/*.txt' $guide c$cyc");
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




	


