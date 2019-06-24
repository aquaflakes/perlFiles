#!/usr/bin/perl -w
# using new version of relative packages  //2016-01-17
# new version of NuselFT_align_bk_winfunc_p.pl, support multifile and multiprocess/file 

use warnings;
use strict;
use Getopt::Long;
#use File::Spec qw();
use Data::Dumper;
use File::Basename;
use POSIX;

#custom lib
use Mylib1::Parallel::file;
use Mylib1::Hash_Merge;
use Mylib1::Parallel::line;
use Me::Nu_cnt::ADM;
use Me::FT::fft_win;
    
#*** INIT ***
    my $this_programme = basename($0);
    my $cmdInput=$this_programme." ".join(" ", @ARGV);
    my $usage = "\n\$$this_programme -f '*.txt' -b '*.txt' -o AlignFT/ -k 2 -t 15 -st 1 --samebk -phpdir /var/www/kaz/zhu_data/FJ6.1_OriMe1_Fangj201605/AlignFT


  ## general ##
   -f 'str' read files, META CHAR OK
   -b 'str'	bk files, META CHAR OK
   --samebk the same bk file applies to all files
   -o       output folder (opt)
   -h       help
   -t int   total thread No. (1 /file)
   
  ## specific ##
   -st int  sub_thread No. /file
   -l int   proc only first -l lines /file
   -k int	specify kmer length
   --bkADM    bk file is ADM
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob, 'b=s' => \my $bkfileGlob,'o=s' => \my $outputPath, 'phpdir=s' => \my $phpdir, 'h' => \my $helpInfo,'bkADM' => \my $bkisADM, 'samebk' => \my $samebk, 'l=i' => \my $lines, 'k=i' => \my $k,'t=i' => \our $threads, 'st=i' => \my $sub_threads, ); # i:int   f:float   s:string
    if ($helpInfo||(!$fileGlob)||(!$k)) { print "\nusage: $usage\n";exit;}
    $outputPath||="./$this_programme/";
    $lines||=100000000000000000;
    $threads||=2; $sub_threads||=10;
	
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(rm -r $outputPath) if (-d $outputPath);
	qx(mkdir -p $outputPath); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd
	
    my @files = glob($fileGlob); chomp @files;
	my @bkfiles = glob($bkfileGlob||0); chomp @bkfiles;
	
	my $sumbk_s, my $temp_s, my $lineCntbk_s;
	if ($samebk)
	{
		($sumbk_s, $temp_s, $lineCntbk_s)=&countFile($bkfiles[0]);
		@bkfiles= map 1, @files;
	}
	
	my @param; for (my $i=0;$i<@files;$i++) { push @param, [$files[$i],$bkfiles[$i]||0]; } # make 1-1 pairs

# -----------


&file::mProcess($threads,\&singleProc,\@param); # using mThread necessary for debug !!
#&file::mThread($threads,\&singleProc,\@param);

chdir $outputPath;

#plot individual well, link for display
system ('Rscript /wrk/data/fangjie/lib/Rlib/FJ4.4/R_NuselFT_align_PHP.R ./ ./forPHP/');
if ($phpdir) {
	qx(mkdir -p $phpdir) if (!(-d $phpdir));
	use File::Spec qw(rel2abs); my $link_source= File::Spec->rel2abs("./forPHP");
	qx(ln -sf $link_source $phpdir); # link to PHP dir for disp
}
exit;

	
# *** multi proc ***
sub singleProc  # each file
{
    my $param_Aref = shift;
	my ($file,$bkfile)= @$param_Aref;
	my ($sum,$liglength,$lineCnt)= &countFile($file); &outputAbscnt($sum,$file,$outputPath,$lineCnt);
	if ($bkfile)
	{
		my ($sumbk,$lineCntbk);
		if ($bkisADM) {
			$sumbk= rdADM_asIs($bkfile);
			foreach (keys $sumbk) {$lineCntbk+=$sumbk->{$_}->[0];} }
		else{
				if($samebk){($sumbk,my $temp,$lineCntbk)=($sumbk_s, $temp_s, $lineCntbk_s);} # do not have to calc if bkfile are the same for all files
				else{ ($sumbk,my $temp,$lineCntbk)= &countFile($bkfile)};
			}
		
		print "liglength in bk file not same as sig" && return if @{$sum->{(keys $sum)[0]}} != @{$sumbk->{(keys $sum)[0]}} ; #check if $liglengthbk != $liglength
		my $bkAmp= &extractBkAmp($sumbk,$lineCnt,$lineCntbk);
		$sum= &Hash_Merge::minus($sum,$bkAmp); # $sum still abs amplitude(counts)
	}
	&FFTandOutput($sum,$file,$outputPath,$lineCnt);
}

sub singleSubProc	# an aliquot of reads
{
    my ($reads_Aref,$param_Aref)=@_;
	my ($subProc,$liglength,$k)= @$param_Aref;	my %subProc= %$subProc;
	foreach my $seq (@$reads_Aref)
	{
		map {my $key= substr($seq, $_, $k); (defined $subProc{$key}) && $subProc{$key}->[$_]++;} (0..($liglength-$k));
	}
	return \%subProc;
}
# ----------------


#-----------------  modules ------------------
sub kmerHashGen
{
	my ($k,$liglength)=@_; #my $k_ori=$k;
	my @kmers= ("A","C","G","T");
	while ($k-1>0)
	{
		@kmers= map {$_."A", $_."C", $_."G", $_."T"} @kmers;
		$k--;
	}
	#push (@kmers, ("A", "T", "C", "G")) if ($k_ori==2);  # calc also single nt when dint
	
	my %kmers= map {$_, [split //, (0 x ($liglength-length($_)+1))]} @kmers;
	return \%kmers;
}

sub countFile
{
	my $file = shift;
	my @allreads;
	#open(my $reads, "<", $file) or die "non exist file $file";
    #my @allreads=<$reads>;  close $reads;   chomp @allreads;
	use Mylib1::seqFile qw(); &seqFile::getallSeq(\@allreads,$file);
	
	
    @allreads=splice(@allreads,0,$lines) if @allreads>$lines; # take specified amount of lines
	
	my $lineCnt= scalar(@allreads);
    my $liglength=length($allreads[0]); my %subProc= %{&kmerHashGen($k,$liglength)};
    print "\nligand length for $file is $liglength\n";
    
	my $temp_Data_Aref= &line::mProcess_l($sub_threads,\&singleSubProc,\@allreads,[\%subProc,$liglength,$k]); # upper lv using mThread necessary for debug
	my $sum= \my %sum;
	($sum= &Hash_Merge::add($sum,$_)) for @$temp_Data_Aref;
	return ($sum,$liglength,$lineCnt);
}

sub mean
{
    my $sum = eval join("+",@_);
	return $sum/@_;
}

sub extractBkAmp 
{
	my ($Href,$lineCnt,$lineCntbk)=@_;
	foreach my $keynt (keys($Href))
    {
		my $meanOfKey= &mean(@{$Href->{$keynt}});
        for (my $i = 0; $i < scalar(@{$Href->{$keynt}}); ++$i)
		{
			$Href->{$keynt}->[$i] = ($Href->{$keynt}->[$i]-$meanOfKey) * ($lineCnt/$lineCntbk);
		}
	}
	return $Href;
}

sub FFTandOutput
{
	my ($Summed_ADM_af_bkCorr,$filename_final_anal,$final_output_path,$lineCnt)=@_;
	use Mylib1::fileIO; $filename_final_anal=&fileIO::getfile($filename_final_anal);
	$filename_final_anal=~s/\.[^\.]*$/\.txt/i; # replace with .txt for subsequent R analysis
	
	
	my $Summed_ADM_freq = $Summed_ADM_af_bkCorr; 
	my ($fftPower,$peakPercent,$phaseAngle,$peakArea) = &fft_win::ftADM($Summed_ADM_freq);
	&ADM::wtADM_asIs_capLine($fftPower,$filename_final_anal,$final_output_path."/FT","freq");#final param is key for cap row
	&ADM::wtADM_asIs_capLine($peakPercent,$filename_final_anal,$final_output_path."/FTS","keys");
	&ADM::wtADM_asIs_capLine($phaseAngle,$filename_final_anal,$final_output_path."/phaseAng","keys");
	&ADM::wtADM_asIs_capLine($peakArea,$filename_final_anal,$final_output_path."/Area","keys");# used to plot in R
	&ADM::wtADM_asIs($Summed_ADM_freq,$filename_final_anal,$final_output_path."/ADM_after_corr_not_freq");
	
	&normalize($Summed_ADM_freq); #  minus average and divide by average !!!!!!!!! to relative amplitude here by ref
	$Summed_ADM_freq->{pos}= \my @emptyArr; push $Summed_ADM_freq->{pos},(1..scalar(@{$Summed_ADM_freq->{(keys $Summed_ADM_freq)[0]}}));
	&ADM::wtADM_asIs_capLine($Summed_ADM_freq,$filename_final_anal,$final_output_path."/ADM","pos");#,"--prefill","--fillLong"); #fill diNt with 0 in the end bit 
}

sub outputAbscnt
{
	my ($Summed_ADM_freq,$filename_final_anal,$final_output_path,$lineCnt)=@_;
	use Mylib1::fileIO; $filename_final_anal=&fileIO::getfile($filename_final_anal);
	$filename_final_anal=~s/\.[^\.]*$/\.txt/i; # replace with .txt for subsequent R analysis
	&ADM::wtADM_asIs($Summed_ADM_freq,$filename_final_anal,$final_output_path."/ADM_absCnt");
}

sub normalize #(hash without cap line)
{
	my $ntFreqHash=shift;
	foreach my $keynt (keys($ntFreqHash))
    {
		my $meanOfKey= &mean(@{$ntFreqHash->{$keynt}});
		next if $meanOfKey==0;
        for (my $i = 0; $i < scalar(@{$ntFreqHash->{$keynt}}); ++$i)
		{
			$ntFreqHash->{$keynt}->[$i] = ($ntFreqHash->{$keynt}->[$i]-$meanOfKey)/$meanOfKey;
		}
	}
	return $ntFreqHash;
}

sub rdADM_asIs #(input_file_name, ligand_length) # no triming
{
	my $admfile=shift; #the adm file name
	my @data;
	open(DATA, "$admfile") || die "sub rdADMv1: Can't open nt_frequency file $admfile: $!\n";
	while (<DATA>) {
		chomp;
	  push @data, [split /\t/];
	}
	close DATA;
	my %admFreq;
	$admFreq{shift(@$_)}=$_ foreach (@data);
	#my $linecnt=0; foreach (@data) {$linecnt+=$_->[1];}
	return \%admFreq; 
}