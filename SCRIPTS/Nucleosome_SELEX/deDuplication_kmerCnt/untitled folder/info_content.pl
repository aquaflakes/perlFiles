#!/usr/bin/perl -w
# 2015-11-11  give log (base2) freq change of k-mers between files
use warnings;
use strict;
use Getopt::Long;
use File::Basename;

my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme -s '*/*_bound_cyc*_*mer.txt' -p 20
           
find c0 files automatically in the same folder

  ## general ##
   -s 'str' input signal files, can be 'xxx/*.txt' or 'xxx/*'
   -r 'str' input ref files, can be './xxx/{file1.txt, file2.txt, ...}' or './xxx/file{1..4}.txt' or 'xxx/*.txt' or 'xxx/*'
   -o       output folder
   -h       get help
   -t int   number of threads to run (1 thread per file)
  ## specific ## 
   -p int   add freq according to pseudo count to denominator (freq in ref file), to prevent fluctuations for kmers with low counts
            Default:20
   ";
	##### options
	GetOptions ('s=s' => \my $fileGlob_sig,'r=s' => \my $fileGlob_ref,'o=s' => \my $outputPath, 'h' => \my $helpInfo, 't=i' => \my $threads, 'p=i' => \my $pseudoCnt); # i:int   f:float   s:string
    # initialize
    if ($helpInfo||(!$fileGlob_sig)) { print "\nusage: $usage\n";exit;}
    $threads||($threads=10); 
    $pseudoCnt||($pseudoCnt=20); 
    $outputPath||($outputPath="./$this_programme"."_pseudo$pseudoCnt");
    
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd


#### get all files to process ######
my @files_sig = glob($fileGlob_sig);chomp @files_sig; my @files_sig_c0 = map s/cyc./cyc0/ir, @files_sig; # specify c0 files for bk corr
#my @files_ref = glob($fileGlob_ref);chomp @files_ref; my @files_ref_c0 = map s/cyc./cyc0/ir, @files_ref;
if ($#files_sig!=$#files_sig_c0) { print "\n c0 file not found for all reads \n usage: $usage\n";exit;} # ensure identical number of sig and ref files
my @params= map [$files_sig[$_],$files_sig_c0[$_]], 0..$#files_sig; # pack parameters into units for each thread

use Mylib1::Parallel::file; file::mThread($threads,\&singleThread,\@params);

sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
  my $param_unit = shift;
  my $file_sig=$param_unit->[0]; my $file_sig_c0=$param_unit->[1]; 
  
  print "\nsig $file_sig / $file_sig_c0 is calced for info content with pseudo count $pseudoCnt\n";
  my %kmerFreq;
  my ($totalCntSig,$not_used)=readFromFile(\%kmerFreq,$file_sig,"sigFreq","sigCnt");
  my ($totalCntC0,$pseudoFreq)=readFromFile(\%kmerFreq,$file_sig_c0,"sigFreqC0","sigCntC0");

  
  my $infoContent=0;
  
  foreach (keys %kmerFreq)
  {
    (delete $kmerFreq{$_} && next) if $_=~m/N/; #delete N containing k mers 
    my $sigFreq=\$kmerFreq{$_}->{sigFreq}; $$sigFreq||($$sigFreq=0); # fill undef with 0s
    #my $refFreq=\$kmerFreq{$_}->{refFreq}; $$refFreq||($$refFreq=0);
    
    my $sigFreqC0=\$kmerFreq{$_}->{sigFreqC0}; $$sigFreqC0||($$sigFreqC0=0); 
    #my $refFreqC0=\$kmerFreq{$_}->{refFreqC0}; $$refFreqC0||($$refFreqC0=0); 

    my $multiply_factor= 1; #(($$sigCnt<10)&&($$refCnt<10))? 5:1;# increase pseudo by 5 when both counts are low
      $infoContent+= $$sigFreq*(log(($$sigFreq-$$sigFreqC0)/($$sigFreqC0+$pseudoFreq*$multiply_factor)+1)/log(2));
  }
  
  my $capLine=$file_sig."\t".$infoContent;
  &prnResult(\%kmerFreq,$file_sig,$capLine);
  return(0);
}

sub readFromFile
{
  my ($href,$file,$sigOrRefFreq,$sigOrRefCnt)=@_;
  open(my $FH, "<", $file) or die "cannot open file handle for read";
  my $capLine=<$FH>; chomp $capLine;
  (split /\t/, $capLine)[3]=~ m/.*=(.*)/; my $totalCnt= $1; my $pseudoFreq=$pseudoCnt/$totalCnt;
  while (<$FH>)
  {
    chomp;
    my ($kmer,$cnt,$freq)=split /\t/;
    $href->{$kmer}->{$sigOrRefFreq}=$freq;
    $href->{$kmer}->{$sigOrRefCnt}=$cnt;
  }
  close $FH;
  return ($totalCnt,$pseudoFreq);
}

sub prnResult {
  my ($khref,$file_sig,$capLine)=@_;

  my ($dir_sig,$pre_sig,$suf_sig) = ($file_sig =~ /(^.+\/|^)(.+)\.(.+$)/);
  #my ($dir_ref,$pre_ref,$suf_ref) = ($file_ref =~ /(^.+\/|^)(.+)\.(.+$)/);
  `mkdir -p $outputPath`;
  open (my $OUT, ">","$outputPath/$pre_sig"."_infocnt".".txt") || die "Can't create file for $pre_sig \n";
  print $OUT $capLine."\n";
  close $OUT;
}


