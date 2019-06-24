#!/usr/bin/perl -w
# 2015-11-11  give log (base2) freq change of k-mers between files
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
#use File::Spec qw();
#use Data::Dumper;
#use POSIX;

my $this_programme = basename($0);
my $usage = "\n\$$this_programme -s '*/*_bound_cyc*_*mer.txt' -r '*/*_unbound_cyc*_*mer.txt' -p 20
background corrected with c0
formula:   logFoldChnSig= (sigFreq - sigFreqC0) / (sigFreqC0 + pseudoFreq) + 1;
           logFoldChnRef= (refFreq - refFreqC0) / (refFreqC0 + pseudoFreq) + 1;
           logFoldChn= (logFoldChnSig / logFoldChnRef) / log(2);
           
NB!! *should be identical numbers of sigFiles and refFiles. 1st sig accords to 1st ref, 2st accords to 2st ...
     *pseudofreq multiplied by 5 if both sigCnt and refCnt < 10, kmer contains 'N' deleted

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
    if ($helpInfo||(!$fileGlob_sig)||(!$fileGlob_ref)) { print "\nusage: $usage\n";exit;}
    $threads||($threads=10); 
    $pseudoCnt||($pseudoCnt=20); 
    $outputPath||($outputPath="./$this_programme"."_pseudo$pseudoCnt"); #system("rm -r $outputPath") if (-d $outputPath);



#### get all files to process ######
my @files_sig = glob($fileGlob_sig);chomp @files_sig; my @files_sig_c0 = map s/cyc./cyc0/ir, @files_sig; # specify c0 files for bk corr
my @files_ref = glob($fileGlob_ref);chomp @files_ref; my @files_ref_c0 = map s/cyc./cyc0/ir, @files_ref;
if ($#files_sig!=$#files_ref) { print "\nusage: $usage\n";exit;} # ensure identical number of sig and ref files
my @params= map [$files_sig[$_],$files_ref[$_],$files_sig_c0[$_],$files_ref_c0[$_]], 0..$#files_sig; # pack parameters into units for each thread

&multithreads(\@params);

sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
  my $param_unit = shift;
  my $file_sig=$param_unit->[0]; my $file_ref=$param_unit->[1]; my $file_sig_c0=$param_unit->[2]; my $file_ref_c0=$param_unit->[3];
  
  print "\nsig $file_sig / $file_sig_c0 is compared to ref $file_ref / $file_ref_c0 with pseudo count $pseudoCnt\n";
  my %kmerFreq;
  my ($totalCntSig,$not_used)=readFromFile(\%kmerFreq,$file_sig,"sigFreq","sigCnt");
  my ($totalCntRef,$pseudoFreq)=readFromFile(\%kmerFreq,$file_ref,"refFreq","refCnt");
  readFromFile(\%kmerFreq,$file_sig_c0,"sigFreqC0","sigCntC0");
  readFromFile(\%kmerFreq,$file_ref_c0,"refFreqC0","refCntC0");
  
  foreach (keys %kmerFreq)
  {
    (delete $kmerFreq{$_} && next) if $_=~m/N/; #delete N containing k mers 
    my $sigFreq=\$kmerFreq{$_}->{sigFreq}; $$sigFreq||($$sigFreq=0); # fill undef with 0s
    my $refFreq=\$kmerFreq{$_}->{refFreq}; $$refFreq||($$refFreq=0);
    my $sigCnt=\$kmerFreq{$_}->{sigCnt}; $$sigCnt||($$sigCnt=0);
    my $refCnt=\$kmerFreq{$_}->{refCnt}; $$refCnt||($$refCnt=0);
    
    my $sigFreqC0=\$kmerFreq{$_}->{sigFreqC0}; $$sigFreqC0||($$sigFreqC0=0); 
    my $refFreqC0=\$kmerFreq{$_}->{refFreqC0}; $$refFreqC0||($$refFreqC0=0); 

    my $multiply_factor= 1; #(($$sigCnt<10)&&($$refCnt<10))? 5:1;# increase pseudo by 5 when both counts are low
    $kmerFreq{$_}->{logFoldChnSig}= ($$sigFreq-$$sigFreqC0)/($$sigFreqC0+$pseudoFreq*$multiply_factor)+1;
    $kmerFreq{$_}->{logFoldChnRef}= ($$refFreq-$$refFreqC0)/($$refFreqC0+$pseudoFreq*$multiply_factor)+1;
    $kmerFreq{$_}->{logFoldChn}= log($kmerFreq{$_}->{logFoldChnSig}/$kmerFreq{$_}->{logFoldChnRef})/log(2);
  }
  
  my $capLine="kmer\tlogFoldChn\tsigCnt\trefCnt\tsigFreq\trefFreq\t"."totalCntSig=$totalCntSig\t"."totalCntRef=$totalCntRef\t"."pseudoCnt=$pseudoCnt\t"."pseudoFreq=$pseudoFreq";
  &prnResult(\%kmerFreq,$file_sig,$file_ref,$capLine);
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
  my ($khref,$file_sig,$file_ref,$capLine)=@_;

  my ($dir_sig,$pre_sig,$suf_sig) = ($file_sig =~ /(^.+\/|^)(.+)\.(.+$)/);
  my ($dir_ref,$pre_ref,$suf_ref) = ($file_ref =~ /(^.+\/|^)(.+)\.(.+$)/);
  `mkdir -p $outputPath`;
  open (my $OUT, ">","$outputPath/$pre_sig"."_cmp_"."$pre_ref.txt") || die "Can't create file for $pre_sig cmped to $pre_ref\n";
  print $OUT $capLine."\n";
  
  print $OUT join("\t", $_, $khref->{$_}->{logFoldChn}, $khref->{$_}->{sigCnt}, $khref->{$_}->{refCnt}, $khref->{$_}->{sigFreq}, $khref->{$_}->{refFreq}), "\n"
  for (reverse sort {$khref->{$a}->{logFoldChn} <=> $khref->{$b}->{logFoldChn}} keys $khref); # sort by fold change
  
  close $OUT;
}



######  multithread start #######
sub multithreads
{
  my $filelist=shift; my @filelist=@{$filelist};
  use threads;
  use Thread::Queue;
 
  my $q = Thread::Queue->new;
  $q->enqueue(@filelist); #needs a file list

  my $num_workers = @filelist < $threads ? @filelist : $threads; # no need to be wasteful :)

  for (1 .. $num_workers) {
    threads->new(\&worker, $q);
  }

  $_->join for threads->list;
}

sub worker
{
  my $queue = shift;
  while (my $filename = $queue->dequeue_nb) {
    singleThread($filename);
  }
  return(0);
}
######  multithread end #######
