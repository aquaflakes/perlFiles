#!/usr/bin/perl -w
# 2015-11-09  kmer counter and converted to frequency with pseudocount, preparation for freq_cmp.pl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Mylib1::Parallel::file;
#use File::Spec qw();
#use Data::Dumper;
#use POSIX;

my $this_programme = basename($0);
my $usage = "\n\$$this_programme -f '*.txt' -o kmer_cnt -k 6 -m -t 10

  ## general ##
   -f 'str' input reads files, can be 'xxx/*.txt' or 'xxx/*'
   -o       output folder
   -h       get help
   -t int   number of threads to run (1 thread per file)
  ## specific ## 
   -k int   size of the kmer to analyze. Default 6
   -m       will count all possible kmer per sequences
            Default: only one kmer is counted per sequence entries
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob,'o=s' => \my $outputPath, 'h' => \my $helpInfo, 'k=i' => \my $k_length, 'm' => \my $multiCntperSeq, 't=i' => \my $threads); # i:int   f:float   s:string
    # initialize
    if ($helpInfo||(!$fileGlob)||(!$k_length)) { print "\nusage: $usage\n";exit;}
    $outputPath||($outputPath="./$this_programme"); $outputPath.="_".$k_length."mer"; system("rm -r $outputPath") if (-d $outputPath);
    $threads||($threads=6); 
    #$pseudoCnt||($pseudoCnt=50); 


#### get all files to process ######
my @files = glob($fileGlob); chomp @files;

&file::mThread($threads,\&singleThread,\@files);

sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
  my $file = shift;
  my %kmer_p; my $totalCnt=0;
  print "Counting the number of $k_length nt long kmers in $file\n";
  SeqAnal($file,\%kmer_p,\$totalCnt);
  calcfreq(\%kmer_p,\$totalCnt);
  my $capLine= "kmer\tcount\tfreq\t"."totalCounts=$totalCnt\t"; ######
  printHits(\%kmer_p,$file,$capLine);
  return(0)
}

sub SeqAnal {
  my ($file,$kmer_p,$totalCntRef) = @_;

  open FH, "$file" || die "Can't open file $file\n";
  my $seq;
  while(<FH>){
   chomp;
   next if "";
   countKmer(\$_,$kmer_p,$totalCntRef);
  }
  close FH;
  return(0);
}

sub countKmer {
  my ($seq_p, $kmer_p, $totalCntRef)=@_;
  my $k = $k_length;
  my %beenThere;

  for (my $i=0;$i <= length(${$seq_p})-$k;$i++){
    my $w = substr(${$seq_p},$i,$k);
    unless ($multiCntperSeq){
      #Count only one occurrence of a kmer per sequence
      $kmer_p->{$w}->{cnt}++ && $$totalCntRef++ if !exists $beenThere{$w};
      $beenThere{$w}=1;
    }else{
      #Count all instances of a kmer per sequence
      $kmer_p->{$w}->{cnt}++; $$totalCntRef++;
    }
  }
  return(0);
}

sub calcfreq
{
  my $kmer_p = shift; my @keys=(keys $kmer_p);
  my $totalCntRef = shift; #my $totalwithPseudo= $$totalCntRef+$pseudoCnt*$#keys;
  foreach (@keys)
  {
    $kmer_p->{$_}->{freq} = ($kmer_p->{$_}->{cnt})/$$totalCntRef;
    #$kmer_p->{$_}->{freqPseudo} = ($kmer_p->{$_}->{cnt}+$pseudoCnt)/$totalwithPseudo;
  }
  return(0)
}

sub printHits {
  my ($kmer_p,$file,$capLine)=@_;

  ##print out the hits
  my ($dir,$pre,$suf) = ($file =~ /(^.+\/|^)(.+)\.(.+$)/);
  `mkdir -p $outputPath`;
  open (OUT, ">","$outputPath/".$pre."_".$k_length."mer.txt") || die "Can't create file $pre.hits\n";
  print OUT $capLine."\n";
  print OUT join("\t", $_, $kmer_p->{$_}->{cnt}, $kmer_p->{$_}->{freq}), "\n" for (reverse sort {$kmer_p->{$a}->{freq} <=> $kmer_p->{$b}->{freq}} keys $kmer_p); # sort by freq
  close OUT;

  return(0);
}



######  multithread start #######
#sub multithreads #only file list
#{
#  my $filelist=shift; my @filelist=@{$filelist};
#  use threads;
#  use Thread::Queue;
# 
#  my $q = Thread::Queue->new;
#  $q->enqueue(@filelist); #needs a file list
#
#  my $num_workers = @filelist < $threads ? @filelist : $threads; # no need to be wasteful :)
#
#  for (1 .. $num_workers) {
#    threads->new(\&worker, $q);
#  }
#
#  $_->join for threads->list;
#}
#
#sub worker
#{
#  my $queue = shift;
#  while (my $filename = $queue->dequeue_nb) {
#    singleThread($filename);
#  }
#  return(0);
#}
#######  multithread end #######
