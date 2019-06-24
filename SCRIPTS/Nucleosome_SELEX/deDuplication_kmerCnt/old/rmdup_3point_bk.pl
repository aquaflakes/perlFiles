#!/usr/bin/perl -w
# remove duplicated reads and closely related reads (grep -v 20mer in duplicated reads)
use warnings;
use strict;
use Getopt::Long;
#use File::Spec qw();
use Data::Dumper;
use File::Basename;
use POSIX;

my $this_programme = basename($0);
my $usage = "\n\$$this_programme -f '*.txt' -k 20 -o dedup20

   -f 'str' input reads files, can be 'xxx/*.txt' or 'xxx/*'
   -o       output folder
   -k int   use kmer window on every duplicated reads to scan and de-dup
   -h       get help
   -l int   head n lines (for test), blank for all lines
   -t int   number of threads to run (1 thread per file)
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob,'o=s' => \my $outputPath, 'h' => \my $helpInfo, 'k=i' => \my $win_size, 'l=i' => \my $lines, 't=i' => \my $threads); # i:int   f:float   s:string
    if ($helpInfo||(!$fileGlob)) { print "\nusage: $usage\n";exit;}
    $outputPath||($outputPath="./de_duped"); $lines||($lines=100000000000000000); $threads||($threads=6); # initialize
    (-d $outputPath)? system("rm -r $outputPath"):1;


#### get all files to process ######
my @files = glob($fileGlob); chomp @files;

&multithreads(\@files);

sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
    my $file = shift;
    open(my $reads, "<", $file) or die "non exist file $file";
    my %allReadsCnt;
    my $linecnt=0;
    while (<$reads>) # first remove identical duplicates
    {
        chomp $_;
        $allReadsCnt{$_}++;
        $linecnt++;
        $linecnt>$lines? last:1 ;
    }
    close $reads;
    my @de_duped; my %existed;    # try to remove sequences with close hamming distances
    
    #check length of lig
    open(my $readsT, "<", $file) or die "non exist file $file";
    my $testseq=<$readsT>;  close $readsT; chomp $testseq;
    my $liglength= length($testseq);
    print "\nliglength for $file is $liglength bp\n";
    ($win_size>$liglength)? ($win_size=$liglength):1;
    # my $segmentNum=ceil($liglength / $win_size);
    my $segmentNum=$liglength - $win_size +1;
    
    LoopThroughReads:
    foreach my $seq (keys %allReadsCnt)
    {
        # extract substrings of (the key to examine) to check if defined in %existed
        my @patterns_curr_seq;
        for (0..$segmentNum-1) {push @patterns_curr_seq, substr($seq,$_,$win_size);} # extract all $win_size-mer from current seq    
     
        foreach (@patterns_curr_seq)
        {
            if (defined($existed{$_})) {next LoopThroughReads; for (0,int($segmentNum/2),$segmentNum-1) {$existed{substr($seq,$_,$win_size)}++;} }#only add the head, middle, end win_size mer into %existed to reduce memory usage
            #$existed{$_}++;# add all win_size mers into %existed, better use >25 mer
        }
        for (0,int($segmentNum/2),$segmentNum-1) {$existed{substr($seq,$_,$win_size)}++;} #only add the head, middle, end win_size mer into %existed to reduce memory usage
        
        push @de_duped, $seq;
    }
       
    if (!(-d $outputPath)) {system("mkdir -p $outputPath");}
    
    my $outputPath_fileName=($outputPath)."/".basename($file);
    open(my $out, ">", $outputPath_fileName) or die "cannot wirte output file for $outputPath_fileName";
    print $out join("\n", @de_duped);
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
