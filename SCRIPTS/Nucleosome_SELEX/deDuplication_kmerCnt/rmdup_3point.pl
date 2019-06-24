#!/usr/bin/perl -w
# remove duplicated reads and closely related reads (grep -v 20mer in duplicated reads)
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Mylib1::Parallel::file;
use POSIX;

my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme -f '*.txt' -k 20 -o 3_dedup20 --gz

   -f 'str' input reads files, can be 'xxx/*.txt' or 'xxx/*'
   -o       output folder
   -k int   use kmer window on every duplicated reads to scan and de-dup
   -h       get help
   -l int   head n lines (for test), blank for all lines
   -t int   number of threads to run (1 thread per file)
   --gz     output gz file
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob,'o=s' => \my $outputPath, 'h' => \my $helpInfo, 'gz' => \my $isgz, 'k=i' => \my $win_size, 'l=i' => \my $lines, 't=i' => \my $threads); # i:int   f:float   s:string
    if ($helpInfo||(!$fileGlob)) { print "\nusage: $usage\n";exit;}
    $outputPath||="./de_duped";  $threads||=10; # initialize
    
    use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    #qx(rm -r $outputPath) if (-d $outputPath);
    qx(mkdir -p $outputPath);
    qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd
    

#### get all files to process ######
my @files = glob($fileGlob); chomp @files;

&file::mProcess($threads,\&singleProc,\@files); # using mThread necessary for debug !!

sub singleProc    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
    my $file = shift;
    print "processing $file \n";
    my @allreads; use Mylib1::seqFile; &seqFile::getallSeq(\@allreads,$file,$lines);
    &seqFile::rmdup_3point(\@allreads,$win_size);
    &seqFile::writeFile(\@allreads,$outputPath, $file, "_u", $isgz, "autoDetect_gz");
}

