#!/usr/bin/perl -w
# to parallel commands

use warnings;
use strict;
use Mylib1::Parallel::file;

#*** INIT ***
    use File::Basename; my $this_programme = basename($0);
    my $cmdInput=$this_programme." ".join(" ", @ARGV);
    my $usage = "\n\$$this_programme -f '*.gz' -o ../Trimmed_Reads/ -t 10 -cmd 'trim_galore -q 10 --illumina --stringency 4 --length 55 --gzip --suppress_warn -o ../Trimmed_Reads/ SingleFiles'

  ## general ##
   -f 'str' all files to be processed, META CHAR OK
   -t int   total thread No. (1 /file)
   
  ## specific ##
   -cmd 'str'   the command to be issued for each file, 'SingleFile' serve as a placeholder to be replaced
   ";
	##### options
	use Getopt::Long;
	GetOptions ('f=s' => \my $fileGlob, 'o=s' => \my $outputPath, 't=i' => \our $threads, 'cmd=s' => \my $cmd, ); # i:int   f:float   s:string
    if (!$fileGlob) { print "\nusage: $usage\n";exit;}
    $threads||=10;
	
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd
	
	
    my @files = glob($fileGlob); chomp @files;
	my @param; for (my $i=0;$i<@files;$i++) { push @param, $files[$i]; } 

# -----------

&file::mProcess($threads,\&singleProc,\@param); # using mThread necessary for debug !!
#&file::mThread($threads,\&singleProc,\@param);

sub singleProc
{
	my $file=shift;
	$cmd=~ s/SingleFiles/$file/ig;
	print "processing file $file \n";
	system ("$cmd"); 	
}








