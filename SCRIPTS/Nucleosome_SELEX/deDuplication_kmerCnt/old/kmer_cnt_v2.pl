#!/usr/bin/perl -w
# 2015-11-09
# kmer counter and converted to frequency,
# may count the reverse strand and bin kmer with its rc
# deduplication function included

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Mylib1::Parallel::file;


my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme -f '*.txt' -o kmer_cnt -k 6 --m --u (--ori --uOut --gz)

  ## general ##
   -f 'str' input reads files, can use wildcards like 'xxx/*.txt' or 'xxx/*'
   -o       output folder
   -h       get help
   -t int   number of threads to run (1 thread per file)
  ## specific ## 
   -k int   size of the kmer to analyze. Default 6
   -w int   window size for deduplication, default 20
   --m      will count all possible kmer per sequences
            Default: only one kmer is counted per sequence entries
   --rc     count also the other strand
   --bin    bin the counts of k mer with its rc, x2 for palindromic k mer
   --u      also count after removing duplicate ligand
   --ori    also output the original read No. of the k mer
   --uOut   output uniq reads
   --gz     output gz file
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob,'o=s' => \my $outputPath, 'gz' => \my $isgz,'h' => \my $helpInfo,'u' => \my $rmdup, 'k=i' => \my $k_length, 'w=i' => \my $win_size,'m' => \my $multiCntperSeq,
				'bin' => \my $bin, 'rc' => \my $rc,'t=i' => \my $threads, 'ori' => \my $ori, 'uOut' => \my $uOut); # i:int   f:float   s:string
    # initialize
    if ($helpInfo||(!$fileGlob)||(!$k_length)) { print "\nusage: $usage\n";exit;}
    my ($name,$dir,$suffix)= fileparse($fileGlob); $outputPath||="$dir/$this_programme"; $outputPath.="_".$k_length."mer"; # to sepcified path || to path with seq files
    $threads||=10; $win_size||=20;
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath/perlcmd);
    qx(echo ">> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd


#### get all files to process ######
my @file_ori= glob($fileGlob);
my @files = grep {(-f $_)&&($_=~ m/(txt$)|(gz$)|(fastq$)|(seq$)/)} glob($fileGlob); chomp @files;

&file::mProcess($threads,\&singleThread,\@files);
#&file::mThread($threads,\&singleThread,\@files);

sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
  my $file = shift;
  my %kmer; my $totalCnt=0; my %kmer_rmdup; my $totalCnt_rmdup=0;
  print "Counting the number of $k_length nt long kmers in $file\n";
  SeqAnal($file,\%kmer,\$totalCnt,\%kmer_rmdup,\$totalCnt_rmdup); # count k mer
  
  binRc(\%kmer, \$totalCnt) if $bin; # combine rc to the same bin if --bin
  calcfreq(\%kmer,\$totalCnt); 
  my $capLine= "kmer\tcount\tfreq\t"."totalCounts=$totalCnt\t";
  printHits(\%kmer,$file,$capLine, $k_length."mer");
  
  if ($rmdup) {  
  binRc(\%kmer_rmdup, \$totalCnt_rmdup) if $bin;# combine rc to the same bin if --bin
  calcfreq(\%kmer_rmdup, \$totalCnt_rmdup);
  my $capLine_rmdup= "kmer\tcount\tfreq\t"."totalCounts=$totalCnt_rmdup\t";
  printHits(\%kmer_rmdup,$file,$capLine_rmdup,$k_length."mer_u");
  }
  return(0)
}

sub SeqAnal {
  my ($file,$kmer_Href,$totalCntRef,$kmer_rmdup_Href,$totalCnt_rmdup_Ref) = @_;
  my @allreads; 
  use Mylib1::seqFile; &seqFile::getallSeq(\@allreads,$file);
  my $cnt=0;
  foreach(@allreads) { $cnt++; countKmer(\$_,$kmer_Href,$totalCntRef,$cnt); } # count kmer of all seq
  foreach (keys $kmer_Href) { delete $kmer_Href->{$_} if ($_=~ m/N/);} # throw away keys with N
  
  if ($rmdup) {  
  &seqFile::rmdup_3point(\@allreads,$win_size);
  # output uniq reads if necessary
  if ($uOut)  {	&seqFile::writeFile(\@allreads, "$outputPath/uniq_reads/", $file, "_uniq", $isgz); }
  
  my $cnt=0;
  foreach(@allreads) { $cnt++; countKmer(\$_,$kmer_rmdup_Href,$totalCnt_rmdup_Ref, $cnt); } # count kmer in deduplicated seq
  foreach (keys $kmer_rmdup_Href) { delete $kmer_rmdup_Href->{$_} if ($_=~ m/N/);} # throw away keys with N
  }
}

sub countKmer {
  my ($seq_p, $kmer_Href, $totalCntRef,$cnt)=@_;
  my $k = $k_length;
  my %beenThere;

  for (my $i=0;$i <= length(${$seq_p})-$k;$i++)
  {
    my $w = substr(${$seq_p},$i,$k); # take k mer
    unless ($multiCntperSeq)
    {
      #Count only one occurrence of a kmer per sequence
      $kmer_Href->{$w}->{cnt}++ && $$totalCntRef++ if !exists $beenThere{$w};
      $beenThere{$w}=1;
    }
    else
    {
      #Count all instances of a kmer per sequence
      $kmer_Href->{$w}->{cnt}++; $$totalCntRef++;
	  $kmer_Href->{$w}->{ori}.="$cnt;" if $ori;
      #Count rc sequences
      if ($rc)
      {
        $kmer_Href->{make_rc($w)}->{cnt}++;
		$kmer_Href->{$w}->{ori}.="$cnt;" if $ori;
		$$totalCntRef++;
	  }
	}
  }
  return(0);
}

sub make_rc
{  my $rc_k = reverse(shift); $rc_k =~ tr/ACGTacgt/TGCAtgca/; $rc_k;}

sub binRc
{
  my ($kmer_Href,$totalCntRef)= @_;
  foreach (keys $kmer_Href)
  {
    next if !(defined $kmer_Href->{$_}); # skip if already deleted due to rc of previous key
	my $rc_k= make_rc($_); 
    if ($_ eq $rc_k)  # double for palindromic kmer
    { $kmer_Href->{$_}->{cnt}*=2; $$totalCntRef += $kmer_Href->{$_}->{cnt};}
	elsif (defined ($kmer_Href->{$rc_k})) # exist rc count
    {  
      if (($_ cmp $rc_k)== -1) { $kmer_Href->{$_}->{cnt} += $kmer_Href->{$rc_k}->{cnt}; delete $kmer_Href->{$rc_k};}
      elsif (($_ cmp $rc_k)== 1) { $kmer_Href->{$rc_k}->{cnt} += $kmer_Href->{$_}->{cnt}; delete $kmer_Href->{$_};}
      else {die "error: not captured palindromic";}
    }
  }
}


sub calcfreq
{
  my $kmer_Href = shift; my @keys=(keys $kmer_Href);
  my $totalCntRef = shift; #my $totalwithPseudo= $$totalCntRef+$pseudoCnt*$#keys;
  foreach (@keys)
  {
    $kmer_Href->{$_}->{freq} = ($kmer_Href->{$_}->{cnt})/$$totalCntRef;
    #$kmer_Href->{$_}->{freqPseudo} = ($kmer_Href->{$_}->{cnt}+$pseudoCnt)/$totalwithPseudo;
  }
  return(0)
}



sub printHits {
  my ($kmer_Href,$file,$capLine,$nameAdd)=@_;
  
  map { $kmer_Href->{$_}->{ori} && chop($kmer_Href->{$_}->{ori}); } keys $kmer_Href; # delete the final ";"

  ##print out the hits
  my ($dir,$pre,$suf) = ($file =~ /(^.+\/|^)(.+)\.(.+$)/);
  open (OUT, ">","$outputPath/".$pre."_".$nameAdd.".txt") || die "Can't create file $pre.hits\n";
  print OUT $capLine."\n";
  print OUT join("\t", $_, $kmer_Href->{$_}->{cnt}, $kmer_Href->{$_}->{freq}, $kmer_Href->{$_}->{ori}||""), "\n" for (reverse sort {$kmer_Href->{$a}->{freq} <=> $kmer_Href->{$b}->{freq}} keys $kmer_Href); # sort by freq
  close OUT;

  return(0);
}

