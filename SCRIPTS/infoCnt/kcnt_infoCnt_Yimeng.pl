#!/usr/bin/perl -w
# 2015-05-06
# calc info content for same cyc0, '-A1-' like string should be in file name, '-A01-' OK

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Mylib1::Parallel::file;
use Mylib1::seqFile;
use Mylib1::fileIO;

my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme -f '*.txt' -b 'c0.txt' -g 'guide.txt' -o infoCnt_bkNum20k/147c2 -k 6 -p 1 -l 100000 -app c4_100k --m --u -bkNum 20000  (--ori --uOut)
  calc info content for same cyc0, '-A1-' like string should be in file name, '-A01-' OK

  ## general ##
   -f 'str'    input reads files, can use wildcards like 'xxx/*.txt' or 'xxx/*'
   -b 'str'    background file (cyc0)
   -o          output folder
   -h          get help
   -t int      number of threads to run (1 thread per file)
  ## specific ##
   -g str      the guidefile to append colum to
   -k int      size of the kmer to analyze. Default 6
   -w int      window size for deduplication, default 20
   --m         will count all possible kmer per sequences
               Default: only one kmer is counted per sequence entries
   --rc        count also the other strand
   --bin       bin the counts of k mer with its rc, x2 for palindromic k mer
   --u         count after removing duplicate ligand
   --ori       also output the original read No. of the k mer
   --uOut      output uniq reads
   -l int      lines to count for each file
   -p int      pseudo count for each kmer. only added to bk counts, df=1 (for good bkfiles)
   -bkNum int  add counts to sigHash as if having n reads from bkfile, df=20000    
   -app int    add additional char to col name of guidefile, like 'c4_100k'
   ";
	##### options
	GetOptions ('f=s' => \my $fileGlob, 'g=s' => \my $guidefile, 'o=s' => \my $outputPath, 'h' => \my $helpInfo,'u' => \my $rmdup, 'k=i' => \my $k_length, 'w=i' => \my $win_size,'m' => \my $multiCntperSeq,
				'bin' => \my $bin, 'app=s' => \my $append, 'rc' => \my $rc,'t=i' => \my $threads, 'bkNum=i' => \my $bkNum, 'l=i' => \my $lines, 'ori' => \my $ori, 'uOut' => \my $uOut, 'p=i' => \my $pseudoCnt); # i:int   f:float   s:string
    # initialize
    if ($helpInfo||(!$fileGlob)||(!$k_length)) { print "\nusage: $usage\n";exit;}
    my ($name,$dir,$suffix)= fileparse($fileGlob); $outputPath||="$dir/$this_programme"; $outputPath.="_".$k_length."mer"; # to sepcified path || to path with seq files
    $threads||=10; $win_size||=20; $pseudoCnt||=1; $bkNum||=20000;
	
	#output cmd info
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(mkdir -p $outputPath); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd

$rmdup="remove";
#### get all files to process ######
my @file_ori= glob($fileGlob);
my @files = glob($fileGlob); chomp @files;

#-----------count bk file------
my ($kmerbk_Href, $totalCntbk)=&kmerHashGen($k_length,$pseudoCnt); my ($kmer_rmdupbk_Href, $totalCnt_rmdupbk) =&kmerHashGen($k_length,$pseudoCnt);
my $bkfile= $files[0]; $bkfile=~ s/YB\d/0/; 
print "Counting the number of $k_length nt long kmers in **bk** $bkfile\n";
SeqAnal($bkfile,$kmerbk_Href,\$totalCntbk, $kmer_rmdupbk_Href,\$totalCnt_rmdupbk); # count k mer

#-----------count sig file, add infoContent to %info
&file::mProcess($threads,\&singleThread,\@files);
#&file::mThread($threads,\&singleThread,\@files);


#my @info = `cat $outputPath/infoCnt.txt`; chomp @info;
#my %info = map {my @line=split /\t/; $line[0] => $line[1]} @info;
#use Mylib1::guideFile; &guideFile::addColbyPos($guidefile,\%info,"info_".$append."_$k_length",$outputPath);
exit;


sub singleThread    #!!!!!!!  the thing done by each thread, to be modified!!!!!!!
{
  my $file = shift;
  
  my $liglength = &seqFile::liglength($file);
  my $kmerPerLig = $liglength-$k_length+1; my $totalForIni=$kmerPerLig*$bkNum;
  my ($totalCnt,$totalCnt_rmdup)=($totalForIni,$totalForIni);
  my ($kmer_Href, $kmer_rmdup_Href) = (\my %Hash1, \my %Hash2);
    #my ($kmer_Href, $totalCnt) =&kmerHashGen($k_length,$pseudoCnt); my ($kmer_rmdup_Href,$totalCnt_rmdup) =&kmerHashGen($k_length,$pseudoCnt);

  if ($rmdup) { &fillWithBk($kmer_rmdup_Href, $kmer_rmdupbk_Href, $totalForIni/$totalCnt_rmdupbk); }
  else { &fillWithBk($kmer_Href, $kmerbk_Href, $totalForIni/$totalCntbk); }
  
  print "Counting the number of $k_length nt long kmers in ".&fileIO::getfile($file)."\n";
  SeqAnal($file,$kmer_Href,\$totalCnt,$kmer_rmdup_Href,\$totalCnt_rmdup); # count k mer
  
  #-------- calc info cnt------
  my $kmerSig_Href= $rmdup? &calcfreq ($kmer_rmdup_Href, $totalCnt_rmdup):  &calcfreq($kmer_Href, $totalCnt);
  my $kmerBk_Href= $rmdup? &calcfreq ($kmer_rmdupbk_Href, $totalCnt_rmdupbk): &calcfreq($kmerbk_Href, $totalCntbk);
  my $infoContent=0;
  foreach (keys $kmerBk_Href)
  {
	my $sigFreq= $kmerSig_Href->{$_}->{freq}; my $bkFreq= $kmerBk_Href->{$_}->{freq};
	$infoContent+= $sigFreq*(log($sigFreq/$bkFreq)/log(2));
  }
  
  $infoContent= sprintf("%.5f", $infoContent);
  
  #print "\n".$infoContent."\n"; die;
  
  #$file =~ m/-([A-P]{1}\d{1,2})-/; my $well=$1;   #(?=[^\d]{1}) to make extration like TF105-NKX2-4IIIc correct
  #$well=~ s/([A-P])0(\d)/$1$2/;  # to change A01 into A1 for SELEX
  my ($plate,$TF,$well,$barcode,$cyc)= $file =~ m/([^-]*)-([^-]*)-([^-]*)-([^_YB]*)[YB]*(.)_/;
  qx(echo "$plate\t$barcode\t$well\t$TF\tcyc$cyc\t$infoContent\t\t$bkfile" >>$outputPath/infoCnt.txt);
  
  
#  if ($rmdup)
#  {  
#	binRc(\%kmer_rmdup, \$totalCnt_rmdup) if $bin;# combine rc to the same bin if --bin
#	calcfreq(\%kmer_rmdup, \$totalCnt_rmdup);
#	my $capLine_rmdup= "kmer\tcount\tfreq\t"."totalCounts=$totalCnt_rmdup\t";
#	printHits(\%kmer_rmdup,$file,$capLine_rmdup,$k_length."mer_u");
#  }
#  else
#  {
#	binRc(\%kmer, \$totalCnt) if $bin; # combine rc to the same bin if --bin
#	calcfreq(\%kmer,\$totalCnt); 
#	my $capLine= "kmer\tcount\tfreq\t"."totalCounts=$totalCnt\t";
#	printHits(\%kmer,$file,$capLine, $k_length."mer");
#  }
#  
  return(0)
}

sub fillWithBk {my ($Href,$bkHref,$frac)=@_; map { $Href->{$_}->{cnt}= $bkHref->{$_}->{cnt}*$frac;} keys $bkHref;}

sub SeqAnal {
  my ($file,$kmer_Href,$totalCntRef,$kmer_rmdup_Href,$totalCnt_rmdup_Ref) = @_;
  my @allreads; 
  &seqFile::getallSeq(\@allreads,$file,$lines);

  if ($rmdup)
  {  
	&seqFile::rmdup_3point(\@allreads,$win_size);
	# output uniq reads if necessary
	if ($uOut)
	{
	  system ("mkdir -p $outputPath/uniq_reads") if !(-d $outputPath."/uniq_reads");
	  my ($dir,$pre,$suf) = ($file =~ /(^.+\/|^)(.+)\.(.+$)/);
	  open(my $out, ">", "$outputPath/uniq_reads/$pre"."_uniq.txt"); print $out $_."\n" foreach @allreads;
	  close $out;
	}  
	my $cnt=0;
	foreach(@allreads) { $cnt++; countKmer(\$_,$kmer_rmdup_Href,$totalCnt_rmdup_Ref, $cnt); } # count kmer in deduplicated seq
	foreach (keys $kmer_rmdup_Href) { delete $kmer_rmdup_Href->{$_} if ($_=~ m/N/);} # throw away keys with N
  }
  else
  {
	my $cnt=0;
	foreach(@allreads) { $cnt++; countKmer(\$_,$kmer_Href,$totalCntRef,$cnt); } # count kmer of all seq
	foreach (keys $kmer_Href) { delete $kmer_Href->{$_} if ($_=~ m/N/);} # throw away keys with N
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
        $kmer_Href->{make_rc($w)}->{cnt}++; $$totalCntRef++;
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
  my $totalCnt = shift; 
  foreach (@keys)
  {
    $kmer_Href->{$_}->{freq} = ($kmer_Href->{$_}->{cnt})/$totalCnt;
  }
  return $kmer_Href;
}



sub printHits {
  my ($kmer_Href,$file,$capLine,$nameAdd)=@_;

  ##print out the hits
  my ($dir,$pre,$suf) = ($file =~ /(^.+\/|^)(.+)\.(.+$)/);
  open (OUT, ">","$outputPath/".$pre."_".$nameAdd.".txt") || die "Can't create file $pre.hits\n";
  print OUT $capLine."\n";
  print OUT join("\t", $_, $kmer_Href->{$_}->{cnt}, $kmer_Href->{$_}->{freq}, $kmer_Href->{$_}->{ori}||""), "\n" for (reverse sort {$kmer_Href->{$a}->{freq} <=> $kmer_Href->{$b}->{freq}} keys $kmer_Href); # sort by freq
  close OUT;

  return(0);
}

sub kmerHashGen
{
	my ($k,$pseudoCnt)=@_; #my $k_ori=$k;
	my @kmers= ("A","C","G","T"); my $totalCnt = 4**$k * $pseudoCnt;
	while ($k-1>0)
	{
		@kmers= map { $_."A", $_."C", $_."G", $_."T"} @kmers;
		$k--;
	}

	my %kmers= map {$_, {"cnt"=> $pseudoCnt}} @kmers;
	return (\%kmers,$totalCnt);
}