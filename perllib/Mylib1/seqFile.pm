package seqFile;
# read, write of seq files
# Array with all seqs

use strict;

#sub getallSeq_old # compatible with fastq and .gz file, read all reads into $$allSeqAref
#{
#  my ($allSeqAref,$file,$lines) =@_;
#  my $fh;
#  use 5.010;
#  given($file)
#  {
#	when (m/\.gz$/) {open $fh, "-|", "gunzip -c $file" || die "Can't open reads file $file\n";}
#	default {open $fh, "<", "$file" || die "Can't open reads file $file\n";}
#  }
#  @$allSeqAref=<$fh>; chomp @$allSeqAref;  close $fh;
#  if ($$allSeqAref[0] =~ m/@/) {	@$allSeqAref= map {$$allSeqAref[$_*4+1]} (0..scalar(@$allSeqAref)/4-1); } # take only the sequence if is fastq file
#  splice(@$allSeqAref,$lines,-1) if ($lines && $lines < @$allSeqAref);
#}

sub getallSeq # compatible with fastq and .gz file, read all reads into $$allSeqAref, line by line so no waste if just read a few lines
{
  my ($allSeqAref,$file,$lines) =@_;
  my $fh;
  use 5.010;
  given($file)
  {
	when (m/\.gz$/) {open $fh, "-|", "gunzip -c $file" || die "Can't open reads file $file\n";}
	default {open $fh, "<", "$file" || die "Can't open reads file $file\n";}
  }
  my $isFastqFlag=0; my $currline=0; my $linesinAref=0;
  while (<$fh>)
  {
	if ($lines) {last if $linesinAref >= $lines}; # take all lines by default
	$currline++;
	if ($currline==1) {$isFastqFlag=1 if $_ =~ m/@/;}
	if ($isFastqFlag)
	{
	  if ($currline % 4 ==2) {chomp; push @$allSeqAref, $_; $linesinAref++;};
	}else
	{
	  chomp; push @$allSeqAref, $_; $linesinAref++;
	}	
  }
  close $fh;
}

sub liglength # compatible with fastq and .gz file, read all reads into $$allSeqAref, line by line so no waste if just read a few lines
{
  my ($file) =@_;
  my @array;
  &getallSeq(\@array,$file,1);
  return length $array[0];
}




sub rmdup_3point    # grep -v kmer at the beginning, middle, end to remove closely resembling kmer
{
    my ($allReads_Aref, $win_size) = @_;
	$win_size||=20;
    my %allReadsCnt;
    foreach (@$allReads_Aref)  { $allReadsCnt{$_}++; } # first remove identical duplicates
	#check length of lig
    my $liglength= length($$allReads_Aref[0]);
    undef @$allReads_Aref; # free memory
    

    ($win_size>$liglength)? ($win_size=$liglength):1;
    my $segmentNum=$liglength - $win_size +1;
    
    my %existed;    # try to remove sequences with close hamming distances
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
        
        push @$allReads_Aref, $seq;
    }
    $allReads_Aref;  
}

sub writeFile # write out put from array
{
   my ($allReads_Aref, $outPath, $file, $outFileNameAdd, $isgz, $autoDetect_gz) = @_;
   
   use Mylib1::fileIO; my $outFileName=&fileIO::getfile($file); # get filename from path
   $outFileName=~m/(\.[^\.]*)$/;
   if ($autoDetect_gz && $1=~m/\.gz/) { $isgz=1; }
   
   my $suffix= $isgz? ".gz": ($1=~m/\.gz/? ".txt":$1); my $replace= $outFileNameAdd.$suffix;
   $outFileName=~s/(\.[^\.]*)$/$replace/;
   
   my $out; `mkdir -p $outPath` if !(-d $outPath);
   if ($isgz)
   {
	  use IO::Compress::Gzip qw(gzip $GzipError) ; $out = new IO::Compress::Gzip("$outPath/$outFileName") or die "gzip failed: $GzipError\n";
   } else
   {
	  open($out, ">", "$outPath/$outFileName") or die "cannot open out file $outPath/$outFileName to write !\n";
   }
   print $out join "\n", @$allReads_Aref;
   close $out;   
}

sub kmerArrGen
{
	my ($k)=@_;
	my @kmers= ("A","C","G","T");
	while ($k-1>0)
	{
		@kmers= map {$_."A", $_."C", $_."G", $_."T"} @kmers;
		$k--;
	}
	return \@kmers;
}
    
1;    