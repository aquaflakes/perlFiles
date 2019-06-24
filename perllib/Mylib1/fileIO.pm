package fileIO;

use strict;
use File::Basename qw(fileparse);

sub getdir
{
  my($path)=@_;
  my($filename, $dir, $suffix) = fileparse($path);
  return $dir;
}

sub getfile
{
  my($path)=@_;
  my($filename, $dir, $suffix) = fileparse($path);
  return $filename;
}

sub getsuffix
{
  my($path)=@_;
  my($filename, $dir, $suffix) = fileparse($path);
  return $suffix;
}

sub readtoHash  # read two vals in everyline to form a key-value pair
{
  my ($filename,$keycol,$valcol,$skipCaption)=@_;
  $keycol||=0;
  my %hash;
  open(my $fh, "<", $filename) or die "cannot open file fileIO::readtoHash";
  <$fh> if $skipCaption;
  while (my $line = <$fh>)
  {
	chomp $line;
    my @allval= split "\t", $line;
    $hash{$allval[$keycol]}=$allval[$valcol];    
  }
  return \%hash;  
}

sub logcmd
{
  my ($outputPath)=@_;
  use File::Basename;  my $this_programme = basename($0);
  my $cmdInput=$this_programme." ".join(" ", @ARGV);
  qx(mkdir -p $outputPath/perlcmd_log);
  qx(echo ">> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd_log/perlcmd.txt); # memo input cmd
}
    
1;    