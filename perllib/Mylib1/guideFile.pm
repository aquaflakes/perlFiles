package guideFile;
# read, write of guide files


use strict;

# add multiple col
sub addColFromFile
{
  my ($file,$guidefile,$outputPath,$cols_Aref,$colCaps_Aref) =@_; #["0","1","2","3"..], ["TF","class","family",".."] well should in 1st col
  my @colInfo= `cat $file`; chomp @colInfo;
  for (0..@$cols_Aref-1)
  {
	my $currind=$_;
	my %info= map {my @line=split /\t/; $line[0] => $line[$cols_Aref->[$currind]]} @colInfo;
	&addColbyPos($guidefile,\%info,$colCaps_Aref->[$currind],$outputPath);
  }
}

# add single col
sub addColbyPos # Href contains $Href->{A1}="xx"
{
  my ($guide,$Href,$newCol_name,$tmpguide_outputPath) =@_;
  
  use Mylib1::fileIO;
  my $guide_bk= $tmpguide_outputPath."/".&fileIO::getfile($guide);
  system("cp -a $guide $guide_bk");
  
  my $guide_temp= $tmpguide_outputPath."/".&fileIO::getfile($guide)."tmp.txt";

  open(my $fh, "<", $guide) or die "cannot open guide file to read";
  open(my $out, ">", $guide_temp) or die "cannot open guide file to write";

  my $posCol_index = 0; my $locate_pos_flag=0;
  while (<$fh>)
  {
	unless ($locate_pos_flag)
    {
		$_= $_=~ s/\r\n//gr;
		my @capline= split /\t/, $_; chomp @capline; my $found_flag;
        while ($posCol_index < $#capline)
		{
		  if ($capline[$posCol_index] eq "pos")
		  {
			$found_flag++;
			last;			
		  }else{$posCol_index++;} 
		}
        die "no 'pos' column in guide" unless $found_flag;
        $locate_pos_flag++;
        chomp $_; print $out $_."\t$newCol_name\n";
	}
    else
    {
      $_= $_=~ s/\r\n//gr;
	  chomp $_; my @line= split /\t/, $_;
      $_.= "\t".($Href->{$line[$posCol_index]}||"NA")."\n";
      print $out $_;
    }
  }
  close $fh;
  close $out;
  system("chmod 666 $guide_temp");
  system("cp -a $guide_temp $guide");
}

sub colIndex
{
  my ($guide,$colCap) =@_;
  my $capline= `head -n 1 $guide`; chomp $capline;
  my @capline= split /\t/, $capline; chomp @capline;
  @capline=map {s/^\s+|\s+$//gr} @capline;
  my $currColIndex=0;
  while ( $currColIndex++ < $#capline)
  {
	return $currColIndex if (uc($capline[$currColIndex]) eq uc($colCap));
  }
  die "$colCap not found in target $guide";	  
}

1;    