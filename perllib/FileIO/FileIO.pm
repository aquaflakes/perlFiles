
#20150818
package FileIO;
use warnings;
use strict;
use File::Basename;


sub trim { my $s = shift; $s =~ s/^[\s\t]+|[\s\t]+$//g; return $s }; #delete begining or ending tab/space

sub readTable #(file seperated with tab,  --capline) !!!!! use two arrRef to accept the returned, use array will cause two ref all into the first 
{
    my ($filename, $capline) = @_;
    my @dataArr; # 2-dim array to store data
    my @caps; # 1-dim array to store caption of each column
    open(my $fh, "<", $filename) or die "cannot open datafile to read (sub readTable)";
    if ($capline) {@caps = split("\t",<$fh>);chomp @caps; foreach (@caps){ $_=&trim ($_);} }
    
    my $i= 0;
    while (my $line=<$fh>) {chomp $line; my @arr_ini=split("\t", $line);
                   $dataArr[$i]=\@arr_ini;
                   chomp(@{$dataArr[$i]});
                   foreach (@{$dataArr[$i]}){ $_= &trim ($_);} $i++;} # no space for each element
    return (\@dataArr, \@caps);    
}

sub readTable_1_lvHash #(input_file_name, ligand_length) # no triming
{
    	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	my @data;
	open(DATA, "$adm") || die "sub rdADMv1: Can't open table file to Hash $adm: $!\n";
    my $linecount;
	while (<DATA>)
    {
		chomp;	  
      m/\w/ && (push @data, [split /\t/] )&& $linecount++;
      
	}
	#print Dumper(@data);die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<$linecount; $i++) {
		$admFreq{$data[$i][0]}= my $array_ref;
		for (my $j=0; $j<($liglength||scalar(@{$data[$i]}))-1; $j++)
		{
			if (!defined($data[$i][$j+1])) {print $i."i".$j."j\n";}
            
			$data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
			$data[$i][$j+1]=~tr/\t//d; # delete /t
			$data[$i][$j+1]=~tr/,/./; # transform , to .
			$admFreq{$data[$i][0]}->[$j]=$data[$i][$j+1];
		}
	}
	close DATA;
	return \%admFreq;
}





sub writeTable # outFilepath data_HR, subkeys_to_print(ordered)_AR, capline_AR
{
    my ($outPath, $data_HR, $keys_print_AR, $caps_AR) = @_;
    my $dirname = dirname($outPath); my $filename= basename($outPath);
    my $dir_no_ext = $dirname =~ s/\.[^\.]*$//ir; # do not know why includes .pl the outdir cannot be created
    if (-d $dir_no_ext) {system ("rm -r $dir_no_ext");}
    system ("mkdir -p $dir_no_ext");# make path first unless cannot write
    
    open(my $fh, ">", $dir_no_ext."/".$filename) or die "cannot open file handle in 'writeTable'";
    if($caps_AR) {print $fh join("\t", @{$caps_AR});   print $fh "\n"}; # the capline
    foreach my $key (keys($data_HR))  # print the content
    {
        my @PrintTmpA;
        foreach my $subkey (@{$keys_print_AR}) { push @PrintTmpA, $data_HR->{$key}->{$subkey}; }
        print $fh join("\t",@PrintTmpA);
        print $fh "\n";
    }
    close $fh;
}

sub writeTable_1_LvHR # outFilepath data_HR, subkeys_to_print(ordered)_AR, capline_AR
{
    my ($outPath, $data_HR, $keys_print_AR, $caps_AR) = @_;
    my $dirname = dirname($outPath); my $filename= basename($outPath);
    my $dir_no_ext = $dirname; #=~ s/\.[^\.]*$//ir; # do not know why includes .pl the outdir cannot be created
    if (!-d $dir_no_ext) #{system ("rm -r $dir_no_ext");}
    {system ("mkdir -p $dir_no_ext")};# make path first unless cannot write
    
    open(my $fh, ">", $dir_no_ext."/".$filename) or die "cannot open file handle in 'writeTable'";
    if($caps_AR) {print $fh join("\t", @{$caps_AR});   print $fh "\n"}; # the capline
    foreach my $key (@$keys_print_AR)  # print the content
    {
        print $fh $key."\t";
        print $fh join("\t",$data_HR->{$key});
        print $fh "\n";
    }
    close $fh;
}

1;