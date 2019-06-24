#!/usr/bin/perl -w
#20150818

use warnings;
use strict;
use FindBin;  use lib "$FindBin::Bin/../..";  # use this directory
use Data::Dumper;
use perllib::FileIO::FileIO;
use File::Basename;

 ###########################!!!!!!!!!!!!! param to adjust
my $indexCol= 0; # field used as key for the storage hash
my $proteinNameCol=1;
#--------------------
my $proteinNameCol_inTarget=0;
my $consensusCol_inTarget=4;


my $this_programme = basename($0);
unless ($ARGV[0] && $ARGV[1]) {print "$this_programme  <file_containing_protein_info> <guide_with_consensus_info1> <guide_with_consensus_info2> <guide_with_consensus_info3>...\n prority of motifs 1..3"; exit;}

my ($allBarData, $capArr)= &FileIO::readTable($ARGV[0], "--cap"); # read in barcodes to search and corresponding info #!!!!! use two arrRef to accept the returned, use array will cause two ref all into the first
my $proteinNameKey=$capArr->[$proteinNameCol];
my %matchKey; #store keys to match the file list
foreach my $Data (@{$allBarData})
{
    my $index= $Data->[$indexCol]; 
    foreach (0..(scalar(@$capArr)-1))
    {
        $matchKey{$index}->{$capArr->[$_]}=$Data->[$_];   
        $matchKey{$index}->{consensus}= \my @arrini1;
    }
}

for (1..(@ARGV-1))
{
    my $cmd= "cat ".$ARGV[$_]." | grep -v 'Methylated' | grep -v 'GROUP'";
    my @allConsensus=`$cmd`;
    chomp @allConsensus;
    my %allConsensus; # storing all consensus info into hash for search
    
    foreach (@allConsensus)
    {
        my $proteinCurr= &FileIO::trim((split("\t",$_))[$proteinNameCol_inTarget]); my $consensusCurr=&FileIO::trim((split("\t",$_))[$consensusCol_inTarget]);
        if (defined($allConsensus{$proteinCurr})) {push $allConsensus{$proteinCurr}, $consensusCurr;} else {$allConsensus{$proteinCurr}= [$consensusCurr];}# multiple entries may have the same name
    }
    
    map {if(defined($allConsensus{$matchKey{$_}->{$proteinNameKey}})&& (scalar(@{$matchKey{$_}->{consensus}})<=5)&& ($allConsensus{$matchKey{$_}->{$proteinNameKey}}->[0] !~m/^(fail)|(AA)$/i) ) #
         { push $matchKey{$_}->{consensus}, @{$allConsensus{$matchKey{$_}->{$proteinNameKey}}}; }} keys(%matchKey);
    # add consensus info into original array only if the consensus is defined in current @allConsensus and less than 5 consensus in the result hash
}

push @$capArr, "consensus";
&writeTable("./$this_programme/foundConsensus.txt",\%matchKey,$capArr,$capArr); #[] is ref already so no \ needed
# outFilepath data_HR, subkeys_to_print(ordered)_AR, capline_AR

exit;






#-------------------------------------
sub writeTable # outFilepath data_HR, subkeys_to_print(ordered)_AR, capline_AR
{
    my ($outPath, $data_HR, $keys_print_AR, $caps_AR) = @_;
    my $dirname = dirname($outPath); my $filename= basename($outPath);
    my $dir_no_ext = $dirname =~ s/\.[^\.]*$//ir; # do not know why includes .pl the outdir cannot be created
    if (-d $dir_no_ext) {system ("rm -r $dir_no_ext");}
    system ("mkdir -p $dir_no_ext");# make path first unless cannot write
    
    open(my $fh, ">", $dir_no_ext."/".$filename) or die "cannot open file handle in 'writeTable'";
    print $fh join("\t", @{$caps_AR});   print $fh "\n"; # the capline
    
    my @keys=keys($data_HR); # my @sortedKeys=sort { (substr($data_HR->{$a}->{well},0,1) cmp substr($data_HR->{$b}->{well},0,1))||substr($data_HR->{$a}->{well},1,2) <=> substr($data_HR->{$b}->{well},1,2)  } @keys; # arrange result in Hash , first by alphabetical order, then by well No.
    my @sortedKeys=sort { $data_HR->{$a}->{no} <=> $data_HR->{$b}->{no} } @keys; # arrange result in Hash by well No.
    
    foreach my $key (@sortedKeys)  # print the content
    {
        my @PrintTmpA;
        foreach my $subkey (@{$keys_print_AR})
        {
            if (($subkey eq "fileLinks")||($subkey eq "consensus"))
            { push @PrintTmpA, @{$data_HR->{$key}->{$subkey}};}
            else
            { push @PrintTmpA, $data_HR->{$key}->{$subkey};}
        }
        print $fh join("\t",@PrintTmpA); # do not care about the undef
        print $fh "\n";
        
    }
    close $fh;
}

sub uniq
{
    my %seen;
    grep !$seen{$_}++, @_;
}

