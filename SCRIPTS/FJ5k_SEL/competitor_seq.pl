#!/usr/bin/perl -w
#20150818

use warnings;
use strict;
use FindBin;  use lib "$FindBin::Bin/../..";  # use this directory
use Data::Dumper;
use perllib::FileIO::FileIO;
use File::Basename;



my $this_programme = basename($0);
unless ($ARGV[0] && $ARGV[1] && $ARGV[2]) {print "$this_programme  <file_containing_barcodes_info> <batch name> <html file> \n"; exit;}

my ($allBarData, $capArr)= &FileIO::readTable($ARGV[0], "--cap"); # read in barcodes to search and corresponding info #!!!!! use two arrRef to accept the returned, use array will cause two ref all into the first
my $batchName= $ARGV[1];
my $html = $ARGV[2];
my %matchKey; #store keys to match the file list
foreach (@{$allBarData})
{
    $matchKey{$_->[0]."_$batchName"}->{barcode}=$_->[0];     #!!!!!!!!!!!!! param to adjust : which col in the input barcode file contains corrsponding info
    $matchKey{$_->[0]."_$batchName"}->{protein}=$_->[5];    #!!!!!!!!!!!!! param to adjust 
    $matchKey{$_->[0]."_$batchName"}->{consensus}= \my @arr_for_ini_1;
    $matchKey{$_->[0]."_$batchName"}->{fileLinks}= \my @arr_for_ini_2;
    $matchKey{$_->[0]."_$batchName"}->{well}=$_->[4];   #!!!!!!!!!!!!! param to adjust 
    
}
#my @succ_svg= `ls  /var/www/kaz/data/svg/svg/`; # svg only exists for successful enriched, check if the file name pattern matches the patten constructed from input
#chomp @succ_svg;
my $ScreenWrd="_$batchName"."_";

# open http://nutcase.biosci.ki.se/kaz/index.php?base=BatchKW and find svg links (consensus seq contained in svg filename)
open(my $php, "<", $ARGV[2]) or die "fail to open SELEX html file";
my $phpContent=join("\t", <$php>);
my @svgs = &uniq ($phpContent =~ m/\/svg\/svg\/(.*?\.svg)/gi);# take uniq records from page


foreach (@svgs)
{
    unless ($_=~ m/$ScreenWrd/) {next;} # first screen out the files with wrong batch name
    my $keyExt=$_=~s/(.*?_.*?)_.*((c4|c3b0)[^\.]*.svg)/$1/ir; # only use cyc4 and c3b0 data #!!!!!!!! prarm to adjust
    # -r: not modifying but just return the modified value, ? after * means non-greedy
        
    if (defined($matchKey{$keyExt})) # if the key is in key list
    {
        push $matchKey{$keyExt}->{consensus}, $_=~s/.*_(?<consensus>[^_]*)_m.{1,3}_c.{1,5}\.svg/$+{consensus}/ir; # push the consensus seq suggested by svg file name
        push $matchKey{$keyExt}->{fileLinks}, "/var/www/kaz/data/svg/svg/".$_;
    }
}

&writeTable("./$this_programme/foundLinks.txt",\%matchKey,["well","barcode","protein","consensus"],["well","barcode","protein","consensus(competitor)"]); #[] is ref already so no \ needed
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
    
    my @keys=keys($data_HR); my @sortedKeys=sort { $data_HR->{$a}->{well} cmp $data_HR->{$b}->{well} } @keys; # arrange result in Hash by well No.
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
        print $fh join("\t",@PrintTmpA);
        print $fh "\n";
        
    }
    close $fh;
}

sub uniq
{
    my %seen;
    grep !$seen{$_}++, @_;
}

