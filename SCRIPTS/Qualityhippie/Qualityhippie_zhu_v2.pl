#!/usr/bin/perl
# operates gz file
# 2015-10-12, used for FJKSEL5.1test, updated 2015-10-25

use warnings;
use strict;

 #initialize stuff
##################
my $title = "Qualityhippie V2.3";
my $usage = "\$ Qualityhippie <fastq file> <row number | all> <filter file> -exact|0 <Batchname> <cyc No.> <-small>|0
-small: move all files smaller than 3M into a seperate folder
* run under the folder to contain 'fastq' and 'JustRead' folder";

print "$title\n\n";
if(!defined($ARGV[0]) or !defined($ARGV[1]) or !defined($ARGV[2]) or !defined($ARGV[3]) or !defined($ARGV[4]) or !defined($ARGV[5])){
	print "$usage\n\n";
	exit;
}

my $exact = "No";	#Exact match of the Fastq length.
if($ARGV[3] eq "-exact"){$exact = "Yes";}
my $BatchName = $ARGV[4]."_";
my $CycNo = "_".$ARGV[5];

my $tiimaalku =localtime();
print "Program started $tiimaalku\n";

my $fileename = $ARGV[0];
open(FILEENAME, $fileename) or die "Cannot open the input file: $!\n ";
my $row_amount = $ARGV[1];
my $filterfile = $ARGV[2];
my $targetfolder = "./" ;  #combines filename to make the folder
#my $targetfolder = $fileename . "_barcodes" ;  #combines filename to make the folder

my $outputfile = $row_amount . "_first_rows_output.txt"; #makes the generic output name (of all nucleotide seq.)
open (BARCODEFILTER, "$filterfile") or die "Can't open filterfile, $!\n";
my @filterarray = <BARCODEFILTER>;
my @names;
my @searchwords;
my @offset1;
my @offset2;
my $randomizedregion;
my $parsedquality;

shift(@filterarray);
foreach (@filterarray){
	my ($rowow1, $rowow2, $rowow3, $rowow4) = split(/\t/);
	if($rowow2 =~ /[a-z]/){next;}
	push (@names, $rowow1);
	push (@searchwords, $rowow2);
	push (@offset1 , $rowow3);
	push (@offset2, $rowow4);
}
my %barcodehash;
my $failedreadcounter = 0;
my $hashnumber = 0;
my @front_length;
my @random_length;
my @rear_length;
my %filter_checker;
my $max_length = 0;
foreach my $barcodehashish (@searchwords){
	#$hashnumber+=1;
	chomp $barcodehashish;
	$barcodehashish =~ /\d*N/;
	my $temp1 = $`;
	my $temp2 = $&;
	my $temp3 = $';
	$temp2 =~ s/N//;
	if(!defined($temp1)){$temp1 = "";}
	if(!defined($temp3)){$temp3 = "";}
	my $temp = length($temp1).$temp2.length($temp3);
	if($max_length < $temp){$max_length = $temp;}
	if(!defined($filter_checker{$temp})){
		push(@front_length, length($temp1));
		push(@random_length, $temp2);
		push(@rear_length, length($temp3));
	}
	$filter_checker{$temp} = 1;
	$barcodehash{$barcodehashish} = "$hashnumber";
 }
my @listofahshkeysh = sort keys %barcodehash;
foreach my $blaaaaargh (@listofahshkeysh){
					print "key $blaaaaargh\thas value \t$barcodehash{$blaaaaargh}\n";
					}

 
 if ($row_amount eq "all"){
	$row_amount = 666666666666666666;
 }

## do not need, count when encountering the end of file
#my $readrows;
#$readrows = $row_amount;
#if ($row_amount > 10000000){#if the row amount is smaller than puny ten millions the program won't bother with calculating the total rows
#	my $f = $ARGV[0];#to check whether there is less reads than what is asked which takes a lot of time
#	my $linecount = "0";
#
#	$linecount = `wc -l $ARGV[0]`;
#	chomp($linecount);
#	$linecount =~ s/\s$ARGV[0]//;
#	print "\nInput fastq file contains $linecount lines.\n\n";
##	close (TXT);
#	if ($row_amount > $linecount){# sees whether the total linenumber is actually smaller than what was asked for
#	$readrows = $linecount;
#	}
#}	

print "$filterfile contains the barcode length below.\n";
for(my $i = 0; $i < $#front_length+1; $i++){
	print("$front_length[$i] $random_length[$i]N $rear_length[$i]\n");
}

#########################################################
unless ( open(OUTPUTROWS, "+>$targetfolder/$outputfile") ) {
       #checks whether the generic outputfile works and opens the path to it
    print "Cannot open file \"$outputfile\" to write to!!\n\n";
    exit;        #checks the ooutput for generic seq.
   }
close OUTPUTROWS;
system("rm $outputfile");

my $subscalar1;
my $subscalar2;
my $subscalar3;
my $subscalar4;
my $subscalar5;
my $roughmatch;
# system ("rm $targetfolder/*TEMPORARY.txt");#destroys TEMPORARY files made by earlier runs on same dataset
 #needed as the filehandle can only append to the formed files, as it would otherwise delete data by putting reads on top of each other
 
open (TRASHBIN, "+>$targetfolder/Trashreads.txt");# makes trashbin for reads that do not match any of the given regular expressions
open (TRASHBINFQ, "+>$targetfolder/Trashreads_fastq.txt");# makes trashbin for fastq that do not match any of the given regular expressions
#open (CRAPPYREADS, "+>$targetfolder/Crapreads.txt");# makes trashbin for reads that have N basecalls
my $definement;# used for checking if the roughmatch is defined as a hash value, runs much faster (han monster) this way 

my $linecounter=0;
while ((my $reads1 = <FILEENAME>) && $linecounter<$row_amount/4) { #first is Illuminatus header
	$definement = 0;
    	my $reads2 = <FILEENAME>;# is the actual sequence
    	  	my $reads3 = <FILEENAME>;#another header
    	  	  	my $reads4 = <FILEENAME>;#fastq quality string

 #if(not $reads1=~ m/^\@HISE/) {
 #  print"File might be corrupt, header doesn't start with HISE on read $i\n";
 #  }
 #if(not $reads3=~ m/^\+/) {
 #  print"File might be corrupt, a quality line header doesn't start with + on read $i \n";
 #  }

	my @matchpatterns;
	my $orig_reads2 = $reads2;
	for(my $i = 0; $i < $#front_length+1; $i++){
		if($exact eq "Yes" and length($reads2) != $front_length[$i]+$random_length[$i]+$rear_length[$i]+1){next;}	#+1 for \n.
		if(length($reads2) < $max_length){chomp($reads2); $reads2 = $reads2."NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";}
		push(@matchpatterns, (substr $reads2, 0, $front_length[$i]).$random_length[$i]."N".(substr $reads2, $front_length[$i]+$random_length[$i], $rear_length[$i])); 
	}
	
	foreach my $roughmatch (@matchpatterns){
		if(defined $barcodehash{$roughmatch}) {#this checks whether the roughmatch is a value on the hash, makes the thing run much faster
			open (OUT1, ">>$targetfolder/$roughmatch.TEMPORARY.txt");	 #opens a filehandle to the temporary filename
			open (OUT2, ">>$targetfolder/$roughmatch.TEMPORARYjustreads.txt");	 #opens a filehandle to the temporary filename
			$definement = 1;
			$roughmatch =~ /\d*N/;
			my $temp1 = $`;
			my $temp2 = $&;
			my $temp3 = $';
			$temp2 =~ s/N//;
			my $char_num = 0;
			if($temp2 < 100){$char_num = 2;}else{$char_num = 3;}
			$randomizedregion = (substr $reads2, length($temp1), $temp2);
			$parsedquality = (substr $reads4, length($temp1), $temp2);

			my $Justreadfilename = "$roughmatch.Justreads";
			print OUT1 "$reads1";
			print OUT1 "$randomizedregion"."\n"; 	#writes the file into it
 			print OUT1 "$reads3";
 			print OUT1 "$parsedquality"."\n";
 			print OUT2 "$randomizedregion"."\n";
 		
 			close OUT1;	#closes the filehandle, PERL was auto-closing them sporadically without any reasonable reason without this
 			close OUT2;
 		} 
 		
 	}
	if ($definement != 1){
		my $reads2 = $orig_reads2;
		print TRASHBIN "$reads2";
		$failedreadcounter+=1;
		print TRASHBINFQ "$reads1";
		print TRASHBINFQ "$reads2";
		print TRASHBINFQ "$reads3";
		print TRASHBINFQ "$reads4";
	}
	$linecounter++;
}
my $numberofsequences=int($./4);
print "\nInput fastq file contains ".$numberofsequences." reads.\n\n";

if (-d "$targetfolder/Fastq") {system ("rm -r $targetfolder/Fastq");}
if (-d "$targetfolder/JustReads") {system ("rm -r $targetfolder/JustReads");}
system ("mkdir $targetfolder/Fastq");
system ("mkdir $targetfolder/JustReads");
my $readnamecounter = 0;


foreach my $readnamestochange (@names){
	if(-e "$targetfolder/$searchwords[$readnamecounter].TEMPORARY.txt"){
	#########################
	$names[$readnamecounter]=$BatchName.$names[$readnamecounter].$CycNo;
	#########################
		system ("cp $targetfolder/$searchwords[$readnamecounter].TEMPORARY.txt $targetfolder/Fastq/$names[$readnamecounter].txt");
		system ("cp $targetfolder/$searchwords[$readnamecounter].TEMPORARYjustreads.txt $targetfolder/JustReads/$names[$readnamecounter].txt");
	}
$readnamecounter+=1;
										}


#everything should be now parsed into cute little files


close FILEENAME;# closes original filename, as it is not needed anymore
								
								
system ("rm $targetfolder/*TEMPORARY*.txt");

########################
#print join(",", @ARGV);
if ($ARGV[6]) {
	chdir "$targetfolder/Fastq/";
	 system ("mkdir -p small_files");
	 system ('find *.* -size -3085760c -exec mv {} ./small_files \;');
	chdir "../JustReads/";
	 system ("mkdir -p small_files");
	 system ('find *.* -size -3085760c -exec mv {} ./small_files \;');	
}
########################

my $tiima =localtime();
print "\nProgram finished $tiima\n";
my $fail_per = $failedreadcounter/$numberofsequences * 100;
print "Number of failed reads are $failedreadcounter in $numberofsequences reads ($fail_per %).\n";
 exit;


