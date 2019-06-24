package Qualityhippie;

#use FindBin;                 # locate this script
#use lib "$FindBin::Bin/../..";  # use the 2 lv up parent directory
use warnings;
use strict;


sub barcodeHash_genFromFile
{
	my $filterfile=shift;
	
	open (BARCODEFILTER, "$filterfile") or die "Can't open filterfile, $!\n";
	my @filterarray = <BARCODEFILTER>;
	my @names;
	my @searchwords;
	my @offset1;
	my @offset2;

	
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
	# print values %barcodehash;
	 #print keys %barcodehash;
	print "$filterfile contains the barcode length below.\n";
	for(my $i = 0; $i < $#front_length+1; $i++)
	{
		print("$front_length[$i] $random_length[$i]N $rear_length[$i]\n");
	}

	return 	[\@names, \@searchwords, \@front_length, \@random_length, \@rear_length, $max_length, \%barcodehash];
}




sub BarcodeSeperation
{
	my ($tmpFile_path, $GroupedReadsArrRef, $filterArrRef, $exact, $lines_per_core, $core_index)=@_;
	
	my ($names_ref, $searchwords_ref, $front_length_ref, $random_length_ref, $rear_length_ref, $max_length, $barcodehash_ref)= @{$filterArrRef};
	#my (@names, @searchwords, @front_length, @random_length, @rear_length, %barcodehash) = (@{$names_ref}, @{$searchwords_ref}, @{$front_length_ref}, @{$random_length_ref}, @{$rear_length_ref}, %{$barcodehash_ref});
	my (@names, @searchwords, @front_length, @random_length, @rear_length, %barcodehash);
	@names=@{$names_ref}; @searchwords=@{$searchwords_ref}; @front_length=@{$front_length_ref};  @random_length=@{$random_length_ref}; @rear_length=@{$rear_length_ref};  %barcodehash=%{$barcodehash_ref};
	
	
	my ($reads1,$reads2,$reads3,$reads4) = @{$GroupedReadsArrRef};
	#my @names;
	#my @searchwords;
	#my @front_length;
	#my @random_length;
	#my @rear_length;
	#my $max_length;
	#my %barcodehash;

	my $randomizedregion;
	my $failedreadcounter;
	my $parsedquality;

	my $targetfolder = $tmpFile_path ."/core$core_index". "_barcodes_tmp" ;  #combines filename to make the folder
	
	#if (!(-d "$targetfolder")) {system ("mkdir -p $targetfolder");} # if nt_count folder not exist create it first
	#mkdir -p $targetfolder;   #makes folder for fastq
	#my $outputfile = $targetfolder .'/' . "core$core_index first rows output.txt"; #makes the generic output name (of all nucleotide seq.)
	
	#my $targetfolder = $fileename . "_barcodes" ;  #combines filename to make the folder
	#print $targetfolder;
	if(-d $targetfolder){
		`rm -r $targetfolder`;
	}
	mkdir $targetfolder;   #makes folder for fastq
	my $roughmatch;
 
	open (TRASHBIN, ">>$targetfolder/Trashreads.txt");# makes trashbin for reads that do not match any of the given regular expressions
	open (TRASHBINFQ, ">>$targetfolder/Trashreads_fastq.txt");# makes trashbin for fastq that do not match any of the given regular expressions
	my $definement;# used for checking if the roughmatch is defined as a hash value, runs much faster (han monster) this way 
	$definement = 0;

		my @matchpatterns;
		my $orig_reads2 = $reads2;
		for(my $i = 0; $i < $#front_length+1; $i++){
			if($exact eq "Yes" and length($reads2) != $front_length[$i]+$random_length[$i]+$rear_length[$i]+1){next;}	#+1 for \n.
			if(length($reads2) < $max_length){chomp($reads2); $reads2 = $reads2."NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";}
			push(@matchpatterns, (substr $reads2, 0, $front_length[$i]).$random_length[$i]."N".(substr $reads2, $front_length[$i]+$random_length[$i], $rear_length[$i])); 
		}
		
		foreach my $roughmatch (@matchpatterns){
			if(defined $barcodehash{$roughmatch}) {#this checks whether the roughmatch is a value on the hash, makes the thing run much faster
				open (OUT1, ">>$targetfolder/$roughmatch._c$core_index.TEMPORARY.txt");	 #opens a filehandle to the temporary filename
				open (OUT2, ">>$targetfolder/$roughmatch._c$core_index.TEMPORARYjustreads.txt");	 #opens a filehandle to the temporary filename
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
	
				#my $Justreadfilename = "$roughmatch.Justreads";
				print OUT1 "$reads1";
				print OUT1 "$randomizedregion"."\n"; 	#writes the file into it
				print OUT1 "$reads3";
				print OUT1 "$parsedquality"."\n";
				print OUT2 "$randomizedregion"."\n";
			
				close OUT1;	#closes the filehandle, PERL was auto-closing them sporadically without any reasonable reason without this
				close OUT2;die;
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
	
	
	#system ("mkdir $targetfolder/Fastq");
	#system ("mkdir $targetfolder/JustReads");
	#my $readnamecounter = 0;
	
	
	
	#foreach my $readnamestochange (@names){
	#	if(-e "$targetfolder/$searchwords[$readnamecounter].TEMPORARY.txt"){
	#		system ("cp $targetfolder/$searchwords[$readnamecounter].TEMPORARY.txt $targetfolder/Fastq/$names[$readnamecounter].txt");
	#		system ("cp $targetfolder/$searchwords[$readnamecounter].TEMPORARYjustreads.txt $targetfolder/JustReads/$names[$readnamecounter].txt");
	#	}
	#$readnamecounter+=1;}
	#
	
	#everything should be now parsed into cute little files
	
									
									
	#system ("rm $targetfolder/*TEMPORARY*.txt");
	
	#my $tiima =localtime();
	#print "\nProgram finished $tiima\n";
	#my $fail_per = $failedreadcounter/$lines_per_core * 100;
	#print "Number of failed reads are $failedreadcounter in $lines_per_core ($fail_per %).\n";
}

1;

