package ADM;
use warnings;
use strict;
#use Data::Dumper;


sub addSR_simple_align #(tmp_storage_hash_ref, readseq_for_process)
{       my $nucleotide_cntTmp_Ref=shift;
        my $readseq=shift;
   
        for (my $i = 0; $i < length($readseq); ++$i)
		{
                my $subreads1 = substr($readseq,$i,1); #extract single nucleotide signal for a specific position
				my $subreads2 = substr($readseq,$i,2); #extract dinucleotide signal for a specific position
				
				if (!($subreads1 eq "N")) {$nucleotide_cntTmp_Ref->{$subreads1}->[$i]++;}
				
				# manipulate dinucleotide counts in %nucleotide_cnt, rule out the final digit where length ($subreads2) ==1 
                if (length ($subreads2) ==2 && !($subreads2=~ m/N/))
                { 
						$nucleotide_cntTmp_Ref->{$subreads2}->[$i]++;
				}
		}
        return $nucleotide_cntTmp_Ref;
}

sub wtADM #(Input_ntFreq_hash, output_file_name, ligand_length, output_path)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $liglength=shift;
		my $output_path=shift; $output_path.="/ADM";
     	foreach my $keynt (keys($nucleotide_cnt))
		{	# replace all the undefined value with 0 (there should not be any if reads number is large enouth)
			for (my $i = 0; $i < $liglength; ++$i)
			{
			if (!defined ($nucleotide_cnt->{$keynt}->[$i])) {$nucleotide_cnt->{$keynt}->[$i]=0;}
			}
		}
		#print "$output_path \n";

        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, ">$output_path/$fileename") || die "unable to open write file handle in ADM::wtADM";
		foreach my $keynt (sort keys($nucleotide_cnt)) {	
			print OUT1 $keynt."\t";
			for (my $i = 0; $i < $liglength; ++$i) {
				print OUT1 $nucleotide_cnt->{$keynt}->[$i]."\t";
                }
			print OUT1 "\n";
		}
		close OUT1;

}

sub rdADM #(input_file_name, ligand_length, input_path)
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift;
	my $path=shift; $path.="/ADM";
	# read in adm file
	my @data;
   # print "here\n"; print ">$path/$adm"."\n";
	open(DATA, "$path/$adm") || die "sub rdADM: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
	  push @data, [split /\t/];
	}
	#print $liglength."\n";die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<1000000; $i++) {
		if (!defined($data[$i][0])) {last;} #quit when no more lines available
        
		$admFreq{$data[$i][0]}= my $array_ref;
		#my $j_length;
		#$j_length=$liglength;
		#if (length($data[$i][0])==2) {
		#	$j_length=$liglength-1; # for dinucleotide freq
		#}else{$j_length=$liglength;} # for single nucleotide freq
		
		for (my $j=0; $j<$liglength; $j++) {
			$data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
			$data[$i][$j+1]=~tr/\t//d; # delete /t
			$data[$i][$j+1]=~tr/,/./; # transform , to .
			$admFreq{$data[$i][0]}->[$j]=$data[$i][$j+1];
			#if ($data[$i][$j+1]==0) {
			#	$data[$i][$j+1]=0.00000000000000000001; #give a very small numer to avoid log0
			#}
			#$admFreq{$data[$i][0]}->[$j]=(log($data[$i][$j+1])/log(10));# transform into log form, admFreq start from 0
		}
	}
	close DATA;
	return %admFreq;
	
	##print adm hash for check
	#if (!(-d "adm_check")) {system ("mkdir adm_check");} # if nt_count folder not exist create it first
	#open (OUT1, ">adm_check/"."_adm_hash.txt");
	#	foreach my $keynt (sort keys(%admFreq)) {	
	#		print OUT1 $keynt."\t";
	#		my $dim_length=@{$admFreq{$keynt}}; #length of the array of each key
	#		for (my $i = 0; $i < $dim_length; ++$i) {
	#			print OUT1 10**$admFreq{$keynt}->[$i]."\t";
	#		}
	#		print OUT1 "\n";
	#	}
	#	close OUT1;
	#######transfer data from adm file into hash######################

}

sub WriteSummedADM
{
		(my $tmpFile_path, my $tmpfiles_arrRef, my $liglength, my $final_output_path, my $filename_final_ADM)=@_;
		#add together counts for all processes from tmp files
		my %nucleotide_cnt; # store the addup from all files
		foreach my $tmpfile (@{$tmpfiles_arrRef})
		{
			my %nt_cnt_temp = &rdADM($tmpfile,$liglength,$tmpFile_path); #param1=file name for input; #param2=ligand length; #param3=input path
			foreach my $keynt (keys(%nt_cnt_temp))
			{	
				for (my $i = 0; $i < $liglength; ++$i)
				{
					$nucleotide_cnt{$keynt}->[$i] += $nt_cnt_temp{$keynt}->[$i];
				}
			}
			
		}

		#undef values already substituted by 0 in wtADM
		&wtADM(\%nucleotide_cnt,$filename_final_ADM,$liglength,$final_output_path);#param1=hashref for Input; #param2=file name for out put; #param3=ligand length; #param4=output path
		if (-d "$tmpFile_path") {system ("rm -r $tmpFile_path");} # if the temp folder exist delete it
}




1;	