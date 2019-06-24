package ADM;
use warnings;
use strict;
use Data::Dumper;


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

sub addSR_simple_align3mer #(tmp_storage_hash_ref, readseq_for_process)
{       my $nucleotide_cntTmp_Ref=shift;
        my $readseq=shift;
   
        for (my $i = 0; $i < length($readseq)-2; ++$i)
		{
               	my $subreads3 = substr($readseq,$i,3); #extract dinucleotide signal for a specific position
				
				if (!($subreads3=~ m/N/)) {$nucleotide_cntTmp_Ref->{$subreads3}->[$i]++;}
		}
        return $nucleotide_cntTmp_Ref;
}


sub wtADM_asIs_capLine #(Input_ntFreq_hash, output_file_name,output_path, captionLine_key, --prefill0, --fillLong)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $output_path=shift;
		my $capKey=shift;
		my $prefill=shift; my $fillLong=shift;
		
		my $nucleotide_cnt_cap->{$capKey}=$nucleotide_cnt->{$capKey};
		&wtADM_asIs($nucleotide_cnt_cap,$fileename,$output_path);
		delete $nucleotide_cnt->{$capKey};
		&wtADM_asIs($nucleotide_cnt,$fileename,$output_path,$prefill,$fillLong);
}

sub wtADM_asIs #(Input_ntFreq_hash, output_file_name,output_path, --prefill0, --fillLong)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $output_path=shift; #$output_path.="/ADM";
		my $prefill=shift; #only fill 0 according to the length of each nt/multint, fill longer add --fillLong
		my $fillLong=shift; #fill all nt with 0 according to the length of single nt 
		
		if ($prefill)
		{
			my $liglength=scalar(@{$nucleotide_cnt->{A}});
            foreach my $keynt (keys($nucleotide_cnt))
			{	# replace all the undefined value with 0 (there should not be any if reads number is large enough)
				my $splice_l;
				if ($fillLong) {$splice_l=0;} else {$splice_l = (length $keynt)-1;}
				for (my $i = 0; $i < ($liglength - $splice_l); ++$i)
				{
				if (!defined ($nucleotide_cnt->{$keynt}->[$i])) {$nucleotide_cnt->{$keynt}->[$i]=0;}
				}
			}
        }
        	
        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, "+>>$output_path/$fileename") || die "unable to open write file handle in ADM::wtADM";
		foreach my $keynt (sort keys($nucleotide_cnt)) {	
			print OUT1 $keynt."\t";
			for (my $i = 0; $i < scalar(@{$nucleotide_cnt->{$keynt}}); ++$i)
			{
				if ($i==0) { print OUT1 $nucleotide_cnt->{$keynt}->[$i]; }
				else{ print OUT1 "\t".$nucleotide_cnt->{$keynt}->[$i];}
            }
			print OUT1 "\n";
		}
		close OUT1;
}

sub wtADM_asIs3mer #(Input_ntFreq_hash, output_file_name,output_path, ligand length, --prefill0, --fillLong)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $output_path=shift; #$output_path.="/ADM";
		my $liglength=shift;
		my $prefill=shift; #only fill 0 according to the length of each nt/multint, fill longer add --fillLong
		my $fillLong=shift; #fill all nt with 0 according to the length of single nt 
		
		if ($prefill)
		{
			foreach my $keynt (keys($nucleotide_cnt))
			{	# replace all the undefined value with 0 (there should not be any if reads number is large enough)
				my $splice_l;
				if ($fillLong) {$splice_l=0;} else {$splice_l = (length $keynt)-1;}
				for (my $i = 0; $i < ($liglength - $splice_l); ++$i)
				{
				if (!defined ($nucleotide_cnt->{$keynt}->[$i])) {$nucleotide_cnt->{$keynt}->[$i]=0;}
				}
			}
        }
        	
        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, "+>>$output_path/$fileename") || die "unable to open write file handle in ADM::wtADM";
		foreach my $keynt (sort keys($nucleotide_cnt)) {	
			print OUT1 $keynt."\t";
			for (my $i = 0; $i < scalar(@{$nucleotide_cnt->{$keynt}}); ++$i)
			{
				if ($i==0) { print OUT1 $nucleotide_cnt->{$keynt}->[$i]; }
				else{ print OUT1 "\t".$nucleotide_cnt->{$keynt}->[$i];}
            }
			print OUT1 "\n";
		}
		close OUT1;
}

sub wtADM_asIs3mer_capLine #(Input_ntFreq_hash, output_file_name,output_path, captionLine_key, liglength, --prefill0, --fillLong)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $output_path=shift;
		my $capKey=shift;
		my $liglength=shift;
		my $prefill=shift; my $fillLong=shift;
		
		my $nucleotide_cnt_cap->{$capKey}=$nucleotide_cnt->{$capKey};
		&wtADM_asIs3mer($nucleotide_cnt_cap,$fileename,$output_path,0,$liglength);
		delete $nucleotide_cnt->{$capKey};
		&wtADM_asIs3mer($nucleotide_cnt,$fileename,$output_path,0,$liglength,$prefill,$fillLong);
}

sub wtADM #(Input_ntFreq_hash, output_file_name, ligand_length, output_path, --trim)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $liglength=shift;
		my $output_path=shift; $output_path.="/ADM";
		my $trim=shift;
		
     	foreach my $keynt (keys($nucleotide_cnt))
		{	# replace all the undefined value with 0 (there should not be any if reads number is large enough)
			for (my $i = 0; $i < $liglength; ++$i)
			{
			if (!defined ($nucleotide_cnt->{$keynt}->[$i])) {$nucleotide_cnt->{$keynt}->[$i]=0;}
			}
		}
		
		if ($trim)
		{
            $nucleotide_cnt= &trim_multint($nucleotide_cnt, $liglength); # delete final 0 for output
			#print Dumper($nucleotide_cnt);die;
        }
        

        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, ">$output_path/$fileename") || die "unable to open write file handle in ADM::wtADM";
		foreach my $keynt (sort keys($nucleotide_cnt)) {	
			print OUT1 $keynt."\t";
			for (my $i = 0; $i < scalar(@{$nucleotide_cnt->{$keynt}}); ++$i)
			{
				print OUT1 $nucleotide_cnt->{$keynt}->[$i]."\t";
            }
			print OUT1 "\n";
		}
		close OUT1;

}

sub trim_multint #(inputhash_ref,#liglength) # delete the final 0s appended to multiNt counts
{
	my %nucleotide_cnt=%{(shift)};
	my $liglength=shift;
	foreach my $keynt (keys(%nucleotide_cnt))
	{
		my $splice_l = (length $keynt)-1;
		{
            my $splice_start = ($liglength||scalar(@{$nucleotide_cnt{A}}))-$splice_l;
			splice $nucleotide_cnt{$keynt}, $splice_start, $splice_l;
        }       
	}
	return \%nucleotide_cnt;
}

sub wtADM_200outADM #(Input_ntFreq_hash, output_file_name, ligand_length, output_path)
{
		my $nucleotide_cnt=shift;
        my $fileename=shift;
        my $liglength=shift;
		############
		$liglength+=53;
		############
		
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
	return \%admFreq;
	
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


sub rdADMv1 #(input_file_name, ligand_length)
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	#my $path=shift; #$path.="/ADM";
	# read in adm file
	my @data;
   # print "here\n"; print ">$path/$adm"."\n";
	open(DATA, "$adm") || die "sub rdADMv1: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
	  push @data, [split /\t/];
	}
	#print $liglength."\n";die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<20; $i++) {
		$admFreq{$data[$i][0]}= my $array_ref;
		my $j_length;
		if (length($data[$i][0])==2) {
			$j_length=93; # for dinucleotide freq
		}else{$j_length=94;} # for single nucleotide freq
		
		for (my $j=0; $j<$j_length; $j++) {
			$data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
			$data[$i][$j+1]=~tr/\t//d; # delete /t
			$data[$i][$j+1]=~tr/,/./; # transform , to .
			if ($data[$i][$j+1]==0) {
				$data[$i][$j+1]=0.00000000000000000001; #give a very small numer to avoid log0
			}
			$admFreq{$data[$i][0]}->[$j]=(log($data[$i][$j+1])/log(10));# transform into log form, admFreq start from 0
		}
	}
	close DATA;
	return \%admFreq;
	
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

sub rdADMv1_notlg #(input_file_name, ligand_length)
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	#my $path=shift; #$path.="/ADM";
	# read in adm file
	my @data;
   # print "here\n"; print ">$path/$adm"."\n";
	open(DATA, "$adm") || die "sub rdADMv1: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
	  push @data, [split /\t/];
	}
	#print $liglength."\n";die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<20; $i++) {
		$admFreq{$data[$i][0]}= my $array_ref;
		my $j_length;
		if (length($data[$i][0])==2) {
			$j_length=93; # for dinucleotide freq
		}else{$j_length=94;} # for single nucleotide freq
		
		for (my $j=0; $j<$j_length; $j++) {
			$data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
			$data[$i][$j+1]=~tr/\t//d; # delete /t
			$data[$i][$j+1]=~tr/,/./; # transform , to .
			$admFreq{$data[$i][0]}->[$j]=$data[$i][$j+1];
		}
	}
	close DATA;
	return \%admFreq;

}

sub rdADM_asIs #(input_file_name, ligand_length) # no triming
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	my @data;
	open(DATA, "$adm") || die "sub rdADMv1: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
		chomp;
	  push @data, [split /\t/];
	}
	#print Dumper(@data);die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<20; $i++) {
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

sub rdADMv200_notlg #(input_file_name, ligand_length)
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	#my $path=shift; #$path.="/ADM";
	# read in adm file
	my @data;
   # print "here\n"; print ">$path/$adm"."\n";
	open(DATA, "$adm") || die "sub rdADMv1: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
	  push @data, [split /\t/];
	}
	#print $liglength."\n";die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<20; $i++) {
		$admFreq{$data[$i][0]}= my $array_ref;
		my $j_length;
		if (length($data[$i][0])==2) {
			$j_length=146; # for dinucleotide freq
		}else{$j_length=147;} # for single nucleotide freq
		
		for (my $j=0; $j<$j_length; $j++) {
			$data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
			$data[$i][$j+1]=~tr/\t//d; # delete /t
			$data[$i][$j+1]=~tr/,/./; # transform , to .
			$admFreq{$data[$i][0]}->[$j]=$data[$i][$j+1];
		}
	}
	close DATA;
	return \%admFreq;

}

sub rdADM_asIs_3mer #(input_file_name, ligand_length) # no triming
{
	#######transfer data from adm file into hash######################
	my $adm=shift; #the adm file name
	my $liglength=shift; ##not necessary
	my @data;
	open(DATA, "$adm") || die "sub rdADM_asIs_3mer: Can't open nt_frequency file $adm: $!\n";
	while (<DATA>) {
		chomp;
	  push @data, [split /\t/];
	}
	#print Dumper(@data);die;
	my %admFreq; # store adm data for match
	for (my $i=0; $i<@data; $i++) {
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

sub WriteSummedADM
{
		(my $tmpFile_path, my $tmpfiles_arrRef, my $liglength, my $final_output_path, my $filename_final_ADM)=@_;
		#add together counts for all processes from tmp files
		my %nucleotide_cnt; # store the addup from all files
		foreach my $tmpfile (@{$tmpfiles_arrRef})
		{
			my $nt_cnt_temp = &rdADM($tmpfile,$liglength,$tmpFile_path); #param1=file name for input; #param2=ligand length; #param3=input path
			my %nt_cnt_temp = %{$nt_cnt_temp};
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

sub WriteSummedADM3mer
{
		(my $tmpFile_path, my $tmpfiles_arrRef, my $liglength, my $final_output_path, my $filename_final_ADM)=@_;
		#add together counts for all processes from tmp files
		my %nucleotide_cnt; # store the addup from all files
		foreach my $tmpfile (@{$tmpfiles_arrRef})
		{
			#print $tmpFile_path.$tmpfile."\n"; die;
			my $nt_cnt_temp = &rdADM_asIs_3mer($tmpFile_path."/ADM/".$tmpfile,$liglength); #param1=file name for input; #param2=ligand length; #param3=input path
			my %nt_cnt_temp = %{$nt_cnt_temp};
			foreach my $keynt (keys(%nt_cnt_temp))
			{	
				for (my $i = 0; $i < $liglength-2; ++$i)
				{
					$nucleotide_cnt{$keynt}->[$i] += $nt_cnt_temp{$keynt}->[$i];
				}
			}
			
		}

		#undef values already substituted by 0 in wtADM
		&wtADM(\%nucleotide_cnt,$filename_final_ADM,$liglength,$final_output_path);#param1=hashref for Input; #param2=file name for out put; #param3=ligand length; #param4=output path
		if (-d "$tmpFile_path") {system ("rm -r $tmpFile_path");} # if the temp folder exist delete it
}

sub WriteSummedADMv1
{
		(my $tmpFile_path, my $tmpfiles_arrRef, my $liglength, my $final_output_path, my $filename_final_ADM)=@_;
		#add together counts for all processes from tmp files
		my %nucleotide_cnt; # store the addup from all files
		foreach my $tmpfile (@{$tmpfiles_arrRef})
		{
			my $nt_cnt_temp = &rdADMv1($tmpfile,$liglength,$tmpFile_path); #param1=file name for input; #param2=ligand length; #param3=input path
			my %nt_cnt_temp = %{$nt_cnt_temp};
			foreach my $keynt (keys(%nt_cnt_temp))
			{	
				for (my $i = 0; $i < $liglength; ++$i)
				{
					$nucleotide_cnt{$keynt}->[$i] += $nt_cnt_temp{$keynt}->[$i];
				}
			}
			
		}

		&wtADM(\%nucleotide_cnt,$filename_final_ADM,$liglength,$final_output_path);#param1=hashref for Input; #param2=file name for out put; #param3=ligand length; #param4=output path
		if (-d "$tmpFile_path") {system ("rm -r $tmpFile_path");} # if the temp folder exist delete it
}


sub cntToFreq #(inputhashref)
{
	my $nucleotide_cnt=shift;
	#print Dumper($nucleotide_cnt);die;
	my $lineCnt = $nucleotide_cnt->{"A"}->[4]+$nucleotide_cnt->{"T"}->[4]+
				$nucleotide_cnt->{"C"}->[4]+$nucleotide_cnt->{"G"}->[4]; #use the sum of ATCG of the 5th position as the line count
	foreach my $keynt (keys(%{$nucleotide_cnt}))
			{	
				for (my $i = 0; $i < scalar(@{$nucleotide_cnt->{$keynt}}); ++$i)
				{
					$nucleotide_cnt->{$keynt}->[$i] = $nucleotide_cnt->{$keynt}->[$i]/$lineCnt;
				}
			}
	return $nucleotide_cnt;
}

sub cntToFreqv1 #(inputhashref)
{
	my $nucleotide_cnt=shift;
	my $lineCnt = linecntFromADM($nucleotide_cnt);
	foreach my $keynt (keys(%{$nucleotide_cnt}))
			{	
				for (my $i = 0; $i < scalar(@{$nucleotide_cnt->{$keynt}}); ++$i)
				{
					$nucleotide_cnt->{$keynt}->[$i] = $nucleotide_cnt->{$keynt}->[$i]/$lineCnt;
				}
			}
	return $nucleotide_cnt;
}

sub linecntFromADM #(inputhashref)
{
	my $nucleotide_cnt=shift;
	my $linecnt;
	map {$linecnt+=$nucleotide_cnt->{$_}->[4]} keys($nucleotide_cnt); #use the sum of the 5th position as the line count, for 3mer only
	return $linecnt;
}


1;	