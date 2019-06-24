package parallel;

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../..";  # use the 2 lv up parent directory
use warnings;
use strict;
use Parallel::ForkManager;
use perllib::Me::Nu_cnt::checkfile;

									
sub calc 
{
		my $filename = shift;	my $filename_no_ext=$filename; $filename_no_ext =~ tr/\.txt//d; 
		my $liglength = shift;
		my $row_param = shift;
		my $cores = shift;
		my $tmpFile_path=shift||"./parallel_tmp";
		my $isFastq = shift;
		my $reads_anal_methodRef=shift; 
		my $tmpHash_output_methodRef=shift; 
		my $reads_output_methodRef=shift;
		my $filterDataRef=shift;
		my $exact =shift;
		
		my @tmp_hash_files; # list of tmp files generated from hash data, for further proc
		my @tmp_seq_files; # list of tmp seq files, for further proc

		
		
		my $readrows = &checkfile::SetRowNum ($filename,$row_param); # ensure readrows in the scope

		#########################################################


		my $lines_per_core = int($readrows/($cores));
		if ($isFastq) {$lines_per_core = $lines_per_core-($lines_per_core%4);} #each core handles 4x* lines, adjusted for only fastq
		my $pm = Parallel::ForkManager->new($cores);
			
	Process:	
	for(my $core_index = 0; $core_index < $cores; ++$core_index)
	{
		my $start_line=$core_index*$lines_per_core;
		my $end_line=($core_index+1)*$lines_per_core;
		if ($core_index==($cores-1)) {$end_line=$readrows;} #final core with more reads
		
		$pm->start and next Process;
		
		my %temp_storage = (A=>[],C=>[],T=>[],G=>[],);	#important for cases only di_nt is calculated, no s_nt information give err when ploting in R(designed for 20 lines)
		my $temp_storage_ref = \%temp_storage;	# store temp data
		
		open(SEQ, $filename) or die "failed to open fh for core[$core_index]"; # crashes unless open a fh for every process
				#print "core[$core_index]"."start:$start_line"."\t"."end:$end_line"."\n";
		#####################
		my $loopcnt1=0;
		for (my $i = 0; $i < $readrows; ++$i)
		{
			my $GroupedReadsArrRef; my $reads2;
			if ($isFastq)
			{
               #print "is Fastq\n";
			    my $reads1 = <SEQ>;#first is Illuminatus header
				$reads2 = <SEQ>;# is the actual sequence
				chomp $reads2; 	
				my $reads3 = <SEQ>;#another header
				my $reads4 = <SEQ>;#fastq quality string#code
				
				if(not $reads1=~ m/^\@HISE/) {
				  print"File might be corrupt, header doesn't start with HISE on read $i\n";
				  }
				if(not $reads3=~ m/^\+/) {
				  print"File might be corrupt, a quality line header doesn't start with + on read $i \n";
				  }
				$GroupedReadsArrRef=[$reads1,$reads2,$reads3,$reads4];
            }
			else{ $reads2 = <SEQ>; chomp $reads2; $GroupedReadsArrRef=["",$reads2];
				 }
            

			if( $.<= $start_line) {next;} if($.>$end_line) {last;} #only analyse the lines in scope
						#if (!defined($reads2)) {
						#	#print $.."line number of core[$core_index] error \n";
						#	last;
						#}else{print $.."line of core[$core_index] is ".$reads2."\n";
						#	  }
			
			################## !!!process of the reads #######################
		    if(defined($reads_anal_methodRef)) {$temp_storage_ref=$reads_anal_methodRef->($temp_storage_ref,$reads2,$liglength, $filterDataRef);}
			
			if(defined($reads_output_methodRef))
			{
				#my $filename_tmpSeq = $filename_no_ext."_core[$core_index]". "_seq_tmp.txt";
				$reads_output_methodRef->($tmpFile_path, $GroupedReadsArrRef, $filterDataRef, $exact, $lines_per_core, $core_index);
			} 
			################## !!!process of the reads #######################
			$loopcnt1++;
		}
		############# !!!!!!!!!!!!!!!! print "core[$core_index] looped $loopcnt1 cycs \n";

		
		#####################
		if(defined($tmpHash_output_methodRef))
		{
			my $filename_tmpHash = $filename_no_ext."_core[$core_index]". "_hash_tmp.txt";
			$tmpHash_output_methodRef->($temp_storage_ref,$filename_tmpHash,$liglength,$tmpFile_path);#param1=hashref for Input; #param2=file name for out put; #param3=ligand length; #param4=output path
		}
		#####################
		close SEQ;
		$pm->finish;


	}
		
	$pm->wait_all_children;
	
	#prepare a tmp file list for final processing
	if(defined($tmpHash_output_methodRef))
	{
		for(my $core_index = 0; $core_index < $cores; ++$core_index)
		{
			my $filename_tmp = $filename_no_ext."_core[$core_index]". "_hash_tmp.txt";
			push (@tmp_hash_files,$filename_tmp);
		}
	}
	
	if(defined($reads_output_methodRef))
	{
		for(my $core_index = 0; $core_index < $cores; ++$core_index)
		{
			my $filename_tmp = $filename_no_ext."_core[$core_index]". "_seq_tmp.txt";
			push (@tmp_seq_files, $filename_tmp);
		}
	}
	return ($tmpFile_path, \@tmp_hash_files, \@tmp_seq_files)	


}

1;	