package parallel_v1;

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../..";  # use the 2 lv up parent directory
use warnings;
use strict;
use Parallel::ForkManager;
use perllib::Me::Nu_cnt::checkfile;

									
sub calc 
{
		my $filename = shift;	$filename =~ m{\/([^/]*)\.\w+$}|| $filename =~m/(.*)...\w+$/; my $filename_no_ext = $1; 
		my $liglength = shift;
		my $row_param = shift;
		my $cores = shift;
		my $tmpFile_path=shift||"./parallel_tmp";
		my $isFastq = shift;
		my $tmp_storage_Ref=shift; 
		my $reads_anal_methodRef=shift; 
		my $tmpStorage_output_methodRef=shift;
		my $matrix_NamesRef=shift;
		my $cut_values_Ref=shift;


		
		my @tmp_storage_files; # list of tmp files generated from hash data, for further proc
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
		my $combinedReads="";
				
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
				
				########### original position for reads processing
								$combinedReads.=$reads2;
								#########################
			$loopcnt1++;
		}
		#print "core[$core_index] looped $loopcnt1 cycs \ncombined Reads:\n$combinedReads\n\n";
		
		################## !!!process of the reads #######################
		$tmp_storage_Ref=$reads_anal_methodRef->($tmp_storage_Ref,$combinedReads,$liglength,$cut_values_Ref);			
		################## !!!process of the reads #######################

		
		#####################
		if(defined($tmpStorage_output_methodRef))
		{
			my $filename_tmpStorage = $filename_no_ext."_core[$core_index]". "_storage_tmp.txt";
			$tmpStorage_output_methodRef->($tmp_storage_Ref,$matrix_NamesRef,$tmpFile_path,$filename_tmpStorage);#param1=hashref for Input; #param2=file name for out put; #param3=ligand length; #param4=output path
		}
		#####################
		close SEQ;
		$pm->finish;


	}
		
	$pm->wait_all_children;
	
	#prepare a tmp file list for final processing
	if(defined($tmpStorage_output_methodRef))
	{
		for(my $core_index = 0; $core_index < $cores; ++$core_index)
		{
			my $filename_tmp = $tmpFile_path.$filename_no_ext."_core[$core_index]". "_storage_tmp.txt";
			push (@tmp_storage_files,$filename_tmp);
		}
	}
	

	return (\@tmp_storage_files,$readrows);	


}

1;	