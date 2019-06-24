package checkfile;
use warnings;
use strict;


sub SetRowNum #param1=filename; #param2=row amount param from cmd
{
		my $fileename = shift;
        my $row_amount = shift;
		if ($row_amount eq "all"){
			$row_amount = 666666666666666666;
		 }		
		 #tests if the file has less reads than what have been asked for, if thats the case, it'll
		 #take only the asked amount
		my $readrows;
		$readrows = $row_amount;

			my $linecount = "0"; 
			$linecount = `wc -l $fileename`;
			chomp($linecount);
			$linecount =~ s/\s$fileename//;
			if ($row_amount > $linecount){# sees whether the total linenumber is actually smaller than what was asked for
			$readrows = $linecount;   #!!!!!!!!all files has less lines than input should reset
			}
            print "Input seq file contains $readrows lines.\n";
        return $readrows;

}

sub formatCheck_and_liglength #($filename,(yes to check if only fastq file required))
{
		my $filename = shift;
		my $fastq_check = shift || 0;
		
		if (!($filename =~ m/\.txt$/i) && !($filename =~ m/\.fastq$/i)) { print "!! $filename is not reads"; return;} #skip non-reads file
		#open the filehandle and count the length of read, then close
		open(FILENAME, $filename) or die "failed to open filehandle in &formatCheck_and_liglength()";
			if($fastq_check) #if true, check whether the file is fastq file
			{
				my $readstmp1 = <FILENAME>;#first is Illuminatus header
		 		if(not $readstmp1=~ m/^\@HISE/)	{ print "$filename is not a fastq file"; return; }
			}
				my $readstmp2 = <FILENAME>;# is the actual sequence
				chomp $readstmp2; #get rid of return
				my $liglength= length($readstmp2);
				
										 
		close FILENAME;
		return $liglength;
		 #####################################################################
}
		 
		 
		 
		 
1;	