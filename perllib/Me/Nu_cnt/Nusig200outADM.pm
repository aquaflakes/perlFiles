package Nusig200outADM;


use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../..";  # use the 2 lv up parent directory
use warnings;
use strict;
use perllib::Me::Nu_cnt::ADM;

sub Nusig200outADM
{	
	my ($temp_storage_ref,$reads2,$liglength,$filterDataRef)=@_;
	
	my %nucleotide_cnt=%{$temp_storage_ref}; #initialized ref with single nucleotide keys
	my %admFreq = %{$filterDataRef};
	
	
				my $p_value=-10000000000; my $start_pos;
				for (my $i = 0; $i < ($liglength-94+1); ++$i) #$i has info about the starting position of match
					{
						my $sublig_94 = substr($reads2,$i,94); #extract single nucleotide signal for a specific position
						if ($sublig_94=~ m/N/) {next;} #no "N" containing reads
						my $p_value_tmp=$admFreq{substr($sublig_94,0,1)}->[0]; #entering probability, single nucleotide freq
												
						for (my $j = 0; $j < 93; ++$j) {
							
							
							
							
							my $di_nt = substr($sublig_94, $j, 2);
							$p_value_tmp += $admFreq{$di_nt}->[$j];
						}					
						if ($p_value_tmp>$p_value) { $p_value=$p_value_tmp; $start_pos=$i;} # choose the most possible position
					}
					
							#if (!defined($start_pos)) { 
							#	print $reads2."   reads2\n";
							#	print $p_value."   p_value\n";
							#}
				if (!defined ($start_pos)) {return \%nucleotide_cnt;}								
				my $prior = 53-$start_pos;
				my $reads_shift=$reads2;
				$reads_shift=~s/^(.*)/'N' x $prior . $1/e; # add N to shift reads for alignment
				$reads_shift=~s/(.*)$/$1 .'N' x $start_pos/e; # add N to shift reads for alignment
				
					
	my $nucleotide_cnt_ref = &ADM::addSR_simple_align (\%nucleotide_cnt,$reads_shift);	
				
		return $nucleotide_cnt_ref;
				##di_nt in adm position -1
				#if ($start_pos_win[0]!=0)
				#{
				#	my $win_di_nt_head = substr($reads2,$start_pos_win[0]-1,2); #extract dinucleotide signal for a specific position
				#	if (!($win_di_nt_head=~ m/N/)) {$nucleotide_cnt{$win_di_nt_head}->[0]++;}
				#}
				##di_nt in adm position 94
				#if ($start_pos_win[$win_num-1]!=53)
				#{
				#	my $win_di_nt_tail = substr($reads2,$start_pos_win[$win_num-1]+93,2); #extract dinucleotide signal for a specific position
				#	if (!($win_di_nt_tail=~ m/N/)) {$nucleotide_cnt{$win_di_nt_tail}->[94]++;}
				#}
				
	
			
			
			
			
			
			
			## replace all the undefined value with 0 (there should not be any if reads number is large enouth)
			#foreach my $keynt (keys(%nucleotide_cnt))
			#{	
			#	for (my $i = 0; $i < 95; ++$i)
			#	{
			#	if (!defined ($nucleotide_cnt{$keynt}->[$i])) {$nucleotide_cnt{$keynt}->[$i]=0;}
			#	}
			#}
			#
			## out put di_nt_cnt in adm to files
			#foreach my $keynt (sort keys(%nucleotide_cnt))
			#{	
			#print out1 $keynt."\t";
			#for (my $i = 0; $i < 95; ++$i) {
			#	print out1 $nucleotide_cnt{$keynt}->[$i]."\t";
			#}
			#print out1 "\n";
			#}
			#
			#close out1;

		
	
}

	
	
	
	
	
	
	
1;
	
	