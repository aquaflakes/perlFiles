package Nusig200_full_largewin;

use warnings;
use strict;

sub Nusig200inADMlargewin
{
	my ($temp_storage_ref,$reads2,$liglength,$filterDataRef)=@_;
	
	my %nucleotide_cnt=%{$temp_storage_ref};
	my %admFreq = %{$filterDataRef};
	
	
	
	
	
		#initialize stuff
	   ##################
	   my $N_in_w = 40; #number of "N" in window, basically =2
	   my $win_size = $N_in_w - 1; # since containg N_in_w -1 di_nt in window
	   my $win_num = int(94/$N_in_w)+1; if (int(94/$N_in_w) == (94/$N_in_w)) {$win_num--;} #number of windows needed
	   

			  
				
			
				   my $p_value=-10000000000; my $start_pos=-100; #give a initial prob and pos for 94mer (no negtive window)
				   
				   
				   my $amdSeq_temp; #to store seq inside adm
				   my (@p_value_win, @start_pos_win); #give initial prob and pos for windows (negtive window, only count outside) inside the adm
				   
	   
				   for (my $i1 = 0; $i1 < $win_num; ++$i1)
				   {
					   $p_value_win[$i1] = -10000000000; $start_pos_win[$i1]=-100;
				   }
				   #print Dumper(@p_value_win); print Dumper(@start_pos_win); die;
				   
				   
				   
				   for (my $i = 0; $i < ($liglength-94+1); ++$i) #$i has info about the starting position of match, 94mer window go through 147mer
				   {
						   my $sublig_94 = substr($reads2,$i,94); #extract single nucleotide signal for a specific position
						   if ($sublig_94=~ m/N/) {next;} #no "N" containing reads, may skip too many for calc with windows
						   my $p_value_tmp=$admFreq{substr($sublig_94,0,1)}->[0]; #entering probability, single nucleotide freq
												   
						   for (my $j = 0; $j < 93; ++$j) {		
							   my $di_nt = substr($sublig_94, $j, 2);
							   $p_value_tmp += $admFreq{$di_nt}->[$j];
						   }					
						   if ($p_value_tmp>$p_value) { $p_value=$p_value_tmp; $start_pos=$i;} # choose the most possible position
						   
						   ############### done for nt count outside the adm model, begin with inside the adm ###########
						   my @p_tmp_win;
						   
						   # calc win[0]				
						   $p_tmp_win[0] = $p_value_tmp - $admFreq{substr($sublig_94,0,1)}->[0] + $admFreq{substr($sublig_94,$N_in_w,1)}->[$N_in_w]; # substract p of first s_nt, add p of s_nt after "N"s
						   for (my $i2=0; $i2<$N_in_w; $i2++) # p value for window 0
						   {
							   my $di_nt = substr($sublig_94, $i2, 2);
							   $p_tmp_win[0] -= $admFreq{$di_nt}->[$i2];
							   #print 10**$admFreq{$di_nt}->[$i2]."\n";							
						   }
						   
						   ## calc last win
						   #$p_tmp_win[$win_num-1] = $p_value_tmp;
						   #for (my $i2=0; $i2<$N_in_w; $i2++) # p value for window 0
						   #{
						   #	my $di_nt = substr($sublig_94, $i2+93-$N_in_w, 2);
						   #	$p_tmp_win[$win_num-1] -= $admFreq{$di_nt}->[$i2+93-$N_in_w];
						   #	#print 10**$admFreq{$di_nt}->[$i2+94-$N_in_w-1]."\n";							
						   #}
						   
						   # calc middle wins
						   for (my $iw=1; $iw<$win_num-1; $iw++)
						   {
							   $p_tmp_win[$iw] = $p_value_tmp + $admFreq{substr($sublig_94,($iw+1)*$N_in_w,1)}->[($iw+1)*$N_in_w]; # add p of s_nt after "N"s
							   for (my $i2=-1; $i2<$N_in_w; $i2++) # p value for window 0
							   {
								   my $di_nt = substr($sublig_94, $i2+$iw*$N_in_w, 2);
								   $p_tmp_win[$iw] -= $admFreq{$di_nt}->[$i2+$iw*$N_in_w];
								   #print 10**$admFreq{$di_nt}->[$i2]."\n";							
							   }				
						   }
						   
						   # compare p value, choose the most possible position for each window
						   for (my $nw=0; $nw<$win_num-1; $nw++)
						   {
							   if ($p_tmp_win[$nw]>$p_value_win[$nw]) { $p_value_win[$nw]=$p_tmp_win[$nw]; $start_pos_win[$nw]=$i;} 
						   }
																	   
				   }
						   
						   
						   
				   if ($start_pos==-100) {return \%nucleotide_cnt;}		# some ligands contains too many N so no sublig available, skip print
				   ###### print seq outside adm ######## 									
				   my $prior = 53-$start_pos; # change to liglength-94
				   my $reads_shift=$reads2;
				   $reads_shift=~s/^(.*)/'N' x $prior . $1/e; # add N to shift reads for alignment
				   $reads_shift=~s/(.*)$/$1 .'N' x $start_pos/e; # add N to shift reads for alignment
				   ###################################
				 #  print out2 $reads_shift."\n";
					   
						   #print "start_pos: ". $start_pos ."\n";
						   #	for (my $nw=0; $nw<$win_num-1; $nw++)
						   #	{
						   #		print "start_pos of win $nw is".$start_pos_win[$nw] ."\n";
						   #	}
							   
				   
				   
				   # write di_nt freq of each pos into hash %nucleotide_cnt
				   for (my $iw=0; $iw<$win_num-1; $iw++) #win_num=93 for N=2
				   {
					   for (my $i=0; $i<$N_in_w; $i++)
					   {
					   my $win_di_nt = substr($reads2,$start_pos_win[$iw]+$iw*$N_in_w+$i,2); #extract dinucleotide signal for a specific position
					   $nucleotide_cnt{$win_di_nt}->[$iw*$N_in_w+$i]++;
					   }
				   }
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
				   
	   
			   	return \%nucleotide_cnt;
			
			   
			   
	
			   


}	
	
	
	
1;
	
	