#!/usr/bin/perl -w
# to do motif mapping for FJ4.4CAP selex with DBD
# count at the center (odd) or center 2 bp (even) of motif

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
#use Me::Nu_cnt::checkfile;
#use POSIX qw(ceil strftime);
use File::Spec qw(rel2abs);
use Mylib1::IUPAC_code;
use POSIX qw/ceil/;

use Bio::Perl;
use Bio::Seq;
use Me::MotifMap::CntArray;
use Mylib1::seqFile;
use MOODS; use MOODS::Tools qw(printResults numResults readMatrix);

#my $cutoff=30; # only take ori or rc motif <$cutoff from end

		
######### global variables to adjust #############
our @allGuides= (
				 "/var/www/kaz/guide/all-human-version4.txt" ,
				 "/var/www/kaz/guide/Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper.txt",
				 "/var/www/kaz/guide/MethylSELEX_v0.7_FinalCall_without_fail.txt",
				 "/var/www/kaz/guide/Methylation_SELEX_Version_0.5.txt",
				 #"/var/www/kaz/guide/FJ_shortest_motif.txt"
			  );
######### usage
my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme -f 'file' -u -mono -k 4 -t 20 -o dimerFFT

	***general***
	-f      files to process
	-o      output path
	-t      threads
	-h      help
	
	***specific***
	-p      p value for Moods
	-l      lines to process per file
	-u      remove duplicated reads
	-k      kmer length
	-mono   prefer using monomer PWM
	-lenCut only consider dimer gaps > lenCut
	";

########## Info from cmd line #################
	GetOptions ('f=s' => \my $files, 't=i' => \my $threads, 'o=s' => \my $outputPath, 'p=f' => \my $p_value, 'l=i' => \my $lines_to_process, 'k=i' => \my $k, 'h' => \my $helpInfo,
				'u' => \my $rmdup, 'test' => \my $test, 'mono' => \my $mono, 'lenCut=i' => \my $lenCut ); # i:int   f:float   s:string    // 'm=s' => \my $matrices
	if ($helpInfo||(!$files)) { print "\nusage: $usage\n";exit;}

	# set default value
	$outputPath||="./$this_programme";
	$threads||=15;
	$p_value||=0.0001; 

	#output cmd info
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath/data/); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd
		
	# packing param
	my @files=glob($files);	chomp @files;
	my @param; for (0..$#files){ push @param, [$files[$_]];} 
	

# -----------
use Mylib1::Parallel::file;
&file::mProcess($threads,\&singleProc,\@param); # using mThread necessary for debug !!
#&file::mThread($threads,\&singleProc,\@param);

chdir $outputPath;
#system ('Rscript /wrk/data/fangjie/lib/Rlib/motifMap/nucleosome_score.R');
#system ("convert -quality 95 ./PNG/* ./TFMotifMap_p_$p_value.pdf");

exit;

# *** multi proc ***
sub singleProc  
{
	my $param_Aref = shift;
	my ($file)= @$param_Aref;
	
	my @allreads; use Mylib1::seqFile; &seqFile::getallSeq(\@allreads,$file,$lines_to_process);
	&seqFile::rmdup_3point(\@allreads) if ($rmdup && (!$test));
	my $combinedReads=join "",@allreads; 
	
	    # figure out ligand length and cut_value for each PWM
    my $liglength= length $allreads[0]; print "length of the reads is $liglength bp for $file \n";
	my $totalReads=@allreads; @allreads=[];
	
	#----------- parse the file name to find out TF name and wells -------------
	$file =~ m/TF.*-(?=[^\d]{1})(.{1,15})IIIc/; my $TFname=$1; my $TFnoNum=$TFname=~s/[0-9\-L]*$//ir;
	$file =~ m/-([A-P]{1}\d{1,2})-/; my $well=$1;   #(?=[^\d]{1}) to make extration like TF105-NKX2-4IIIc correct
	$well=~ s/([A-P])0(\d)/$1$2/;  # to change A01 into A1 for SELEX
	
	#----------- find out pfm file to use from TF name -----
	my @matrix_files; my @matrix_Names;
	#if (!$test)
	#{	
		if ($mono) {
			&pfmFromGuide($_,$TFname, \@matrix_files, \@matrix_Names,"monomeric") for @allGuides; # first the monomeric
			&pfmFromGuide($_,$TFnoNum, \@matrix_files, \@matrix_Names,"monomeric") for @allGuides; 
		}
		&pfmFromGuide($_,$TFname, \@matrix_files, \@matrix_Names) for @allGuides; # first the exact hit
		&pfmFromGuide($_,$TFnoNum, \@matrix_files, \@matrix_Names) for @allGuides; # then everything related
		
		if ($#matrix_files>19) {splice(@matrix_files,19,-1) ;splice(@matrix_Names,19,-1); } # cannot plot too many
		if ($#matrix_files== -1) {print "\n!!!!!!!!!!!!!!\n!!!!!!!!!\nno PWM found for $file\n\n"; return;}; # do not calc if no matrix found
	#}
	
	my @matrix_Names_ori= @matrix_Names;
	for my $nameAdd (qw(r_ rc_ c_)) {push @matrix_Names, $nameAdd.$matrix_Names_ori[$_] for (0..$#matrix_Names_ori)}; #construct names for reverse, cmp, rev_cmp      
	my $Num_of_Matrix=scalar(@matrix_files);
	my (@matrices, @rc_matrices, @r_matrices, @c_matrices, @all_matrices,          @Motif_L_of_matrices, @thresholds);
	for (my $i=0;$i<$Num_of_Matrix;$i++)
	{
		$matrices[$i] = readMatrix($matrix_files[$i]);
		$c_matrices[$i] = [reverse (@{$matrices[$i]})]; # reverse fun result in cmp
	
		$rc_matrices[$i] = MOODS::Tools::reverseComplement($matrices[$i]);
		$r_matrices[$i] = [reverse (@{$rc_matrices[$i]})];
	
		$Motif_L_of_matrices[$i] = scalar(@{$matrices[$i][0]});
	}
                    #signal       #control
	@all_matrices = (@matrices, @r_matrices, @rc_matrices, @c_matrices); #join all matrices together to get faster mapping,   //@rc_matrices, @c_matrices not used here, pay attention to MapSR indices of @{$results[$i+$Num_of_Matrix]}
	@thresholds= (($p_value) x scalar(@all_matrices)); #()surronding p_value is important
	
    my @cut_values; # account for blank positions
    for (my $i=0; $i<$Num_of_Matrix;$i++) {$cut_values[$i] = $liglength-$Motif_L_of_matrices[$i]+1;}    
    

	my %results;
	&MapSR($combinedReads,$liglength,\@cut_values,\@all_matrices,\@thresholds,$Num_of_Matrix,\@matrix_Names,\%results, \@Motif_L_of_matrices);
	
	#use Mylib1::fileIO; my $outfilename=$well."_".&fileIO::getfile($file);
	my $outfilename=$TFname;	
	
	&writeDimerGapSeqs(\%results, $outputPath, $outfilename.".xls") ;
	map {&writeDimerGapSeqs_dir(\%results, $outputPath, $_."_".$outfilename.".txt",$_) } qw(HT TT HH);
	&writeDimerGapSeqs_dir_all(\%results, $outputPath, "all_".$outfilename.".txt")
	
	#print"";
	# gap pos rel to most flanking positions (lig rotated) !!!!!!!!!!!!
	#&write_dimer_gap_rel_1st_flanking(\%results, $outputPath, $outfilename.".xls");
	
	## gap pos to all abs positions (lig rotated) !!!!!!!!!!
	#&write_dimer_gap_abs_all_pos(\%results, $outputPath, $outfilename, $liglength);
	#
	#&write_kcnt_adjusted_align(\%results, $outputPath, $outfilename,$liglength);
	#
	#my $Rcmd="~/R/R320/bin/Rscript ~/lib/Rlib/FJ4.4/dimerGapNusig.R $outputPath/data/dimer_spacing_rel/$outfilename.xls '~/seqFiles/FAC_PE_EMSA_FJ003_0,25/FAC_Cyc4/lane".
	#			"8_NoIndex_L008_Assembly.fastq.assembled.fastq_barcodes/JustReads/FAC_14704GGTC94N_4.txt' $outputPath/";
	##print $Rcmd."\n\n";
	#system ($Rcmd);
}

sub writeDimerGapSeqs_dir_all
{
	my ($results_Href, $outputPath, $outfilename, $direction)=@_;
	my $dataPath="$outputPath/data/dimerGasSeqs/";
	qx(mkdir -p $dataPath) if (!(-d $dataPath));
	open(my $fh, ">", $dataPath.$outfilename) or die "cannot open output file";

	foreach (keys $results_Href) #PWM name
	{
		my $temp=$results_Href->{$_};
		foreach (keys $temp->{dimerGaps})
		{
			my $temp = $temp->{dimerGaps}->{$_};
			foreach (@$temp) #dimerGaps
			{
							next unless $_;
							print $fh $_."\n";	
			}
		}
	}
	close $fh;
}

sub writeDimerGapSeqs_dir
{
	my ($results_Href, $outputPath, $outfilename, $direction)=@_;
	my $dataPath="$outputPath/data/dimerGasSeqs/";
	qx(mkdir -p $dataPath) if (!(-d $dataPath));
	open(my $fh, ">", $dataPath.$outfilename) or die "cannot open output file";

	foreach (keys $results_Href) #PWM name
	{
		my $temp=$results_Href->{$_};
		foreach (@{$temp->{dimerGaps}->{$direction}}) #dimerGaps
		{
						next unless $_;
						print $fh $_."\n";	
		}
	}
	close $fh;
}

sub writeDimerGapSeqs
{
	my ($results_Href, $outputPath, $outfilename)=@_;
	my $dataPath="$outputPath/data/dimerGasSeqs/";
	qx(mkdir -p $dataPath) if (!(-d $dataPath));
	open(my $fh, ">", $dataPath.$outfilename) or die "cannot open output file";
	print $fh "PWM	Direction	seq\n";

	foreach (keys $results_Href) #PWM name
	{
		my $temp= $results_Href->{$_}; my $prefix= ""; $prefix.=$_."\t";
		foreach (keys $temp->{dimerGaps}) #dimerGaps
		{
			my $temp= $temp->{dimerGaps}->{$_}; my $prefix=$prefix.$_."\t";
			foreach (@$temp) # direction
			{
						my $temp= $_; next unless $temp;
						#my $prefix=$prefix.$_."\t";
						print $fh $prefix.$temp."\n";
			}		
		}
	}
	close $fh;
}



sub write_kcnt_adjusted_align
{
	my ($results_Href, $outputPath, $outfilename, $liglength)=@_;
	qx(mkdir -p $outputPath/data/kcnt_adj_align/) if (!(-d "$outputPath/data/kcnt_adj_align/"));
	open(my $fh, ">", $outputPath."/data/kcnt_adj_align/".$outfilename.".xls") or die "cannot open output file";
	my $posInfo= join "\t", ((-$liglength+1)..($liglength-1)); print $fh "\t\t\t".$posInfo."\n";
	
	foreach my $PWM_name (keys $results_Href)
	{

		map { print $fh $PWM_name."	ori_1st	".$_."\t".(join "\t", @{$results_Href->{$PWM_name}->{k_cnt}->{ori_1st}->{$_}})."\n"; } keys $results_Href->{$PWM_name}->{k_cnt}->{ori_1st};

		map { print $fh $PWM_name."	rc_1st	".$_."\t".(join "\t", @{$results_Href->{$PWM_name}->{k_cnt}->{rc_1st}->{$_}})."\n"; } keys $results_Href->{$PWM_name}->{k_cnt}->{rc_1st};

		map { print $fh $PWM_name."	r_c_1st	".$_."\t".(join "\t", @{$results_Href->{$PWM_name}->{k_cnt}->{r_c_1st}->{$_}})."\n"; } keys $results_Href->{$PWM_name}->{k_cnt}->{r_c_1st};

	}
	close $fh;
}

sub write_dimer_gap_rel_1st_flanking
{
	my ($results_Href, $outputPath, $outfilename)=@_;
	qx(mkdir -p $outputPath/data/dimer_spacing_rel/) if (!(-d "$outputPath/data/dimer_spacing_rel/"));
	open(my $fh, ">", $outputPath."/data/dimer_spacing_rel/".$outfilename) or die "cannot open output file";
	foreach my $PWM_name (keys $results_Href)
	{
		print $fh "FF_".$PWM_name."\t".(join "\t", @{$results_Href->{$PWM_name}->{dimer_spacing}->{ori_1st}->{oriRel}})."\n";
		print $fh "FR_".$PWM_name."\t".(join "\t", @{$results_Href->{$PWM_name}->{dimer_spacing}->{ori_1st}->{rcRel}})."\n";
		print $fh "RF_".$PWM_name."\t".(join "\t", @{$results_Href->{$PWM_name}->{dimer_spacing}->{rc_1st}->{oriRel}})."\n";
		print $fh "RR_".$PWM_name."\t".(join "\t", @{$results_Href->{$PWM_name}->{dimer_spacing}->{rc_1st}->{rcRel}})."\n";
	}	
	close $fh;
}



sub MapSR
{
    my ($combinedReads,$liglength,$cut_values_Ref,$all_matrices_Aref,$thresholds_Aref,$Num_of_Matrix,$matrix_Names_Aref,$results_Href,$Motif_L_of_matrices)=@_;  #liglength, @matrices, @Motif_L_of_matrices should be passed if not in the same file
	my @all_matrices = @$all_matrices_Aref;
	my @thresholds= @$thresholds_Aref;
	my $seq = Bio::Seq->new(-seq  => $combinedReads, -alphabet => 'dna' );

        ## We convert the count matrix into log-odds scores in base 2 logarithm and find hits scoring
        my @results = MOODS::search(-seq => $seq, -matrices => \@all_matrices,  -thresholds => \@thresholds, -threshold_from_p => 1);

        #add SR count into countArrays of each PWM
        for (my $i=0; $i<$Num_of_Matrix;$i++)
        {
            my $currMotifLen=$Motif_L_of_matrices->[$i];
			
			#collect mapping hits
			my ($positions,$positions_rc,$positions_c,$positions_r); 
            {my $j = 0; $positions= [grep {!($j++ % 2)} @{$results[$i]}];                     $positions = &CntArray::getRealPosition($positions,$cut_values_Ref->[$i],$liglength);};
            {my $j = 0; $positions_rc =[grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix*2]}]; $positions_rc = &CntArray::getRealPosition($positions_rc,$cut_values_Ref->[$i],$liglength);};            
            {my $j = 0;  $positions_c =[grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix*3]}]; $positions_c = &CntArray::getRealPosition($positions_c,$cut_values_Ref->[$i],$liglength);}; 
			{my $j = 0;  $positions_r =[grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix]}];   $positions_r = &CntArray::getRealPosition($positions_r,$cut_values_Ref->[$i],$liglength);};
			
			use POSIX "fmod"; use List::Util qw(min max);
			my %matches;
			
			parsePosSeq(\%matches,$positions,$liglength,"ori",$combinedReads);
			parsePosSeq(\%matches,$positions_rc,$liglength,"rc",$combinedReads);
			parsePosSeq(\%matches,$positions_r,$liglength,"r",$combinedReads);
			parsePosSeq(\%matches,$positions_c,$liglength,"c",$combinedReads);
			
			my @dimers= grep { @{$matches{$_}->{ori}||[]}+@{$matches{$_}->{rc}||[]} >1 } keys %matches; # !!!!!!!!! maybe still monomer if palindromic
			
			my $dimerGapSeqs= {HH=>[], TT =>[], HT=>[]}; # array for storing sequences
			extDimerGap(\%matches, \@dimers, $dimerGapSeqs, $currMotifLen,$liglength,$lenCut);
			fillN($dimerGapSeqs, $liglength);
			
			

			#@dimers=@matches{@dimers}; # collect ref from @matches ({ ori => ["49.5"], rc => ["51.5"] }, ...,  ...)
			$results_Href->{$matrix_Names_Aref->[$i]}={dimerGaps => $dimerGapSeqs};
			
        }
}

sub fillN   # fill end with N to make all ligands with the same length for FFT
{
	my ($dimerGapSeqs, $liglength)= @_;
	foreach (keys $dimerGapSeqs)
	{
		foreach (@{$dimerGapSeqs->{$_}})
		{
			if (length($_)<$liglength)
			{
				$_.= "N" x ($liglength-length($_));
			}
		}
	}
}

sub extDimerGap
{
	my ($matches, $dimers, $dimerGapSeqs, $motifL,$liglength, $lenCut)=@_;
	foreach (@$dimers)
	{
		my $temp= $matches->{$_};
		my $oriHitNo=scalar(@{$temp->{ori}||[]});
		my $rcHitNo=scalar(@{$temp->{rc}||[]});
		
		if ($rcHitNo==0 || $oriHitNo==0)
		{
			rotate($temp, $liglength) if $oriHitNo==0;			
			my $start= $temp->{ori}->[0]+($motifL-1)/2;
			my $len= $temp->{ori}->[1]-$temp->{ori}->[0]-1;
			next if $len < $lenCut;
			push @{$dimerGapSeqs->{HT}}, substr($temp->{seq}, $start, $len); 
		}
		
		elsif($temp->{ori}->[0] < $temp->{rc}->[0])
		{
			my $start= $temp->{ori}->[0]+($motifL-1)/2;
			my $len= $temp->{rc}->[0] - $temp->{ori}->[0]-1;
			next if $len < $lenCut;
			push @{$dimerGapSeqs->{HH}}, substr($temp->{seq}, $start, $len);
		}
		
		elsif($temp->{ori}->[0] > $temp->{rc}->[0])
		{
			my $start= $temp->{rc}->[0]+($motifL-1)/2;
			my $len= $temp->{ori}->[0] - $temp->{rc}->[0]-1;
			next if $len < $lenCut;
			push @{$dimerGapSeqs->{TT}}, substr($temp->{seq}, $start, $len);
		}
	}
}

sub cnt_gap_absPos_1_lig 
{
	my ($curr,$gap_absPos,$motif_l)=@_;
	foreach my $dir1 (qw(ori rc))
	{ foreach my $dir2 (qw(ori rc))
	 {	foreach my $pos1 (@{$curr->{$dir1}})
		{ foreach my $pos2 (@{$curr->{$dir2}})
		 {
			next if $pos1==$pos2;
			my $posAbs= ($motif_l % 2 ? $pos1 : $pos1-0.5); # if even motif_l then -0.5bp in starting pos
			my $gapL=$pos2-$pos1;
			my $dir1t= ($dir1 eq "ori")? "ori_1st": "rc_1st";
			my $dir2t= ($dir2 eq "ori")? "oriRel": "rcRel";
			$gap_absPos->{$dir1t}->{$dir2t}->{$posAbs}->{$gapL}++; 
		  }
		}
	  }
	}
}
sub cnt_adj
{
	my ($matches,$adj_k_cnt,$liglength)= @_;
	foreach (keys $matches)
	{
		
		my $temp= $matches->{$_}; my $kmerNo=length($temp->{seq})-($k-1);
		my @kmers= map {substr($temp->{seq},$_,$k)} (0..$kmerNo-1);
		my $dir= $temp->{is_r_c}? "r_c_1st": ($temp->{minDir} eq "rc"? "rc_1st": "ori_1st");
		map { $adj_k_cnt->{$dir}->{$kmers[$_]}->[ ceil( ($_+($k+1)/2) - $temp->{minDist} )+ $liglength-1]++ if $kmers[$_] !~ m/N/; } (0..$kmerNo-1); # No. $liglength of the arr is rel zero
	}
}

sub parsePosSeq
{
	my ($matches,$positions,$liglength,$direction,$combinedReads)=@_;
	foreach (@$positions)
	{
		my $ligNo=int($_/$liglength);
		my $actualPos=fmod($_, $liglength)+1; # start pos is 1 now, not 0
		#my $minDistToEnd= min($actualPos-1,$liglength-$actualPos)
		$matches->{$ligNo}->{$direction}||=[];
		push @{$matches->{$ligNo}->{$direction}}, $actualPos;
		if (!$matches->{$ligNo}->{seq}) 
		{
			$matches->{$ligNo}->{seq}=substr($combinedReads,$liglength*$ligNo,$liglength)
		}	
	}
}

sub autoRotate
{	# rotate if the dist to the right is smaller than dist to the left (both ori and rc)
	my ($matches,$liglength)=@_;
	foreach (keys $matches)
	{
		my $minLeft; my $minRight;
		if ((!$matches->{$_}->{ori}) && (!$matches->{$_}->{rc}))
		{
			$matches->{$_}->{is_r_c}=1; # mark if only rev or cmp hit
			$minLeft= (min @{$matches->{$_}->{r}}, @{$matches->{$_}->{c}})-1; # pos 1 and 101 is of dist 0
			$minRight= $liglength-(max @{$matches->{$_}->{r}}, @{$matches->{$_}->{c}});
		} else
		{		
			$minLeft= (min @{$matches->{$_}->{ori}}, @{$matches->{$_}->{rc}})-1; # pos 1 and 101 is of dist 0
			$minRight= $liglength-(max @{$matches->{$_}->{ori}}, @{$matches->{$_}->{rc}});
		}
		$matches->{$_}->{minDist}= $minRight<$minLeft? $minRight+1: $minLeft+1; # change 0 to 1
		if ($minRight<$minLeft)
		{
			rotate($matches->{$_},$liglength);
			$matches->{$_}->{rotated}=1;
		}
		if (!$matches->{$_}->{is_r_c}) # check whether the smallest dist is ori or rc
		{
			my $min_ori= min (@{$matches->{$_}->{ori}})||(100000); #should not omit () from the arr otherwise it will be a scalar
			my $min_rc= min (@{$matches->{$_}->{rc}})||(100000);
			if ($min_ori<=$min_rc) {$matches->{$_}->{minDir}="ori";} else {$matches->{$_}->{minDir}="rc";}
		}
	}
}

sub rotate
{
	my ($ligHref,$liglength)=@_;
	my $tmp= $ligHref->{ori}; $ligHref->{ori}= $ligHref->{rc}; $ligHref->{rc}= $tmp;
	   $tmp= $ligHref->{r}; $ligHref->{r}= $ligHref->{c}; $ligHref->{c}= $tmp;
	foreach my $dir ("ori","rc","r","c") {foreach (@{$ligHref->{$dir}}) {$_= $liglength-$_+1;}} #calc the new pos
	$ligHref->{seq}= &IUPAC_code::reverse_complement_IUPAC($ligHref->{seq});
}

sub addCnt
{
    my ($countArrayRef, $matches,$cut_value,$liglength,$direction)=@_;

    foreach my $key (keys $matches)
    {
		next if (!$matches->{$key}->{$direction});
		foreach (@{$matches->{$key}->{$direction}})	{ $countArrayRef->{$_}++; }
    }  
}

sub pfmFromGuide #extract pfm path from guide
{
	my ($guidefile, $search_string, $matrix_files_Aref,$matrix_Names_Aref,$monomeric)=@_;
	my $string= "cat $guidefile | grep -P -i '^$search_string\[\\d-L]*\\t' | "; #some TF name like HOXB2L2 NKD2-4
	$string.= "grep 'mono' | " if $monomeric;
	$string.=q{perl -lne '@all=split /\t/; $pfm_key="$all[2]_$all[3]_$all[4]"; @allpfm= glob("/var/www/kaz/data/pfm/shortspacek/$pfm_key*.pfm"); chomp @allpfm; print "$all[0]_$all[3]_$all[4]_$all[7]_$all[8]"."\t".$allpfm[0] if @allpfm;'}; # some like TAGTTA40NGTG_XEPAQ_NNNNTGCTGAC
	my @matrix_info_tmp=`$string`; chomp @matrix_info_tmp;
	map {my @info=split/\t/; push @$matrix_Names_Aref, $info[0]; push @$matrix_files_Aref, $info[1];} @matrix_info_tmp;
}

sub normalize
{
	my ($Aref,$totalReads)=@_;
	foreach my $Href (@$Aref)
	{
		#if (ref($_)) {&normalize($_,$totalReads) } else{$_/=$totalReads};
		foreach (keys $Href) {$Href->{$_}/=$totalReads;}
	}
}