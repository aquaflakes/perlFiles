package CntArray;
use warnings;
use strict;
use Data::Dumper;


sub addCnt_to_CntSt
{
    my ($countArrayRef, $positions)=@_;
    foreach my $position (@{$positions})
    {
        $countArrayRef->[$position]++;
    }
    return $countArrayRef;    
}

sub addCnt_to_CntSt_combinedRead #count at the center
{
    my ($countArrayRef, $positions,$cut_value,$liglength)=@_;
    my @positions = @{$positions};

    for (@positions) { $_= ($_) % $liglength;} # take mod, figure out the position of motif, assume the pos start from 0 
    @positions = grep { !($_ > ($cut_value-1)) } @positions; # delete mapping hit on the border of ligand joints, for 100bp lig, 6bp motif, mapping hit on >95bp should be cut
    
	my $motif_l= $liglength+1-$cut_value; 
    foreach my $position (@positions)
    {
			$countArrayRef->{$position+($motif_l-1)/2}++; # count the position in the middle
			#print (($position+($motif_l-1)/2)."\n");
			#$countArrayRef->[$position+($motif_l-1)/2]++; # count the position in the middle
    }
}

sub getRealPosition #count at the center
{
    my ($positions,$cut_value,$liglength)=@_;
    my @positions = @{$positions};
	@positions = grep { !( ($_% $liglength) > ($cut_value-1)) } @positions;

    #for (@positions) { $_= ($_) % $liglength;} # take mod, figure out the position of motif, assume the pos start from 0 
	my $motif_l= $liglength+1-$cut_value; 
    foreach (@positions)
    {
		$_+=($motif_l-1)/2; # count the position in the middle
    }
	return \@positions;
}

sub addCnt_Test
{
    my ($countArrayRef, $positions,$cut_value,$liglength)=@_;
    #my @positions = @{$positions};
	use POSIX "fmod";

    for (@$positions) { $_= fmod($_, $liglength);} # take mod, figure out the position of motif, assume the pos start from 0 
    foreach (@$positions)
    {
			$countArrayRef->{$_}++; 
    }  
}

#sub addCnt_to_CntSt_combinedRead
#{
#    my ($countArrayRef, $positions,$cut_value,$liglength)=@_;
#    my @positions = @{$positions};
#
#    for (@positions) { $_= $_ % $liglength;} # take mod, figure out the position of motif
#    @positions = grep { !($_ > $cut_value) } @positions; # delete mapping hit on the border of ligand joints, for 100bp lig, 6bp motif, mapping hit on >95bp should be cut
#    
#    foreach my $position (@positions)
#    {
#        $countArrayRef->[$position]++;
#    }
#    return $countArrayRef;    
#}

sub dist_Cnt_KatjaE2F #(distance between two motifs on the same ligand)
{
    my ($distCntArrayRef, $positions,$cut_value,$liglength)=@_;
    my @positions = @{$positions};

    for (@positions) { $_= $_ % $liglength;} # take mod, figure out the position of motif
    @positions = grep { !($_ > $cut_value) } @positions; # delete mapping hit on the border of ligand joints, for 100bp lig, 6bp motif, mapping hit on >95bp should be cut
    #print join ("\t",@positions);print "before shift\n\n";    
    while (scalar(@positions) > 2)
    {
        my $this_pos = shift @positions;
        my $next_pos = shift @positions; # consumes the next element
        
        #print join ("\t",@positions);print "before unshift\n\n";
        if ($next_pos>$this_pos)
        {$distCntArrayRef->[$next_pos-$this_pos]++;}
        unshift (@positions, $next_pos); # retruns the next element
        #print join ("\t",@positions);print "after unshift\n\n";die;
    }
    return $distCntArrayRef;    
}

sub addTwoCntSt
{
     my ($CntSt1_Ref,$CntSt2_Ref)=@_;
     my %CntSt1=%{$CntSt1_Ref}; my %CntSt2=%{$CntSt2_Ref};
     
     my $size_dim1_St1 = keys %CntSt1;  my $size_dim1_St2 = keys %CntSt2;
     my %St_as_std = %CntSt1; if ($size_dim1_St1<$size_dim1_St2) { %St_as_std = %CntSt2; }
     
     foreach my $key (keys(%St_as_std))
     {
        for (my $i=0; $i<scalar(@{$St_as_std{$key}}); $i++)
        {
                if (!defined($CntSt1{$key}->[$i])) {$CntSt1{$key}->[$i]=0};
                if (!defined($CntSt2{$key}->[$i])) {$CntSt2{$key}->[$i]=0};
                $CntSt2{$key}->[$i] += $CntSt1{$key}->[$i];
        }
      }
        return \%CntSt2;      
}



sub wtCountFile 
{
        my ($countArrayRef,$matrix_NamesRef,$output_path,$outFilename)=@_;
        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, ">$output_path/$outFilename") || die "unable to open write file handle in CntArray::wtCountFile";
		for (my $i=0;$i<scalar(@{$countArrayRef});$i++)
        {
                print OUT1 $matrix_NamesRef->[$i]."\t";
                print OUT1 join("\t",@{$countArrayRef->[$i]});
                print OUT1 "\n";
        }
        close OUT1;
}

sub wtCountFile_Hash
{
        my ($countArrayRef,$matrix_NamesRef,$output_path,$outFilename,$isOri)=@_;
        if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open (OUT1, ">$output_path/$outFilename") || die "unable to open write file handle in CntArray::wtCountFile";
		my $motifNo = @{$countArrayRef}/4;
		for (my $i=0;$i<(@{$countArrayRef});$i++)
        {
                next if ($isOri eq "ori" && $i+1>$motifNo*2);
				next if ($isOri eq "rc" && $i+1<=$motifNo*2);
				my @keys= sort { $a <=> $b } keys $countArrayRef->[$i];
				#print OUT1 "\t".(join "\t", @keys)."\n" if $i==0;
				print OUT1 $matrix_NamesRef->[$i];
				foreach (@keys) {print OUT1 "\t".$countArrayRef->[$i]->{$_};}
                print OUT1 "\n";
        }
        close OUT1;
}



sub rdCountFile
{
	my $CntFile=shift;
	
	my @data;
   # print "here\n"; print ">$path/$adm"."\n";
	open(DATA, "$CntFile") || die "sub  rdCountFile: Can't open nt_frequency file $CntFile: $!\n";
	while (<DATA>) {
	  push @data, [split /\t/];
	}
    #print Dumper(@data); die;
	my %Cnt;
	for (my $i=0; $i<scalar(@data); $i++)
    {
        for (my $j=0; $j<(scalar(@{$data[0]})-1);$j++)
        {
                $data[$i][$j+1]=~tr/\n//d; # delete /n, do not know why chomp does not work
                $data[$i][$j+1]=~tr/\t//d; # delete /t
                $data[$i][$j+1]=~tr/,/./; # transform , to .
                $Cnt{$data[$i][0]}->[$j]= $data[$i][$j+1];
				#if ($data[$i][$j+1]==0) {$data[$i][$j+1]=0.00000000000000000001;} #give a very small numer to avoid log0
        }

	}
	close DATA;
    return \%Cnt;
}


sub sumCountFiles
{
     my ($tmp_storage_files)=@_;
     my %summedCnt;
     foreach my $file (@{$tmp_storage_files})
     {
        my %tmpStorageCnt= %{rdCountFile($file)};
        %summedCnt= %{&addTwoCntSt(\%summedCnt,\%tmpStorageCnt)}; 
     }
     return \%summedCnt;
}


1;	