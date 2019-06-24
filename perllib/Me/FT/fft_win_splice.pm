package fft_win_splice;
use warnings;
use strict;
use Data::Dumper;
use List::Util qw(sum);
use Storable 'dclone';
use Math::Complex;



sub ftADM #(ADM containing frequencies)
{
    my %ADM_freq=%{dclone shift};
    my (%fftPower,%peakPercent,%phaseAngle,%peakArea);
    foreach my $keynt (keys(%ADM_freq))
    {
        $fftPower{freq}=\my @emptyarr; #storage for freq
        $fftPower{$keynt}=\my @emptyarr1; # storage for amplitude
        #!!!!!remove the first 5 point and the end 5 point with high variation in bk
        ##########################
        my $splice_l = (length $keynt)-1;
        my $frontSpl=0; my $endSpl = -(0 + $splice_l); if($endSpl==0) {$endSpl=1000;} #adjust the input length to be the same for multi nt
        ##endSpl =0 will cause problem so if endSpl=0 then take the whole seq using 1000
        
        my @spliced= splice(@{$ADM_freq{$keynt}}, 0,-1);#$frontSpl, $endSpl; #return subarray from 4 to (end-4)
#        #print Dumper(@spliced);
#		print "\n".scalar(@{$ADM_freq{$keynt}})."\n";
#		print "\n".scalar(@spliced)."\n";
#		print "\n".$frontSpl."\n";
#		print "\n".$endSpl."\n"; #die;
        ##########################
        # substract mean freq before calc
        my $meanOfKey= &mean(@spliced);
        for (my $i = 0; $i < scalar(@spliced); ++$i)
				{
					$spliced[$i] = $spliced[$i]-$meanOfKey;
                    $spliced[$i] = $spliced[$i]* &win_Welch(scalar(@spliced),$i); # total_points, value_position
                }
                
      
        ####### padding with 0 ############
        #my @x1=@spliced;
        #my $n=@x1; #n input
        #my $m=log($n)/log(2); #FFT padding
        #my $k= int($m) + ($m != int($m)); #round up for power
        #my $np=2**$k;   #n after padding (for FFT)
        ### fill undef data in array with 0 , padding for FFT ##
        #my @f_x;
        # for (my $i=0; $i< $np; $i++)
        #{
        #    push @f_x,((defined $x1[$i])? 1*$x1[$i] : 0); 
        #}
        ### fill undef data in array with 0 ##
        
        my @f_x=@spliced; #FT
        my $nx=@f_x;

        ### C(h) = sigma(x=0..nx-1) f_x* Exp(i*(2(pi)h/c)*x)
        ###      = sigma(x=0..nx-1) f_x*cos((2(pi)h/c)*x) + i sin((2(pi)h/c)*x)
        my $pi=3.14159265;
        my (@xre,@xim,@Power,$peakArea,$peakBk);
        for ( my $h = 0; $h < 600; $h++ ) # should be 150 to 251
        {
            #real part
            $xre[$h]=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xre[$h]+=$f_x[$x]*cos(2*$pi*$h*$x/2000); #frequency: $h/2000
            }
            
            #imaginary part
            $xim[$h]=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xim[$h]-=$f_x[$x]*sin(2*$pi*$h*$x/2000);
            }
             #power spectra
            $Power[$h]=sqrt((($xre[$h]*$xre[$h])+($xim[$h]*$xim[$h])));
            #push @Px, $P[$h];
            
            my $freq = $h/2000;
            my $normalized_power=$Power[$h]/($meanOfKey*5.1*scalar(@f_x)/10.2); #!!!!!!!!!! devide by mean freq to normalize
            push $fftPower{freq}, $freq;
            push $fftPower{$keynt}, $normalized_power;
            
            if ($freq <= 0.1020408 && $freq >= 0.09434){$peakArea+= $normalized_power;} # take peak area from 9.8 to 10.6 bp
            if ($freq <= 0.125 && $freq >= 0.074074){$peakBk+= $normalized_power;} # take peak area bk from 8 to 13.5 bp
           
            
        }
        $peakArea{$keynt}=[$peakArea];     $peakArea{keys}=["abs area of 10.2bp peak"]; # add cap line
        $peakPercent{$keynt}=[$peakArea/$peakBk];    $peakPercent{keys}=["area% of 10.2bp peak"]; # add cap line
        my $cplx= $xre[196] + $xim[196]*i;
        my $phaseAng_h196=arg($cplx)/$pi;# in units of pi (-1...1 output)
        my $phaseAngle_h196_corr=$phaseAng_h196-2*(28/10.204); # correct for the phase of adaptor at the start
        $phaseAngle{$keynt}= [$phaseAngle_h196_corr]; 
    } 

    $phaseAngle{keys}=["phase angle of 10.2bp peak (-1..1)π"]; # add cap line
#print Dumper (%phaseAngle);die;
    return (\%fftPower,\%peakPercent,\%phaseAngle,\%peakArea);
}


#sub max {
#    splice(@_, ($_[0] > $_[1]) ? 1 : 0, 1);
#    return ($#_ == 0) ? $_[0] : max(@_);
#}

sub mean
{
    return sum(@_)/@_;
}

sub win_flatTop # total points, position of value
{
    my $tot_p=shift; my $pos=shift; my $pi=3.14159265;
    my $normFactor= 1-1.93*cos(2*$pi*$pos/($tot_p-1))+ 1.29*cos(4*$pi*$pos/($tot_p-1)) -0.388*cos(6*$pi*$pos/($tot_p-1)) + 0.028*cos(8*$pi*$pos/($tot_p-1));
    return $normFactor;
}

sub win_Blackman_Harris # total points, position of value
{
    my $tot_p=shift; my $pos=shift; my $pi=3.14159265;
    my $normFactor= 0.35875-0.48829*cos(2*$pi*$pos/($tot_p-1))+ 0.14128*cos(4*$pi*$pos/($tot_p-1)) -0.01168*cos(6*$pi*$pos/($tot_p-1)) ;
    return $normFactor;
}
sub win_Welch # total points, position of value
{
    my $tot_p=shift; my $pos=shift; my $pi=3.14159265;
    my $normFactor= 1-( ($pos-($tot_p-1)/2) / (($tot_p-1)/2) )**2;   
    return $normFactor;
}


sub ftSRfreq #(ADM containing frequencies)
{
    my %ADM_freq=%{$_[0]};
    my (%peakArea, %ntCount);
    foreach my $keynt (keys(%ADM_freq))
    {
       

        my @spliced= @{$ADM_freq{$keynt}}; #return subarray from 4 to (end-4)
        $ntCount{$keynt}=sum(@spliced);
        #print Dumper(@spliced);die;
        ##########################
        # substract mean freq before calc
        my $meanOfKey= &mean(@spliced);
        for (my $i = 0; $i < scalar(@spliced); ++$i)
				{
					$spliced[$i] = $spliced[$i]-$meanOfKey;
                    #$spliced[$i] = $spliced[$i]* &win_Welch(scalar(@spliced),$i); # total_points, value_position
                }
                
      
        ####### padding with 0 ############
        #my @x1=@spliced;
        #my $n=@x1; #n input
        #my $m=log($n)/log(2); #FFT padding
        #my $k= int($m) + ($m != int($m)); #round up for power
        #my $np=2**$k;   #n after padding (for FFT)
        ### fill undef data in array with 0 , padding for FFT ##
        #my @f_x;
        # for (my $i=0; $i< $np; $i++)
        #{
        #    push @f_x,((defined $x1[$i])? 1*$x1[$i] : 0); 
        #}
        ### fill undef data in array with 0 ##
        
        my @f_x=@spliced; #FT
        my $nx=@f_x;

        ### C(h) = sigma(x=0..nx-1) f_x* Exp(i*(2(pi)h/c)*x)
        ###      = sigma(x=0..nx-1) f_x*cos((2(pi)h/c)*x) + i sin((2(pi)h/c)*x)
        my $pi=3.14159265;
        my ($xre,$xim,$Power,$peakArea,$peakBk);
       
        
            #real part
            $xre=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xre+=$f_x[$x]*cos(2*$pi*196*$x/2000); #frequency: 196/2000
            }
            
            #imaginary part
            $xim=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xim-=$f_x[$x]*sin(2*$pi*196*$x/2000);
            }
             #power spectra
            $Power=sqrt((($xre*$xre)+($xim*$xim)));
            #push @Px, $P;
            
            #my $freq = 196/2000;


            
            #if ($freq <= 0.1020408 && $freq >= 0.09434){$peakArea+= $normalized_power;} # take peak area from 9.8 to 10.6 bp
            #if ($freq <= 0.125 && $freq >= 0.074074){$peakBk+= $normalized_power;} # take peak area bk from 8 to 13.5 bp
           
            
        
        $peakArea{$keynt}=$Power;    # $peakArea{keys}=["abs area of 10.2bp peak"]; # add cap line
        #$peakPercent{$keynt}=[$peakArea/$peakBk];    $peakPercent{keys}=["area% of 10.2bp peak"]; # add cap line
        #my $cplx= $xre + $xim*i;
        #my $phaseAng_h196=arg($cplx)/$pi;# in units of pi (-1...1 output)
        ##my $phaseAngle_h196_corr=$phaseAng_h196-2*($frontSpl/10.204); # correct for the phase of splicing at the start
        #$phaseAngle{$keynt}= [$phaseAng_h196]; 
    } 

   # $phaseAngle{keys}=["phase angle of 10.2bp peak (-1..1)π"]; # add cap line
#print Dumper (%phaseAngle);die;
    return (\%peakArea,\%ntCount);#\%peakPercent,
}

#use PDL;
#use PDL::FFT qw(realfft);
sub ftSRfreqPDL #(ADM containing frequencies)
{
    my %ADM_freq=%{$_[0]};
    my (%peakHeight, %ntCount, %ntCntSqr);
    foreach my $keynt (keys(%ADM_freq))
    {
       

        my @spliced= @{$ADM_freq{$keynt}}; #return subarray from 4 to (end-4)
        $ntCount{$keynt}=sum(@spliced); #$ntCntSqr{$keynt}=sum(@spliced**2);
        #print Dumper(@spliced);die;
        ##########################
        # substract mean freq before calc
        # my $meanOfKey= &mean(@spliced);
#        for (my $i = 0; $i < scalar(@spliced); ++$i)
#				{
#					$spliced[$i] = $spliced[$i]-$meanOfKey;
#                    #$spliced[$i] = $spliced[$i]* &win_Welch(scalar(@spliced),$i); # total_points, value_position
#                }
                
      
        ####### padding with 0 ############
        #my @x1=@spliced;
        #my $n=@x1; #n input
        #my $m=log($n)/log(2); #FFT padding
        #my $k= int($m) + ($m != int($m)); #round up for power
        #my $np=2**$k;   #n after padding (for FFT)
        ### fill undef data in array with 0 , padding for FFT ##
        #my @f_x;
        # for (my $i=0; $i< $np; $i++)
        #{
        #    push @f_x,((defined $x1[$i])? 1*$x1[$i] : 0); 
        #}
        ### fill undef data in array with 0 ##
        
        my @f_x=@spliced; #FT
        my $nx=@f_x;

        ### C(h) = sigma(x=0..nx-1) f_x* Exp(i*(2(pi)h/c)*x)
        ###      = sigma(x=0..nx-1) f_x*cos((2(pi)h/c)*x) + i sin((2(pi)h/c)*x)
        my $pi=3.14159265;
        my ($xre,$xim,$Power,$peakArea,$peakBk);
       
        
            #real part
            $xre=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xre+=$f_x[$x]*cos(2*$pi*196*$x/2000); #frequency: 196/2000
            }
            
            #imaginary part
            $xim=0;
            for ( my $x = 0; $x < $nx; $x++ )
            {
                $xim-=$f_x[$x]*sin(2*$pi*196*$x/2000);
            }
             #power spectra
            $Power=sqrt((($xre*$xre)+($xim*$xim)));
            #push @Px, $P;
            
            #my $freq = 196/2000;


            
            #if ($freq <= 0.1020408 && $freq >= 0.09434){$peakArea+= $normalized_power;} # take peak area from 9.8 to 10.6 bp
            #if ($freq <= 0.125 && $freq >= 0.074074){$peakBk+= $normalized_power;} # take peak area bk from 8 to 13.5 bp
           
            
        
        $peakHeight{$keynt}=$Power;    # $peakArea{keys}=["abs area of 10.2bp peak"]; # add cap line
        #$peakPercent{$keynt}=[$peakArea/$peakBk];    $peakPercent{keys}=["area% of 10.2bp peak"]; # add cap line
        #my $cplx= $xre + $xim*i;
        #my $phaseAng_h196=arg($cplx)/$pi;# in units of pi (-1...1 output)
        ##my $phaseAngle_h196_corr=$phaseAng_h196-2*($frontSpl/10.204); # correct for the phase of splicing at the start
        #$phaseAngle{$keynt}= [$phaseAng_h196]; 
    } 

   # $phaseAngle{keys}=["phase angle of 10.2bp peak (-1..1)π"]; # add cap line
#print Dumper (%phaseAngle);die;
    return (\%peakHeight,\%ntCount,\%ntCntSqr);#\%peakPercent,
}


1;
 
