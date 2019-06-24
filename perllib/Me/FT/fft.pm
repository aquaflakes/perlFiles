package fft;
use warnings;
use strict;
use Data::Dumper;
use List::Util qw(sum);
use Storable 'dclone';
use Math::Complex;

sub fftADM #(ADM containing frequencies, sampling rate)
{
    my %ADM_freq=%{dclone shift}; my $SR=shift;
    my %fftPower;
    foreach my $keynt (keys(%ADM_freq))
    {
        #!!!!!remove the first 5 point and the end 5 point with high variation in bk
        ##########################
        my $splice_l = (length $keynt)-1;
        my $frontSpl=5; my $endSpl =-5 + $splice_l; #adjust the input length to be the same for multi nt
        my @spliced= splice @{$ADM_freq{$keynt}}, $frontSpl, $endSpl; #return subarray from 4 to (end-4)
        ##########################
        # substract mean freq before calc
        my $meanOfKey= &mean(@spliced);
        for (my $i = 0; $i < scalar(@spliced); ++$i)
				{
					$spliced[$i] = $spliced[$i]-$meanOfKey;
				}
                
        my @x1=@spliced;
        
        #print Dumper($ADM_freq{$keynt});die;
        $fftPower{freq}=\my @emptyarr;
        $fftPower{$keynt}=\my @emptyarr1;
         
        my $n=@x1; #n input
        my $m=log($n)/log(2); #FFT padding
        my $k= int($m) + ($m != int($m)); #round up for power
        my $np=2**$k;   #n after padding (for FFT)
        #my $np=$n;
        
        my @x;
         for (my $i=0; $i< $np; $i++)
        {
            push @x,((defined $x1[$i])? 1*$x1[$i] : 0); 
        }
         
        my $nx=@x; 

        my $pi=3.14159;
        my (@xre,@xim,@P,@Px);
        for ( my $k = 0; $k < $nx; $k++ )
        {
        #real part
        $xre[$k]=0;
        for ( $n = 0; $n < $nx; $n++ )
        {
        $xre[$k]+=$x[$n]*cos(((2*$pi*$k*$n)/$nx));
        }
        
        #imaginary part
        $xim[$k]=0;
        for ( $n = 0; $n < $nx; $n++ )
        {
        $xim[$k]-=$x[$n]*sin(((2*$pi*$k*$n)/$nx));
        }
         #power spectra
        $P[$k]=sqrt((($xre[$k]*$xre[$k])+($xim[$k]*$xim[$k])));
        push @Px, $P[$k];
         }
         
       # open OUT," > outputfft.txt" or die "$!\n";
         #my $maxP=max(@Px);
         for ( my $k = 0; $k <= (@Px/2); $k++ )
        {
         my $f=($k/(@Px/2))*(1/(2*$SR));
         push $fftPower{freq}, $f;
         push $fftPower{$keynt}, $Px[$k]/($meanOfKey*5.1*scalar(@x1)/10.2);
        }
    }
    
    #print Dumper(\%fftPower); die;
    return \%fftPower;
}


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
        my $frontSpl=5; my $endSpl =-5 + $splice_l; #adjust the input length to be the same for multi nt
        my @spliced= splice @{$ADM_freq{$keynt}}, $frontSpl, $endSpl; #return subarray from 4 to (end-4)
        ##########################
        # substract mean freq before calc
        my $meanOfKey= &mean(@spliced);
        for (my $i = 0; $i < scalar(@spliced); ++$i)
				{
					$spliced[$i] = $spliced[$i]-$meanOfKey;
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
        my $phaseAngle_h196_corr=$phaseAng_h196-2*($frontSpl/10.204); # correct for the phase of splicing at the start
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

1;
 
########### PLOT#############################################
#open my $GP, '|-', 'gnuplot -persist';
#print {$GP}  <<'__GNUPLOT__';
#
#          set multiplot;      # get into multiplot mode
#          set size 1,0.5;  
#          set origin 0.0,0.5;  
#    set xlabel 'Time (Second)'
#    set ylabel 'Amplitude'
#    set title 'TIME DOMAIN' 
#    set xtics font "Verdana,7" 
#    set ytics font "Verdana,7"
#          plot 'inputfft.txt' using 1:2 title '' with lines;
#
#
#          set origin 0.0,0.0;  
#    set xlabel 'Frequency (Hz)'
#    set ylabel 'Amplitude'
#    set title 'FREQUENCY DOMAIN' 
#    set xtics font "Verdana,7" 
#    set ytics font "Verdana,7"
#    plot 'outputfft.txt' using 1:2 title '' with lines
#          unset multiplot    # exit multiplot mode
#
#__GNUPLOT__
#
#close $GP;
