package PWM;
# handling PWM and pfm

use strict;

sub new
{
  my ($class,$PWM_Href)=@_;
  bless $PWM_Href, $class;
  return $PWM_Href;
}

sub readPWM
{
  my ($PWM,$sep,$order)=@_;
  my %PWM;
  open(my $fh, "<", $PWM) or die "cannot open PWM file $PWM\n";
  foreach (split //,$order)
  {
	my $line=<$fh>; chomp $line;
	$PWM{$_}=[split($sep,$line)];
  } 
  return new PWM(\%PWM);  
}

sub revComp
{
  my ($PWM_Href)=@_;
  my %PWM_rc= map {my $cmp= $_; $cmp =~ tr/ACGTacgt/TGCAtgca/; my @val= reverse @{$PWM_Href->{$_}}; ($cmp,\@val);} keys %$PWM_Href;
  return new PWM(\%PWM_rc);
}

sub output
{
  my ($PWM_Href,$path)=@_;
  use Mylib1::fileIO;
  my $dir = &fileIO::getdir($path);
  qx(mkdir -p $dir);
  open(my $out, ">", $path) or die "cannot open file $path to write pwm";
  foreach (qw(A C G T)) {print $out ($_, "\t", join("\t",@{$PWM_Href->{$_}}),"\n")};
}

sub tofreq
{
  my ($PWM_Href,$bkfreq)=@_;
  #$PWM_Href->output("./pfmtest/ori.pfm");
  my $sumPos1; $sumPos1 += $PWM_Href->{$_}->[0] for keys(%$PWM_Href);
  if ($sumPos1>1.2 || ($sumPos1<0.8 && $sumPos1>0)) # if eq 1 then is already freq pfm
  {
	my @sum_eachPos= map {my $sum_curr; for my $key (qw(A T C G)) {$sum_curr+= $PWM_Href->{$key}->[$_];} $sum_curr;} (0..@{$PWM_Href->{A}}-1);
	map { for my $pos (0..@{$PWM_Href->{$_}}-1) {$PWM_Href->{$_}->[$pos]/= $sum_eachPos[$pos]; $PWM_Href->{$_}->[$pos] ||=$bkfreq; } } (qw(A T C G));
  }
  #$PWM_Href->output("./pfmtest/freq.pfm");
}

sub tolog
{
  my ($PWM_Href,$bkfreq)=@_;
  $PWM_Href->tofreq($bkfreq);
  foreach (keys %$PWM_Href) {foreach (@{$PWM_Href->{$_}}) {$_=log($_)/log(10);} }
  #$PWM_Href->output("./pfmtest/log.pfm");
}






sub genbk_even
{
  my ($PWM_Href,$evenfreq)=@_;
  use Storable;  my $bkPWM_Href= &Storable::dclone($PWM_Href);
  map {my $curr=$evenfreq->{$_}; foreach (@{$bkPWM_Href->{$_}}) {$_=$curr;} } keys %$bkPWM_Href;
  $bkPWM_Href->tolog(0.000001);
  return $bkPWM_Href;
}
    
1;    