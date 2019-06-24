package Hash_Merge;

use Hash::Merge qw( merge );
use strict;

Hash::Merge::specify_behavior(  ## addition for all equal positions
    {
                'SCALAR' => {
                        'SCALAR' => sub { $_[0]+$_[1] }, # scalar with scalar
                        'ARRAY'  => sub { #[ $_[0], @{$_[1]} ];
										 "type diffs for this pos"; }, # scalar with array
                        'HASH'   => sub { #$_[1] ;
										 "type diffs for this pos";}, # scalar with hash
                },
                'ARRAY' => {
                        'SCALAR' => sub { #$_[1];
										 "type diffs for this pos"; },
                        'ARRAY'  => sub { return "two array length differs" if @{$_[0]} != @{$_[1]}; my $ArrL= $#{$_[0]};  my @Arr1=@{$_[0]}; my @Arr2=@{$_[1]};  my @array=map $Arr1[$_]+$Arr2[$_], (0..$ArrL);  \@array; },
                        'HASH'   => sub { #$_[1] ;
										 "type diffs for this pos";}, 
                },
                'HASH' => {
                        'SCALAR' => sub { #$_[1] ;
										 "type diffs for this pos";},
                        'ARRAY'  => sub { #[ values %{$_[0]}, @{$_[1]} ] ;
										 "type diffs for this pos";},
                        'HASH'   => sub { Hash::Merge::_merge_hashes( $_[0], $_[1] ) }, 
                },
        }, 
        'addition', 
);

Hash::Merge::specify_behavior(  ## addition for all equal positions
    {
                'SCALAR' => {
                        'SCALAR' => sub { $_[0]-$_[1] }, # scalar with scalar
                        'ARRAY'  => sub { #[ $_[0], @{$_[1]} ];
										 "type diffs for this pos"; }, # scalar with array
                        'HASH'   => sub { #$_[1] ;
										 "type diffs for this pos";}, # scalar with hash
                },
                'ARRAY' => {
                        'SCALAR' => sub { #$_[1];
										 "type diffs for this pos"; },
                        'ARRAY'  => sub { return "two array length differs" if @{$_[0]} != @{$_[1]}; my $ArrL= $#{$_[0]};  my @Arr1=@{$_[0]}; my @Arr2=@{$_[1]};  my @array=map $Arr1[$_]-$Arr2[$_], (0..$ArrL);  \@array; },
                        'HASH'   => sub { #$_[1] ;
										 "type diffs for this pos";}, 
                },
                'HASH' => {
                        'SCALAR' => sub { #$_[1] ;
										 "type diffs for this pos";},
                        'ARRAY'  => sub { #[ values %{$_[0]}, @{$_[1]} ] ;
										 "type diffs for this pos";},
                        'HASH'   => sub { Hash::Merge::_merge_hashes( $_[0], $_[1] ) }, 
                },
        }, 
        'substraction', 
);

Hash::Merge::specify_behavior(  ## addition for all equal positions
    {
                'SCALAR' => {
                        'SCALAR' => sub { $_[0]+$_[1] }, # scalar with scalar
                        'ARRAY'  => sub { $_[1] }, # scalar with array
                        'HASH'   => sub { $_[1] }, # scalar with hash
                },
                'ARRAY' => {
                        'SCALAR' => sub { #$_[1];
										 "type diffs for this pos"; },
                        'ARRAY'  => sub { return "two array length differs" if @{$_[0]} != @{$_[1]}; my $ArrL= $#{$_[0]};  my @Arr1=@{$_[0]}; my @Arr2=@{$_[1]};  my @array=map $Arr1[$_]+$Arr2[$_], (0..$ArrL);  \@array; },
                        'HASH'   => sub { #$_[1] ;
										 "type diffs for this pos";}, 
                },
                'HASH' => {
                        'SCALAR' => sub { #$_[1] ;
										 "type diffs for this pos";},
                        'ARRAY'  => sub { #[ values %{$_[0]}, @{$_[1]} ] ;
										 "type diffs for this pos";},
                        'HASH'   => sub { Hash::Merge::_merge_hashes( $_[0], $_[1] ) }, 
                },
        }, 
        'addition_scalar_replacement', 
);

sub add # return Href
{
    my ($Href1,$Href2)=@_;
	my $merge = Hash::Merge->new( 'addition' );
    return $merge->merge($Href1,$Href2);
}

sub add_scalar_replacement # return Href
{
    my ($Href1,$Href2)=@_;
	my $merge = Hash::Merge->new( 'addition_scalar_replacement' );
    return $merge->merge($Href1,$Href2);
}

sub minus # return Href
{
    my ($Href1,$Href2)=@_;
	my $merge = Hash::Merge->new( 'substraction' );
    return $merge->merge($Href1,$Href2);
}
    
#    
#my %a = ( 
#        'AA'    => [1,12,3,5,6],
#        'CC'    => [ 1,2,3,4],
#        'GG' => { 'GGC' => [1,2,3,4,5], 'GGNN' => 1 },
#        'TTAA' => 10,
#    );
#my %b = ( 
#        'AA'    => [1,12,3,5,6],
#        'CC'    => [ 1,2,3,4],
#        'GG' => { 'GGC' => [1,2,3,4,5], 'GGNN' => 1 },
#        'TTAA' => 10,
#    );
#
#       #my %c = %{ merge( \%a, \%b ) };
#	   my $c = add( \%a, \%b );
#
#	   
#    print "";
    
1;    