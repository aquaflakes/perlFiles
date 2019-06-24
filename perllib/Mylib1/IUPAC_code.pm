package IUPAC_code;
# handling oligos with IUPAC_code

use strict;

sub decode_v2 # take like "AATN AAAN" or "AAAN", return a list of expanded oligos of all elements
{
	my $oligo = uc(shift); chomp $oligo;
	my @oligos=split(/\s/,$oligo);
	my $decode_ref={
	R=> [qw(A G)],
	Y=> [qw(C T)],
	N=> [qw(A T C G)],
	W=> [qw(A T)],
	S=> [qw(G C)],
	M=> [qw(A C)],
	K=> [qw(G T)],
	B=> [qw(G C T)],
	H=> [qw(A C T)],
	D=> [qw(A G T)],
	V=> [qw(A G C)]
	};
	my $allAbbr=join("", keys $decode_ref);
	
	foreach (@oligos)
	{
		if ($_=~m/[$allAbbr]/) # only replace the first match
		{
			my @subset=@{$decode_ref->{$&}}; my $bf=$`; my $af=$';
			my @decoded= map {decode_v2($bf.$_.$af)} @subset; # use of $' here cause problem
			$_= join (" ", @decoded);
		}
	}
	return split (/\s/, join (" ", @oligos));
}

sub decode # take single degenerated oligo, return a list of expanded oligos
{
	my $oligo=shift;
	use Bio::Tools::IUPAC;
	use Bio::PrimarySeq;
	my $ambiseq = Bio::PrimarySeq->new(-seq => $oligo, -alphabet => 'dna');
	
	# Create all possible non-degenerate sequences
	my $iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq);
	my @oligos;
	while (my $uniqSeq= $iupac->next_seq())
	{
		# process the unique Bio::Seq object.
		push @oligos, $uniqSeq->{seq};
	}
	return @oligos;
}

sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}
    
1;    