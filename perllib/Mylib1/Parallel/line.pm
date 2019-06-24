package line;


use warnings;
use strict;
use POSIX;


sub mProcess_l
{
    my($processNum,$singleproc_Cref,$allreads_Aref,$param_Aref)=@_;
    
    use Parallel::ForkManager; # v1.17 required for Href return
    
	# aliquot the reads for each proc
	my $lines_each=ceil(scalar(@$allreads_Aref)/$processNum);
	my @spliced;
	while (@$allreads_Aref)
	{
		my @nextArr = @$allreads_Aref>$lines_each? splice(@$allreads_Aref,0,$lines_each) : splice(@$allreads_Aref,0,scalar(@$allreads_Aref));
		push (@spliced,\@nextArr);
	}
	
	my @temp_Data;
	my $pm = Parallel::ForkManager->new($processNum); my $cnt;
	$pm->run_on_finish(sub{ 
							my($pid,$exit_code,$ident,$exit_signal,$core_dump,$v)=@_;
							push @temp_Data,$v;
						   });  # run at the finish of every child, $v is Href returned in each $pm->finish
	Process:
	foreach (@spliced)
	{
		#$cnt++; print "proc ".$cnt." ".scalar(@$_)."\n";
		$pm->start and next Process;
		my $Data_Href=$singleproc_Cref->($_, $param_Aref);
		$pm->finish(0,$Data_Href);
	}
	$pm->wait_all_children;

	return \@temp_Data;
}


sub mThread_l
{
	my($threadNum,$singleproc_Cref,$allreads_Aref,$param_Aref)=@_;
	
	use threads;
    use Thread::Queue;
	
	# aliquot the reads for each thread
	my $lines_each=ceil(scalar(@$allreads_Aref)/$threadNum);
	my @spliced;
	while (@$allreads_Aref)
	{
		my @nextArr = @$allreads_Aref>$lines_each? splice(@$allreads_Aref,0,$lines_each) : splice(@$allreads_Aref,0,scalar(@$allreads_Aref));
		push (@spliced,\@nextArr);
	}
	
	 
    my $q = Thread::Queue->new;
    $q->enqueue(@spliced); #needs a file list
	
	my @threads; my @temp_Data;
    while (my $reads = $q->dequeue_nb)
	{
		#my $thread=threads->new(sub{ my $Data_Href= $singleproc_Cref->($reads,$param_Aref);	return $Data_Href;});
		my $thread=threads->new($singleproc_Cref,$reads,$param_Aref);
		push @threads, $thread;
    }
	
	#do { sleep 10; print $_."\n";  foreach (threads->list) {push (@temp_Data, $_->join) if $_->is_joinable() }; } while threads->list;	
	
	foreach (@threads)	{ push @temp_Data, $_->join(); }
	die "not joint for all threads" if threads->list;
	return \@temp_Data;
}


1;	