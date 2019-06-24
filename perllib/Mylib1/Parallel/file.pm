package file;
# multi thread is still buggy  2015/12/01

use warnings;
use strict;


sub mProcess
{
    my($processMaxNum,$singleproc_Cref,$params_Aref)=@_;
    
    use Parallel::ForkManager;
    
	my $pm = Parallel::ForkManager->new($processMaxNum); my $cnt=0;
	FileLoop:
	foreach my $param (@$params_Aref)
	{	
		$cnt++;
        $pm->start and next FileLoop;
        $singleproc_Cref->($param);	
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub mThread
{
	my($threadMaxNum,$singleproc_Cref,$params_Aref)=@_;
	
	use threads;
    use Thread::Queue;
  
    my $q = Thread::Queue->new;
    $q->enqueue(@$params_Aref); #needs a file list

    my $num_workers = @$params_Aref < $threadMaxNum ? @$params_Aref : $threadMaxNum; # no need to be wasteful :)
	
	my $worker = sub { my $queue = shift; while (my $param = $queue->dequeue_nb) {$singleproc_Cref->($param);} return 1};
    FileLoop:
	for (1 .. $num_workers)
	{
		threads->new($worker, $q);
    }
	#while (threads->list(threads::running)) { sleep 10;}
	
    $_->join for threads->list;
}

sub packParams
{
	my @allparamArr=@_;
	for(my $i=0; $i<=$#allparamArr;$i++)	# check if all paramters have equal num of elements
	{
		my $lenArr_curr=@{$allparamArr[$i]};
		for (my $j=$i+1; $j<=$#allparamArr;$j++)
		{
			die "input do not contain equa number of elements for each params\n" if $lenArr_curr!=@{$allparamArr[$j]};
		}
	}
	
	my @params = map { my @unit; foreach my $param (@allparamArr) {push @unit, $param->[$_]}; \@unit;} (0..@{$_[0]}-1);
	\@params;
}



1;	