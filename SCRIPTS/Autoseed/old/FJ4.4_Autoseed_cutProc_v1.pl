#!/usr/bin/perl -w
##perl programme #1
#used to cut long reads into 40bp shorter ones for autoseed analysis

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use autodie;
use Me::Nu_cnt::checkfile;
use POSIX qw(ceil strftime);
use File::Spec qw(abs2rel);
use Parallel::ForkManager;

		
######### global variables to adjust #############
my $wkPATH= "/wrk/data/fangjie/"; # set the rootpath of working

######### time and usage
my $tiimaalku =localtime(); my $start_time=time();
my $this_programme = basename($0);
my $cmdInput=$this_programme." ".join(" ", @ARGV);
my $usage = "\n\$$this_programme [JustReads file | folder]\n
	Processes all files undr JustRead subfolders
	Example: $this_programme -f 'file' -g 'FACcut_t' (-b 0) (--jr -l 51) --copy (--test first and check guidefile) --phponly
	cyc0 use: $this_programme <folder> (-l 51) --copy (without -g!!)\n
	***options***
	-l	set length of ligand manually (for shorter reads containing NNNN... tails induced by PEAR)
	--jr	specify the target folder as JustRead folder (should apply for link folder)
	--guide	generate guidefiles for non-c0
	--copy	copy seqfiles and guidefile to Autoseed dir and run
	--test	only process 50 lines for each file
	-c	num of parallel processes (default 6)
	-g ''	generate guidefiles for non-c0 and add string '' to name of the guide file
	-b	specify the cyc to be used as background, if -g exists, usually 0
	--phponly  only run the php file, no cut and copy (do not forget -b0)
	";


########## Info from cmd line #################
	##### options
	GetOptions ('l=i' => \my $lengthDef, 'guide' => \my $guideFile, 'phponly' => \my $phponly,'test' => \my $isTest, 'copy' => \my $toCopy, 'c=i' => \my $cores, 'g=s' => \my $guideName_add,
				'jr' => \my $isJustReadDir, 'b=s' => \my $bk_cyc, 'o=s' => \my $outputPath, 'overwrite' => \my $overwrite, 'f=s' => \my $files);#, 'c=i' => \$cores, 'n=s' => \$appended_foldername, 'p=f' => \$p_value ); # i:int   f:float   s:string
	
	##### remaining param

	my $specified_path_absolute = "./"; #File::Spec->rel2abs($ARGV[0]);
	### extract the name list of reads files to be processed
	my @JustReads_files=glob($files);#$specified_path_absolute."/*.{gz,fastq,txt}");
	chomp @JustReads_files;
		#---print Dumper(@JustReads_files);die;
	### construct output path (extract only the path info from $ARGV[0])
	my $outputPathAdd = "/40bp_for_autoseed";
		#folder name of the output folder under the current
	$outputPath||="."; $outputPath .= $outputPathAdd;
		# if it is a file, get its path or use ".", if input is directory than use it for concatanation directly
		

	unless($phponly) {
	use File::Spec qw(rel2abs); my $currpath=File::Spec->rel2abs("./");
    qx(rm -r $outputPath) if (-d $outputPath);    qx(mkdir -p $outputPath); qx(mkdir -p $outputPath/perlcmd/);
    qx(echo ">> launching folder\n""$currpath""\n\n>> ORIGINAL CMD\n""$cmdInput""\n\n>> META CHAR EXPANDED\n"$cmdInput >$outputPath/perlcmd/perlcmd.txt); # memo input cmd
	}
		
	##### generate guideFile if required
	my $guideFileName; my $guideFileNameOnly;
	if ($guideName_add)
	{
		if (!(-d "$outputPath")) {system ("mkdir -p $outputPath");} # if nt_count folder not exist create it first
		system("cp -a $wkPATH/DefaultFilesforProc/guidefiles/template/blank_guide.txt $outputPath");
		$guideFileName="$outputPath/$guideName_add"."_Fangj".strftime('%Y%m',localtime).".txt";
		system("mv $outputPath/blank_guide.txt $guideFileName") unless $phponly;
		$guideFileNameOnly="$guideName_add"."_Fangj".strftime('%Y%m',localtime);
	}
	
	if ($phponly) {  # only run php
		$outputPath=File::Spec->rel2abs($outputPath); 
		chdir("/var/www/kaz");
		my $guidefile=glob("$outputPath/*Fangj*.txt"); $guidefile=~m{(?P<nameonly>[^//]*)\.[^\.]*$}; $guideFileNameOnly=$+{nameonly};
		system("ls $outputPath/*Fangj*.txt"); print "copying guide\n";
		system("cp -a $outputPath/*Fangj*.txt /var/www/kaz/guide");
		# run all cycs
		#system("nohup php asks_gen.php $guideFileNameOnly all autoseed 8 thread=10 >$outputPath/nohup_ask_gen.txt&"); # "&" to run at background
		
		# run also the bk cyc if defined b0
		if (defined($bk_cyc))
		{	# cyc like 4b0 will not run unless delete "all"
			system("nohup php asks_gen.php $guideFileNameOnly autoseed 8 thread=4 >$outputPath/nohup_ask_gen_b0_phponly.txt&"); #ksow
		}
		print("nohup log to $outputPath/nohup_ask_gen_b0_phponly.txt\n");
		exit;
	}

######### Loop through all JustReads files gen cut files
my $temp_cnt=1;

my $pm = Parallel::ForkManager->new($cores||10); #PPPPPP
Process:#PPPPPP
foreach my $JustReadsFile (@JustReads_files)
{   
    #print "\n########### file $temp_cnt ##############\n";
	$pm->start and next Process;#PPPPPPPP
	
	use Mylib1::seqFile; my @allseq;
	&seqFile::getallSeq(\@allseq,$JustReadsFile); chomp @allseq;
	
    ####### figure out ligand length
    my $liglength= length ($allseq[0]);
	print "length of the reads set to $liglength bp for $JustReadsFile \n";

	
	####### calc the dividing pos (into 40bp) for the long reads
	my $division = ceil($liglength/40); #round up to decide how many cuts (output file numbers) of 40bp
		# set position of each division
	my @subSeq40; for (my $i=0; $i<($division-1);$i++) {$subSeq40[$i]->{start}=40*$i+1; $subSeq40[$i]->{end}=40*($i+1);}
		# set position of final division
	$subSeq40[scalar(@subSeq40)]->{start}=$liglength-39; $subSeq40[scalar(@subSeq40)-1]->{end}=$liglength;
		# !! note the position should -1 before use
	
	####### read the original seq file, extract subseq

	
	my $counter=0;
	foreach (@allseq)
	{
		if (!($_=~ m/.*\w.*/)) {next;}
			#skip blank row		
		for (my $i=0; $i<$division;$i++) {push(@{$subSeq40[$i]->{seq}}, substr($_, $subSeq40[$i]->{start}-1, 40));}
		if($isTest){$counter++; if ($counter>49) {last;}} # process only 50 lines when --test
	}
	
	##### extract info form the name of JustReadfile
	$JustReadsFile =~ s/\+/S40N\+/; #!!!!!S40N is important for Autoseed to run
	$JustReadsFile =~ m{\/([^/_]*)[^/]*\.\w+$}; 
	my $outFilename = $1;
	$1 =~ m{(?P<barcode>.*)(?P<batch>\+.*\+c)(?P<cyc>.)(?=$)};
	$outFilename =~ tr/\+-/../;
	$outFilename.="_sig.seq";
	
	
	#my $aa="+"; my $bb="alfkj+dsf+sfd";
	#print $bb=~ s/\Q$aa\E/**/rg;
	
	
	my $bar_in_guide=($+{barcode}); $bar_in_guide=~ tr/\+-/../;
	my $batch = $+{batch}; $batch=~ tr/\+-/../;
	my $cyc = $+{cyc}; $cyc=~ tr/\+-/../;
	
	 
	##### write subseq to files
	for (my $i=0; $i<$division; $i++)
	{
		if ($i=$division-1) {  # only take the last division		
		##### construct output Filename
		#my $outFilename = ($+{barcode}).$subSeq40[$i]->{start}."S"."40N".($+{batch}).($+{cyc})."_sig.seq";		
		if (!(-d "$outputPath")) {system ("mkdir -p $outputPath");} # if nt_count folder not exist create it first
		open(my $outSeq, ">", "$outputPath/$outFilename") or die "cannot open write fh for subseq output";
		print $outSeq join ("\n", @{$subSeq40[$i]->{seq}});
		close $outSeq;
		}
	}
	
	##### write guideFile if required
	if ($guideName_add && (defined($cyc)) && !($+{cyc}==0)) # not writing guide file when cyc0
	{
		open(my $guide_fh, ">>", $guideFileName) or die "cannot open guidefile to append";
			print $guide_fh "NU\t"."Full\t"."$bar_in_guide\t".$batch."\t"."AA\t"."1\t".$cyc."\t-\t-\t-;\n";
		close $guide_fh;
	}
	

    #------------- fixed message for ending a JustReadfile proc
    #print $JustReadsFile." processing finished";
	#print "\n";
	#print "############ file $temp_cnt #############\n\n";
	$pm->finish;#PPPPPP
    $temp_cnt++;
}
$pm->wait_all_children;	#PPPPPP

##### sort and deduplicate the guide file ####
if ($guideName_add)
{
    #no warnings;
	#open (my $data , '<', $guideFileName)|| die "could not open $guideFileName to add bk_cyc:\n$!"; # write and read, if file not exist return error
	#my @array=<$data>;
    #my @sorted=sort {(split(/\t/,$a))[6]<=>(split(/\t/,$b))[6]} @array; # sort according to col7 (cyc)
    #open (my $data1 , '>', $guideFileName)|| die "could not open $guideFileName:\n$!"; 
    #print $data1 @sorted; close $data1;
	#system ("sort -k 4 $guideFileName| tac | awk -F\"\\t\" '!_[\$4]++'"." |tac >$guideFileName"."1"); # awk extract only uniq records according to filed 3, 4b0
	system ("sort -k 4 $guideFileName| awk -F\"\\t\" '!_[\$4]++'"." >$guideFileName"."1"); # awk extract only uniq records according to filed 3, 3b0
	system ("mv "."$guideFileName"."1"." $guideFileName"); #should save another file first otherwise will return blank file
	
	#if bk cyc defined then modify the cyc in guidefile like "4b0"
	#print $bk_cyc."\n";die;
	if (defined($bk_cyc))
	{
		# add b0 to cyc numbers, use cyc0 as background in auto seed
		my $replace_cmd = "awk -v OFS=\"\\t\" '\$7=\$7\"b$bk_cyc\"' $guideFileName"." | sed 's/cycleb0/cycle/' "." >$guideFileName"."1"; #first replace all col7 than replace the first row back
		system ($replace_cmd);
		system ("mv "."$guideFileName"."1"." $guideFileName"); 
	}
	system ("chmod 666 "."$guideFileName"); #important

}

########## copy for autoseed ##########
if ($toCopy)
{
	$outputPath=File::Spec->rel2abs($outputPath);
	system("ln -sf $outputPath/*.seq /var/www/kaz/data/sequence");
	if ($guideName_add)
	{
		print "\n\n"; system("ls $outputPath/*Fangj*.txt"); print "copying guide\n";
		system("cp -a $outputPath/*Fangj*.txt /var/www/kaz/guide");
		chdir("/var/www/kaz");
		
		my $ow= $overwrite? " overwrite":"";
		# run all cycs
		system("nohup php asks_gen.php $guideFileNameOnly all order=0123 8$ow thread=10 >$outputPath/nohup_ask_gen.txt&"); # "&" to run at background
		
		# run also the bk cyc if defined b0
		if (defined($bk_cyc))
		{	# cyc like 4b0 will not run unless delete "all"
			system("nohup php asks_gen.php $guideFileNameOnly order=0123 8$ow thread=10 >$outputPath/nohup_ask_gen_b0.txt&"); #ksow
		}
		print("nohup log to $outputPath/nohup_ask_gen.txt\n");
		system("nohup php asks_gen.php $guideFileNameOnly all order=0123 8$ow thread=10 >$outputPath/nohup_ask_gen.txt&"); # "&" to run at background
		system("nohup php sequence_checker.php >$outputPath/nohup_seq_checker.txt&");
		
	}
}




#------------- fixed ending of program 
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
exit;







