#!/usr/bin/perl -w
##perl programme #1
#used to cut long reads into 40bp shorter ones for autoseed analysis

use FindBin;               # locate this script
use lib "$FindBin::Bin/../.."; 
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use autodie;
use perllib::Me::Nu_cnt::checkfile;
use POSIX qw(ceil strftime);
use File::Spec qw(abs2rel);
use Parallel::ForkManager;
#use Cwd; # get current dir
# use File::Temp qw(tempdir);
		
######### global variables to adjust #############
my $wkPATH= "/wrk/data/fangjie/"; # set the rootpath of working

######### time and usage
my $tiimaalku =localtime(); my $start_time=time();
my $this_programme = basename($0);
my $usage = "\n\$$this_programme [JustReads file | folder]\n
	Processes all files undr JustRead subfolders
	Example: $this_programme <folder> -g 'FACcut_t' (-b 0) (--jr -l 51) --copy (--test first and check guidefile)
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
	";


########## Info from cmd line #################
	##### options
	GetOptions ('l=i' => \my $lengthDef, 'guide' => \my $guideFile,'test' => \my $isTest, 'copy' => \my $toCopy, 'c=i' => \my $cores, 'g=s' => \my $guideName_add,
				'jr' => \my $isJustReadDir, 'b=s' => \my $bk_cyc,);#, 'c=i' => \$cores, 'n=s' => \$appended_foldername, 'p=f' => \$p_value ); # i:int   f:float   s:string
	
	##### remaining param
	if (!$ARGV[0]) { print "usage: $usage\n\n";exit;}
	my $specified_path_absolute = File::Spec->rel2abs($ARGV[0]);
	### extract the name list of reads files to be processed
	my @JustReads_files;
		if ( -f $specified_path_absolute ) {@JustReads_files = $specified_path_absolute;}
			# pick the file name if filename input
		else{
				# if the current path is "JustRead then take the files included"
				if(($specified_path_absolute =~ m/justread/i)||$isJustReadDir) {@JustReads_files = <$ARGV[0]/*.txt>;}  
				# if current path is not JustRead then search for all JustRead subfolder
				else{ 
						my @JustReads_file_folders =`find $specified_path_absolute -iname "*Justread*" -type d 2>/dev/null | grep -v "copied"`;
						foreach (@JustReads_file_folders)
						{
							chomp; # eliminate the \n
							my @filestmp=glob("$_/*.txt");
							push(@JustReads_files,@filestmp);
						}
					}	
			}  
		#---print Dumper(@JustReads_files);die;
	### construct output path (extract only the path info from $ARGV[0])
	my $output_pathAdd = "/40bp_for_autoseed";
		#folder name of the output folder under the current
	my $output_path;    if ( -f $specified_path_absolute ) {$specified_path_absolute =~ m/^(.*\/)/;  my $matched = $1 || "."; $output_path = $matched.$output_pathAdd;} else { $output_path = $specified_path_absolute.$output_pathAdd;}
		# if it is a file, get its path or use ".", if input is directory than use it for concatanation directly
	if (-d "$output_path") {system ("rm -r $output_path");}
		# if the output folder already exists delete it
		
	##### generate guideFile if required
	my $guideFileName; my $guideFileNameOnly;
	if ($guideName_add)
	{
		if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		system("cp -a $wkPATH/DefaultFilesforProc/guidefiles/template/blank_guide.txt $output_path");
		$guideFileName="$output_path/$guideName_add"."_Fangj".strftime('%Y%m',localtime).".txt";
		system("mv $output_path/blank_guide.txt $guideFileName");
		$guideFileNameOnly="$guideName_add"."_Fangj".strftime('%Y%m',localtime);
	}

######### Loop through all JustReads files gen cut files
my $temp_cnt=1;

my $pm = Parallel::ForkManager->new($cores||6); #PPPPPP
Process:#PPPPPP
foreach my $JustReadsFile (@JustReads_files)
{   
    #print "\n########### file $temp_cnt ##############\n";
	$pm->start and next Process;#PPPPPPPP
	
    ####### figure out ligand length
    my $liglength= $lengthDef || &checkfile::formatCheck_and_liglength($JustReadsFile);
	print "length of the reads set to $liglength bp for $JustReadsFile \n";

	
	####### calc the dividing pos (into 40bp) for the long reads
	my $division = ceil($liglength/40); #round up to decide how many cuts (output file numbers) of 40bp
		# set position of each division
	my @subSeq40; for (my $i=0; $i<($division-1);$i++) {$subSeq40[$i]->{start}=40*$i+1; $subSeq40[$i]->{end}=40*($i+1);}
		# set position of final division
	$subSeq40[scalar(@subSeq40)]->{start}=$liglength-39; $subSeq40[scalar(@subSeq40)-1]->{end}=$liglength;
		# !! note the position should -1 before use
	
	####### read the original seq file, extract subseq
	open(my $singleSeq, "<", $JustReadsFile) or die "cannot open fh for reads file $JustReadsFile";
	my $counter=0;
	while (<$singleSeq>)
	{
		if (!($_=~ m/.*\w.*/)) {next;}
			#skip blank row		
		for (my $i=0; $i<$division;$i++) {push(@{$subSeq40[$i]->{seq}}, substr($_, $subSeq40[$i]->{start}-1, 40));}
		if($isTest){$counter++; if ($counter>49) {last;}} # process only 50 lines when --test
	}
	close $singleSeq;
	
	##### extract info form the name of JustReadfile			
	$JustReadsFile =~ m{\/([^/]*)\.\w+$} ||  $JustReadsFile =~ m/(.*)...\w+$/;
	$1 =~ m/(?<batch>^[\w\.]*)_(?<barcode>[0-9A-Z]*[A-Z])[0-9]*N_(?<cyc>[0-9]$)/;
	 
	##### write subseq to files
	for (my $i=0; $i<$division; $i++)
	{
		##### construct output Filename
		my $outFilename = ($+{barcode}).$subSeq40[$i]->{start}."S"."40N".($+{batch}).($+{cyc})."_sig.seq";		
		if (!(-d "$output_path")) {system ("mkdir -p $output_path");} # if nt_count folder not exist create it first
		open(my $outSeq, ">", "$output_path/$outFilename") or die "cannot open write fh for subseq output";
		print $outSeq join ("\n", @{$subSeq40[$i]->{seq}});
		close $outSeq;
	}
	
	##### write guideFile if required
	if ($guideName_add && !($+{cyc}==0)) # not writing guide file when cyc0
	{
		open(my $guide_fh, ">>", $guideFileName) or die "cannot open guidefile to append";
		for (my $i=0; $i<$division; $i++)
		{
			my $bar_in_guide=($+{barcode}).$subSeq40[$i]->{start}."S"."40N";
			print $guide_fh "NU\t"."Full\t"."$bar_in_guide\t".($+{batch})."\t"."AA\t"."1\t".($+{cyc})."\t-\t-\t-;\n";
		}
		close $guide_fh;
	}
	

    #------------- fixed message for ending a JustReadfile proc
    print $JustReadsFile." processing finished";
	print "\n";
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
	system ("sort -k 7 $guideFileName| tac | awk -F\"\\t\" '!_[\$3]++'".">$guideFileName"."1"); # awk extract only uniq records according to filed 3
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
	system ("chmod 666 "."$guideFileName");

}

########## copy for autoseed ##########
if ($toCopy)
{
	system("cp -a $output_path/*.seq /var/www/kaz/data/sequence");
	if ($guideName_add)
	{
		system("cp -a $output_path/*.txt /var/www/kaz/guide");
		chdir("/var/www/kaz");
		system("nohup php asks_gen.php $guideFileNameOnly all kmerspace 8 thread=10 >$output_path/nohup_ask_gen.txt&"); # "&" to run at background
		if (defined($bk_cyc))
		{	# cyc like 4b0 will not run unless delete "all"
			system("nohup php asks_gen.php $guideFileNameOnly kmerspace 8 thread=10 >$output_path/nohup_ask_gen_bk.txt&"); 
		}
		print("nohup log to $output_path/nohup_ask_gen.txt\n");
	}
}


#------------- fixed ending of program 
my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
exit;







