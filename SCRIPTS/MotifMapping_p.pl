#!/usr/bin/perl -w

use FindBin;               # locate this script
use lib "$FindBin::Bin/..";  # use this directory
use warnings;
use strict;
use perllib::Me::parallel_v1;
use Bio::Perl;
use Bio::Seq;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use perllib::Me::MotifMap::CntArray;
use perllib::Me::Nu_cnt::checkfile;
use MOODS; use MOODS::Tools qw(printResults numResults readMatrix);
my $tiimaalku =localtime(); my $start_time=time();
my $this_programme = basename($0);
my $usage = "\$$this_programme [JustReads file | folder] [PWM file | folder] \n ***options***
-l lines to process (default:all lines) 
-c number of parallel porcess
-o output folder
-p pValue for Moods' motif match";


########## Info from cmd line #################
#options
my $row_to_process=1000000000000; my $cores=15;  #use default values if not overwrited
my $p_value=0.0001;
GetOptions ('l=i' => \$row_to_process, 'c=i' => \$cores, 'o=s' => \my $outputFolder, 'p=f' => \$p_value ); # i:int   f:float   s:string

#remaining param
if ((!$ARGV[0]) || (!$ARGV[1])) { print "usage: $usage\n\n";exit;}
my @JustReads_files; if ( -f $ARGV[0] ) {@JustReads_files = $ARGV[0];} else{@JustReads_files = <$ARGV[0]/*.txt>;}
my @matrix_files; if ( -f $ARGV[1] ) {@matrix_files = $ARGV[1];} else{@matrix_files = <$ARGV[1]/*.pfm>;}
#construct output path (extract only the path info from $ARGV[0])
my $output_path;    if ( -f $ARGV[0] ) {$ARGV[0] =~ m/^(.*\/)/;  my $matched = $1 || "."; $output_path = $matched."/analysis/$this_programme";} else { $output_path = $ARGV[0]."/analysis/$this_programme";}
$output_path=$outputFolder if $outputFolder;
$output_path.= "_pValue$p_value";
if (-d "$output_path") {system ("rm -r $output_path");}
my $tmpFile_path="$output_path"."/_tmp/";                 
######### variables to adjust #############
my $reads_anal_methodRef = \&MapSR; # in this file, process reads and store result in tmp hash/arr, return a tmp hash/arr
my $tmp_storage_output_methodRef = \&CntArray::wtCountFile; # write hash in the process of each core into a tmp file
my $summing_tmpStorage_methodRef = \&CntArray::sumCountFiles;
my $sum_storage_output_methodRef = \&CntArray::wtCountFile_Hash;#\&CntArray::wtCountFile;
		                                         
my $tmp_storage_ref;  #= \%tmp_storage;	# store temp data  ####my %tmp_storage = (A=>[],C=>[],T=>[],G=>[],);	#important for cases only di_nt is calculated, no s_nt information give err when ploting in R(designed for 20 lines)
my $isFastq =0; # 1 if input file is fastq format, 0 if is JustReads
################################################

my $Num_of_Matrix=scalar(@matrix_files);
my (@matrices, @rc_matrices, @r_matrices, @c_matrices, @all_matrices,         @matrix_Names, @Motif_L_of_matrices, @thresholds);
for (my $i=0;$i<$Num_of_Matrix;$i++)
{
    $matrices[$i] = readMatrix($matrix_files[$i]);
    my @c_tmp = reverse (@{$matrices[$i]});
    $c_matrices[$i] = \@c_tmp;

    $rc_matrices[$i] = MOODS::Tools::reverseComplement($matrices[$i]);
    my @r_tmp = reverse (@{$rc_matrices[$i]});
    $r_matrices[$i] = \@r_tmp;

    $matrix_Names[$i] = $matrix_files[$i]; $matrix_Names[$i] =~ s/\.pfm$//; $matrix_Names[$i]=~ m{[\/\\]([^\/\\]+)$}; $matrix_Names[$i]= $1; $matrix_Names[$i+$Num_of_Matrix]= "rv&c_".$1;#construt names for plot labeling
    $Motif_L_of_matrices[$i] = scalar(@{$matrices[$i][0]});
}
                   #signal                        # control
@all_matrices = (@matrices, @rc_matrices,    @r_matrices, @c_matrices); #join all matrices together to get faster mapping
@thresholds= (($p_value) x scalar(@all_matrices)); #()surronding p_value is important, ######## specify p value

## Loop through all JustReads files
my $temp_cnt=1;
foreach my $JustReadsFile (@JustReads_files)
{   
    print "\n########### file $temp_cnt ##############\n";
	 
    # construct output Filename
    my $row_to_process_tmp= $row_to_process;
    $JustReadsFile =~ m{\/([^/]*)\.\w+$} ||  $JustReadsFile =~ m/(.*)...\w+$/; my $outFilename = $1."_MotifMapCnt.txt"; #construct output filename for each JustReads file
    
    # figure out ligand length and cut_value for each PWM
    my $liglength= &checkfile::formatCheck_and_liglength($JustReadsFile,$isFastq); print "length of the reads set to $liglength bp for $JustReadsFile \n"; 
    my @cut_values;
    for (my $i=0; $i<$Num_of_Matrix;$i++) {$cut_values[$i] = $liglength-$Motif_L_of_matrices[$i]+1;}    
    
    
    #creat N countArrayRef and initialize with 0 for different PWMs
    my $countArrayRef;
    for (my $i=0; $i<$Num_of_Matrix;$i++)  { for (my $j=0; $j<$liglength; $j++) {$countArrayRef->[$i]->[$j]=0;$countArrayRef->[$i+$Num_of_Matrix]->[$j]=0};}
    $tmp_storage_ref = $countArrayRef;
    
    
    my ($tmp_storage_files,$readrows) =
    &parallel_v1::calc($JustReadsFile, $liglength, $row_to_process, $cores, $tmpFile_path, $isFastq,      $tmp_storage_ref, $reads_anal_methodRef, $tmp_storage_output_methodRef, \@matrix_Names, \@cut_values); #process for the files <<<<<<<<<<<<<<<
   
	my $countHashRef=$summing_tmpStorage_methodRef->($tmp_storage_files);
    
    #normalization against lines in reads file
    foreach my $key (keys($countHashRef)) {foreach (@{$countHashRef->{$key}}) {$_=$_/$readrows;} }
    
    $sum_storage_output_methodRef->($countHashRef,$output_path,$outFilename); #($countArrayRef,$matrix_NamesRef,$output_path,$outFilename)
    
    print $JustReadsFile." processing finished";
	print "\n";
	print "############ file $temp_cnt #############\n\n";
    $temp_cnt++;
}
if (-d "$tmpFile_path") {system ("rm -r $tmpFile_path");} # if the temp folder exist delete it
########## plot ##########
chdir $output_path;
system ('Rscript /wrk/data/fangjie/lib/Rlib/MotifMapPlot.R');
system ("convert -quality 95 ./PNG/* ./TFMotifMap_p_$p_value.pdf");


my $tiima =localtime(); my $end_time = time(); my $time_used = -$start_time+$end_time;
print "Program finished $tiima, processing time is $time_used\n\n";
exit;


sub MapSR
{
    my ($tmp_storage_ref,$combinedReads,$liglength,$cut_values_Ref)=@_;  #liglength, @matrices, @Motif_L_of_matrices should be passed if not in the same file
    my $countArrayRef =$tmp_storage_ref;
    my $seq = Bio::Seq->new(-seq  => $combinedReads, -alphabet => 'dna' );

        ## We convert the count matrix into log-odds scores in base 2 logarithm and find hits scoring
        my @results = MOODS::search(-seq => $seq, -matrices => \@all_matrices,  -thresholds => \@thresholds, -threshold_from_p => 1 );

        #add SR count into countArrays of each PWM
        my $Num_of_Matrix=scalar(@matrix_files);
        for (my $i=0; $i<$Num_of_Matrix;$i++)
        {
            #count mapping hits
            my @positions = do {my $j = 0; grep {!($j++ % 2)} @{$results[$i]}}; #extract only the position info but not scores for matrix[i]
            my @positions_rc = do {my $j = 0; grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix]}}; 
            $countArrayRef->[$i]= &CntArray::addCnt_to_CntSt_combinedRead($countArrayRef->[$i],\@positions,$cut_values_Ref->[$i],$liglength);            
            $countArrayRef->[$i]= &CntArray::addCnt_to_CntSt_combinedRead($countArrayRef->[$i],\@positions_rc,$cut_values_Ref->[$i],$liglength);
            #print join("\t",@{$countArrayRef->[$i]}); print "\n\n";die;
            
            #count ctrl hits
            my @positions_c = do {my $j = 0; grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix*2]}}; #extract only the position info but not scores for matrix[i]
            my @positions_r = do {my $j = 0; grep {!($j++ % 2)} @{$results[$i+$Num_of_Matrix*3]}}; 
            $countArrayRef->[$i+$Num_of_Matrix]= &CntArray::addCnt_to_CntSt_combinedRead($countArrayRef->[$i+$Num_of_Matrix],\@positions_c,$cut_values_Ref->[$i],$liglength);            
            $countArrayRef->[$i+$Num_of_Matrix]= &CntArray::addCnt_to_CntSt_combinedRead($countArrayRef->[$i+$Num_of_Matrix],\@positions_r,$cut_values_Ref->[$i],$liglength);
        }

    return $countArrayRef;
}







