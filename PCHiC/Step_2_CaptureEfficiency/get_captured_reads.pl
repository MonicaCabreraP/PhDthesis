#!/usr/bin/perl

#############################################
#Reads a SAM/BAM file and a baits file
#and separates captured / non captured reads
#Summarises capture efficiency
#############################################

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Math::Round;

use Data::Dumper;

my $version = '0.1.1';

##############################################
#Option variables
my %config = (
	display_version => '',
	baits => '',
	header => '',
	help => '',
	interactions => ''
);
	
my $config_result = GetOptions(
	"baits=s"	=> \$config{baits},
	"header=i" => \$config{header},
    "help"     	=> \$config{help},
	"interactions" => \$config{interactions},
	"version"	=> \$config{display_version},
);
die "Could not parse options\n" unless ($config_result);

if ( $config{help} ) {    #Print help and exit
    print while (<DATA>);
    exit(0);
}

if ( $config{display_version} ) {    #Print help and exit
    print "get_captured_reads v$version\n";
    exit(0);
}


print "Running get_captured_reads.pl v$version\n";

##############################################
#Check input
unless(@ARGV){
warn "Please specify SAM/BAM files to process.\n";
    print while (<DATA>);
    die "\n";
}
die "Please specify a --baits regions file (tab-delimited: chr  start  end) (see --help for more details)\n" if $config{baits} eq '';
die "Please adjust configuration.\n" unless ( check_files_exist(\@ARGV, 'EXISTS') );
die "Capture file '$config{baits}' does not exits, please adjust configuration.\n" unless ( -e $config{baits} );


##############################################
#Check output files don't already exist
checkNoOutfiles(@ARGV);

##############################################
#Read in capture file

print "Reading baits file '$config{baits}'\n";


if ( $config{baits} =~ /\.gz$/ ) {
	open( BAITS, "zcat $config{baits} |" ) or die "Couldn't read file '$config{baits}' : $!";
} else {
	open(BAITS, '<', $config{baits}) or die "Could not open '$config{baits}' : $!";
}

if($config{header}){    #Number of lines to skip
	$config{header} = abs($config{header});    #Ensure +ve
	for (1..$config{header}){
		scalar <BAITS>;
	}
}

my %baits;    #%{chromosome}->{10kb region}->{Start} = End
my %bait_ids;    #%{chr\tstart\tend} = id

while(<BAITS>){
	my $line = $_;
	chomp $line;
	next if $line =~ /^\s*$/;    #Ignore blank lines
	my ($csome, $start, $end) = split(/\t/, $line);
	add_bait($csome, $start, $end, \%baits);
	
	#If calculating interactions
	if($config{interactions}){
		my $id = (split(/\t/, $line))[3];	
		$bait_ids{"$csome\t$start\t$end"} = $id;
		die "May not name baits '0' or 'NOT_CAPTURED' or '', please edit '$config{baits}'.\n" if( ($id eq '0') or ($id eq 'NOT_CAPTURED') or ($id eq '') );
	}	
}
close BAITS or die "Could not close '$config{baits}' : $!";


#Create a hash of all possible bait-bait interactions (ordered readF and readR)

# foreach my $pos1 (keys %bait_ids){
# 	my $id1 = $bait_ids{$pos1};	
# 	foreach my $pos2 (keys %bait_ids){
# 		my $id2 = $bait_ids{$pos2};		
# 		$bait_bait_interactions{"$id1\t$id2"} = 0;
# 		$bait_bait_interactions{"$id1\tNOT_CAPTURED"} = 0;
# 		$bait_bait_interactions{"NOT_CAPTURED\t$id2"} = 0;
# 		$bait_bait_interactions{"NOT_CAPTURED\tNOT_CAPTURED"} = 0;
# 	}	
# }


###########################################################
#Create R script for producing charts of results
createRscript();

##############################################
#Read in SAM/BAM files and identify captured/uncaptured ditags
foreach my $file (@ARGV){

	print "Processing $file\n";

	my %bait_bait_interactions;    #%{bait1_chr\tbait1_start\tbait1_end\tbait2_chr\tbait2_start\tbait2_end} = count

	#Categorise and record the data
	#initialise_hash(\%bait_bait_interactions);    #Sets all terms to 0
	
	my %results_counter;
	my @results_categories = ('Ditags_Processed', 'Both_Captured', 'Forward_Captured',
								'Reverse_Captured', 'Uncaptured', 'Captured_Cis',
								'Captured_Trans', 'Percent_Both_Captured', 
								'Percent_Forward_Captured', 'Percent_Reverse_Captured', 
								'Percent_Uncaptured', 'Captured_Cis', 'Captured_Trans',
								'Total_Captured', 'Percent_Total_Captured'
							);
	
	foreach my $category (@results_categories){
		$results_counter{$category} = 0;
	}
			
	#Open input file and write to output files
	if( $file =~ /\.bam$/ ) {
		open( IN, "samtools view -h $file |" ) or die "Couldn't read '$file' : $!";
	}elsif( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }
	
	my $captured_filename = determineOutputFilenames($file, 'CAPTURED');
	my $uncaptured_filename = determineOutputFilenames($file, 'UNCAPTURED');
	my $summary_filename = determineOutputFilenames($file, 'SUMMARY');
	
	open(CAPTURED, "| samtools view -bSh - > $captured_filename"  ) or die "Could not write to '$captured_filename' : $!";
	open(UNCAPTURED, "| samtools view -bSh - > $uncaptured_filename"  ) or die "Could not write to '$uncaptured_filename' : $!";	
	open(SUMMARY, '>', $summary_filename) or die "Could not open '$summary_filename' : $!";
	
	{    #Print header to summary file
		my $summary_header = "File\t";
		foreach my $category (@results_categories){
			$summary_header = $summary_header . "$category\t";
		}
		$summary_header =~ s/\t$/\n/;
		print SUMMARY $summary_header;
	}

	while(<IN>){
		my $readF = $_;    #First in sequenced pair (R1 file)
		chomp $readF; 
		if( substr($readF, 0, 1) eq '@'){
			print CAPTURED "$readF\n";
			print UNCAPTURED "$readF\n";
			next;
		}
		
		my $readR = scalar <IN>;   #Second in sequenced pair (R2, R3, R4 file)
		chomp $readR;
		
		my ($csomeF, $midposF) = samMidPoint($readF);    #Convert to csome and read mid-point position
		my ($csomeR, $midposR) = samMidPoint($readR);
		
		
		#Are either of the reads captured?
		my $readF_maps_to_bait_pos = coord2bin($csomeF, $midposF, \%baits);
		my $readR_maps_to_bait_pos = coord2bin($csomeR, $midposR, \%baits);
		my $is_capturedF = 0; 
		my $is_capturedR = 0; 
		$is_capturedF = 1 unless ( $readF_maps_to_bait_pos eq '0' );
		$is_capturedR = 1 unless ( $readR_maps_to_bait_pos eq '0' );
		my $bait_bait_interaction;
		
		$results_counter{Ditags_Processed}++;
			
		if($is_capturedF and $is_capturedR){    #Both reads captured				
			print CAPTURED "$readF\n$readR\n";
			$results_counter{Both_Captured}++;
			$bait_bait_interaction = "$bait_ids{$readF_maps_to_bait_pos}\t$bait_ids{$readR_maps_to_bait_pos}" if($config{interactions});
				
		}elsif($is_capturedF){	
			print CAPTURED "$readF\n$readR\n";
			$results_counter{Forward_Captured}++;
			$bait_bait_interaction = "$bait_ids{$readF_maps_to_bait_pos}\tNOT_CAPTURED" if($config{interactions});
		
		}elsif($is_capturedR){
			print CAPTURED "$readF\n$readR\n";
			$results_counter{Reverse_Captured}++;
			$bait_bait_interaction = "NOT_CAPTURED\t$bait_ids{$readR_maps_to_bait_pos}" if($config{interactions});
				
		}else{
			print UNCAPTURED "$readF\n$readR\n";
			$results_counter{Uncaptured}++;	
			$bait_bait_interaction = "NOT_CAPTURED\tNOT_CAPTURED" if($config{interactions});
			$bait_bait_interactions{$bait_bait_interaction}++ if($config{interactions});    #Need to count here for uncaptured di-tags
			next;
		}
		
		#Is this a cis or trans ditag? (Will not get to this part of script if not captured)
		if( (split(/\t/, $readF))[6] eq '=' ){
			$results_counter{Captured_Cis}++;
		}else{
			$results_counter{Captured_Trans}++;
		}

		$bait_bait_interactions{$bait_bait_interaction}++ if($config{interactions});	
	}

	#Determine results
	$results_counter{Total_Captured} = $results_counter{Both_Captured} + $results_counter{Forward_Captured} + $results_counter{Reverse_Captured};
	$results_counter{Percent_Total_Captured} = calc_perc( $results_counter{Total_Captured}, $results_counter{Ditags_Processed} );
	$results_counter{Percent_Both_Captured} = calc_perc( $results_counter{Both_Captured}, $results_counter{Ditags_Processed} );
	$results_counter{Percent_Forward_Captured} = calc_perc( $results_counter{Forward_Captured}, $results_counter{Ditags_Processed} );
	$results_counter{Percent_Reverse_Captured} = calc_perc( $results_counter{Reverse_Captured}, $results_counter{Ditags_Processed} );
	$results_counter{Percent_Uncaptured} = calc_perc( $results_counter{Uncaptured}, $results_counter{Ditags_Processed} );
	$results_counter{Percent_Captured_Cis} = calc_perc( $results_counter{Captured_Cis}, $results_counter{Total_Captured} );	
	$results_counter{Percent_Captured_Trans} = calc_perc( $results_counter{Captured_Trans}, $results_counter{Total_Captured} );	
		
	{    #Print results to the summary file
		my $summary_results = "$file\t";
		foreach my $category (@results_categories){
			$summary_results = $summary_results . "$results_counter{$category}\t";
		}
		$summary_results =~ s/\t$/\n/;
		print SUMMARY $summary_results;
	}
	
	close IN or warn "Could not close '$file' : $!";
	close CAPTURED or warn "Could not close '$captured_filename' : $!";
	close UNCAPTURED or warn "Could not close '$uncaptured_filename' : $!";
	close SUMMARY or warn "Could not close '$summary_filename' : $!";
	
	
	#Create summary chart
	my $chart_filename = determineOutputFilenames($file, 'CHART');
	!system("Rscript capture_chart_generator.r $summary_filename $file $chart_filename >/dev/null") or warn "Could not execute command 'Rscript capture_chart_generator.r $summary_filename $file $chart_filename'\n";


	#Print out interactions results
	if($config{interactions}){
		
		my $interactions_ordered_filename = determineOutputFilenames($file, 'INTERACTIONS_ORDERED');
		open( INTERACTIONS_ORDERED, '>', $interactions_ordered_filename ) or die "Could not write to '$interactions_ordered_filename' : $!";
		print INTERACTIONS_ORDERED "Bait1\tBait2\tCount\n";
		foreach my $interaction (keys %bait_bait_interactions){
			my $count = $bait_bait_interactions{$interaction};
			print INTERACTIONS_ORDERED "$interaction\t$count\n";
		}
		close INTERACTIONS_ORDERED or die "Could not close filehandle on '$interactions_ordered_filename' : $!";



		my %unordered_bait_bait_interactions = unorder_interactions(%bait_bait_interactions);
		my $interactions_not_ordered_filename = determineOutputFilenames($file, 'INTERACTIONS_NOT_ORDERED');
		open( INTERACTIONS_NOT_ORDERED, '>', $interactions_not_ordered_filename ) or die "Could not write to '$interactions_not_ordered_filename' : $!";
		print INTERACTIONS_NOT_ORDERED "Bait1\tBait2\tCount\n";
		foreach my $interaction (keys %unordered_bait_bait_interactions){
			my $count = $unordered_bait_bait_interactions{$interaction};
			print INTERACTIONS_NOT_ORDERED "$interaction\t$count\n";
		}
		close INTERACTIONS_NOT_ORDERED or die "Could not close filehandle on '$interactions_ordered_filename' : $!";
	}
}

unlink 'capture_chart_generator.r' or warn "Could not delete capture_chart_generator.r\n";

print "Processing complete.\n";

exit (0);




##########################################################################
#Subroutines
##########################################################################


#Subroutine unorder_interactions
#Takes the %bait_bait_interactions hash and returns a hash
#in which the order of the reads constituting the di-tag is unimportant
sub unorder_interactions{
	my %ordered = @_;
	my %unorderd;

	foreach my $interaction (keys %ordered){
		my $count = $ordered{$interaction};
		my ($bait1, $bait2) = split(/\t/, $interaction);
		my $interaction;

		if( ($bait1 cmp $bait2) == 1){
			$interaction = "$bait2\t$bait1";
		}else{
			$interaction = "$bait1\t$bait2";
		}
		$unorderd{$interaction} += $count;
	}
	return %unorderd;
}



#Subroutine: initialise_hash
#Takes a reference to a hash and sets all its values to 0,
# sub initialise_hash{
# 	my $hash_ref = $_[0];	
# 	foreach my $key ( keys %{ $hash_ref } ){
# 		$$hash_ref{$key} = 0;
# 	}
# }


#Subroutine: calc_perc
#Receives a number and a total and returns the perentage value
#Optional argument: decimal places of the output
sub calc_perc {
    my ( $n, $total, $dp ) = @_;
	
	if(defined $dp){
		$dp = abs( int($dp) );
	}else{
		$dp = 2;
	}
	$dp = 1 / (10**$dp);   #'Nearest' function needs 1, 10, 100 etc. 
	
	return 'NA' unless(defined $n and defined $total);   #Avoid initialisation error
	return 'NA' if($total == 0);    #Avoid division by zero error
	
	my $pc = 100 * $n / $total;
	$pc = nearest($dp, $pc);
	
	return $pc;
}


##########################################
#Subroutine: determineOutputFilenames
#Receives i) input filename and ii) category of the file to return
#(can be either 'SUMMARY', 'CAPTURED', 'UNCAPTURED' 'CHART')
#and returns the output filename
sub determineOutputFilenames{
	my $filename = $_[0];
	my $type_required = $_[1];
	
	$filename =~ s/\.bam$|\.sam$//g;
	
	if($type_required eq 'SUMMARY'){
		return "$filename.capture_summary.txt";
	}elsif($type_required eq 'CAPTURED'){
		return "$filename.captured.bam";	
	}elsif($type_required eq 'UNCAPTURED'){
		return "$filename.uncaptured.bam";	
	}elsif($type_required eq 'CHART'){
		return "$filename.capture_charts.pdf";
	}elsif($type_required eq 'INTERACTIONS_ORDERED'){
		return "$filename.interactions_ordered.txt";	
	}elsif($type_required eq 'INTERACTIONS_NOT_ORDERED'){
		return "$filename.interactions_not_ordered.txt";	
	}else{
		die "Subroutine 'determineOutputFilenames' passed invalid argument '$type_required'.\n";
	}	
}


##########################################
#Subroutine: deduplicate_array
#Takes and array and returns the array with duplicates removed
#(keeping 1 copy of each unique entry).
sub deduplicate_array {
    my @array = @_;
    my %uniques;

    foreach (@array) {
		$uniques{$_} = '';
    }
    my @uniques_array = keys %uniques;
    return @uniques_array;
}


#Subroutine: checkNoOutfiles
#Takes and array of all the input files and checks that
#none of the outputfiles already exits or are duplicates
sub checkNoOutfiles{
	my @inputfiles = @_;
	my @outputfiles;
	
	my $parameters_ok = 1;
	foreach my $inputfile (@inputfiles){
		push( @outputfiles, determineOutputFilenames($inputfile, 'CAPTURED') );
		push( @outputfiles, determineOutputFilenames($inputfile, 'UNCAPTURED') );
		push( @outputfiles, determineOutputFilenames($inputfile, 'CHART') );
		push( @outputfiles, determineOutputFilenames($inputfile, 'SUMMARY') );
		push( @outputfiles, determineOutputFilenames($inputfile, 'INTERACTIONS_ORDERED') );
		push( @outputfiles, determineOutputFilenames($inputfile, 'INTERACTIONS_NOT_ORDERED') );
	}
	
	die "Outputfile(s) already exist, please delete.\n" unless check_files_exist(\@outputfiles, 'NOT_EXISTS');

	#Check all the outputfiles have unique names
	my @uniques = deduplicate_array(@outputfiles);
	die "Some outputfile(s) will have identical names, rename inputfiles\n." if(scalar @outputfiles != scalar @uniques);
}


##########################################
#Subroutine: check_files_exist
#Takes a reference to an array containing paths to filenames and verifies they exist
#Warns of files that do no exit. Returns 1 if all files exist but 0 if this is not
#the case.
#
#Also, takes a second argument:
#$_[1] should be 'EXISTS' or 'NOT_EXISTS'
#If 'NOT_EXISTS' warns if file already exists.  Returns '1' if none of the
#files exists and '0' if one or multiple files already exist
sub check_files_exist {
    my $files      = $_[0];    #Reference to array
    my $check_for  = $_[1];
    my $all_exist  = 1;
    my $not_exists = 1;

    if ( $check_for eq 'EXISTS' ) {
        foreach my $file (@$files) {
            unless ( -e $file ) {
                warn "File '$file' does not exist\n";
                $all_exist = 0;
            }
        }
    } elsif ( $check_for eq 'NOT_EXISTS' ) {
        foreach my $file (@$files) {
            if ( -e $file ) {
                warn "File '$file' already exists\n";
                $not_exists = 0;
            }
        }
    } else {
        die "Subroutine 'check_files_exist' requires argument 'EXISTS' or 'NOT_EXISTS'.\n";
    }

    if ( $check_for eq 'EXISTS' ) {
        return $all_exist;
    } else {
        return $not_exists;
    }
}



##########################################
#Subroutine: add_bait
#Takes the bait chromosome/start/end
#and populates the passed hash accordingly:
#%{chromosome}->{10kb region}->{Start} = End
#Note: if the bin/fragment spans more than one 10kb region,
#then multiple 10 regions will be populated
sub add_bait {
	my ($csome, $start, $end, $hash_ref) = @_;
	
	my $ten_kb_start = ceil($start / 10_000);
	my $ten_kb_end = ceil($end/ 10_000);
	
	for (my $ten_kb_region = $ten_kb_start; $ten_kb_region <= $ten_kb_end; $ten_kb_region++){
		${$hash_ref}{$csome}->{$ten_kb_region}->{$start} = $end;
	}
}


##########################################
#Subroutine: samMidPoint
#Receives a SAM format line and returns
#chromosome, midpoint of the position of then
#read
sub samMidPoint {
	my $read = $_[0];

	my $csome = ( split( /\t/, $read ) )[2];
	my $start_genome_perspective = ( split( /\t/, $read ) )[3];
	my $seq   = ( split( /\t/, $read ) )[9];
	
	my $length = length($seq);	
	my $midpoint = $start_genome_perspective + ceil($length / 2);

	return ( $csome, $midpoint );
}


##########################################
#Subroutine: coord2bin
#Receives a chromosome name and a position and reference to the baits hash
#and returns the bait co-ordinates where this location is found (else returns 0)
#%lookup_hash{chromosome}->{10kb region}->{Start} = End
sub coord2bin{
	my ($csome, $pos, $hash_ref) = @_;
	my $ten_kb_region = ceil($pos / 10_000);

	foreach my $start ( keys %{ $baits{$csome}->{$ten_kb_region} }){
		my $end = ${ $hash_ref }{$csome}->{$ten_kb_region}->{$start};
		if ( ($start <= $pos) and ($end >= $pos) ){
			return ("$csome\t$start\t$end");
		}
	}
	return 0;    #Not captured
}


##########################################
#Subroutine: createRscript
#Produces the R script that is used to produce graphics of the results
sub createRscript{

	my $r_script_text = <<'R_SCRIPT';
args <- commandArgs(TRUE)

summaryFile <- args[1]
sampleName <- args[2]
outputfilename <- args[3]

pdf(file=outputfilename, paper="a4")

data <- read.delim(summaryFile, header=FALSE, skip=1)
ditags_processed <- data[,2]
both_captured <- data[,3]
forward_captured <- data[,4]
reverse_captured <- data[,5]
uncaptured <- data[,6]
total_captured <- data[,15]
percent_total_captured <-data[,16]
percent_captured_trans <- round( (100 * data[,14] / total_captured), 1 )

pieTitle <- paste( "Proportion Captured Di-tags\n", "Percent Captured: ", percent_total_captured, "%",
                   "\nPercent Captured Trans: ", percent_captured_trans, "%",  sep = "")

pcData <- c(both_captured, forward_captured, reverse_captured, uncaptured)
percLabels <- round(pcData / ditags_processed * 100, 1)
percLabels<- paste(percLabels, "%", sep="")

pie  (pcData, 
      labels=percLabels,
      main=pieTitle,
      cex.main = 1,
      col = rainbow(4),
)

legend("bottom", c("Both Captured", "Forward Captured", "Reverse Captured", "Uncaptured"), 
       ncol=2, cex=0.8, fill=rainbow(4))

dev.off()
R_SCRIPT

	open( RSCRIPT, '>', 'capture_chart_generator.r') or warn "Could not write to 'capture_chart_generator.r' : $!";
	print RSCRIPT $r_script_text;
	close RSCRIPT or warn "Could not close filehandle on 'capture_chart_generator.r' : $!";
}



__DATA__

SYNOPSIS

get_captured_reads.pl

get_captured_reads.pl [OPTIONS] --baits [baits file] [BAM/SAM files]
get_captured_reads.pl [OPTIONS]

FUNCTION

Takes a baits file and BAM/SAM Hi-C files (output from HiCUP) and separates 'captured'
di-tags from 'uncaptured' di-tags, writing into two different BAM files. 
Reports summary statistics on the results.

The baits file should be a tab-delimited file of format:
Chromosome    Start    End
 
COMMAND LINE OPTIONS

--baits            Baits format file
--header           Specify number of header lines in the baits file (i.e. skip these)
--help             Print help message and exit
--interactions     Calculate interaction frequecies between baits
--version          Print the program version and exit

Steven Wingett, Babraham Institute, Cambridge, UK (steven.wingett@babraham.ac.uk)
