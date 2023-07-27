#!/usr/bin/perl
# Analyzes coverage for each promoter from a de novo promoter bed6 file. Output is a bed6 file with coordinates of start site and then FC in coverage in the 7th column.
# NOTES: 
    # 1. $gene_annotation_file is my bed file that has a weird extra column for transcript (col 4) then gene (col 5).
    # 2. $start_site_bed_file is a bed 6 file with coordinates of de novo promoters
    # 3. $coverage_distance is the number of bases downstream from start site to use to quantify differential coverage
    # 4. $min_distance_of_start_from_exon is minimum distance that an annotated exon must be downstream from the de novo start site in order to be included in the output. Also NOTE that if there is an annotated start site within the $coverage_distance, differential coverage will be analyzed from the de novo promoter start site to the start of that exon.
    # 5. IMPORTANT: @bws must be in the order of cntl_bw_strand_1, cntl_bw_strand_2, test_bw_strand_1, test_bw_strand_2 for each set of samples to be analyzed (multiple experiments can be analyzed at a time as long as they are in this order for each set).

# USAGE:
	# perl denovo_promoter_individual_coverage_values.pl $gene_annotation_file, $start_site_bed_file, $coverage_distance, $min_distance_of_start_from_exon, $experiment_name, @bws

# EKF (Flemington lab) 12/10/2020.

use warnings;
use strict;  
#use diagnostics;
use File::Basename;
use List::MoreUtils qw(uniq);

my ($gene_annotation_file, $start_site_bed_file, $coverage_distance, $min_distance_of_start_from_exon, $experiment_name, @bws) = @ARGV;
my $main_directory = dirname(__FILE__);
my $output_directory = $main_directory."\/".$experiment_name."_output";
my $feature_coverage_from_wig_bed12 = dirname(__FILE__)."\/1_1_bin\/feature_coverage_from_wig_bed12";
my $start_site_bed_basename = basename($start_site_bed_file);
my $bed_file_for_quantification = $output_directory."/".$start_site_bed_basename."_for_quantification.bed";
my @experiment_basename_array = ();
my $number_of_denovo_promoters = 0;
my @total_annotated_coverage_array = ();
my @denovo_annotation_information = ();
my @denovo_coverage_data_across_samples = ();
my $first_sample_count = 0;

sub make_directories {
    if (-e "$output_directory") {
        print "\n- directory exists...skipping mkdir (output directory)\n";
    }
    else {
        `mkdir $output_directory`;
    }
    if (-e "$output_directory/1_annotated_transcript_files") {
        print "\n- directory exists...skipping mkdir (annotated transcript files)\n";
    }
    else {
        `mkdir $output_directory/1_annotated_transcript_files`;
    }
    if (-e "$output_directory/2_annotated_genes_coverage_files") {
        print "- directory exists...skipping mkdir (annotated genes coverage files)\n";
    }
    else {
        `mkdir $output_directory/2_annotated_genes_coverage_files`;
    }
    if (-e "$output_directory/output_files") {
        print "- directory exists...skipping mkdir (output files)\n";
    }
    else {
        `mkdir $output_directory/output_files`;
    }
}
sub print_output_parameters_file {
    my $parameters_file = $output_directory."/1_2_".$experiment_name. "_run_paramaters.txt";
    open (OUT, ">$parameters_file") or die "couldn't open file";
    print OUT "Experiment name: ", $experiment_name, "\n\nGene Annotation file: ", basename($gene_annotation_file), "\nStart site bed file: ", basename($start_site_bed_file), "\nLength of downstream sequence for coverage analysis: ", $coverage_distance, "\nMinimum distance to annotated exon to be included in output: ", $min_distance_of_start_from_exon, "\n\nInput bw files:\n", join("\n", basename(@bws));
}
sub process_annotation_file {
    print "\n- Processing annotated gene file\n";
    my $file = $_[0];
    open (INF, "<$file") or die "couldn't open input file";
    open (OUT1, ">$output_directory/1_annotated_transcript_files/annotated_gene_starts.bed") or die "couldn't open file";
    open (OUT2, ">$output_directory/1_annotated_transcript_files/annotated_genes_bed12.bed") or die "couldn't open file";
    open (OUT3, ">$output_directory/1_annotated_transcript_files/annotated_exons.bed") or die "couldn't open file";
    print OUT1 "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n";
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        next if ($. == 1);
        next if ($split_line[4] =~ m/^SNORA|^MIR/g); # Delete SNORA and microRNA genes from annotation

        if ($split_line[5] eq "+") {
        print OUT1 $split_line[0], "\t", $split_line[1], "\t", $split_line[1], "\t", $split_line[4]."\_".$split_line[3], "\t1000\t", $split_line[5], "\n";
        }
        else {
            print OUT1 $split_line[0], "\t", $split_line[2], "\t", $split_line[2], "\t", $split_line[4]."\_".$split_line[3], "\t1000\t", $split_line[5], "\n";
        }

        my $exon_size_array = $split_line[10];
        my $exon_start_array = $split_line[11];
        $exon_size_array =~ s/\,$//;
        $exon_start_array =~ s/\,$//;

        print OUT2 $split_line[0], "\t", $split_line[1], "\t", $split_line[2], "\t", $split_line[4]."\_".$split_line[3], "\t1000\t", $split_line[5], "\t", $split_line[6], "\t", $split_line[7], "\t", $split_line[8], "\t", $split_line[9], "\t", $exon_size_array, "\t", $exon_start_array, "\n";

        my @split_exon_lengths = split("\,", $split_line[10]);
        my @split_exon_starts = split("\,", $split_line[11]);

        for(my $i = 0; $i < $split_line[9]; $i++) {
            print OUT3 $split_line[0], "\t", $split_line[1]+$split_exon_starts[$i], "\t", $split_line[1]+$split_exon_starts[$i]+$split_exon_lengths[$i], "\t", $split_line[4]."\_".$split_line[3]."\_exon".($i+1), "\t1000\t", $split_line[5], "\n";
        }
    }
    close(INF);
    close(OUT1);
    close(OUT2);
    close(OUT3); 
}
sub process_input_de_novo_promoter_file {
    print "\n- Processing de novo promoter file\n";
    my $file = $_[0];
    my $sorted_start_site_file = $file.".sorted.bed";
    `sort -k1,1 -k2,2n -k3,3n $file > $sorted_start_site_file`;

    my $annotated_exons_file = $output_directory."/1_annotated_transcript_files/annotated_exons.bed";
    my $sorted_annotated_exon_file = $output_directory."/1_annotated_transcript_files/annotated_exons.bed.sorted.bed";
    `sort -k1,1 -k2,2n -k3,3n $annotated_exons_file > $sorted_annotated_exon_file`;

    my $subtracted_denovo_promoter_file = $output_directory."/".$start_site_bed_basename."_bedtools_subtract_start_site_file.tsv";
    `bedtools subtract -a $sorted_start_site_file -b $sorted_annotated_exon_file -s -A > $subtracted_denovo_promoter_file`;

    my $closest_output_file = $output_directory."/".$start_site_bed_basename."_bedtools_closest_exon.tsv";
    `bedtools closest -a $subtracted_denovo_promoter_file -b $sorted_annotated_exon_file -fd -D a -iu -s -t first> $closest_output_file`;

    open (INF, "<$closest_output_file") or die "couldn't open input file";
    open (OUT, ">$bed_file_for_quantification") or die "couldn't open file"; #### Defined in beginning
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        next if ($split_line[12] < $min_distance_of_start_from_exon and $split_line[12] != -1);
        $number_of_denovo_promoters++;
        my $midpoint = ($split_line[1]+$split_line[2])/2;
        my $midpoint_int = int($midpoint);
        if ($split_line[12] > $coverage_distance || $split_line[12] == -1) {
            if ($split_line[5] eq "+") {
                print OUT $split_line[0], "\t", $midpoint_int, "\t", $midpoint_int+$coverage_distance, "\t", join("\t", @split_line[3..5]), "\t", $midpoint_int, "\t", $midpoint_int+$coverage_distance, "\t0\t1\t", ($midpoint_int+$coverage_distance)-$midpoint_int, "\t0\n";
            } 
            elsif ($split_line[5] eq "-") {
                print OUT $split_line[0], "\t", $midpoint_int-$coverage_distance, "\t", $midpoint_int, "\t", join("\t", @split_line[3..5]), "\t", $midpoint_int-$coverage_distance, "\t", $midpoint_int, "\t0\t1\t", $midpoint_int-($midpoint_int-$coverage_distance), "\t0\n";
            } 
        }
        else {
            if ($split_line[5] eq "+") {
                print OUT $split_line[0], "\t", $midpoint_int, "\t", $midpoint_int+$split_line[12], "\t", join("\t", @split_line[3..5]), "\t", $midpoint_int, "\t", $midpoint_int+$split_line[12], "\t0\t1\t", ($midpoint_int+$split_line[12])-$midpoint_int, "\t0\n";
            } 
            if ($split_line[5] eq "-") {
                print OUT $split_line[0], "\t", $midpoint_int-$split_line[12], "\t", $midpoint_int, "\t", join("\t", @split_line[3..5]), "\t", $midpoint_int-$split_line[12], "\t", $midpoint_int, "\t0\t1\t", $midpoint_int-($midpoint_int-$split_line[12]), "\t0\n";
            } 
        }
    }
    close(INF);
    close(OUT);
    `rm $sorted_start_site_file`;
    `rm $sorted_annotated_exon_file`;
}
sub annotated_gene_coverage {
    print "     - Calculating coverage of annotated exons\n";
    my ($bw_cntl_file_str1, $bw_cntl_file_str2) = @_;
    my $sample_basename = basename($bw_cntl_file_str1);
    $sample_basename =~ s/Signal.*$//g;

    print "          - annotated gene coverage for sample: ", $sample_basename, "\n";
    if (-e "$output_directory/2_annotated_genes_coverage_files/$sample_basename.annotated_gene_coverage") {
        print "               - control coverage for annotated genes file already exists...skipping step\n";
        open (INF, "<$output_directory/2_annotated_genes_coverage_files/$sample_basename.annotated_gene_coverage") or die "couldn't open gene annotation coverage file";
    }
    else {
        `python $feature_coverage_from_wig_bed12 $output_directory/1_annotated_transcript_files/annotated_genes_bed12.bed $bw_cntl_file_str1 $bw_cntl_file_str2 > $output_directory/2_annotated_genes_coverage_files/$sample_basename.annotated_gene_coverage`;
        open (INF, "<$output_directory/2_annotated_genes_coverage_files/$sample_basename.annotated_gene_coverage") or die "couldn't open gene annotation coverage file";
    }

    my @annotated_genes_coverage_array_minus_EBV = ();
    while (my $line = <INF>) {
	    chomp($line);
        next if ($. == 1);
        my @split_line = split("\t", $line);
        if ($split_line[0] ne "chrEBV_Akata_inverted") {
            push(@annotated_genes_coverage_array_minus_EBV, $split_line[12]);
        }
    }
    close(INF);

    my $annotated_genes_coverage_minus_EBV_sum = 0;
    for my $value(@annotated_genes_coverage_array_minus_EBV) {
        $annotated_genes_coverage_minus_EBV_sum = $value + $annotated_genes_coverage_minus_EBV_sum;
    }
    push(@total_annotated_coverage_array, $annotated_genes_coverage_minus_EBV_sum);
    print "          Total coverage of annotated genes (minus EBV): ", $annotated_genes_coverage_minus_EBV_sum, "\n";
}
sub de_novo_promoter_coverage {
    print "\n     - Calculating de novo promoter coverage\n";
    my ($bw_cntl_file_str1, $bw_cntl_file_str2) = @_;
    my $experiment_basename = basename($bw_cntl_file_str1);
    $experiment_basename =~ s/Signal.*$//g;
    push(@experiment_basename_array, $experiment_basename);

    print "          - calculation\n";
    `python $feature_coverage_from_wig_bed12 $bed_file_for_quantification $bw_cntl_file_str1 $bw_cntl_file_str2 > $output_directory/output_files/$experiment_basename.de_novo_promoter_specific_coverage`;

    open (INF, "<$output_directory/output_files/$experiment_basename.de_novo_promoter_specific_coverage") or die "couldn't open gene annotation coverage file";
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        next if ($. == 1); 
        push(@denovo_coverage_data_across_samples, $split_line[12]);
        if ($first_sample_count == 0) {
            push(@denovo_annotation_information, join("\t", @split_line[0..5]));
        }
    }
    close(INF);
    $first_sample_count++;
}
sub coverage {
    my $bw_array_count = 0;
    my $strand_1_bw;
    my $strand_2_bw;

    foreach my $bw_name(@bws) {
        $bw_array_count++;
        if($bw_array_count == 1) {
            print "\nStrand 1: ", basename($bw_name);
            $strand_1_bw = $bw_name;
        }
        elsif($bw_array_count == 2) {
            print "\nStrand 2: ", basename($bw_name);
            $strand_2_bw = $bw_name;
            print "\n- Calculating coverages (annotated exons and de novo promoters) for:\n";
            print "  ", basename($strand_1_bw), "\n  ", basename($strand_2_bw), "\n";
            annotated_gene_coverage ($strand_1_bw, $strand_2_bw);
            de_novo_promoter_coverage ($strand_1_bw, $strand_2_bw);
            $bw_array_count = 0;
        }        
    }
    my $number_of_samples = @total_annotated_coverage_array;
    my $total_annotated_coverage_all_samples = 0;
    foreach my $each_annotated_coverage(@total_annotated_coverage_array) {
        $total_annotated_coverage_all_samples = $each_annotated_coverage + $total_annotated_coverage_all_samples;
    }
    my $average_annotated_coverage = $total_annotated_coverage_all_samples/$number_of_samples;
    print "     - Average total annotated coverage across all samples: ", $average_annotated_coverage, "\n";

    my $final_output_file = $output_directory."/output_files/".$experiment_name."_normalized_output.tsv";
    my $final_output_file_test = $output_directory."/output_files/".$experiment_name."_non_normalized_output.tsv";
    open (OUTfinal, ">$final_output_file") or die "couldn't open file";
    open (OUTfinal_test, ">$final_output_file_test") or die "couldn't open file";
    print OUTfinal "chrom\tstart\tend\tname\ttotal CAGE test reads\tstrand\t", join("\t", @experiment_basename_array), "\n";
    print OUTfinal_test "chrom\tstart\tend\tname\ttotal CAGE test reads\tstrand\t", join("\t", @experiment_basename_array), "\n";

    for (my $j = 0; $j < $number_of_denovo_promoters; $j++) {
        print OUTfinal $denovo_annotation_information[$j];
        print OUTfinal_test $denovo_annotation_information[$j];
        for (my $k = 0; $k < $number_of_samples; $k++) {
            print OUTfinal "\t", sprintf("%.2f", ($denovo_coverage_data_across_samples[$k*$number_of_denovo_promoters+$j]*$average_annotated_coverage/$total_annotated_coverage_array[$k]));
            print OUTfinal_test "\t",  sprintf("%2f", $denovo_coverage_data_across_samples[$k*$number_of_denovo_promoters+$j]);
        }
        print OUTfinal "\n";
        print OUTfinal_test "\n";
    }
    close(OUTfinal);
    close(OUTfinal_test);
}

make_directories;
print_output_parameters_file;
process_annotation_file($gene_annotation_file);
process_input_de_novo_promoter_file($start_site_bed_file);
coverage;

`rm -r $output_directory/1_annotated_transcript_files`;
`rm -r $output_directory/2_annotated_genes_coverage_files`;
