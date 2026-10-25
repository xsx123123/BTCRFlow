#!/usr/bin/env perl
# Vgene_Freq.pl
# Bin V-gene frequency percentages into divergence categories.
# Usage: perl Vgene_Freq.pl <input_file> <output_prefix>
# Input: tab-separated with V-gene in col1, frequency% in col2

use strict;
use warnings;

my ($input_file, $output_prefix) = @ARGV;
die "Usage: perl Vgene_Freq.pl <input_file> <output_prefix>\n" unless $input_file && $output_prefix;

open(my $in, "<", $input_file) or die "Cannot open $input_file: $!\n";
my $outfile = "${output_prefix}_Vgene_Freq.xls";
open(my $out, ">", $outfile) or die "Cannot open $outfile: $!\n";

print $out "Vgene\tFrequency\tDivergence_class\n";

my @classes = (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27);

while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^#/ || $line =~ /^Vgene/;
    my @fields = split(/\t/, $line);
    next unless @fields >= 2;

    my ($vgene, $freq) = @fields[0, 1];
    $freq = 0 unless defined $freq && $freq =~ /^[\d.]+$/;

    my $class = 27;
    for my $c (@classes) {
        if ($freq <= $c) {
            $class = $c;
            last;
        }
    }

    print $out "$vgene\t$freq\t$class\n";
}

close($in);
close($out);
print "V-gene frequency table written: $outfile\n";
