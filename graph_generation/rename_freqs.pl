#!/usr/bin/perl
##############################################################################
# SCRIPT NAME:	rename_freqs.pl
# DESCRIPTION:	rename freqs.POP.csv.gz to POPS.freqs.gz
#               run this in the frequencies directory if you have frequencies
#               that need to be renamed
#
# DATE WRITTEN: 2021-04-06
# WRITTEN BY:   mmaiers
#
##############################################################################
use strict;    # always
use warnings;  # or else

foreach my $file (`/bin/ls -1 *.csv.gz`) {
  chomp $file;
  my ($f, $pop, $csv, $gz)  = split /\./, $file;
  my $newfile = join ('.', $pop, $f, $gz);
  my $cmd = "/bin/mv $file $newfile";
  print $cmd, "\n";
}
exit 0;
