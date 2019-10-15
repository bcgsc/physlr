#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Handle;
use IO::File;

my $groundTruthFile = "";
my $thresholdEnd    = 1;
my @minOvlLengths   = qw(0 1000 2000 5000 10000);
my $result = GetOptions( 'g=s' => \$groundTruthFile, 't=s' => \$thresholdEnd );

my %groundTruth;

main();

sub main {

	my %groundTruth;
	my $fh = IO::File->new( $groundTruthFile, "r" );
	my $line = $fh->getline();

	#skip header
	$line = $fh->getline();

	my $count = 0;

	#load in ground truth
	while ($line) {
		chomp($line);
		my @split  = split( /\t/, $line );
		my $U      = $split[0];
		my $V      = $split[1];
		my $length = $split[2];

		#convert overlaps to cannoical form
		if ( ( $U cmp $V ) > 0 ) {
			my $temp = $U;
			$U = $V;
			$V = $temp;
		}

		$groundTruth{ $U . "\t" . $V } += $length;
		$line = $fh->getline();
		$count++;
	}

	#print $count . " " . (scalar(keys %groundTruth)) . "\n";
	$fh->close();

	#compare ground truth overlaps to list of files (ARGV)
	foreach my $minLen (@minOvlLengths) {
		for ( my $i = 1 ; $i <= $thresholdEnd ; $i++ ) {
			foreach my $file (@ARGV) {

				#summary table: filename, threshold, TP, FP, TN, FN, Sens, FDR
				my $TP   = 0;
				my $FP   = 0;
				my $fh1  = IO::File->new( $file, "r" );
				my $line = $fh1->getline();

				while ($line) {
					my @split = split( /\t/, $line );
					my $U     = $split[0];
					my $V     = $split[1];
					if ( scalar(@split) > 2 ) {
						last;
					}
					$line = $fh1->getline();
				}

				#skip header
				$line = $fh1->getline();

				#				my %dup;

				while ($line) {
					chomp($line);
					my @split = split( /\t/, $line );
					my $U     = $split[0];
					my $V     = $split[1];
					my $n     = $split[2];

					if ( $n >= $i ) {

						#convert overlaps to cannoical form
						if ( ( $U cmp $V ) > 0 ) {
							my $temp = $U;
							$U = $V;
							$V = $temp;
						}
						if ( exists( $groundTruth{ $U . "\t" . $V } ) ) {
							if ( $groundTruth{ $U . "\t" . $V } > $minLen ) {
								$TP++;
							}
						}
						else {
							print $line . "\n";
							$FP++;
						}
					}
					$line = $fh1->getline();
				}
				if ( $FP + $TP != 0 ) {
					my $totalEdges = 0;
					foreach my $key ( keys %groundTruth ) {
						if ( $groundTruth{$key} > $minLen ) {
							$totalEdges++;
						}
					}
					my $FN   = $totalEdges - $TP;
					my $FDR  = $FP / ( $FP + $TP );
					my $sens = $TP / $totalEdges;
					print "$file\t$TP\t$FP\t$FN\t$totalEdges\t$FDR\t$sens\t$i\t$minLen\n";
				}
				$fh1->close();
			}
		}
	}
}

#repeat for every threshold level of the overlaps

