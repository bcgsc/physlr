#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Handle;
use IO::File;

my $whiteList = "";
my $result = GetOptions( 'w=s' => \$whiteList );
main();

#load in whitelsit file

sub main {
	my %whiteListSet;
	my $fh1 = IO::File->new($whiteList, "r");
	my $line = $fh1->getline();
	while ($line) {
		chomp($line);
		$whiteListSet{$line} = 1;
		$line = $fh1->getline();
	}
	$fh1->close();

	my $fh = IO::Handle->new();
	$fh->fdopen( fileno(STDIN), "r" );
	$line = $fh->getline();

	#print header
	print $line;
	$line = $fh->getline();

	while ($line) {
		chomp($line);
		my @split = split( /\t/, $line );
		my $U     = $split[0];
		my $V     = $split[1];
		if (   exists( $whiteListSet{$U} )
			&& exists( $whiteListSet{$V} ) )
		{
			if ( ($U cmp $V) > 0 ) {
				my $temp = $U;
				$U = $V;
				$V = $temp;
			}
			print $U . "\t" . $V . "\t". $split[2] . "\n";
		}
		$line = $fh->getline();
	}
	$fh->close();
}
