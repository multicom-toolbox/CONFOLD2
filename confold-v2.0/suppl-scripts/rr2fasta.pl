#!/usr/bin/perl
# 12/14/2016, Badri Adhikari

use strict;
use warnings;
use Carp;
use File::Basename;
use File::Temp qw(tempfile);

my $rr_file    = shift;

if (not $rr_file){
	print STDERR "Usage: $0 <RR>\n";
	exit 1;
}

if (not -f $rr_file){
	print STDERR "ERROR! Input RR file does not exist!\n";
	print STDERR "RR: $rr_file\n" if $rr_file;
	print STDERR "Usage: $0 <RR>\n";
	exit 1;
}

my $id = basename($rr_file);
my $seq = seq_rr($rr_file);

print ">".$id."\n";
print $seq."\n";

sub seq_rr{
	my $rr_file = shift;
	confess "ERROR! Input file $rr_file does not exist!" if not -f $rr_file;
	my $seq;
	open RR, $rr_file or confess "ERROR! Could not open $rr_file! $!";
	while(<RR>){
		chomp $_;
		$_ =~ s/\r//g; # chomp does not remove \r
		$_ =~ s/^\s+//;
		$_ =~ s/\s+//g;
		next if ($_ =~ /^>/);
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	confess "Input RR file does not have sequence row!\nRR-file : $rr_file\nPlease make sure that all input RR files have sequence headers!" if not defined $seq;
	return $seq;
}
