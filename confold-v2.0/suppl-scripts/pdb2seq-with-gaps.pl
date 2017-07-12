#!/usr/bin/perl -w
# Badri Adhikari, 3/26/2016
# PDB to Sequence - with gaps

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;

my $pdb  = shift;

if (not $pdb){
	print STDERR "Usage  : $0 <pdb>\n";
	print STDERR "Example: $0 ./abc.pdb\n";
	exit 1;
}
if (not -f $pdb){
	print STDERR "Usage  : $0 <pdb>\n";
	print STDERR "Example: $0 ./abc.pdb\n";
	exit 1;
}

our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;
our @AA3;
our @AA1;

print "".seqChainWithGaps($pdb)."\n";

################################################################################
sub seqChainWithGaps{
	my $chain = shift;

	# 1.find end residue number
	my $start; my $end;
	open CHAIN, $chain or die "ERROR! Could not open $chain";
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		$end = parse_pdb_row($_,"rnum");
	}
	close CHAIN;

	# 2.initialize
	my $seq = "";
	for (my $i = 1; $i <= $end; $i++){
		$seq .= "-";
	}
 
	# 3.replace with residues
	open CHAIN, $chain or die "ERROR! Could not open $chain";
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		substr $seq, (parse_pdb_row($_,"rnum") - 1), 1, $AA3TO1{parse_pdb_row($_,"rname")}; 
	}
	close CHAIN;

	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}

################################################################################
sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}
