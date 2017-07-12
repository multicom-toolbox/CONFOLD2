#!/usr/bin/perl -w
# Badri Adhikari, 4/15/2016

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;

my $inpdb  = shift;
my $start  = shift;
my $end    = shift;
my $outpdb = shift;

if (not $inpdb){
	print "Purpose: Trim a PDB chain starting from [start,end] residue.\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <indexStart> <indexEnd> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 2 200 abc-new.pdb\n";
	exit 1;
}
if (not $start){
	print "ERROR! start Index not supplied!\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <indexStart> <indexEnd> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 2 200 abc-new.pdb\n";
	exit 1;
}
if (not $end){
	print "ERROR! end Index not supplied!\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <indexStart> <indexEnd> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 2 200 abc-new.pdb\n";
	exit 1;
}
if (not $outpdb){
	print "ERROR! out-pdb not supplied!\n";
	print "Purpose: Trim a PDB chain starting from [start,end] residue.\n";
	exit 1;
}

my @pdb_list = ();
if (-f $inpdb){
	push @pdb_list, $inpdb;
}
else{
	@pdb_list = load_pdb($inpdb);
	mkdir $outpdb or confess $! if not -d $outpdb;
}

foreach my $p (@pdb_list){
	if (-d $outpdb){
		my $id = basename($p);
		trim_pdb($p, $start, $end, "$outpdb/$id");
		print "Trimmed $p to $outpdb/$id\n";
	}
	else{
		trim_pdb($p, $start, $end, $outpdb);
		print "Trimmed $p to $outpdb\n";
	}
}

sub trim_pdb{
	my $inPDB  = shift;
	my $start  = shift;
	my $end    = shift;
	my $outPDB = shift;

	open PDBFILE, $inPDB or confess $!;
	my @pdbLines = <PDBFILE>;
	close PDBFILE;

	open OUTPDB, ">$outPDB" or confess $!;
	foreach (@pdbLines) {
		next if $_ !~ m/^ATOM/;
		my $row = $_;
		my $this_rnum = parse_pdb_row($_,"rnum");
		next if $this_rnum < $start;
		next if $this_rnum > $end;
		print OUTPDB $row;
	}
	close OUTPDB;
}

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

sub load_pdb{
	my $dir_chains = shift;
	confess "\n:( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdb_list = <$dir_chains/*.pdb>;
	if(not (@pdb_list)){
		@pdb_list = <$dir_chains/*.ent>;
	}
	confess "\nERROR! Directory $dir_chains has no pdb files!\n" unless(@pdb_list);
	return @pdb_list;
}

