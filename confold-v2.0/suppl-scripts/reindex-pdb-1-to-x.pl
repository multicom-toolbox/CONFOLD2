#!/usr/bin/perl -w
# Badri Adhikari, 4/15/2016

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;

my $inpdb  = shift;
my $index  = shift;
my $outpdb = shift;

if (not $inpdb){
	print "Purpose: Reindex a PDB chain starting from X residue.\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <index> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 24 abc-new.pdb\n";
	exit 1;
}
if (not $index){
	print "ERROR! Index not supplied!\n";
	print "Purpose: Reindex a PDB chain starting from X residue.\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <index> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 24 abc-new.pdb\n";
	exit 1;
}
if (not $outpdb){
	print "ERROR! out-pdb not supplied!\n";
	print "Purpose: Reindex a PDB chain starting from X residue.\n";
	print "Usage  : [$0] <in-pdb/in-pdb-dir> <index> <out-pdb-dir>\n";
	print "Example: $0 abc.pdb 24 abc-new.pdb\n";
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
		reindexChain($p, $index, "$outpdb/$id");
		print "Reindexed $p to $outpdb/$id\n";
	}
	else{
		reindexChain($p, $index, $outpdb);
		print "Reindexed $p to $outpdb\n";
	}
}

sub reindexChain{
	my $inPDB = shift;
	my $index = shift;
	my $outPDB = shift;

	open PDBFILE, $inPDB or confess $!;
	my @pdbLines = <PDBFILE>;
	close PDBFILE;

	# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
	my $resCounter = $index - 1;
	my $atomCounter = 0;
	my $prevrNum = "XX";
	open OUTPDB, ">$outPDB" or confess $!;
	foreach (@pdbLines) {
		next if $_ !~ m/^ATOM/;
		my $this_rnum = parse_pdb_row($_,"rnum");
		if ($prevrNum ne $this_rnum) {
			$prevrNum = $this_rnum;
			$resCounter++;
		}
		$atomCounter++;
		my $rnum_string = sprintf("%4s", $resCounter);
		my $anum_string = sprintf("%5s", $atomCounter);
		my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
		print OUTPDB $row;
	}
	print OUTPDB "END\n";
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

