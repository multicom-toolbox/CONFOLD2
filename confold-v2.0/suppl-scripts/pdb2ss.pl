#!/usr/bin/perl -w
# Badri Adhikari, 3/26/2016
# Use DSSP to extract SS information from a PDB file

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;

my $dir_base = dirname(abs_path($0));
my $program_dssp   = "$dir_base/../third-party-programs/dssp-2.0.4-linux-amd64";
confess "DSSP program not found at $program_dssp!" if not -f $program_dssp;

my $file_pdb = shift;
confess "Usage: $0 <in>\n" if not $file_pdb;

confess "ERROR! Fasta file $file_pdb does not exist!" if not -f $file_pdb;
my %ss = dssp_result($file_pdb, "ss");
my $ssrow = "";
my $max = (sort {$b <=> $a} keys %ss)[0];
for(my $i = 1; $i <= $max; $i++){
	$ssrow .= "-" if not defined $ss{$i};
	$ssrow .= $ss{$i} if defined $ss{$i};
}
confess "looks like DSSP failed!" if length($ssrow) < 2;
print $ssrow."\n";

sub dssp_result{
	my $file_pdb = shift;
	my $selection = shift;
	confess "ERROR! $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! selection not defined!" if not defined $selection;
	my %RESIDUE = ();
	my %SS  = ();
	my %PHI = ();
	my %PSI = ();
	my $command = "$program_dssp $file_pdb | grep -C 0 -A 1000 \'  #  RESIDUE\' | tail -n +2";
	my @dssp_rows = `$command`;
	foreach(@dssp_rows){
		my $rnum = substr($_,  5, 6);
		my $res  = substr($_, 13, 1);
		my $sstr = substr($_, 16, 1);
		my $phia = substr($_,103, 6);
		my $psia = substr($_,109, 6);
		$rnum =~ s/\s+//g;
		$res  =~ s/\s+//g;
		$sstr =~ s/\s+//g;
		$phia =~ s/\s+//g;
		$psia =~ s/\s+//g;
		# alternate residue locations may have alphabets
		#$rnum =~ s/[^0-9]//g;
		$res  =~ s/[^A-Z]//g;
		$sstr =~ s/[^A-Z]//g;
		next if length($rnum) < 1;
		confess ":( residue not defined for $rnum" if length($res) < 1;
		confess ":( phi not defined for $rnum" if length($phia) < 1;
		confess ":( psi not defined for $rnum" if length($psia) < 1;
		$sstr = "C" if length($sstr) < 1;
		$sstr =~ s/\./C/g;
		$sstr =~ s/I/C/g;
		$sstr =~ s/S/C/g;
		$sstr =~ s/T/C/g;
		$sstr =~ s/B/C/g;
		$sstr =~ s/G/C/g;
		$RESIDUE{$rnum} = $res;
		$SS{$rnum}      = $sstr;
		$PHI{$rnum}     = $phia;
		$PSI{$rnum}     = $psia;
	}
	# special temporary change
	#return if not scalar keys %PHI;
	confess "$file_pdb seems empty!" if not scalar keys %PHI;
	return %SS if ($selection eq "ss");
	return %PHI if ($selection eq "phi");
	return %PSI if ($selection eq "psi");
	confess "ERROR! Invalid selection string!";
}
