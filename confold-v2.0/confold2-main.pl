#!/usr/bin/perl -w
# Badri Adhikari, 07-09-2017
# Main CONFOLD2 script for protein structure prediction

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;

my $DIR_BASE   = dirname(abs_path($0));
my $CONFOLD    = "$DIR_BASE/core.pl";
my $EVALUATE   = "$DIR_BASE/suppl-scripts/eval-using-tmscore.pl";
my $PSIPRED2SS = "$DIR_BASE/suppl-scripts/psipred2ss.pl";
my $MAXCLUSTER = "$DIR_BASE/third-party-programs/maxcluster64bit";

confess "Oops!! $CONFOLD not found!" if not -f $CONFOLD;

my ($file_fasta, $dir_out, $file_rr, $file_ss, $file_native, $mcount);
GetOptions(
	"out=s"     => \$dir_out,
	"rr=s"      => \$file_rr,
	"mcount=i"  => \$mcount,
	"ss=s"      => \$file_ss,
	# Supply native structure for evaluation after prediction (hidden parameter)
	"native=s"  => \$file_native)
or confess "ERROR! Error in command line arguments!";

print_usage() if not defined $file_rr;
print_usage() if not defined $file_ss;
print_usage() if not defined $dir_out;

if (not -f $file_rr){
	print "RR file $file_rr does not exist!\n";
	print_usage();
	exit(1);
}
if (not -f $file_ss){
	print "SS file $file_ss does not exist!\n";
	print_usage();
	exit(1);
}
if (not $dir_out){
	print 'Output directory not defined!\n';
	print_usage();
	exit(1);
}

$mcount       = 20 if not $mcount;
$file_rr      = abs_path($file_rr);
$dir_out      = abs_path($dir_out);
$file_ss      = abs_path($file_ss);

if (defined $file_native){
	$file_native = abs_path($file_native) if -f $file_native;
	confess "Oops!! $file_native doest not exist!" if not -f $file_native;
	confess "Oops!! $EVALUATE not found!" if not -f $EVALUATE;
}

my @xL = ();
for(my $i = 0.1; $i < 4.01; $i = $i + 0.1){
	my $x = sprintf "%.1fL", $i;
	push @xL, $x;
}

print "Started [$0]: ".(localtime)."\n";
my $id = basename($file_rr, ".rr");
system_cmd("mkdir -p $dir_out") if not -d $dir_out;
chdir $dir_out or confess $!;

if (-f "$dir_out/$id-model1.pdb"){
	print "Looks like everything is already done for this prediction.. quitting..\n";
	exit 0;
}

$| = 1;

my $flag_new_efunct = 1;
foreach my $x (@xL){
	if (-f "$dir_out/top-$x/stage2/${id}_model1.pdb" or -f "$dir_out/top-models/top-$x-model-1.pdb"){
		print "\nLooks like CONFOLD job with top $x contacts has already finished.. Not running anything..";
		system_cmd("rm -f $dir_out/top-$x.running");
		next;
	}
	print "\n";
	print "Forking CONFOLD job with top $x contacts\n";
	print "Log file at -> $dir_out/top-$x/\n";
	system_cmd("rm -rf $dir_out/top-$x");
	system_cmd("mkdir -p $dir_out/top-$x");
	system_cmd("touch $dir_out/top-$x.queued");
	chdir "$dir_out/top-$x" or confess $!;
	open SB, ">job-$id-top-$x.sh" or confess $!;
	print SB "#!/bin/bash -l\n";
	print SB "#SBATCH -J $id-$x\n";
	print SB "#SBATCH -o $id-$x.log\n";
	print SB "#SBATCH -p Bonus,Lewis\n";
	print SB "#SBATCH -n 1\n";
	print SB "#SBATCH --mem 2G\n";
	print SB "#SBATCH -t 2-00:00\n";
	print SB "mv $dir_out/top-$x.queued $dir_out/top-$x.running\n";
	if ($flag_new_efunct){
		print SB "$CONFOLD -asymptote 0 -rswitch 1.8 -mcount $mcount -rr $file_rr -ss $file_ss -selectrr $x -o ./\n";
		$flag_new_efunct = 0;
	}
	else{
	print SB "$CONFOLD -mcount $mcount -rr $file_rr -ss $file_ss -selectrr $x -o ./\n";
		$flag_new_efunct = 1;
	}
	print SB "rm $dir_out/top-$x.running\n";
	close SB;
	system "chmod +x job-$id-top-$x.sh";
	# Option 1 for parallelization (slurm managed HPC cluster)
	system "sbatch job-$id-top-$x.sh";
	# Option 2 for parallelization (multiple CPUs in one computer)
	# Remove the "&" for running serially
	#system "./job-$id-top-$x.sh > top-$x.log &";
	sleep 1;
}

print "\n";
print("Wait for all CONFOLD jobs to finish (check individual log files for progress)..\n");

my $running = 1;
while(int($running) > 0){
	sleep 10;
	$running = 0;
	foreach my $x (@xL){
		if (-f "$dir_out/top-$x.queued"){
			$running = 1;
			last;
		}
		if (-f "$dir_out/top-$x.running"){
			$running = 1;
			last;
		}
	}
}

print "\n";
print "Looks like all jobs have completed, checking to see if expected model files are present..\n";

my $c = 0;
system_cmd("mkdir -p $dir_out/top-models");
foreach my $x (@xL){
	my $flag_zip = 0;
	for (my $i = 1; $i <= 5; $i++){
		my $modela = "$dir_out/top-$x/stage2/${id}_model$i.pdb";
		my $modelb = "$dir_out/top-models/top-$x-model-$i.pdb";
		if (-f $modelb){
			print "Expected model $modelb already present! Do nothing!\n";
			$c++;
			next;
		}
		if (-f $modela){
			print "Copying $modela to $modelb\n";
			system_cmd("cp $modela $modelb");
			$c++;
			$flag_zip = 1;
			next;
		}
		warn "Oops!! Expected model $modela not found!";
	}
	next if $flag_zip == 0;
	chdir "$dir_out/top-$x" or confess $!;
	if (not -f "stage1n2.tar.gz"){
		next if not -d "stage1";
		next if not -d "stage2";
		next if not -d "input";
		print "Zipping intermediate files for top-$x ..\n";
		system_cmd("rm -f stage1n2.tar.gz");
		system_cmd("tar -czf stage1n2.tar.gz stage1 stage2 input");
		system_cmd("rm -r stage1");
		system_cmd("rm -r stage2");
		system_cmd("rm -r input");
	}
}

print "\n";
confess "Oops!! Only $c models found in mkdir -p $dir_out/top-models [expected = 200]" if $c < 100;
print "$c models found in mkdir -p $dir_out/top-models [expected = 200]\n";

print "\n";
my $models_to_remove = 50;
print "Removing bottom $models_to_remove models using contact satisfaction scores..\n";

my $seq = seq_rr($file_rr);
my $L   = length($seq);
my %rr = topxLongRange($file_rr, 0.2 * $L);
confess "Cannot calculate precision! Empty selected contacts" if not scalar keys %rr;
system_cmd("mkdir -p $dir_out/clustering");
chdir "$dir_out/clustering" or confess $!;
system_cmd("ls $dir_out/top-models/*.pdb > original.lst");
open LIST, "original.lst" or confess $!;
my @all_pdbs = <LIST>;
chomp @all_pdbs;
close LIST;

my %pdbscore = ();
foreach my $pdb (@all_pdbs){
	my %pdb_rr = native_contacts($pdb);
	my $score = 0.0;
	foreach my $row (keys %rr){
		next if not defined $pdb_rr{$row};
		my @C = split /\s+/, $row;
		$score = $score + $rr{$row};
	}
	$pdbscore{$pdb} = $score;
}

my @good_pdbs = ();
my $i = 0;
for(sort {$pdbscore{$b} <=> $pdbscore{$a}} keys %pdbscore){
	push @good_pdbs, $_;
	$i++;
	last if $i == (200 - $models_to_remove);
}

open LIST, ">updated.lst" or confess $!;
foreach (@good_pdbs){
	print LIST $_."\n";
}
close LIST;

print "\n";
print "Run clustering with the updated list and pick 5 centroids..\n";

my $clustersize = int((200 - $models_to_remove) / 5);
system_cmd("$MAXCLUSTER -l updated.lst -ms 5 -C 5 -is $clustersize > maxcluster.results");
system_cmd("grep -A 8 \"INFO  : Cluster  Centroid  Size        Spread\" maxcluster.results | grep $id > centroids.txt");

open CENT, "centroids.txt" or confess $!;
my @centroids = <CENT>;
chomp @centroids;
close CENT;

my @pdb_list = ();
for (my $i = 0; $i < 5; $i++){
	next if not defined $centroids[$i];
	my @C = split /\s+/, $centroids[$i];
	next if not defined $C[7];
	next if not -f $C[7];
	push @pdb_list, $C[7]; 
	print "Added ".$C[7]." to top 5 list\n";
}

print "\n";
print "Rank the ".($#pdb_list)." models selected [expected = 5] ..\n";
confess "Error! Top models could not be found!" if scalar(@pdb_list) < 1;

my %top5pdbscores = ();
foreach (@pdb_list){
	$top5pdbscores{$_} = $pdbscore{$_};
}

$i = 1;
for (sort {$top5pdbscores{$b} <=> $top5pdbscores{$a}} keys %top5pdbscores){
	print "coping ".$_." as model$i.pdb\n";
	system_cmd("cp ".$_." $dir_out/$id-model$i.pdb");
	$i++;
}

if (defined $file_native){
	print "\n";
	print "Native structure supplied.. Evaluating best-of-200 models ..\n";
	system_cmd("$EVALUATE $file_native $dir_out/top-models best &> $dir_out/full-evaluation.txt");

	print "Native structure supplied.. Evaluating best-of-5 models ..\n";
	system_cmd("$EVALUATE $file_native $dir_out best &> $dir_out/top5-evaluation.txt");

	print "Native structure supplied.. Evaluating top-5 models ..\n";
	system_cmd("$EVALUATE $file_native $dir_out all &> $dir_out/top1-evaluation.txt");
}

print "\nFinished [$0]: ".(localtime)."\n";

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "EXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command > $log 2>&1");
	}
	else{
		#print "[[Executing: $command]]\n\n";
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "ERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}

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

sub native_contacts{
	my $file_pdb = shift;
	my $separation = shift;
	my $d_threshold = shift;
	my $flg_atoms_not_dist = shift;
	confess ":( $file_pdb does not exist!" if not -f $file_pdb;
	my %contacts = ();
	my %xyz = xyz_pdb($file_pdb);
	foreach my $r1(sort {$a <=> $b} keys %xyz){
		foreach my $r2(sort {$a <=> $b} keys %xyz){
			next if abs($r1 - $r2) < 1;
			my $d = calc_dist($xyz{$r1}, $xyz{$r2});
			next if $d >= 8.0;
			$contacts{$r1." ".$r2} = 1;
		}
	}
	return %contacts;
}

sub xyz_pdb{
	my $chain = shift;
	confess "\nERROR! file $chain does not exist!" if not -f $chain;
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "\nERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	my %native_res_list = res_num_res_name($chain);
	my %selected_xyz = ();
	foreach (sort keys %xyz_pdb){
		my @C = split /\s+/, $_;
		my $atom_of_interest = "CB";
		$atom_of_interest = "CA" if $native_res_list{$C[0]} eq "GLY";
		next if $C[1] ne $atom_of_interest;
		$selected_xyz{$C[0]} = $xyz_pdb{$_};
	}
	confess "\nERROR! Empty xyz coordinates in the pdb file!" if not scalar keys %selected_xyz;
	return %selected_xyz;
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
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}

sub res_num_res_name{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}

sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}

sub topxLongRange{
	my $rr_file = shift;
	my $count = shift;
	confess "Input file not defined" if not defined $rr_file;
	confess "Input file $rr_file does not exist!" if not -f $rr_file;
	confess "No contact count!" if not defined $count;
	my %rr = ();
	my $i = 1;
	open RR, $rr_file or error_exit($!);
	while(<RR>){
		my $row = $_;
		chomp $row;
		$row =~ s/\r//g;
		$row =~ s/^\s+//;
		next unless $row =~ /^[0-9]/;
		my @C = split /\s+/, $row ;
		confess "Expecting a pair in row [".$row."]!" if (not defined $C[0] || not defined $C[1]);
		confess "Confidence column not defined in row $row" if not defined $C[4];
		my $d = abs($C[0]-$C[1]);
		next if $d < 24;
		$rr{$C[0]." ".$C[1]} = $C[4];
		$i++;
		last if $i > $count;
	}
	close RR;
	return %rr;
}

sub print_usage{
	my $param_info = <<EOF;
===============================================================================
CONFOLD version 2.0 (CONFOLD2)
===============================================================================

-------------------------------------------------------------------------------
PARAMETER   DESCRIPTION 
rr        : Predicted Contacts in CASP RR format
ss        : SCRATCH predicted secondary structure ('.ss' file in fasta format)
out       : Output directory
mcount    : Number of models for each CONFOLD job (default 20; change to 5 for faster results)

-------------------------------------------------------------------------------
Example Usage:
\$ ./confold2-main.pl -rr ./dry-run/input/1guu.rr -ss ./dry-run/input/1guu.ss -out ./output

-------------------------------------------------------------------------------
REFERENCES:
(A) CONFOLD v2.0:
    "Submitted"

(B) CONFOLD v1.0:
    "CONFOLD: Residue-Residue Contact-guided ab initio Protein Folding",
    Proteins: Structure, Function, and Bioinformatics, 2015.
    B. Adhikari, D. Bhattacharya, R. Cao, J. Cheng. 
-------------------------------------------------------------------------------
EOF
	print $param_info;
	print "\n";
	exit 1;
}
