#!/usr/bin/perl -w
# Badri Adhikari, 8/2/2016
# Evaluate predicted PDB files against native structure

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;

my $dir_base = dirname(abs_path($0));

my $native  = shift;
my $dir_in  = shift;
my $option  = shift;
my $header  = shift;

if (not $dir_in){
	print STDERR "Usage: $0 <native> <dir-pdb> <top/best/all/avg> [header]\n";
	exit 1;
}
if (not (-f $dir_in or -d $dir_in)){
	print STDERR "ERROR! Input directory $dir_in does not exist!\n";
	print STDERR "Usage: $0 <native> <dir-pdb> <top/best/all/avg> [header]\n";
	exit 1;
}
if (not $native){
	print STDERR "ERROR! Native file not defined!\n";
	print STDERR "Usage: $0 <native> <dir-pdb> <top/best/all/avg> [header]\n";
	exit 1;
}
if (not -f $native){
	print STDERR "ERROR! Native file $native does not exist!\n";
	print STDERR "Usage: $0 <native> <dir-pdb> <top/best/all/avg> [header]\n";
	exit 1;
}

if (not $option){
	print STDERR "ERROR! Option (top/best/all) not supplied!\n";
	print STDERR "Usage: $0 <native> <dir-pdb> <top/best/all/avg> [header]\n";
	exit 1;
}
$header = "header" if not defined $header;

my $tmscore = "$dir_base/../third-party-programs/TMscore";
confess "TM-score program not found at $tmscore!" if not -f $tmscore;

my @pdb_list = ();
if (-f $dir_in){
	push @pdb_list, $dir_in;
}
else{
	@pdb_list = load_pdb($dir_in);
}
my %tmscores = ();
my %rmsds = ();
my %gdtts = ();
foreach my $p (@pdb_list){
	next if $p =~ /sub_embed/;
	next if ($p !~ /model1.pdb/ and $option eq "top");
	warn "Reading $p\n";
	my %tm = evaluate_using_TMscore($p, $native);
	$tmscores{$p} = $tm{"tm-score"};
	$rmsds{$p} = $tm{"rmsd"};
	$gdtts{$p} = $tm{"gdt-ts"};
}
confess "ERROR! Empty folder $dir_in!" if not scalar keys %tmscores;

if($option eq "avg"){
	print "TM\tRMSD\tMODEL\n" if $header eq "header";
	my $tm_sum = 0.0;
	my $rmsd_sum = 0.0;
	my $count = 0;
	for(keys %tmscores){
		$count++;
		$tm_sum += $tmscores{$_};
		$rmsd_sum += $rmsds{$_}
	}
	printf "%.4f\t%.2f\t$dir_in (models = $count)\n", $tm_sum/$count, $rmsd_sum/$count;
	exit 0;
}

if($option eq "all"){
	print "TM\tRMSD\tGDTTS\tMODEL\n" if $header eq "header";
	for(sort {$tmscores{$a} <=> $tmscores{$b}} keys %tmscores){
		printf "%.4f\t%.2f\t%.2f\t$_\n", $tmscores{$_}, $rmsds{$_}, 100 * $gdtts{$_};
	}
}

if($option eq "top"){
	print "TM\tRMSD\tGDTTS\tCOUNT\tMODEL\n" if $header eq "header";
	for(sort {$tmscores{$b} <=> $tmscores{$a}} keys %tmscores){
		printf "%.4f\t%.2f\t%.2f\t".(scalar keys %tmscores)."\t$_\n", $tmscores{$_}, $rmsds{$_}, 100 * $gdtts{$_};
	}
}

if($option eq "best"){
	print "TM\tRMSD\tGDTTS\tCOUNT\tMODEL\n" if $header eq "header";
	for(sort {$tmscores{$b} <=> $tmscores{$a}} keys %tmscores){
		printf "%.4f\t%.2f\t%.2f\t".(scalar keys %tmscores)."\t$_\n", $tmscores{$_}, $rmsds{$_}, 100 * $gdtts{$_};
		last;
	}
}

sub load_pdb{
	my $dir_chains = shift;
	confess "\n:( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdb_list = <$dir_chains/*.pdb>;
	if(not (@pdb_list)){
		@pdb_list = `find $dir_chains -name '*.ent'`;
	}
	chomp @pdb_list;
	confess "\nERROR! Directory $dir_chains has no pdb files!\n" unless(@pdb_list);
	return @pdb_list;
}

sub evaluate_using_TMscore{
	my $predicted = shift;
	my $native = shift;
	confess "\nERROR! Predicted pdb $predicted does not exit!" if not -f $predicted;
	confess "\nERROR! Native pdb $native does not exist!" if not -f $native;
	my @results = `$tmscore $predicted $native | grep -e RMSD\\ of -e TM-score\\ \\ \\  -e MaxSub-score -e GDT-TS-score -e GDT-HA-score`;
	if (not defined $results[0]){
		print "Executing: [$tmscore $predicted $native]\n";
		system("$tmscore $predicted $native");
		confess "\nError! TM-score failed!";
	}
	$results[0] =~ s/RMSD of  the common residues//g;
	$results[1] =~ s/TM-score//g;
	$results[2] =~ s/MaxSub-score//g;
	$results[3] =~ s/GDT-TS-score//g;
	$results[4] =~ s/GDT-HA-score//g;
	for(my $i = 0; $i <= 4; $i++ ){
        $results[$i] =~ s/\s+//g;
        $results[$i] =~ s/=//;
	}
	my %result = ();
	$result{"rmsd"}		= substr $results[0], 0, 5;
	$result{"tm-score"}	= substr $results[1], 0, 6;
	$result{"max-sub"}	= substr $results[2], 0, 6;
	$result{"gdt-ts"}	= substr $results[3], 0, 6;
	$result{"gdt-ha"}	= substr $results[4], 0, 6;
	return %result;
}

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "\nEXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "\nERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}
