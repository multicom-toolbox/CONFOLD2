#!/usr/bin/perl -w
# Badri Adhikari, 8/2/2016
# Convert PSIPRED predicted '.ss' file to 3-STATE secondary structure file in FASTA format.

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;

my $file_psipred = shift;
my $out_file_ss = shift;

confess "Usage: $0 <in> <out>\n" if not $file_psipred;
confess "Usage: $0 <in> <out>\n" if not -f $file_psipred;

my $ss  = "";
open INPUT, $file_psipred or confess $!;
while(<INPUT>){
	$_ =~ s/^\s+//;
	next if $_ !~ m/^[0-9]/;
	my @columns = split(/\s+/, $_);
	$ss .= $columns[2];
}
close INPUT;

system_cmd("rm -f $out_file_ss");
print2file($out_file_ss, ">".basename($file_psipred));
print2file($out_file_ss, $ss);

sub print2file{
	my $file = shift;
	my $message = shift;
	my $newline = shift;
	$newline = "\n" if not defined $newline;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
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
