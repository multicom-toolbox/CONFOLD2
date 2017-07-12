#!/usr/bin/perl
# 3/27/2016, Badri Adhikari
# Evaluate the precision of predicted long-range contacts against a PDB file

use strict;
use warnings;
use Carp;
use File::Basename;
use File::Temp qw(tempfile);

my $native_pdb = shift;
my $rr_file    = shift;
my $x_in_xL    = shift;
my $xxL        = "";

if (not $rr_file){
	print STDERR "Usage: $0 <native> <RR> <xL or #ofContacts>\n";
	print STDERR "Example (evaluate top-L/10 contacts): $0 ./native.pdb ./contact.rr 0.1L\n";
	print STDERR "Example (evaluate top-5    contacts): $0 ./native.pdb ./contact.rr 5\n";
	exit 1;
}

if (not -f $rr_file){
	print STDERR "ERROR! Input RR file does not exist!\n";
	print STDERR "RR: $rr_file\n" if $rr_file;
	print STDERR "Usage: $0 <native> <RR> <xL or #ofContacts>\n";
	print STDERR "Example (evaluate top-L/10 contacts): $0 ./native.pdb ./contact.rr 0.1L\n";
	print STDERR "Example (evaluate top-5    contacts): $0 ./native.pdb ./contact.rr 5\n";
	exit 1;
}

if (not $native_pdb){
	print STDERR "Usage: $0 <native> <RR> <xL or #ofContacts>\n";
	print STDERR "Example (evaluate top-L/10 contacts): $0 ./native.pdb ./contact.rr 0.1L\n";
	print STDERR "Example (evaluate top-5    contacts): $0 ./native.pdb ./contact.rr 5\n";
	exit 1;
}

if (not -f $native_pdb){
	print STDERR "ERROR! Native PDB file $native_pdb does not exist!\n";
	print STDERR "Usage: $0 <native> <RR> <xL or #ofContacts>\n";
	print STDERR "Example (evaluate top-L/10 contacts): $0 ./native.pdb ./contact.rr 0.1L\n";
	print STDERR "Example (evaluate top-5    contacts): $0 ./native.pdb ./contact.rr 5\n";
	exit 1;
}

my $seq = seq_rr($rr_file);

my ($rr_sorted_fh, $rr_sorted_file) = tempfile();
my ($rr_sorted_fh_temp, $rr_sorted_file_temp) = tempfile();

system_cmd("rm -f $rr_sorted_file");
system_cmd("cp $rr_file $rr_sorted_file_temp");
# Keep contact rows only
system_cmd("sed -i '/^[A-Z]/d' $rr_sorted_file_temp");
system_cmd("sed -i '/^-]/d' $rr_sorted_file_temp");
# Some files have leading white spaces
system_cmd("sed -i 's/^ *//' $rr_sorted_file_temp");
# Stable sort with -s option, i.e. maintain order in case confidence are equal
# Also using -g option instead of -n because some predictors have exponential values in confidence
system_cmd("echo \"$seq\" > $rr_sorted_file");
system_cmd("sort -gr -s -k5 $rr_sorted_file_temp >> $rr_sorted_file");

my $min_seq_sep = 6;
my $short_range_min  = 6;
my $short_range_max  = 11;
my $medum_range_min  = 12;
my $medium_range_max = 23;
my $long_range_min   = 24;
my $long_range_max   = 10000;

# http://proteopedia.org/wiki/index.php/Standard_Residues and http://prody.csb.pitt.edu/manual/reference/atomic/flags.html
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V UNK -);
our %AA1TO3 = reverse %AA3TO1;

my $d_threshold = 8;
my $param_atom_type = "cb";
my $rr_sequence = seq_rr($rr_sorted_file);
my $rr_pdb = seq_chain($native_pdb);

my $seq_pdb_chk = seqChainWithGaps($native_pdb);
if ($rr_sequence ne $seq_pdb_chk){
	warn "Warning! PDB and RR do not have same sequence!\n";
	warn "$rr_sequence [RR]\n";
	warn "$seq_pdb_chk [PDB]\n";
}

if (not defined $x_in_xL){
	$x_in_xL = int(length($rr_pdb)/5);
	$xxL = "top-L/5";
}
elsif ($x_in_xL =~ m/L/){
	$x_in_xL =~ s/L//g;
	$xxL = "top-${x_in_xL}L";
	$x_in_xL = int(sprintf("%.2f",length($rr_pdb)) * $x_in_xL + 0.5);
}

my %min_true_dist = all_pairs_min_dist($native_pdb, $min_seq_sep, $d_threshold, 0);
my %rr = rr_rows_ordered_in_hash($rr_sorted_file, 10000, "all", "long");

confess "Cannot calculate precision! Empty selected contacts" if not scalar keys %rr;
my %distances = rrhash2dist(\%rr, $x_in_xL);
confess "Distances could not be calculated for selected contacts! Something went wrong!" if not scalar keys %distances;
my $satisfied = 0;
foreach (sort {$distances{$a} <=> $distances{$b}}keys %distances){
#	print $_."\n";
	my @R = split /\s+/, $_;
	$satisfied++ if $distances{$_} <= $R[3];
#	print "->".$_."\n" if $distances{$_} <= $R[3];
}
printf "Precision-for-$xxL=%.2f\n", 100 * ($satisfied/(scalar keys %distances));

sub rrhash2dist{
	# If input is empty output will be empty
	my $rrhash = shift;
	my %rr_hash = %{$rrhash};
	my $xL = shift;
	my %output = ();
	my $count = 0;
	foreach (sort {$a <=> $b} keys %rr_hash){
		my @C = split /\s+/, $rr_hash{$_};
		next if not defined $min_true_dist{$C[0]." ".$C[1]};
		$output{$rr_hash{$_}} = $min_true_dist{$C[0]." ".$C[1]};
		$count++;
		last if $count == $xL;
	}
	return %output;
}

sub all_pairs_min_dist{
	my $file_pdb = shift;
	my $separation = shift;
	my $d_threshold = shift;
	my $flg_atoms_not_dist = shift;
	confess ":( $file_pdb does not exist!" if not -f $file_pdb;
	my %xyz = ();
	my %pairs_dist = ();
	my %pairs_atoms = ();
	if($param_atom_type eq "ca" or $param_atom_type eq "cb"){
		%xyz = xyz_pdb($file_pdb, $param_atom_type);
		foreach my $r1(sort {$a <=> $b} keys %xyz){
			foreach my $r2(sort {$a <=> $b} keys %xyz){
				next if $r1 >= $r2;
				next if abs($r1 - $r2) < $separation;
				my $d = calc_dist($xyz{$r1}, $xyz{$r2});
				$pairs_dist{$r1." ".$r2} = $d;
				if ($param_atom_type eq "cb"){
					$pairs_atoms{$r1." ".$r2} = "".return_cb_or_ca_atom($r1)." ".return_cb_or_ca_atom($r2);
				}
				else{
					$pairs_atoms{$r1." ".$r2} = "CA CA";
				}
			}
		}
		if ($flg_atoms_not_dist){
			return %pairs_atoms;
		}
		#append2log("Returning ".(scalar keys %pairs_dist)." rows for $file_pdb at separation $separation");
		return %pairs_dist;
	}
	else{
		%xyz = xyz_pdb($file_pdb, "ALL");
		foreach my $row1(keys %xyz){
			my @row1 = split /\s+/, $row1;
			my $res1 = $row1[0];
			my $atm1 = $row1[1];
			if ($param_atom_type eq "heavyatoms"){
				next if not ($atm1 eq "N" or $atm1 eq "CA" or $atm1 eq "C" or $atm1 eq "O");  
			}
			foreach my $row2(keys %xyz){
				my @row2 = split /\s+/, $row2;
				my $res2 = $row2[0];
				my $atm2 = $row2[1];
				if ($param_atom_type eq "heavyatoms"){
					next if not ($atm2 eq "N" or $atm2 eq "CA" or $atm2 eq "C" or $atm2 eq "O");  
				}
				next if $res1 >= $res2;
				next if abs($res1 - $res2) < $separation;
				my $d = calc_dist($xyz{$row1}, $xyz{$row2});
				if (not defined $pairs_dist{$res1." ".$res2}){
					$pairs_dist{$res1." ".$res2} = $d;
					$pairs_atoms{$res1." ".$res2} = "$atm1 $atm2" if $flg_atoms_not_dist;
				}
				if ($pairs_dist{$res1." ".$res2} > $d){
					$pairs_dist{$res1." ".$res2} = $d;
					$pairs_atoms{$res1." ".$res2} = "$atm1 $atm2" if $flg_atoms_not_dist;
				}
			}
		}
		if($flg_atoms_not_dist){
			return %pairs_atoms;
		}
		append2log("Returning ".(scalar keys %pairs_dist)." rows for $file_pdb at separation $separation");
		return %pairs_dist;
	}
}

sub return_cb_or_ca_atom{
	my $rnum = shift;
	my $residue = substr $rr_sequence, $rnum - 1, 1;
	confess "rnum not defined!" if not defined $rnum;
	confess "Could not find residue name for $rnum!" if not $residue;
	return "CA" if $residue eq "G";
	return "CB";
}

sub xyz_pdb{
	my $chain = shift;
	my $atom_selection = shift; # ca or cb or all/any
	$atom_selection = "all" if $atom_selection eq "any";
	confess "\nERROR! file $chain does not exist!" if not -f $chain;
	confess "\nERROR! Selection must be ca or cb or all or heavyatoms" if not (uc($atom_selection) eq "CA" or uc($atom_selection) eq "ALL" or uc($atom_selection) eq "CB" or uc($atom_selection) eq "HEAVYATOMS");
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "\nERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	if (uc($atom_selection) eq "ALL"){
		return %xyz_pdb;
	}
	elsif (uc($atom_selection) eq "HEAVYATOMS"){
		foreach (sort keys %xyz_pdb){
			my @C = split /\s+/, $_;
			if (not($C[1] eq "N" or $C[1] eq "CA" or $C[1] eq "C" or $C[1] eq "O")){
				delete $xyz_pdb{$_};
			}
		}
		return %xyz_pdb;
	}
	my %native_res_list = res_num_res_name($chain);
	my %selected_xyz = ();
	foreach (sort keys %xyz_pdb){
		my @C = split /\s+/, $_;
		my $atom_of_interest = uc($atom_selection);
		$atom_of_interest = "CA" if $native_res_list{$C[0]} eq "GLY";
		# Some pdb files have errors. Non-gly residues do not have CB atoms.
		# For instance, http://sysbio.rnet.missouri.edu/conassess/preloaded_data/fragfold/native/1ej8A_reindexed.pdb
		# Need to throw errors in such cases.
		confess("The pdb file $native_pdb does not have CB atom for residue ".$C[0]."! Try assessing CA-contacts instead!") if not defined $xyz_pdb{$C[0]." ".$atom_of_interest};
		next if $C[1] ne $atom_of_interest;
		$selected_xyz{$C[0]} = $xyz_pdb{$_};
	}
	confess "\nERROR! Empty xyz coordinates in the pdb file!" if not scalar keys %selected_xyz;
	return %selected_xyz;
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


sub rr_rows_ordered_in_hash{
	my $rr_file = shift;
	my $count = shift;
	my $source = shift;
	my $type = shift;
	my $confidence_thres = shift;
	my %native_residues = xyz_pdb($native_pdb, "cb");
	confess "Input file not defined" if not defined $rr_file;
	confess "Input file $rr_file does not exist!" if not -f $rr_file;
	confess "No File RR!" if not -f $rr_file;
	confess "No contact count!" if not defined $count;
	confess "No contact source!" if not defined $source;
	confess "No contact type" if not defined $type;
	confess "Invalid type" if not ($type eq "everything" or $type eq "all" or $type eq "short" or $type eq "medium" or $type eq "long");
	my %rr = ();
	my $i = 1;
	my %i_for_each_source = ();
	# Find all Contact sources
	if ($source eq "all"){
		open RR, $rr_file or confess $!;
		while (<RR>){
			next unless $_ =~ /[0-9]/;
			chomp $_;
			$_ =~ s/\r//g;
			$_ =~ s/^\s+//;
			next unless $_ =~ /^[0-9]/;
			my @C = split /\s+/, $_ ;
			last if not defined $C[5];
			$i_for_each_source{$C[5]} = 1;
		}
		close RR;
	}
	open RR, $rr_file or error_exit($!);
	while(<RR>){
		my $row = $_;
		next unless $row =~ /[0-9]/;
		chomp $row;
		$row =~ s/\r//g;
		$row =~ s/^\s+//;
		next unless $row =~ /^[0-9]/;
		my @C = split /\s+/, $row ;
		error_exit("Expecting a pair in row [".$row."]!\n") if (not defined $C[0] || not defined $C[1]);
		error_exit("Confidence column not defined in row [".$row."] in file <b>$rr_file</b>! </br>Please make sure that the input RR file is in 5-column format!") if not defined $C[4];
		# Fix order 
		if ($C[0] > $C[1]){
			$row = $C[1]." ".$C[0]." ".$C[2]." ".$C[3]." ".$C[4];
			$row = $C[1]." ".$C[0]." ".$C[2]." ".$C[3]." ".$C[4]." ".$C[5] if defined $C[5];
		}
		if (defined $confidence_thres){
			next if $C[4] < $confidence_thres;
		}
		# skip all those that are not in native
		next if not defined $native_residues{$C[0]};
		next if not defined $native_residues{$C[1]};
		# Select only LR, MR, SR or all
		my $d = abs($C[0]-$C[1]);
		if ($type eq "long"){
			next if $d < $long_range_min;
		}
		if ($type eq "medium"){
			next if $d < $medum_range_min;
			next if $d > $medium_range_max;
		}
		if ($type eq "short"){
			next if $d < $short_range_min;
			next if $d > $short_range_max;
		}
		next if ($d < $short_range_min and $type eq "all");
		# Select only those matching the source
		$C[5] = "all" if not defined $C[5];
		if ($source eq "all" and $C[5] eq "all"){
			$rr{$i} = $row;
			$i++;
			last if $i > $count;
		}
		elsif($source ne "all" and $C[5] eq "all"){
			confess "Specified some specific source but there is a whose source is not specified";
		}
		elsif($source eq "all" and $C[5] ne "all"){
			next if not defined $i_for_each_source{$C[5]};
			next if $i_for_each_source{$C[5]} > $count;
			$i++;
			last if $i > $count;
			$rr{$i} = $row;
			$i_for_each_source{$C[5]}++;
		}
		else{
			next if $C[5] ne $source;
			$rr{$i} = $row;
			$i++;
			last if $i > $count;
		}
	}
	close RR;
	return %rr;
}

sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}

################################################################################
# for extracting contacts from native structures
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

sub system_cmd{
	my $command = shift;
	confess "EXECUTE [$command]" if (length($command) < 5  and $command =~ m/^rm/);
	system($command);
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "Failed executing [$command]!<br>ERROR: $!";
	}
}