--------------------------------------------------------------------------------
CONFOLD version 2.0 (CONFOLD2) 7/12/2016
--------------------------------------------------------------------------------
CONFOLD2 is a contact-guided protein structure prediction tool. Most scripts are written in Perl.
It takes predicted contacts file (CASP RR format) and 3-state secondary structure (SCRATCH fasta format) as input and 
delivers top 5 models.

--------------------------------------------------------------------------------
Installing CONFOLD
--------------------------------------------------------------------------------
1. Install CNS-suite (see below for instructions).
2. Update paths for the variables '$cns_suite' in 'core.pl'.
3. Configure how CONFOLD2 will be parallelized.
   Current (default) configuration is to use HPC cluster.
   Make appropriate changes at lines 113-117 of 'confold2-main.pl' to run in a local machine.
4. Test the program:
   $ ./confold2-main.pl -rr ./dry-run/input/1guu.rr -ss ./dry-run/input/1guu.ss -out ./output-1guu 

--------------------------------------------------------------------------------
Installing CNS Suite
--------------------------------------------------------------------------------
1. To download CNS suite, provide your academic profile related 
   information at http://cns-online.org/cns_request/. An email
   with (a) link to download, (b) login, and (c) password
   will be sent to you. Follow the link, possibly
   http://cns-online.org/download/, and download 
   CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
2. Unzip
   $ tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
3. Change directory to cns_solve
   $ cd cns_solve_1.3
4. Unhide the file '.cns_solve_env_sh'
   $ mv .cns_solve_env_sh cns_solve_env.sh
5. Edit (a) 'cns_solve_env.sh' and (b) 'cns_solve_env' to replace
   '_CNSsolve_location_' with CNS installation directory.
   For instance, if your CNS installation path is
   '/home/user/programs/cns_solve_1.3' replace
   '_CNSsolve_location_' with this path
6. Test CNS installation
   $ source cns_solve_env.sh
   $ cd test 
   $ ../bin/run_tests -tidy *.inp

--------------------------------------------------------------------------------
Predictions for PSICOV, CASP11, and CASP12 datasets
--------------------------------------------------------------------------------
The contact predictions, secondary structure predictions, and CONFOLD2 reconstructed models for (a) PSICOV dataset using PSICOV contacts, (b) PSICOV dataset using MetaPSICOV contacts, (c) CASP11 dataset using CONSIP2 contacts, and (d) CASP12 dataset using RaptorX contacts are available at http://sysbio.rnet.missouri.edu/bdm_download/confold2/.

--------------------------------------------------------------------------------
Contact
--------------------------------------------------------------------------------
bap54@mail.missouri.edu (developer)
chengji@missouri.edu (PI)

