#!/usr/bin/perl -w
# Badri Adhikari, 4/2/2016

use strict;
use warnings;
use LWP::UserAgent;
 
my $ua = LWP::UserAgent->new;
 
my $server_endpoint = "http://cactus.rnet.missouri.edu/cgi-bin/dncon/freecontact.cgi";
# Can be tested from http://cactus.rnet.missouri.edu/dncon/freecontact.html
# Can be tested from http://cactus.rnet.missouri.edu/dncon/jobs-freecontact

my $id = "badri8";
my $aln = "";
$aln .= "MALSLCVLFTLASVVSGHVAHPSLGRGDGFPFLWDNAASTLDQLNGTDTTI\n";
$aln .= "---------------LAPNPESDTDENDSYPPFWDQINGDIAEFPVQNNKT\n";
$aln .= "---------TVIGCNFASFAVVS-SLSDAHPPLWKDCPDQLSDYKIENGKY\n";
$aln .= "---VLLTAVALSPTTVVALTQEDVCKNDVYPPLWHQAPGSIEDFPVHGNKI\n";
$aln .= "---VLLTAVALSPTTVVALTQEDVCKNDVYPPLWHQAPGSIEDFPVHGNKI\n";
$aln .= "MAFNIFWACVIIGCIFASIAKAS-NISDVYPPLWKESPGQFSDYKIENGKY\n";
$aln .= "----LIVTITLSPVTAAGQSKNDATREDVYPPLWDLAPESLLNFLVKDNKI\n";
$aln .= "---LLVVTVTLSPVTAEAQSENDVTGEDVYPPLWDLAPGNLLDFPVKDNKI\n";
$aln .= "----LIVTVTLSPVTASAQSEQDATREDAYPPLWDLAPENLMDFPVKDSKI\n";
$aln .= "----LVVPVMLSPVTATSWSEKDATTEDAYPPLWDLAPENLLDFLVKDNKI\n";
$aln .= "----LTVAITLSPVTAASWSERDATGEDIYPPLWNLAPENLSDFLVKDNKI\n";
$aln .= "----LIVTVTLSPVIATAQSEKDATGEDVYPPLWDLAPENLLDFLVKDNKS\n";
$aln .= "----LVVIVTFSPVPAATQSENDVTGEDLYPPLWDLAPGNLLDFPVKDNKI\n";
$aln .= "MAFLPPWTCVLVGCFSASLAKES-NFSDLYPPLWKKSASQFSDYRVENGKY\n";
$aln .= "MAFPSLWACLLIGCFSVSLAEAY-NLSDLYPPLWKESPGQFSDYRVKNGKY\n";

# set custom HTTP request header fields
use HTTP::Request::Common qw( POST );
my $req = POST($server_endpoint,
   Content_Type => 'application/x-www-form-urlencoded',
   Content => [ id => "$id", aln => "$aln" ],
);

my $resp = $ua->request($req);
if ($resp->is_success) {
	my $message = $resp->decoded_content;
	print "Received reply: $message\n";
}
else {
    print "HTTP POST error code: ", $resp->code, "\n";
    print "HTTP POST error message: ", $resp->message, "\n";
}

print "Wait 10 seconds for the job to finish in CACTUS server..\n";
sleep 10;

system("wget http://cactus.rnet.missouri.edu/dncon/jobs-freecontact/$id/$id.freecontact");
system("head $id.freecontact");
