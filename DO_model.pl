
#!/usr/bin/perl -wT

use strict;
use warnings;

use CGI qw/:standard/;
use DBI;
use lib '/home/rackham/modules/';
use rackham;

my ( $dbh, $sth );
$dbh = rackham::DBConnect('superfamily');

my $seqid = $ARGV[0];
my $genome = $ARGV[1];
my $marginal = 0.01;
my $TEMPDIR = '/home/luca/rackham/astral/alignments/';
my $s = $dbh->prepare("SELECT ass.evalue, ass.region, ass.model, ass.sf, t1.description, align.alignment, comb_index.comb, family.evalue, family.px, family.fa, t2.description, genome_sequence.length,protein.protein
                      FROM align, ass, des AS t1, des AS t2, comb, family, protein, genome_sequence, comb_index
                      WHERE ass.auto = align.auto AND
                      protein.protein = genome_sequence.protein AND
                            comb.comb_id=comb_index.id AND
                            protein.protein = comb.protein AND
                            protein.protein = ass.protein AND
                            ass.auto = family.auto AND
                            protein.seqid = '$seqid' AND
                            ass.sf = t1.id AND
                            family.fa = t2.id AND
                            protein.genome = '$genome' AND
                            ass.evalue <= $marginal
                      ORDER BY ass.evalue;");
$s->execute;

                   
my $domain_no = 1;
my %domain_details;
while ( my @temp = $s->fetchrow_array ) {
	$domain_details{$domain_no}{'evalue'} = $temp[0];
	my $region = $temp[1];
	my @regions = split(/-/,$region);
	$region =~ s/-/_/g;
	$domain_details{$domain_no}{'start'} = $regions[0];
	$domain_details{$domain_no}{'stop'} = $regions[1];
	$domain_details{$domain_no}{'region'} = $region;
	$domain_details{$domain_no}{'model'} = $temp[2];
	$domain_details{$domain_no}{'sf'} = $temp[3];
	$domain_details{$domain_no}{'description'} = $temp[4];
	$domain_details{$domain_no}{'alignment'} = $temp[5];
	$domain_details{$domain_no}{'comb_index'} = $temp[6];
	$domain_details{$domain_no}{'family_evalue'} = $temp[7];
	$domain_details{$domain_no}{'px'} = $temp[8];
	$domain_details{$domain_no}{'fa'} = $temp[9];
	$domain_details{$domain_no}{'description'} = $temp[10];
	$domain_details{$domain_no}{'length'} = $temp[11]; 
	$domain_details{$domain_no}{'protein'} = $temp[12]; 
    my $sth2=$dbh->prepare( "select name from des where id =$temp[8];" );
	$sth2->execute;
	while ( my @temp2 = $sth2->fetchrow_array ) {
		$domain_details{$domain_no}{'closest_structure'} = $temp2[0];
	}
	$domain_no++;
}

foreach my $domain_number (keys %domain_details){

my $id = substr $domain_details{$domain_number}{'closest_structure'},1 ,4 ;
my $ch = uc(substr $domain_details{$domain_number}{'closest_structure'},5 ,1 );
my $file = "$domain_details{$domain_number}{'closest_structure'}.ent";
my $folder = substr $file, 2,2;
		my $unique ="$domain_details{$domain_number}{'protein'}"."_"."$domain_details{$domain_number}{'region'}";
        #my $unique = int( rand(999999999999999) );
 # unless (-e "$TEMPDIR/align"."$unique".".temp") {      
            open ALI, '>', "$TEMPDIR/$genome/align"."$unique".".temp" or die "Cannot open $TEMPDIR/align"."$unique".".temp".": $!\n";
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:/home/luca/rackham/astral/$folder/$file:   FIRST : $ch : LAST : $ch ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
my $altemp = uc($domain_details{$domain_number}{'alignment'});
my $tl = length($altemp);
print ALI "$altemp"."*\n";

close ALI;
 my @args = ("./model-single.py","align"."$unique".".temp","$id"."$ch","$unique","$genome",">&/dev/null" );

 system( '/usr/bin/python', @args );
 my $cut = 0;
 my @aligns;
 while (length($domain_details{$domain_number}{'alignment'}) > 0){
 	my $inc = 60;
 	if(length($domain_details{$domain_number}{'alignment'}) < 60){
 		$inc = length($domain_details{$domain_number}{'alignment'});
 	}
 	push @aligns,substr($domain_details{$domain_number}{'alignment'},0,$inc,'');
 }
 
 open(FILE,"$TEMPDIR"."/"."$genome"."/"."$unique".".B99990001.pdb");
 open(FILEOUT,">$TEMPDIR"."/"."$genome"."/"."$unique".".pdb");
 my $flag = 0;
 while(<FILE>){
 	if($flag == 1){
 		#REMARK GENOME3D TOTAL TEMPLATES 1
 		print FILEOUT "REMARK GENOME3D TOTAL TEMPLATES 1\n";
 		#REMARK GENOME3D SELECTION X
		#REMARK GENOME3D MODELLING X
		print FILEOUT "REMARK GENOME3D SELECTION SUPERFAMILY\n";
		print FILEOUT "REMARK GENOME3D MODELLING MODELLER\n";
		#REMARK GENOME3D START X
		#REMARK GENOME3D STOP X
		print FILEOUT "REMARK GENOME3D START $domain_details{$domain_number}{'start'}\n";
		print FILEOUT "REMARK GENOME3D STOP $domain_details{$domain_number}{'stop'}\n";
		print FILEOUT "REMARK GENOME3D >$id"."$ch".":$domain_details{$domain_number}{'region'}\n";
		foreach my $al (@aligns){
		print FILEOUT "REMARK GENOME3D $al\n";		
 		}
 		print FILEOUT $_;
 		$flag++;
 	}else{
 		print FILEOUT $_;
 		$flag++;
 	}
 }
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".B99990001.pdb");
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".ini");
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".rsr");
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".D00000001");
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".V99990001");
 unlink("$TEMPDIR"."/"."$genome"."/"."$unique".".sch");
 
#}
}
