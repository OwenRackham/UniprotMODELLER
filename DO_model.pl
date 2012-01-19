
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
my $fo = $ARGV[2];
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
	my @aligned = split(/,/,$temp[1]);
	my @start = split(/-/,$aligned[0]);
	my $start = $start[0];
	my @finish = split(/-/,$aligned[$#aligned]);
	my $finish = $finish[1];
	my $region = "$start"."_"."$finish";
	$domain_details{$domain_no}{'start'} = $start;
	$domain_details{$domain_no}{'stop'} = $finish;
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
if(-e "$TEMPDIR"."$seqid"."_"."$domain_details{$domain_number}{'start'}"."_"."$domain_details{$domain_number}{'stop'}".".pdb"){
	my $cmd = "cp "."$TEMPDIR"."$seqid"."_"."$domain_details{$domain_number}{'start'}"."_"."$domain_details{$domain_number}{'stop'}".".pdb"." .";
	system($cmd);
}else{

my $id = substr $domain_details{$domain_number}{'closest_structure'},1 ,4 ;
my $ch = uc(substr $domain_details{$domain_number}{'closest_structure'},5 ,1 );
my $file = "$domain_details{$domain_number}{'closest_structure'}.ent";
my $folder = substr $file, 2,2;
		my $unique ="$seqid"."_"."$domain_details{$domain_number}{'region'}";
        #my $unique = int( rand(999999999999999) );
 # unless (-e "$TEMPDIR/align"."$unique".".temp") {      
            open ALI, '>', "$TEMPDIR"."$fo/align"."$unique".".temp" or die "Cannot open $TEMPDIR/$fo/align"."$unique".".temp".": $!\n";
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:/home/luca/rackham/astral/$folder/$file:   FIRST : $ch : LAST : $ch ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
my $altemp = uc($domain_details{$domain_number}{'alignment'});
my $tl = length($altemp);
print ALI "$altemp"."*\n";

close ALI;
 my @args = ("./model-single.py","align"."$unique".".temp","$id"."$ch","$unique","$fo",">&/dev/null" );

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
 
 open(FILE,"$TEMPDIR"."/"."$fo"."/"."$unique".".B99990001.pdb");
 open(FILEOUT,">$TEMPDIR"."/"."$fo"."/"."$unique".".pdb");
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
 
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".B99990001.pdb");
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".ini");
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".rsr");
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".D00000001");
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".V99990001");
 unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".sch");
 unlink("$TEMPDIR"."/"."$fo"."/"."align"."$unique".".temp");
 
#}
}
}
