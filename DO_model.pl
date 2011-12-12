
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
my $TEMPDIR = '../alignments/';
my $s = $dbh->prepare("SELECT ass.evalue, ass.region, ass.model, ass.sf, t1.description, align.alignment, comb_index.comb, family.evalue, family.px, family.fa, t2.description, genome_sequence.length
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
	$domain_details{$domain_no}{'region'} = $temp[1];
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
		my $unique = "$seqid"."_"."$genome"."_"."$id"."_"."$ch"."_"."$domain_details{$domain_number}{'model'}"."_"."$domain_details{$domain_number}{'region'}"."_"."$domain_details{$domain_number}{'px'}";
        #my $unique = int( rand(999999999999999) );
            open ALI, '>', "$TEMPDIR/align"."$unique".".temp" or die "Cannot open $TEMPDIR/align"."$unique".".temp".": $!\n";
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:/home/luca/rackham/astral/$folder/$file:   FIRST : $ch : LAST : $ch ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
my $altemp = uc($domain_details{$domain_number}{'alignment'});
my $tl = length($altemp);
print ALI "$altemp"."*\n";

close ALI;
 my @args = ("./model-single.py","align"."$unique".".temp","$id"."$ch","$unique",">&/dev/null" );

 system( '/usr/bin/python', @args );
 
}
