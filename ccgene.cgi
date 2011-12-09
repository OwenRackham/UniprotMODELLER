#!/usr/bin/perl -wT

use strict;
use warnings;

use CGI qw/:standard/;
use DBI;

use lib 'modules';
use Jgough;
use JgoughSQL;
use UtilityBelt;

$| = 1;

#This program is supposed to: show the domains in a single gene
#Arguements:  genome seqid
#Julian Gough 27.4.01 modified 14.10.01

#ASSIGN-VARIABLES---------------------------------------
my $DEBUG;

my $SPIRITEMP = '../temp/spiricoil';
my $TEMPDIR   = $SPIRITEMP;
#my $TEMPDIR   = '../temp/genes';
my $threshold = 0.0001;
my $marginal  = 0.1;
my $arch      = '';
my $passlink  = '';
my @domain_attributes = qw( sf_eval region model sf superfamily alignment arch fa_eval px fa family length domain_no  );
my @domains_AofHRs;
my $genome;
my $seqid;
my %coils;
my $seqid_file;
my $genpass;
my @href;
my $residue_count = 1;
my $chn_store;
my $exclude_link;
my @include_links = ( "#topofpage\tTop of page",
                      "#section_domain_architecture\tDomain architecture",
                      "#section_domain_assignments\tDomain assignment details",
                      "#section_protein_sequence\tProtein sequence",
                    );
my %model_lookup;
#---------------------------------------------------------------------
#READ-IN-CGI-VARIABLES----------------------------------
my $query = new CGI;
my $sunid = $query->param('sunid') || '';
my $real = $query->param('real');
unless ( defined $real ) {
    $real = 0;
}
$genome = $query->param('genome') || '';
if ( $genome =~ /^(\w{2})$/ ) {
    $genome  = $1;
    $TEMPDIR = "$TEMPDIR/$genome";
    mkdir $TEMPDIR if ! -e $TEMPDIR;
}
else {
    print CGI::redirect('../nogenome.html');
    exit;
}

$seqid = $query->param('seqid');
unless ( defined $seqid ) {
    $seqid = '';
}

if ( $seqid =~ /(\n|')/ ) {
    my $title = 'Illegal characters in sequence identifier';
    my $explanation = 'Sequence identifier contains invalid characters: new line or apostrophe';
    UtilityBelt::BombOut( $title, $explanation );
}
elsif ( $seqid =~ /^(\S+)$/ ) {
    $seqid = $1;
}

$genpass = $query->param('password');
if ( defined $genpass ) {
    $passlink = "&password=$genpass";
}
else {
    $genpass = '';
}

#-------------------------------------------------------
# connect to SQL
my $dbh = JgoughSQL::DBConnect;
my $sth;
my @temp;

if ( GenomePassword($genome) == 0 ) {
    $dbh->disconnect;
    print CGI::redirect('../nogenome.html');
    exit;
}

#-------------------------------------------------------

# Untaint the path
$ENV{'PATH'} = '';
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

#----------------
my @models;
# LOAD IN DATA ------------------------------------------
#open the data.txt file
my %locations;
open DATA, "../spiricoil/regions.txt" or die "couldn't open data file";
while (<DATA>) {
	my $line = $_;

	chomp($line);
	if ($line =~ /(\S+)\s+(\S+)/){

		$locations{$1} = $2;
	}
}
#create a hash of arrays of stops and starts for the models


# RETRIEVE FROM SQL -------------------------------------
if ( defined $genome ) {
    $sth =
      $dbh->prepare( "SELECT ass.evalue, ass.region, ass.model, ass.sf, t1.description, align.alignment, comb_index.comb, family.evalue, family.px, family.fa, t2.description, genome_sequence.length
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
                      ORDER BY ass.evalue;" );
}
else {
    $sth =
      $dbh->prepare( "SELECT ass.evalue, ass.region, ass.model, ass.sf, t1.description, align.alignment, comb_index.comb, family.evalue, family.px, family.fa, t2.description, genome_sequence.length
                      FROM align, ass, des AS t1, des AS t2, comb, family, protein, genome_sequence, comb_index
                      WHERE protein.protein = ass.protein AND
                      protein.protein = genome_sequence.protein AND
                            comb.comb_id=comb_index.id AND
                            protein.protein = comb.protein AND
                            align.auto = ass.auto AND
                            ass.auto = family.auto AND
                            protein.seqid = '$seqid' AND
                            protein.genome = '$genome' AND
                            ass.sf = t1.id AND
                            family.fa = t2.id AND
                            ass.evalue <= $marginal
                      ORDER BY ass.evalue;" );
}
$sth->execute;

my $domain_no = 1;
while ( my @domain_details = $sth->fetchrow_array ) {
    push @domain_details, $domain_no;
    push @domains_AofHRs, getHRfor( \@domain_details, \@domain_attributes );
    $arch = $domain_details[6];
    push @models , $domain_details[2];
    $domain_no++;
}
#15705
#| dz     | FBpp0204060                                                  | 0036192 | 24-108          | 4.85e-21 | 46585 | 85259167 | dz_FBpp0204060                                                  |
#get the model numbers for each of the seeds used
 $sth =
      $dbh->prepare( "select model, seed from model where include ='y';" );

$sth->execute;

while ( my @lookups = $sth->fetchrow_array ) {
	#check to see if there is an entry for that seed in the file
    $model_lookup{$lookups[0]} = $lookups[1];
	}


#check whether that model has a coiled coil in it....
my $tempflag = 0;
#foreach my $model (@models){
#	    $sth =
#      $dbh->prepare( "select model.model, des.name from model,des where model.seed = des.id and model.model = '$model';" );

#$sth->execute;

#while ( my @model_des = $sth->fetchrow_array ) {
	#check to see if there is an entry for that seed in the file
	#print STDERR "$model_des[0]\t$model_lookup{$model_des[0]}\t$locations{$model_lookup{$model_des[0]}}\n";
#	if(exists($locations{$model_lookup{$model_des[0]}})){
#my @positions = split(/:/, $locations{$model_lookup{$model_des[0]}});
#foreach (@positions) {
#	my @terminals = split(/:/, $_);
#	my $start = $terminals[0];
#	my $finish = $terminals[1];
	
#	$tempflag = 1;
#}

		
#	}
	#take the start and stop for that file and map it to the allignment
	
	#create a new object for the file for PNG creation that has the coiled coil in it PLUS SOME FLAG FOR CHANGING THE COLOUR		
#}
#}

$sth = $dbh->prepare( "SELECT genome_sequence.sequence, genome_sequence.length, protein.comment, genome.name, protein.protein
                       FROM genome_sequence, genome, protein
                       WHERE protein.seqid = '$seqid' AND
                             protein.genome = '$genome' AND
                             protein.protein = genome_sequence.protein AND
                             protein.genome = genome.genome;" );
$sth->execute;
my ( $aa_seq, $aa_seq_length, $aa_seq_comment, $genome_name , $protein) = $sth->fetchrow_array;


#unless ( exists( $aa_seq ) && $aa_seq =~ /\w/ ) {
unless ( defined $aa_seq && $aa_seq =~ /\w/ ) {
    $sth->finish;
    $dbh->disconnect;
    my $title       = 'Sequence not in database, or no assignment for this sequence';
    my $explanation = "<strong>Sequence</strong>: $seqid not in database, or no assignment for this sequence.  <br/>\n" .
                      'If you think this is an error please report it using the feedback form below.';
    UtilityBelt::BombOut( $title, $explanation );
}


# Get link details for identical sequences
$sth = $dbh->prepare( "SELECT genome.genome, genome.name, protein.seqid
                       FROM genome, protein
                       WHERE protein.genome = genome.genome AND
                             genome.password = '' AND
                             protein.protein = '$protein' AND
                             (protein.genome != '$genome' AND protein.seqid != '$seqid')
                       ORDER BY genome.genome, protein.seqid;" );
$sth->execute;

my $previous_genome = '';
my $previous_genome_name;
my %id_seqs_H;
my @id_seqs;
while ( my ( $id_seq_genome, $id_seq_genome_name, $id_seq ) = $sth->fetchrow_array ) {
    if ( $id_seq_genome ne $previous_genome && $previous_genome ne '' ) {
        # Using anonymous array instead of reference to array because deleting @id_seqs
        # later.  Otherwise probs with identical sequences in _multiple_ genomes.
	    $id_seqs_H{ $previous_genome }->{ 'seqids' }      = [ @id_seqs ];
        $id_seqs_H{ $previous_genome }->{ 'genome_name' } = $previous_genome_name;
        @id_seqs = ();
    }
    push @id_seqs, $id_seq;

    $previous_genome      = $id_seq_genome;
    $previous_genome_name = $id_seq_genome_name;
}

# load in identical sequence link details for _final/only_ genome
if ( $sth->rows > 0 ) {
    $id_seqs_H{ $previous_genome }->{ 'seqids' }      = [ @id_seqs ];
    $id_seqs_H{ $previous_genome }->{ 'genome_name' } = $previous_genome_name;
}


$sth->finish;
$dbh->disconnect;


#------------------------------------------------------------
#Generate-PNG-FILES----------------------------------
my %regions;
my $start;
my $finish;
$seqid_file = $seqid;
$seqid_file =~ s/\#/_/g;
$seqid_file =~ s/\|/_/g;
$seqid_file =~ s/\//_/g;
#unless ( -e "$TEMPDIR/$seqid_file.png" ) {
    my $one  = 0;
my $coilsequence;
my $modelname;
my $olistate;
my $ccfound = 0;
    open PNGTMP, '>', "$TEMPDIR/$seqid_file"."cc.temp" or die "Cannot open $TEMPDIR/$seqid_file"."cc.temp: $!\n";
    foreach my $domains_HR ( @domains_AofHRs ) {
	my %domains_H = %{ $domains_HR };
        if ( $threshold >= $domains_H{'sf_eval'} ) {
        #            print PNGTMP "$length[$count] $unique $sfevalue[$count] ",
        #                 "$region[$count] $model[$count] $sf ",
        #                 "$superfamily\n";
            
 				my %sequence_hash;
         #   print STDERR "MODEL --> $domains_H{'model'}\n";
            	if(exists($locations{$model_lookup{$domains_H{'model'}}})){
            		$ccfound = 1;
					my @positions = split(/,/, $locations{$model_lookup{$domains_H{'model'}}});
					foreach (@positions) {
						$coilsequence = '';
						my @terminals = split(/-/, $_);
						$start = $terminals[0];
						$finish = $terminals[1];
					#	print STDERR "start->$start\tfinish->$finish";
						my @allign = split(//, $domains_H{'alignment'});
						my $seqcounter = 0;
						my $alligncounter = 0;
						my $on = 0; #=0
						my $off = 0;#=0
						my $coil_sequence;
						
						foreach (@allign){
							if($_ ne '-'){
							$seqcounter += 1;
							}
							if (($_ eq uc($_))&&($_ ne '-')){
								$alligncounter +=1;	
							}
							if ($alligncounter == $start){
								$on = $seqcounter;
							}	
							if (($on != 0)&&($off == 0)) {
								$coilsequence .= $_;
								if ($_ eq uc($_)){
								$sequence_hash{$seqcounter} = $_;
								}else{
									$sequence_hash{$seqcounter} = '-';
								}
							}
							if ($alligncounter == $finish){
								$off = $seqcounter;
								my $tempc = $seqcounter + 0.01;
								$sequence_hash{$tempc} = '/';	
							}
							
						}
						
						if ($off == 0){
							$off = $on + ($finish - $start);
							$coilsequence = substr $coilsequence, 0, ($off-$on);
						}
						if ($on == 0){
							$on = $start;	
							$off = $finish;
							$coilsequence = "HMM region not reachable";
						}
						my @tempreg = split(/-/, $domains_H{'region'});
						$on = $on+$tempreg[0];
						$off = $off+$tempreg[0];
						unless(exists($regions{"$start"."-"."$finish"})){
						print PNGTMP "$domains_H{'length'} $seqid 0 ",
                         "$on-$off 0 0 ",
                         "coiled coil\n";
                        
                        push @{$coils{$domains_H{'sf'}}{'region'}},"$on-$off";
                        push @{$coils{$domains_H{'sf'}}{'alignment'}},"$domains_H{'alignment'}";
                        $coils{$domains_H{'sf'}}{'align'} ="$domains_H{'alignment'}";
                        push @{$coils{$domains_H{'sf'}}{'sequence'}},"$coilsequence";
                        push @{$coils{$domains_H{'sf'}}{'strstp'}},"$start-$finish";
                        $coils{$domains_H{'sf'}}{'seqhash'} = \%sequence_hash;
                        my @tempreg = split(/-/, $domains_H{'region'});
#                        	my @seqArray= split (//, $domains_H{'alignment'});
#						my $cutAllign = Array2String(\@seqArray, $tempreg[0], $tempreg[1]);
 my $cutAllign = $domains_H{'alignment'};
						
                        $coils{$domains_H{'sf'}}{'cutallign'} ="$cutAllign";
                       $sth =    $dbh->prepare( "select name from des where id =$domains_H{'px'};" );
						$sth->execute;
						
						while ( my @lookups = $sth->fetchrow_array ) {
   							 $modelname = $lookups[0];
						}
						$sth =    $dbh->prepare( "select code from functional_annotation where sunid in ('$domains_H{'sf'}','$domains_H{'fa'}' )and level != '';" );
						$sth->execute;
						
						while ( my @lookups = $sth->fetchrow_array ) {
   							 $olistate = $lookups[0];
						}
						$coils{$domains_H{'sf'}}{'oli'} ="$olistate";
                        push @{$coils{$domains_H{'sf'}}{'nearest'}},"$modelname";
                        push @{$coils{$domains_H{'sf'}}{'px'}},"$domains_H{'px'}";
						$tempflag = 1;
						$regions{"$start"."-"."$finish"} = 1;
                        }
					}	
				}else{
					if (($domains_H{'sf'} == $sunid)||($domains_H{'fa'} == $sunid)){
					push @{$coils{$domains_H{'sf'}}{'region'}},"N/A";
					$coils{$domains_H{'sf'}}{'alignment'} ="$domains_H{'alignment'}";
                        push @{$coils{$domains_H{'sf'}}{'sequence'}},"There was no SOCKET prediction for this seed it is assumed due to family membership";
                          $sth =    $dbh->prepare( "select name from des where id =$domains_H{'px'};" );
						$sth->execute;
						
						while ( my @lookups = $sth->fetchrow_array ) {
   							 $modelname = $lookups[0];
						}
						$sth =    $dbh->prepare( "select code from functional_annotation where sunid = '$domains_H{'sf'}' and level != '';" );
						$sth->execute;
						
						while ( my @lookups = $sth->fetchrow_array ) {
   							 $olistate = $lookups[0];
						}
						$coils{$domains_H{'sf'}}{'oli'} ="$olistate";
                        
                        push @{$coils{$domains_H{'sf'}}{'nearest'}},"$modelname";

                        push @{$coils{$domains_H{'sf'}}{'px'}},"$domains_H{'px'}";
				}
				}
				print PNGTMP "$domains_H{'length'} $seqid $domains_H{'sf_eval'} ",
                         "$domains_H{'region'} $domains_H{'model'} $domains_H{'sf'} ",
                         "$domains_H{'superfamily'}\n";
                        #print PNGTMP "$domains_H{'alignment'} $seqid $domains_H{'sf_eval'} ",
                        # "$domains_H{'region'} $domains_H{'model'} $domains_H{'sf'} ",
                        # "$domains_H{'superfamily'}\n";
            $one = 1;
        }
    }

    if ( $one == 0 ) {
        print PNGTMP "NoDomain\t$aa_seq_length\t$aa_seq_length\n";
    }
    if ($tempflag == 1){

    }
   
    close PNGTMP;
	
    my @args = ( "$TEMPDIR/$seqid_file"."cc.temp", 'gene' );
    system( './domain_plot.pl', @args );

    if ( $one == 0 ) {
        $seqid_file = $aa_seq_length;
    }
#}


#----------------------------------------------------

#Start-PRINTING--------------------------------------

my $heading = "Domain assignment for $seqid from $genome_name";
my $title   = Jgough->htmlStrip( $heading );
my @breadcrumbs = ( "taxonomic_gen_list.cgi\tGenomes",
                    "gen_list.cgi?genome=$genome\t$genome_name",
                    $seqid,
                  );
my $unformatted_genome_name = Jgough->htmlStrip( $genome_name );
my $meta_desc = 'Domain architecture and assignment details '
            .   "(superfamily, family, region, evalue) for $seqid from "
            .   "$unformatted_genome_name.  Plus protein sequence and "
            .   'external database links.';

print UtilityBelt::spiri_header( $title, \@breadcrumbs, $meta_desc );

#print "$DEBUG";
print "<a name=\"section_domain_architecture\" />\n",
      "<h1 class=\"c2\">$heading</h1>\n",
      "<strong>Domain architecture</strong>\n",
      "<p><a href=\"./allcombs.cgi?comb=$arch$passlink\"><img src=\"$TEMPDIR/$seqid_file",
      "gene.png\" alt = \"allcombs\" border=\"0\"/></a></p>\n",
      "<p>Click on the picture above to see other proteins with the same domain architecture.</p>\n\n";
# print     "$TEMPDIR/$seqid_file"."cc.temp";

if($ccfound == 0){
	print "<h3> This is an exception </h3>";
	print '<p>The closest structure to this sequence is assumed to contain a coiled coil based on the fact that the majority of the members of its superfamily/family contained a structure. For this reasaon there is no SOCKET annotation of the structure and we are unable to identify the coiled coil. Hopefully there is some obvious reason that this has occured, for instance the structure is low resolution or has no side chains, if you would like to enquire further please contact <a href="mailto:owen.rackham@bristol.ac.uk?subject=Query_about_a_missing_SOCKET_output_for_';
print "$seqid";
print "_from_$genome".'">Owen Rackham</a></p>';
	
}else{   
print "<a name=\"loading\"></a>" ;
print "<div id=\"loading\" class=\"loading-invisible\">  <p><img src=\"../spiricoil/images/loading.gif\" alt=\"Please 
wait.....\" /> please be patient the model is being generated, this can take up
to a minute or two.</p> To reload the original page press F5 or refresh at any time. </div>"	;
print "<div id=\"jmolviewer\" >";
foreach my $sfam (keys %coils){
	my $name = $coils{$sfam}{'oli'};
	
print "<strong>Details of the coiled coil assignments for $sfam in the state $name </strong>\n",
      "\n\n";
	print "<table>";
	print "<tr class=\"toprow\"><td class = \"top\">region</td><td class = \"top\">sequence</td></tr><hr/>";
	my @regions = @{$coils{$sfam}{'region'}};	
	my @sequences = @{$coils{$sfam}{'sequence'}};
	my %seqhash = %{$coils{$sfam}{'seqhash'}};
	my $totalseq;
	foreach my $key (sort { $a <=> $b } keys %seqhash) {
     $totalseq .= "$seqhash{$key}";
}





	for my $i ( 0 .. scalar(@regions) - 1 ) {
		print "<tr class =\"otherrow\"><td class = \"other\">$regions[$i]</td><td class = \"other\">$sequences[$i]</td></tr>";

	}
	print "<tr><td class = \"bottom\"></td><td class = \"bottom\"></td></tr>";

print "</table>";
#my @locations = @{$coils{$sfam}{'region'}};
my @locations = @{$coils{$sfam}{'strstp'}};
my $file = "$coils{$sfam}{'nearest'}[0].ent";
my $folder = substr $file, 2,2;
my $res; 
my $chn;
$residue_count = 1;


if ($real == 0){

print "<script type=\"text/javascript\">";
print "  document.getElementById(\"loading\").className = \"loading-visible\"; ";
print "  var hideDiv = function(){document.getElementById(\"loading\").className = \"loading-invisible\";}; ";
print "  var oldLoad = window.onload; ";
print "  var newLoad = oldLoad ? function(){hideDiv.call(this);oldLoad.call(this);} : hideDiv; ";
print "  window.onload = newLoad;  ";
print "</script> ";
print "<p>";
		($res, $chn) = getNearestLocation("../db/astral_pdb/$folder/$file", "$SPIRITEMP/ORIGMOD$coils{$sfam}{'nearest'}[0].ent",\@locations, $coils{$sfam}{'nearest'}[0]);
print "</p>";	

print "This is the section of the closest known structure that contains the coiled coil <br/>";
print "<script type=\"text/javascript\"> jmolInitialize(\"../spiricoil/jmol\"); jmolApplet(400, \"load ../temp/spiricoil/ORIGMOD$coils{$sfam}{'nearest'}[0].ent;cartoon; color cartoon structure\"";

print ");</script>";

print "<br/>";
#print "$totalseq";

#ccsubmit.cgi?fn=../temp/spiricoil/11801356290565.B99990001.pdb
#print "<br/><input type=\"button\" value=\"SOCKET\" onclick=\"location.href='./ccsubmit.cgi?fn=../temp/spiricoil/$SPIRITEMP/ORIGMOD$coils{$sfam}{'nearest'}[0].ent'\">";
print "<br/><input type=\"button\" value=\"MODELLER\" onclick=\"location.href='#loading';location.href='./ccgene.cgi?genome=$genome&seqid=$seqid&sunid=$sunid&real=2';document.getElementById(\'jmolviewer\').className = \'loading-invisible\';document.getElementById(\'loading\').className = \'loading-visible\'; \"/>";
print "(may take a few minutes to load)<br/>";
print "click here to automatically generate a 3D homology model
using MODELLER <link> based on the structure and alignment above. We strongly advise, to understand the meaning of the  
model, that you refer to relevant literature on homology modelling and the <a href = 
\"http://predictioncenter.org/index.cgi?page=proceedings\">CASP</a> competition.";
print "<br/><input type=\"button\" value=\"SOCKET\" onclick=\"location.href='./ccsubmit.cgi?fn=../db/astral_pdb/$folder/$file'\"/>";
print "(may take a few minutes to load)<br/>";
print "click here to submit the original structure to the SOCKET server based at Sussex University. (nb. This is not the same version of SOCKET as we used.)";
print "<p/>";
print "<br/>";


}
############################################HACK FIRST MODEL SECOND#############################
if ($real == 1){

		($res, $chn) = getNearestLocation("../db/astral_pdb/$folder/$file", "$SPIRITEMP/ORIGMOD$coils{$sfam}{'nearest'}[0].ent",\@locations, $coils{$sfam}{'nearest'}[0]);
my $id = substr $coils{$sfam}{'nearest'}[0],1 ,4 ;
my $ch = uc(substr $coils{$sfam}{'nearest'}[0],5 ,1 );

	my $unique = int( rand(999999999999999) );
	    open ALI, '>', "$SPIRITEMP/align"."$unique"."cc.temp" or die "Cannot open $TEMPDIR/$seqid_file"."cc.temp: $!\n";
$totalseq = uc($totalseq);
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:MOD$coils{$sfam}{'nearest'}[0]:   1 : $chn : $res : $chn ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
print ALI "$totalseq"."*\n";
close ALI;
 my @args = ("./model-single.py","align"."$unique"."cc.temp","$id"."$ch","$unique",">&/dev/null" );
 
 system( '/usr/bin/python', @args );


print "<script type=\"text/javascript\"> jmolInitialize(\"../spiricoil/jmol\"); jmolApplet(400, \"load ../temp/spiricoil/$unique.B99990001.pdb;cartoon; color cartoon structure\");</script>";
print "<br/>TOOLS<br/>";
print "<br/><input type=\"button\" value=\"SOCKET\" onclick=\"location.href='./ccsubmit.cgi?fn=../temp/spiricoil/$unique.B99990001.pdb'\"/>";
print "<p/>";
}
############################################################################################
#####################################WHOLE STRUCTURE#########################################
if ($real == 2){


#Ichi"244
#	print("cp  $SPIRITEMP/MOD$coils{$sfam}{'nearest'}[0].ent");
	
my $id = substr $coils{$sfam}{'nearest'}[0],1 ,4 ;
my $ch = uc(substr $coils{$sfam}{'nearest'}[0],5 ,1 );

	my $unique = int( rand(999999999999999) );
	    open ALI, '>', "$SPIRITEMP/align"."$unique"."cc.temp" or die "Cannot open $TEMPDIR/$seqid_file"."cc.temp: $!\n";
$totalseq = uc($totalseq);
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:../../db/astral_pdb/$folder/$file:   FIRST : $ch : LAST : $ch ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
#print ALI "$coils{$sfam}{'cutallign'}"."*\n"; 
#$coils{$sfam}{'align'};
my $altemp = $coils{$sfam}{'cutallign'};
$altemp = uc($altemp);
my $tl = length($altemp);
print ALI "$altemp"."*\n";

close ALI;
 my @args = ("./model-single.py","align"."$unique"."cc.temp","$id"."$ch","$unique",">&/dev/null" );
 
 system( '/usr/bin/python', @args );

my ($res, $chn) = getNearestLocation("$SPIRITEMP/$unique.B99990001.pdb", "$SPIRITEMP/FULLMOD$coils{$sfam}{'nearest'}[0].ent",\@locations, $coils{$sfam}{'nearest'}[0]);
#print "<script type=\"text/javascript\"> jmolInitialize(\"../spiricoil/jmol\"); jmolApplet(400, \"load $SPIRITEMP/FULLMOD$coils{$sfam}{'nearest'}[0].ent;cartoon; color cartoon structure\");</script>";
#px --> $coils{$sfam}{'px'}[0] .. $coils{$sfam}{'px'}[1].. $coils{$sfam}{'px'}[2] .. $coils{$sfam}{'px'}[3] ID --> $id --> chain --> $ch  --> align"."$unique"."cc.temp --> 
print "In order to download this model <a href=\"$SPIRITEMP/$unique.B99990001.pdb\">click here</a> which is a 3D model produced using MODELLER  <br/>";
print "<script type=\"text/javascript\"> jmolInitialize(\"../spiricoil/jmol\"); jmolApplet(400, \"load $SPIRITEMP/$unique.B99990001.pdb;cartoon;cartoon structure;";
foreach my $loc (@locations){
	print "select $loc:$chn; color blue;";
}
print "\");</script>";
#print "LENGTH --> $tl";

print "<br/><input type=\"button\" value=\"SOCKET\" onclick=\"location.href='./ccsubmit.cgi?fn=$SPIRITEMP/$unique.B99990001.pdb'\"/>";
print "(may take a few minutes to load)<br/>";
print "click here to submit the modified structure to the SOCKET server based at Sussex University. (nb. This is not the same version of SOCKET as we used.)";


print "<p/>";
}
##############################################################################################
}



print "<strong>Domain combinations with similar phylogenetic distribution</strong>\n",
      "<p><a href=\"yiduo.cgi?comb=$arch$passlink\">See other domain ",
      "architectures which have a similar genomic distribution.</a></p>\n",
      "<hr/>\n\n";


#-------------------------------------------------------------------------------------------------------
#
# Domain assignments
#


@href = Href( $seqid, $genome, $dbh );
my $strong_hit_count = getHitCount( 'sf_eval', 'Strong', \@domains_AofHRs );
my $weak_hit_count   = getHitCount( 'sf_eval', 'Weak', \@domains_AofHRs );


$exclude_link = '#section_domain_assignments';
print UtilityBelt::internal_page_section_links( \@include_links, $exclude_link );

print "<a name=\"section_domain_assignments\" />\n",
      "<h3>Domain assignment details</h3>\n";
if ( $strong_hit_count >= 1 ) {
    my $hit_type = 'Strong';
    print "<strong>Strong hits</strong>\n";
    printDomainTableHeaders( $genome, $seqid, $href[0], $hit_type );
}
else {
    print "<p><span class=\"red_text\">No Significant Hits.</span></p>\n\n";
}


#Go through printing out all of the domains assigned--------------
my $seen_weak_hits_yet = 0;
foreach my $domains_HR ( @domains_AofHRs ) {
    my %domains_H = %{ $domains_HR };

    if ( $domains_H{'sf_eval'} > $threshold && $seen_weak_hits_yet == 0 ) {
	    $seen_weak_hits_yet = 1;
        my $hit_type = 'Weak';
	    print "</table>\n\n" if ( $strong_hit_count >= 1 );
	    print "<strong>Weak hits</strong>\n";
        printDomainTableHeaders( $genome, $seqid, $href[0], $hit_type );
    }

    printDomainDetails( \%domains_H );

    my $unique = int( rand(999999999999999) );
    printDomainTableFamilyLink( \%domains_H, $genome, $seqid, $seqid_file );
    printDomainTableAlignLink( \%domains_H, $seqid, $unique, $passlink, $genome );
    printDomainTableGenomeLink( \%domains_H, $genome, $passlink );
    printDomainTableAllcombsLink( \%domains_H, $passlink );
}

print "</table>\n" if ( $strong_hit_count >= 1 or $weak_hit_count >= 1 );


print "\n<p class=\"text_justify\">The results are sorted from lowest E-value to highest E-value.  ",
      'Strong classifications have a low E-value.  ',
      "Weak classifications have an E-value greater than $threshold.  ",
      "Weak hits are shown in gray.  Weak hits are not shown on the domain architecture.</p>\n",
      "<p class=\"text_justify\">The family level classification is conditional on the domain being a member of the specified superfamily.  ",
      "There is a possibility that the selected domain is a member of a sub-family for which no structure has yet been solved.</p>\n",
      "<hr/>\n\n";


#-------------------------------------------------------------------------------------------------------
#
# Protein sequence details
#


$exclude_link = '#section_protein_sequence';
print UtilityBelt::internal_page_section_links( \@include_links, $exclude_link );

print "<a name=\"section_protein_sequence\" />\n",
      "<h3>Protein sequence</h3>\n",
      "<table class=\"funcsum\" cellpadding=\"2\" cellspacing=\"2\">\n",
      "  <tr>\n",
      "    <td>External link</td>\n",
      "    <td><a href=\"$href[0]\" target=\"$seqid\">$seqid</a></td>\n",
      "  </tr>\n";

print "  <tr>\n",
      "    <td>Sequence length</td>\n",
      "    <td>$aa_seq_length</td>\n",
      "  </tr>\n";

#if ( $aa_seq_comment !~ /^$/ ) {
if ( $aa_seq_comment ne '' ) {
    print "  <tr>\n",
          "    <td valign=\"top\">Comment</td>\n",
          "    <td wrap=\"wrap\">$aa_seq_comment</td>\n",
          "  </tr>\n";
}

print "  <tr>\n",
      "    <td valign=\"top\">Sequence</td>\n",
      "    <td>\n<pre>\n",
      Jgough->WrapSeq( $aa_seq ),
      "</pre></td>\n",
      "  </tr>\n";

$seqid_file = $seqid;
$seqid_file =~ s/\#/_/g;
$seqid_file =~ s/\|/_/g;
$seqid_file =~ s/\//_/g;
print "  <tr>\n",
      "    <td valign=\"top\">Download sequence</td>\n",
      "    <td><a href=\"save.cgi?var=$genome\_$seqid;type=protein_sequence;filename=$seqid_file.fa\">",
      "<img src=\"../pics/save.png\" alt = \"save.png\" title=\"FASTA formatted sequence download\" border=\"0\"/></a></td></tr>\n";

if ( keys %id_seqs_H ) {
    print "  <tr>\n",
          "    <td valign=\"top\">Identical sequences</td>\n",
          '    <td>';

    # Link to UniProt seqs in sf.o first because we will add
    # the iproclass/biothesaurus data to up seqs Real Soon Now (TM)
    if ( exists $id_seqs_H{ 'up' } ) {
        my $id_genome_name = $id_seqs_H{ 'up' }->{ 'genome_name' };
        $id_genome_name = Jgough->htmlStrip( $id_genome_name );
        #print $id_genome_name . ": ";

        foreach my $id_seqid ( @{ $id_seqs_H{ 'up' }->{ 'seqids' } } ) {
            print "<a href=\"gene.cgi?genome=up;seqid=$id_seqid\" title=\"$id_genome_name\">$id_seqid</a> ";
        }

        print "<br/>\n";
    }

    foreach my $id_gn ( keys %id_seqs_H ) {
   	next if ( $id_gn eq 'up' );

	my $id_genome_name = $id_seqs_H{ $id_gn }->{ 'genome_name' };
        $id_genome_name = Jgough->htmlStrip( $id_genome_name );
        #print $id_genome_name . ": ";

        foreach my $id_seqid ( @{ $id_seqs_H{ $id_gn }->{ 'seqids' } } ) {
            print "<a href=\"gene.cgi?genome=$id_gn;seqid=$id_seqid\" title=\"$id_genome_name\">$id_seqid</a> ";
        }
    }

    print "    </td>\n",
          "  </tr>\n";
}

print "</table>\n\n";

}
print "</div>";
print UtilityBelt::spiri_footer;


#--------------------------------------------------------------------------------------------------------
#
# Subroutines
#
sub getNearestLocation {
    my $file = shift;
    my $output = shift;
    my $location    = shift;
    my $id = shift;
    my @location = @{$location};
#ATOM      1  N   ARG     1      41.059  40.597   3.734  1.00122.74 
#ATOM   4657  CA  MET C 578     482.228 -78.615 155.302  1.00174.00           C  
#ATOM      1  N   LYS A   1      -0.769  25.488  87.359  1.00 48.43           N
open OUT, ">$output" or die $!;
#open OUT2, ">$output"."new" or die $!;
print OUT "HEADER    SCOP/ASTRAL/SPIRICOIL coiled coil structure   0000\n";
my $init;
my $flag = 0;
foreach my $region (@location){
	
	my @regions = split(/-/, $region);
	my $start = $regions[0];
	my $finish = $regions[1];
open FILE, "<", "$file" or die $!;
	my $prevpos = 1;
	my $pos = 99999;
	my $count = 1;
while (<FILE>) { 
	my $linestart;
	my $lineend;
	my $line = $_;
	#ATOM    303  CD1 ILE A  38      21.475  23.101  39.677  1.00  | ATOM    402  NH1AARG A  48       25  24.345  46.083  0.50 37.
	#ATOM      1  N   ARG     1      41.059  40.597   3.734  1.00122.74
	if ($line =~ m|(^.{23}).{3}(.*)|){
	$linestart = $1;
	$lineend = $2;
	}
	if ($line =~ m|^ATOM.{16}(.{3})(.{3})|){
		chomp;
	#if ($line =~ m|^ATOM\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)|){
		if ($count == 1){
			$prevpos = $2;	
			$count += 1;
		}	
		my $chn = $1;
		$pos = $2;
		
		if ($flag == 0){
			$init = $pos;
			#print "INIT $init";
			$flag = 1;	
		}
		
		if (($pos >= ($start+$init))&&($pos<=($finish+$init))){
			unless ($pos == $prevpos){
			$residue_count += 1;
			
			}
			my $num = pad($residue_count);
			print OUT "$linestart"."$num"."$lineend\n";
			$chn_store = $chn;
			#print OUT2 $line;
			
			
		}else{
			$count = 1;
		}
		$prevpos = $pos;		
	}else{
		#print "THERE WAS NO MATCH";
	}
}
}
print ".<br/>";
print OUT "END\n";
return ("$residue_count","$chn_store");
}

sub pad {
	my $number  = shift;
	
	if (length($number) == 1){
		return "  $number";	
	}elsif(length($number) == 2){
		return " $number";	
	}else{
		return $number;	
	} 
}

sub Href {
    my $seqid  = shift;
    my $genome = shift;
    my $dbh    = shift;
    my $href;
    my $pre_parse_seqid = $seqid;

    my $sth = $dbh->prepare( 'SELECT gene_link, parse FROM genome WHERE genome = ?' );
    $sth->execute($genome);
    my @temp = $sth->fetchrow_array;
    $sth->finish;

    # evaluate the parse regular expression
    # if the regex fails seqid is null and errors are sent to the apache log file
    # it is tempting to add, || '', to the eval.
    # don't, monitor apache log file and fix broken regexs in parse field of genome table
    # mysql requires backslashes to be escaped, but this is done in the Info.pm module
    if ( length $temp[1] > 3 ) {
        $_     = $seqid;
        $seqid = eval( $temp[1] );
    }

    if ( length $temp[0] > 0 ) {
        my $url_seqid = Jgough->urlEncode($seqid) if defined $seqid;

        if ( defined $url_seqid && $url_seqid ne '' ) {
            $href = $temp[0] . $url_seqid;
	    }
        else {
            $href = '../nolink.html';
            print STDERR "Check genome.parse problem with: $genome $pre_parse_seqid\n";
        }
    }
    else {
        $href = '../nolink.html';
    }

    return ($href);
}

#cmd = 'rm %s.*' % sys.argv[3]
#---------------------------------------------------------

#GENOME-PASSWORD---------------------------------------
#check password
sub GenomePassword {
    my $theone = $_[0];

    my $go;
    my $sth;

    $sth = $dbh->prepare( "SELECT password FROM genome WHERE genome = '$theone';" );
    $sth->execute;
    my $passwd = $sth->fetchrow;
    $sth->finish;

    unless ( defined $passwd ) {
	$passwd = '';
    }

    if ( $sth->rows == 0 ) {
        $go = 0;
    }
    elsif ( $passwd eq '' ) {
        $go = 1;
    }
    else {
        if ( $genpass eq $passwd ) {
            $go = -1;
        }
        else {
            $go = 0;
        }
    }

    return ( $go, );
}


#-----------------------------------------------------



sub getHRfor {
    my $sql_out_array_ref = shift;
    my $attributes_ref    = shift;
    my @attributes        = @{ $attributes_ref };
    my @domains;

    my %rec;
    @rec{@attributes} = @{ $sql_out_array_ref };
    push @domains, \%rec;

    return @domains;
}


sub getHitCount {
    my $attr               = shift;
    my $hit_type           = shift;
    my $domains_AofHRs_ref = shift;
    my @domains_AofHRs     = @{ $domains_AofHRs_ref };
    my $count = 0;

    foreach my $domains_HR ( @domains_AofHRs ) {
        my %domains_H = %{ $domains_HR };

        if ( $hit_type eq 'Strong' ) {
            $count++ if ( $domains_H{ "$attr" } <= $threshold );
        }
        else {
            $count++ if ( $domains_H{ "$attr" } > $threshold );
        }
    }

    return $count;
}


sub printDomainTableHeaders {
    my $genome   = shift;
    my $seqid    = shift;
    my $out_link = shift;
    my $hit_type = shift;
    my $hit_color = ' class="domaintable_header"';

    if ( $hit_type eq 'Weak' ) {
        $hit_color = ' class="grey_table_text"';
    }

    print "<p></p>\n",
          "<table class=\"domaintable\">\n",
          "  <tr$hit_color style=\"background-color: #CCCCFF;\">\n",
          "     <td><strong>Sequence:&nbsp;</strong></td>\n",
          "     <td colspan=\"2\"><a href=\"$out_link\" target=\"outlink_$genome\">$seqid</a></td>\n",
          "  </tr>\n";

    return 1;
}


sub printDomainDetails {
    my $domains_H_ref = shift;
    my %domains_H = %{ $domains_H_ref };
    my $sf_hit_colour   = '';
    my $fa_hit_colour   = '';
    my $weak_hit_colour = '';

    if ( $domains_H{'sf_eval'} > $threshold ) {
        $sf_hit_colour = ' class="grey_table_text"';
        $domains_H{'domain_no'} = '-';
    }

    if ( $domains_H{'fa_eval'} > $threshold ) {
        $fa_hit_colour = ' class="grey_table_text"';
    }

    if ( $domains_H{'sf_eval'} > $threshold && $domains_H{'fa_eval'} > $threshold ) {
        $weak_hit_colour = ' class="grey_table_text"';
    }

    print "  <tr$weak_hit_colour>\n",
          "    <td class=\"grey_col_bg\"><strong>Domain Number&nbsp;</strong>$domains_H{'domain_no'}</td>\n",
          "    <td colspan=\"2\"><strong>Region:&nbsp;</strong>$domains_H{'region'}</td>\n",
          "  </tr>\n",
          "  <tr$weak_hit_colour>\n",
          "    <td><strong>Classification Level</strong></td>\n",
          "    <td><strong>Classification</strong></td>\n",
          "    <td><strong>E-value</strong></td>\n",
          "  </tr>\n",
          "  <tr$sf_hit_colour>\n",
          "    <td>Superfamily</td>\n",
          "    <td nowrap=\"nowrap\"><a href=\"scop.cgi?sunid=$domains_H{'sf'}\" target=\"$domains_H{'model'}\">$domains_H{'superfamily'}</a></td>\n",
          "    <td>$domains_H{'sf_eval'} </td>\n",
          "  </tr>\n",
          "  <tr$fa_hit_colour>\n",
          "    <td>Family</td>\n",
          "    <td nowrap=\"nowrap\"><a href=\"scop.cgi?sunid=$domains_H{'fa'}\" target=\"$domains_H{'fa'}\">$domains_H{'family'}</a></td>\n",
          "    <td>$domains_H{'fa_eval'} </td>\n",
          "  </tr>\n";

    return 1;
}


sub printDomainTableFamilyLink {
    my $domains_H_ref = shift;
    my $genome        = shift;
    my $seqid         = shift;
    my $seqid_file    = shift;
    my %domains_H = %{ $domains_H_ref };
    my $weak_hit_colour = '';

    if ( $domains_H{'sf_eval'} > $threshold && $domains_H{'fa_eval'} > $threshold ) {
        $weak_hit_colour = ' class="grey_table_text"';
        $domains_H{'domain_no'} = '-';
    }

    print "  <tr$weak_hit_colour>\n",
          "    <td><strong>Further Details: </strong></td>\n",
          "    <td align=\"left\" colspan=\"2\">\n",
          '      <a href="family.cgi?gene=y;',
                                     "unique=$genome;",
                                     "px=$domains_H{'px'};",
                                     "fa=$domains_H{'fa'};",
                                     "familyevalue=$domains_H{'fa_eval'};",
                                     "id=$seqid;",
                                     "evalue=$domains_H{'sf_eval'};",
                                     "model=$domains_H{'model'};",
                                     "domain=$domains_H{'domain_no'};",
                                     "picture=$TEMPDIR/$seqid_file",
                                     'genecc.png;',
                                     "threshold=$threshold\" border=\"0\">";
    print "<img src=\"../pics/family-details.png\" alt = \"familiy details\" border=\"0\" title=\"Family details for this domain\"/></a>&nbsp;\n";

    return 1;
}


sub printDomainTableAlignLink {
    my $domains_H_ref = shift;
    my $seqid         = shift;
    my $unique        = shift;
    my $passlink      = shift;
    my $genome        = shift;
    my %domains_H = %{ $domains_H_ref };

    print "      <a href=\"align.cgi?model=$domains_H{'model'};",
                                    "sf=$domains_H{'sf'};",
                                    "cgi_$seqid\_$domains_H{'region'}=1;",
#                                    "gen_$genome=1;",
                                    "unique=$unique;",
                                    'seed=1;',
                                    'local=Local',
                                    "$passlink\" border=\"0\">";
    print "<img src=\"../pics/alignments.png\" border=\"0\" alt = \"alignments\" title=\"Alignment for this sequence to a model\"/></a>&nbsp;\n";

    return 1;
}


sub printDomainTableGenomeLink {
    my $domains_H_ref = shift;
    my $genome        = shift;
    my $passlink      = shift;
    my %domains_H = %{ $domains_H_ref };

    print "      <a href=\"genome.cgi?model=$domains_H{'model'};",
                                     "cgi_$genome=yes;",
                                     "sf=$domains_H{'sf'}",
                                     "$passlink\" border=\"0\">";
    print "<img src=\"../pics/genome-assignments.png\" alt = \"genome assignments\" border=\"0\" title=\"Genome assignments for this superfamily\"/></a>&nbsp;\n";

    return 1;
}


sub printDomainTableAllcombsLink {
    my $domains_H_ref = shift;
    my $passlink      = shift;
    my %domains_H = %{ $domains_H_ref };

    print "      <a href=\"allcombs.cgi?sf=$domains_H{'sf'}",
                                       "$passlink\" border=\"0\">";
    print "<img src=\"../pics/domain-combinations.png\" alt = \"domain combinations\" border=\"0\" title=\"Domain combinations for this superfamily\"/></a>\n",
          "    </td>\n",
          "  </tr>\n\n";

    print "  <tr>\n",
          "    <td style=\"border: none;\" align=\"center\" colspan=\"4\">&nbsp;</td>\n",
          "  </tr>\n\n";

    return 1;
}

#sub Array2String {
#	my $array = shift;
#	my $start = shift;
#	my $finish = shift;
#	my @array = @{$array};
#	my $string;
#
#	for my $i ( ($start-1) .. ($finish-1) ) {
#		#if ($array[$i] eq uc($array[$i])){
#			$string .= $array[$i];	
#	#	}else{
#	#		$string .= '-';
#	#	}
#	}
#	return $string;
#}


