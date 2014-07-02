#!/usr/bin/perl -w

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use CGI qw/:standard/;
use DBI;
use lib '/home/rackham/modules/';
use rackham;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
#Digest::MD5::md5_hex()
my ( $dbh, $sth );
$dbh = rackham::DBConnect('superfamily');

sub checkCase {
    if ($_[0] =~ /^[[:upper:]]/) {
        return 1;
    }
    else {
        return 0;
    }
}

sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };

sub get_astral_missing {
        my ($id) = @_;
        my $filename = '/home/rackham/workspace/data/astral_map/astral-rapid-access-1.75.raf.tsv';
        open my $fh, $filename or die "Could not open $filename: $!";
        my @lines = grep /\Q$id/, <$fh>;
	my @in_or_out;
        my @cellMessage = split("\t", $lines[0]);
		my $counter = 1;
		my $string_location;
		my $l;
                foreach(@cellMessage){
			$l = $_;
                    	if($counter++ < 6){
				#print "next\n";
				next;
			}elsif($_ eq 'B'){
				#print "reset\n";
				$counter = 4;
			}elsif($_ eq 'M'){
				$string_location = "missing in the middle";
				push(@in_or_out,0);
			}elsif($_ eq 'E'){
				$string_location = "missing at the end";
				push(@in_or_out,0);	
			}elsif(looks_like_number($l)){
				#print "string\n";
				push(@in_or_out,1);
				$string_location = $l;
			}else{
				#print "$string_location --> $l\n";
				#push(@in_or_out,1);
			}

		$counter++;
                }
	return \@in_or_out;
}

#1n9pA 0.02 38 070604 111010   43  370    B .g   B .s   B .k   B .k  43 rr  44 qq  45 rr  46 ff  47 vv  48 dd  49 kk  50 nn  51 gg  52 rr  53 cc  54 nn  55 vv  56 qq  57 hh   M .g   M .n   M .l   M .g   M .s  63 ee 190 rr 191 aa 192 ee 193 tt 194 ll 195 mm 196 ff 197 ss 198 ee 199 hh 200 aa 201 vv 202 ii 203 ss 204 mm 205 rr 206 dd 207 gg 208 kk 209 ll 210 tt 211 ll 212 mm 213 ff 214 rr 215 vv 216 gg 217 nn 218 ll 219 rr 220 nn 221 ss 222 hh 223 mm 224 vv 225 ss 226 aa 227 qq 228 ii 229 rr 230 cc 231 kk 232 ll 233 ll 234 kk 235 ss 236 rr 237 qq 238 tt 239 pp 240 ee 241 gg 242 ee 243 ff 244 ll 245 pp 246 ll 247 dd 248 qq 249 ll 250 ee 251 ll 252 dd 253 vv 254 gg 255 ff 256 ss 257 tt 258 gg 259 aa 260 dd 261 qq 262 ll 263 ff 264 ll 265 vv 266 ss 267 pp 268 ll 269 tt 270 ii 271 cc 272 hh 273 vv 274 ii 275 dd 276 aa 277 kk 278 ss 279 pp 280 ff 281 yy 282 dd 283 ll 284 ss 285 qq 286 rr 287 ss 288 mm 289 qq 290 tt 291 ee 292 qq 293 ff 294 ee 295 vv 296 vv 297 vv 298 ii 299 ll 300 ee 301 gg 302 ii 303 vv 304 ee 305 tt 306 tt 307 gg 308 mm 309 tt 310 cc 311 qq 312 aa 313 rr 314 tt 315 ss 316 yy 317 tt 318 ee 319 dd 320 ee 321 vv 322 ll 323 ww 324 gg 325 hh 326 rr 327 ff 328 ff 329 pp 330 vv 331 ii 332 ss 333 ll 334 ee 335 ee 336 gg 337 ff 338 ff 339 kk 340 vv 341 dd 342 yy 343 ss 344 qq 345 ff 346 hh 347 aa 348 tt 349 ff 350 ee 351 vv 352 pp 353 tt 354 pp 355 pp 356 yy 357 ss 358 vv 359 kk 360 ee 361 qq 362 ee 363 ee 364 mm 365 ll 366 ll 367 mm 368 ss 369 ss 370 pp   E .l

my $seqid = $ARGV[0];
my $genome = $ARGV[1];

my $digest;
my $sth2=$dbh->prepare( "select sequence from genome_sequence,protein where protein.protein = genome_sequence.protein and protein.seqid = '$seqid' and protein.genome='$genome';" );
	$sth2->execute;
	while ( my @temp2 = $sth2->fetchrow_array ) {
		$digest = md5_hex($temp2[0]);
	}

my $fo = $ARGV[2];
my $marginal = 0.01;
#my $TEMPDIR = '/media/scratch/genome3D/pfam_assfiles/models/';
my $TEMPDIR = '../data';
my $s = $dbh->prepare("SELECT ass.evalue, ass.region, ass.model, ass.sf, t1.description, align.alignment, comb_index.comb, family.evalue, family.px, family.fa, t2.description, genome_sequence.length,protein.protein,align.modstart
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
	 $domain_details{$domain_no}{'modstart'} = $temp[13]; 
    my $sth2=$dbh->prepare( "select name from des where id =$temp[8];" );
	$sth2->execute;
	while ( my @temp2 = $sth2->fetchrow_array ) {
		$domain_details{$domain_no}{'closest_structure'} = $temp2[0];
	}
	$domain_no++;
}
#print Dumper \%domain_details;
foreach my $domain_number (keys %domain_details){
if(-e "$TEMPDIR"."$seqid"."_"."$domain_details{$domain_number}{'start'}"."_"."$domain_details{$domain_number}{'stop'}".".pdb"){
	my $cmd = "cp "."$TEMPDIR"."$seqid"."_"."$domain_details{$domain_number}{'start'}"."_"."$domain_details{$domain_number}{'stop'}".".pdb"." .";
	system($cmd);
}else{

my $id = substr $domain_details{$domain_number}{'closest_structure'},1 ,4 ;
my $ch = uc(substr $domain_details{$domain_number}{'closest_structure'},5 ,1 );
my $in_or_out = get_astral_missing($id.$ch);
my $file = "$domain_details{$domain_number}{'closest_structure'}.ent";
my $folder = substr $file, 2,2;
		my $unique ="$seqid"."_"."$domain_details{$domain_number}{'region'}";
        #my $unique = int( rand(999999999999999) );
     
unless (-e "$TEMPDIR"."/"."$fo"."/"."$unique".".pdb") {      
open ALI, '>', "$TEMPDIR"."/"."$fo/align"."$unique".".temp" or die "Cannot open $TEMPDIR/$fo/align"."$unique".".temp".": $!\n";
print ALI ">P1;$id"."$ch\n";
print ALI "structureX:/media/scratch/genome3D/astral/pdbstyle-2.03/$folder/$file:   FIRST : $ch : LAST : $ch ::::\n";
print ALI "*\n";
print ALI ">P1;$unique\n";
print ALI "sequence:"."$unique"."::::::::\n";
print "$domain_details{$domain_number}{'alignment'}\n";
#'-' x 5
my @altemp = split('',$domain_details{$domain_number}{'alignment'});
my $altemp = '-' x ($domain_details{$domain_number}{'modstart'}-1);
my @res_nums;
my $res_num = $domain_details{$domain_number}{'start'};
#print "START --> $domain_details{$domain_number}{'start'}";
my $tla = scalar(@altemp)-1;
    for (my $i=0; $i <= $tla; $i++) {
    	if(@{$in_or_out}[$i]){
		if(checkCase($altemp[$i])){
			$altemp .= $altemp[$i];
			#print "$res_num\n";
			push(@res_nums,$res_num);
		}elsif($altemp[$i] eq '-'){
			$altemp .= $altemp[$i];
			push(@res_nums,$res_num);
		}
	}
	$res_num++;
    } 
#print @res_nums;
#my $altemp = uc($domain_details{$domain_number}{'alignment'});
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
print "/tmp/"."$unique".".pdb\n"; 
open(FILEOUT,">/media/scratch/genome3D/pfam_assfiles/models/tests/"."$unique".".pdb") or die $!;
 my $flag = 0;
my $index = 0;
my $previous = 1;
my $pos_in_file;
my $val2print;
 while(<FILE>){
 	if($flag == 1){
 		#REMARK GENOME3D NAME uniprot_start_stop
 		print FILEOUT "REMARK GENOME3D NAME $unique\n";
 		#REMARK GENOME3D UNIPROT_ID
 		print FILEOUT "REMARK GENOME3D UNIPROT_ID $seqid \n";
		#REMARK GENOME3D UNIPROT_MD5
		print FILEOUT "REMARK GENOME3D UNIPROT_MD5 $digest \n";
		#REMARK GENOME3D TIMESTAMP YYYY-MM-DD
 		print FILEOUT "REMARK GENOME3D TIMESTAMP 2012-02-28\n";
 		#REMARK GENOME3D TOTAL TEMPLATES 1
 		print FILEOUT "REMARK GENOME3D TOTAL TEMPLATES 1\n";
 		#REMARK GENOME3D SELECTION X
		#REMARK GENOME3D MODELLING X
		print FILEOUT "REMARK GENOME3D SELECTION SUPERFAMILY\n";
		print FILEOUT "REMARK GENOME3D MODELLING MODELLER\n";
		
		print FILEOUT "REMARK GENOME3D ALIGNMENT_SOURCE SUPERFAMILY \n";
		
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
 	}elsif($_ =~ /ATOM/){
#ATOM      1  N   GLY     1     -19.234   1.412  78.819  1.00 10.27           N
		my $line = $_;
		$pos_in_file = ltrim(substr $line, 22, 4); 		
 		
		unless($pos_in_file == $previous){
			$index++;
		}
		$previous = $pos_in_file;
		$val2print = sprintf("% 4d", $res_nums[$index]);	
		my $lineout = substr $line, 22, 4, $val2print;
		print FILEOUT $line;
		
 		$flag++;
 	}elsif($_ =~ /TER/){
		my $line = $_;
		my $o = substr $line, 22, 4, $val2print;
		print FILEOUT $line; 
	}else{
		print FILEOUT $_;
		$flag++;
	}
 }
 
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".B99990001.pdb");
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".ini");
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".rsr");
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".D00000001");
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".V99990001");
# unlink("$TEMPDIR"."/"."$fo"."/"."$unique".".sch");
# unlink("$TEMPDIR"."/"."$fo"."/"."align"."$unique".".temp");
 
}
}
}




