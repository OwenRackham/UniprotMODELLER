use File::Slurp;
use lib '/home/rackham/modules';
use rackham;
my $dir = "/home/luca/rackham/astral/alignments/up/";
my @files = read_dir $dir;
my %lu;
$dbh = rackham::DBConnect('superfamily');
$sth =   $dbh->prepare( "select seqid,protein from protein where genome = 'up';" );
        $sth->execute;
        while (my ($seq,$prot)= $sth->fetchrow_array ) {
			$lu{$prot} = $seq;
		}


for my $file ( @files ) {
	my $orig = $file;
 if($file =~ /\.pdb/){
 	unless($file =~ /,/){
 	if($file =~ /^(\S+)(_\S+_\S+)/){
 		if(exists($lu{$1})){
 		my $new = "$lu{$1}"."$2";
 		my $cmd = "mv "."$dir"."$file "."$dir"."$new";
 		system($cmd);
 		
 		}
 	}
 	}
 }
}