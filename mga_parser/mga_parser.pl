#!/usr/bin/perl
# VERSION 0.1.0
########################################################################################
# MetaGeneAnnotator wrapper. Produce gff, tsv, nt and/or aa output                     #
########################################################################################

#########################################################################################
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ##
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A       ##
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  ##
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ##
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ##
## SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                              ##
#########################################################################################

=head1 NAME

    mga_parser.pl

=head1 SYNOPSIS

    perl mga_parser.pl [-i file.fasta] [-aa file.faa] [-nt file.fna] [-tsv file.tsv] [-gff file.gff]

=head1 DESCRIPTION

    This script automatically predict gene/protein from a fasta file withe MetaGeneAnnotator predictor tool
    (Hideki Noguchi, Takeaki Taniguchi and Takehiko Itoh, DNA Research, 2008). 

=head1 OPTIONS

Global : 

	--Help|help|h, produces this help file.    

	--Verbose[no-Verbose]|verbose[no-verbose]|v[no-v]
	boolean option to print out warnings during execution. 
	Warnings and errors are redirected to STDERR.
	Defaults to no verbose (silent mode).

	-f, --force, force the script by ERASED early project. 

Mandatory :

	-i, --input, nucleotid fasta file as input.

Optional :

	--procedure [m|s]	
		Multi or Single species option in MetageneAnnotator. 
		Meta or single parameter for -p parameter in prodigal (default 'm').

	--mga, call existing mga file for nucleotid fasta file (-i).

	--nt, output predicted gene in fasta format.

	--aa, output predicted protein in fasta format.

	--gff, output predicted gene in gff format.

	--tsv, output predicted protein and gene in tabular format. 


=head1 AUTHORS

Written by Corentin Hochart (corentin.hochart@uca.fr)
UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE).
Released under the terms of the GNU General Public License v3.

=head1 VERSION

version 0.1.0

=cut

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# libraries
use warnings;
use strict;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Perl;
use FindBin;
use Getopt::Long;
use Pod::Usage;
use POSIX;


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#scalars
my $help; 		# help flag
my $verbose;	# debugging flag
my $exe = $FindBin::RealScript; # executable path
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>3); # gff version to use 
my $version; # version flag
my $VERSION="0.1.0"; # script version
my $force; # force flag
my $date=date(); # date 

my $procedure="m";
my $GFF;
my $amin;
my $nucleic;
my $gff="";
my $tsv;
my $counter = 0; #counter for file name in file verification fonction
my $input_fasta;

my $VIEWDATASQUARE = 3 ;
my $sid;	#sequence id
my $gc_count;
my $snb=1; 	#sequence number
my $gnb=1 ;	#gene number

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# hashes
my %seq;

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# functions
sub error {
	# management of error messages and help page layout, will stop execution
	# local arguments passed:1st, error message to output
	my $error = shift;
	my $filename = ($0);
	pod2usage(-message => "$filename (error): $error. Execution halted.", -verbose => 1, -noperldoc => 1);
	exit(2);
}

sub warning {
	# management of warnings and execution carry on
	# local arguments passed: 1st, warning message to output
	if ($verbose) {
		my $message = shift;
		my $filename = $0;
		warn("$filename (info): ".$message."\n");
	}
}

sub date { 
  my $time = shift || time;    #$time par defaut vaut le time actuel 
  my ( $seconde, $minute, $heure, $jour, $mois, $annee, $jour_semaine, $jour_annee, $heure_hiver_ou_ete ) 
    = localtime($time); 
  $mois  += 1; 
  $annee += 1900; 
  
  # On rajoute 0 si le chiffre est compris entre 1 et 9 
  foreach ( $seconde, $minute, $heure, $jour, $mois, $annee ) { 
    s/^(\d)$/0$1/; 
  } 
  
  my %ladate = ( 
    "date"         => "$jour-$mois-$annee", 
    "heure"        => "$heure:$minute:$seconde", 
    "jour_semaine" => $jour_semaine, 
    "jour_annee"   => $jour_annee, 
    "hiverOuEte"   => $heure_hiver_ou_ete, 
  ); 
  return \%ladate; 
}

sub fasta_format {
	my $seq=shift;
	$seq=~s/(.{60})/$1\n/g;
	chomp $seq;
    return $seq;
}

sub fileVerification{
	my $file = $_[0] ;
	my $name = "INPUT".$counter ;
	if ($file){
		no strict "refs";
		open ($name,$file) or &error("Unable to read '$file' from '--$_[1]' option !")
	}
	$counter ++ ;
	return
}

sub find_exe {
	my $tool = shift;
	for my $dir (File::Spec->path) {
		my $exe = File::Spec->catfile($dir, $tool);
		return $exe if -x $exe; 
	}
	return;
}

sub version {
  print STDERR "$exe $VERSION\n";
  exit;
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

MAIN: {
	GetOptions(
		"Help|help|h" 			=> \$help,		#help flag
		"Verbose|verbose|v!" 	=> \$verbose,	#verbose flag
		"version!"				=> \$version,	#version flag
		"f|force!"				=> \$force,		#force flag
		
		"i|input=s"	 			=> \$input_fasta,	#input file in fasta format.
		"procedure=s"			=> \$procedure,		#prediction procedure for mga.
		
		"aa=s" 					=> \$amin,		#output for predicted proteins in fasta format.
		"nt=s" 					=> \$nucleic,	#output for predicted genes in fasta format.
		"gff=s" 				=> \$gff,		#output for predicted genes in gff format.
		"tsv=s" 				=> \$tsv,
	);
	
	&version if ($version);
	
	if ($help) {
		pod2usage(-verbose => 2, -noperldoc => 1);
		exit;
	}
	my $user = $ENV{ USER };
	&warning("Hi $user! Let's do some good jobs together.");
	&warning("First. Check tools and options.");
	
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
# check mga tools 
my $fp = find_exe('mga');
&error("Can't find 'mga' in your \$PATH") if !$fp;

# check files and options

unless ($input_fasta){&error("No contigs file")}
&fileVerification ("$input_fasta","input_fasta");

&error("Wrong parameter for '--procedure' option.") unless $procedure eq "m" || $procedure eq "s";  

if($amin){if (-e $amin){
	unless($force){&error("$amin already exists. Use '--force' option to overwrite.")}
	`rm $amin`
}}
if($nucleic){if (-e $nucleic){
	unless($force){&error("$nucleic already exists. Use '--force' option to overwrite.")}
	`rm $nucleic`
}}
if($gff){if (-e $gff){
	unless($force){&error("$gff already exists. Use '--force' option to overwrite.")}
	else{`rm $gff`}
}}
if($tsv){if (-e $tsv){
	unless($force){&error("$tsv already exists. Use '--force' option to overwrite.")}
	else{`rm $tsv`}
}}


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
&warning("Everything seems to be good.");
&warning("Start of script.");

&warning("Load fasta file.");	
my $seqio_object = Bio::SeqIO->new(  -format => 'Fasta',
									-file   => $input_fasta,
						);


my $cmd="mga $input_fasta -$procedure";
&warning("Running: '$cmd'");
open MetaGene, "$cmd|";
while (<MetaGene>){
	chomp;
	if(m/^# gc = (\S+),/){$gc_count=$1}
	elsif(m/^# self/){next}
	elsif(m/^# (\S+)/){$sid = $1}
	else {
		(my $gene_id,my $start,my $end,my $strand,my $frame,my $partial,my $score,my $model,my $rbs_s,my $rbs_e,my $rbs_sc)=split(/\t/);
		push @{$seq{$sid}{FEATURE}}, Bio::SeqFeature::Generic->new( 
			-primary    => "CDS", 
			-seq_id     => $sid,
			-source     => "mga",
			-start      => $start,
			-end        => $end,
			-strand     => $strand,
			-score      => $score,
			-frame      => $frame,
			-tag        => {
				'partial' 	=> $partial,
			},
		);
	}
}
close MetaGene;

&warning("Print output file.");
if ($tsv){
	open (CSV, ">$tsv");
	print CSV "Sequence ID\tStart\tStop\tStrand\tna_sequence\taa_sequence\n";
	close CSV;
}

if($gff){
	open GFF, ">",$gff;
	select GFF;
}
print "##gff-version 3
##source: metagene:$VERSION
##date $date->{date}
##Type DNA\n";
print join ("\t",'##seq_name','method','feature','start','end','score','strand','frame','gene',"\n");
		
while(my $seq_object = $seqio_object->next_seq() ){
	my $sseq=$seq_object->seq;
	my $sid=$seq_object->id;
	$gnb=1;
	for my $feature ( sort { $a->start <=> $b->start } @{ $seq{$seq_object->id}{FEATURE} }) {
		my $start=$feature->start;
		my $end=$feature->end;
		my $frame=$feature->frame;
		my $strand=$feature->strand;
		my $ID = sprintf("%d\_%d",$snb,$gnb);
		$feature->add_tag_value('ID', $ID);
		if($strand eq "-"){
			my $length=$seq_object->length;
			$start=$length-$start+1;
			$end=$length-$end+1;
			($start,$end)=($end,$start);		
		}
		my $CDS=$seq_object->subseq(
			$start,
			$end
		);
		$CDS = Bio::Seq->new(
			'-seq' => $CDS,
		);	
		my $header=$feature->gff_string($gff_factory);
		$header=~m/(\S+\t){8}(\S+)/;
		$header=sprintf(">%s.%d # %d # %d # %d # %s\n",$sid,$gnb,$start,$end,$strand,$2);
		if($strand eq "-1"){
			if($nucleic){
				open (NT, ">>$nucleic");
				print NT $header;
				print NT &fasta_format($CDS->revcom->seq),"\n";
				close NT;
			}
			if($amin){
				open (AA, ">>$amin");
				print AA $header;
				print AA &fasta_format($CDS->revcom->translate(
							'-unknown' => 'X',
							'-frame' => $frame)->seq),"\n";
				close AA;
			}
			if ($tsv){
				open (CSV, ">>$tsv");
				print CSV "$sid.$gnb\t$start\t$end\t$strand\t";
				print CSV $CDS->revcom->seq, "\t";
				print CSV ($CDS->revcom->translate(
							'-unknown' => 'X',
							'-frame' => $frame)->seq),"\n";
				close CSV;
			}
		}
		else{
			if($nucleic){
				open (NT, ">>$nucleic");
				print NT $header;
				print NT &fasta_format($CDS->seq),"\n";
				close NT;
			}
			if($amin){
				open (AA, ">>$amin");
				print AA $header;
				print AA &fasta_format($CDS->translate(
							'-unknown' => 'X',
							'-frame' => $frame)->seq),"\n";
				close AA;
			}	
			if ($tsv){
				open (CSV, ">>$tsv");
				print CSV "$sid\_$gnb\t$start\t$end\t$strand\t";
				print CSV $CDS->seq, "\t";
				print CSV ($CDS->translate(
							'-unknown' => 'X',
							'-frame' => $frame)->seq),"\n";
				close CSV;
			}
		}
		if($gff){
			open GFF, ">>",$gff;
			select GFF;
		}
		print $feature->gff_string($gff_factory),"\n";
		$gnb++;
	}
	$snb++;
}

&warning("End of script.");
}

