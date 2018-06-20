#!/usr/bin/env perl
# VERSION 0.1.0
#########################################################################################
# This script detect circular contigs by looking for exact identical k-mer at the two   #
# ends on a cadre sequence of the sequences prodvide in fasta file.                     #
#########################################################################################

#########################################################################################
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ##
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A       ##
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  ##
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ##
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ##
## SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                              ##  
#########################################################################################

=head1 NAME

detect.circular.seq.pl

=head1 SYNOPSIS

perl detect.circular.seqpl [--fasta] [--kmer] [--cadre] > output.fasta

=head1 DESCRIPTION

This script detect circular contigs by looking for exact identical k-mer at the two ends on a cadre sequence of the sequences prodvide in fasta file. In order to be able to predict genes spanning the orgin of circular contigs, the first 1,000 nucleotides of each circular contigs are dulicated and added at the contig's end. 

=head1 OPTIONS

--Help|help|h, produces this help file.

--Verbose[no-Verbose]|verbose[no-verbose]|v[no-v], boolean option to print out warnings during execution. Warnings and errors are redirected to STDERR. Defaults to no verbose (silent mode).

--fasta|Fasta|f, sequence fasta file

--kmer|Kmer|k, motif length detected to identify circular sequence (default: 10)

--cadre|Cadre|c, length of the inspect fragment at the sequence 5' end for  kmer identity finding. If 0 screen all the sequence. Default: 0.


=head1 AUTHORS

HOCHART Corentin

Inspired by ROUX Simon work for Metavir2 (2014).

=head1 VERSION

V0.1.0

=head1 DATE

Creation : 15/06/2015
Last modification : 30/08/2016

=cut


# libraries
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;

#scalars
my $help; 	# help flag
my $verbose;	# debugging flag
my $counter = 0; #counter for file name in file verification fonction
my $fastaFile="";
my $kmerl=10;
my $CADRE=0;
my $dupnucleotide=1000;
my $seq;
my $header;

#tables

#hashes

#function
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

MAIN: {
    GetOptions(
            "Help|help|h" => \$help,
            "Verbose|verbose|v!" => \$verbose,
            "fasta|Fasta|f=s" => \$fastaFile,
            "kmer|Kmer|k=s" => \$kmerl,
            "cadre|Cadre|c=s" => \$CADRE,
        );

    &warning("Start of script.");
    
    if ($help) {
        pod2usage(-verbose => 2, -noperldoc => 1);
        exit;
    }

unless (-e $fastaFile){&error("Could not find $fastaFile")}
unless ($kmerl =~ m/\d+/){&error("$kmerl : need an appropriate value for kmer option!")};
unless ($CADRE =~ m/\d+/){&error("$CADRE : need an appropriate value for cadre option!")};



open (FASTA, "$fastaFile") or die ("Could not open $fastaFile : $!");
while (my $line = <FASTA>){
    chomp $line ;
    if ($line =~ m/^>/){chomp $line;$header=$line}
    elsif ($line=~ m/^[ATCGatcgnN]/){
        my @line = split(//,$line);
        my $begin;
        my $end;
        my $cadre=$CADRE;
        if ($CADRE > scalar(@line) - $kmerl){
            $cadre = 0;
        }
        if ($CADRE==0){
            $cadre = scalar(@line) - $kmerl;
        }
        for (my $i=0 ; $i<$kmerl ; $i++){
            $begin.=$line[$i];	
        }
        for (my $i=scalar(@line)-1 ; $i>scalar(@line)-$cadre-1 ; $i-- ){
            $end.=$line[$i];
        }
        $end = reverse($end);
        my @end = split(//,$end);
        my $statue="NOP";
        my @endlimit;
        my $endlimit;
        my $k = 0 ;
        while ($end =~ m/$begin/g){
            my $pos = pos($end);
            $statue = "OK";
            $endlimit[$k] = $pos ;
            $k++;
        }
        if ($statue eq "OK"){
            for (my $j=0 ; $j < scalar(@endlimit) ; $j++){
                my $scale = scalar(@line) - $endlimit[$j] ;
                for (my $i = 0 ; $i < $scale; $i++){
                    my $k = $i + $endlimit[$j];
                    unless ($line[$i] eq $line[$k]){
                        $statue = "NOP";
                    }
                    else {$statue = "OK"}
                }
                $endlimit = scalar(@line) - $endlimit[$j];
                if ($statue eq"OK"){last}
            }
        }
        if ($statue eq "OK"){
            my @seq = @line;
            for (my $k=0 ; $k<$endlimit ; $k++ ){
                pop (@line);
            }
            if (scalar(@line)<$dupnucleotide){
                push @line, @line;
            }
            else {
                my $nuc1000;
                for (my $k=0 ; $k<$dupnucleotide ; $k++){
                    $nuc1000.="$seq[$k]";
                }
                push @line, $nuc1000;
            }
            my $seq = join("",@line);
            print "$header-circ\n";
            print $seq,"\n";
        }
        else {
            my $seq = join("",@line);
            print $header,"\n";
            print $seq,"\n";
        }
    }
}
close FASTA;
&warning("End of script.");
}

