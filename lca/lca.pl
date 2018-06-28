#!/usr/bin/env perl
# VERSION 0.1.0
#########################################################################################
# This script produce a lowest common ancestor affiliation (LCA) of each sequence       #
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

lca.affiliation.pl

=head1 SYNOPSIS

    perl affiliation.pl [-t file] [-a abundance] [-lca taxonomy] [-env] [-longer] 
 
=head1 DESCRIPTION
    
    This script produce a lowest common ancestor affiliation (LCA) of each sequence and the abundance of each lineage from a tabular file.

=head1 OPTIONS

    The options for the program as as follows:

Global:

    -h,	--help,	    Print this help file and exit.
    
    -v,	--verbose,  Boolean option to print out warnings during execution. Warnings and errors are redirected to STDERR. 
                    Defaults to no verbose (silent mode).

Input:

    -t,	--taxo,         File with complete taxonomy for each contigs.

    -d,  --delimiter,   Delimiter between the sequence name and features id. Must be dot, pipe or underscore (default: pipe)


Ouptut:

    --lca,           Output file with the result of the LCA.

    -a,             Output file with the abundance of each lineage. 

    --env,          Boolean option to not take account of 'unclultured viruses' into LCA process. 

    --longer,       Keep the longest lineage as possible during LCA process.  

    -f,	--force,    Force the script by ERASE early output file.


=head1 AUTHORS

Written by Corentin Hochart (corentin.hochart@uca.fr), UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE). Released under the terms of the GNU General Public License v3. preprocess version v0.1.0.

=head1 VERSION

V0.1.0

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
my $taxofile;
my $abundance="abundance.tsv";
my $LCA="lca.tsv";
my $env;
my $force;
my $longer;
my $DoOrNot;
my $d = 'pipe' ;

#tables
my @affiliation;

#hashes
my %lca ;
my %lcacount;

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
            "taxo|Taxo|t=s" => \$taxofile,
            "lca=s"=> \$LCA,
            "abundance=s"=> \$abundance,
            "env!"=> \$env,
            "force|f!" =>\$force,
            "longer!" =>\$longer,
            "d=s" => \$d,
        );

    &warning("Start of script.");
    
    if ($help) {
        pod2usage(-verbose => 2, -noperldoc => 1);
        exit;
    }

unless($taxofile){&error("No input file")}
else{
    unless(-e $taxofile){&error("$taxofile do not exist")}
}
if(-e $abundance){
    if($force){`rm $abundance`}	
    else{
        &error("$abundance already exists: use force option to erase early output file")
    }
}
if(-e $LCA){
    if($force){`rm $LCA`}	
    else{
        &error("$LCA already exists: use force option to erase early output file")
    }
}

error("'$d' : Wrong delimiter for '-d' option. Must be 'dot' or 'pipe' or 'underscore'") if ($d ne 'pipe' && $d ne 'dot' && $d ne 'underscore');

&warning("Read $taxofile.");
open (TAXO, $taxofile) or die ("Could not open $taxofile:$!");
while (my $line = <TAXO>){
    chomp $line;
    my @line = split(/\t/,$line);
    my @contig;
    @contig=split(/\|/,$line[0]) if $d eq 'pipe';
    @contig=split(/\./,$line[0]) if $d eq 'dot';
    @contig=split(/_/,$line[0]) if $d eq 'underscore';    
    pop(@contig);
    my $contig = '' ;
    $contig = join('|',@contig) if $d eq 'pipe';
    $contig = join('.',@contig) if $d eq 'dot';
    $contig= join('.',@contig) if $d eq 'underscore';        
    if ($line[1]){
        unless(exists $lca{$contig}){$lca{$contig}=$line[1]; $DoOrNot="do"}
        else{
            if ($env){
                if($line[1]=~m/environmental samples/){
                    next;
                }
                else{
                    if($lca{$contig}=~m/environmental samples/){
                        $lca{$contig}="$line[1]";
                        next;
                    }
                }
            }
            my @taxo1 = split(/;/,$lca{$contig});
            my @taxo2 = split(/;/,$line[1]);
            $lca{$contig}="";
            if (scalar(@taxo2)>scalar(@taxo1)){
                if($longer){
                    for (my$i=0;$i<=$#taxo2;$i++){
                        unless($taxo1[$i]){
                            unless($DoOrNot eq "not"){
                                $lca{$contig}=join(";",@taxo2);	
                            }
                        }
                        else{
                            if ($taxo1[$i] eq $taxo2[$i]){
                                $lca{$contig}.="$taxo1[$i];";
                            }
                            else{
                                if($env){
                                    if($taxo1[$i]=~m/unclassified/){
                                        $lca{$contig}=join(";",@taxo2);	
                                    }
                                    elsif($taxo2[$i]=~m/unclassified/){
                                        $lca{$contig}=join(";",@taxo1);
                                    }	
                                }
                                $DoOrNot="not";
                                last;
                            }
                        }
                    }
                }
                else{
                    for (my$i=0;$i<=$#taxo1;$i++){
                        if ($taxo1[$i] eq $taxo2[$i]){
                            $lca{$contig}.="$taxo1[$i];" 				
                        }
                        else{
                            if($env){
                                if($taxo1[$i]=~m/unclassified/){
                                    $lca{$contig}=join(";",@taxo2);	
                                }
                                elsif($taxo2[$i]=~m/unclassified/){
                                    $lca{$contig}=join(";",@taxo1);
                                }
                                last;	
                            }
                        }
                    }
                }
            }
            else {
                if($longer){
                    for (my$i=0;$i<=$#taxo1;$i++){
                        unless($taxo2[$i]){
                            unless($DoOrNot eq "not"){
                                $lca{$contig}=join(";",@taxo1);
                            }
                        }
                        else{
                            if ($taxo1[$i] eq $taxo2[$i]){
                                $lca{$contig}.="$taxo1[$i];" 				
                            }
                            else{
                                if($env){
                                    if($taxo1[$i]=~m/unclassified/){
                                        $lca{$contig}=join(";",@taxo2);	
                                    }
                                    elsif($taxo2[$i]=~m/unclassified/){
                                        $lca{$contig}=join(";",@taxo1);
                                    }	
                                }
                                $DoOrNot="not";
                                last;	
                                
                            }
                        }
                    }
                }
                else{
                    for (my$i=0;$i<=$#taxo2;$i++){
                        if ($taxo1[$i] eq $taxo2[$i]){
                            $lca{$contig}.="$taxo1[$i];" 				
                        }
                        else{
                            if($env){
                                if($taxo1[$i]=~m/unclassified/){
                                    $lca{$contig}=join(";",@taxo2);	
                                }
                                elsif($taxo2[$i]=~m/unclassified/){
                                    $lca{$contig}=join(";",@taxo1);
                                }	
                            }
                            last;
                        }
                    }
                }
            }
        }
    }
}

&warning("Produce lca affiliation.");
foreach my $keys (sort keys %lca){
    my $lca = $lca{$keys};
    open (LCA, ">>$LCA");
    print LCA "$keys\t$lca\n";
    close LCA;
    if (exists $lcacount{$lca}){
        $lcacount{$lca}++;	
    }
    else {
        $lcacount{$lca}=1;
    }
}

&warning("Produce lineage abundance.");
foreach my $keys (sort keys %lcacount){
    my $lca = $keys;
    $lca =~ s/;/\t/g;
    open (ABUN, ">>$abundance");
    print ABUN "$lcacount{$keys}\t$lca\n";
    close ABUN;
}

&warning("End of script.");

}
