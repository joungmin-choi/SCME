#!/usr/bin/perl -w
# get the position of CpG from genome sequence file

=head1 Description

This script is used to get the position of CpG dyad from
genome file of FASTA format.

=head1 Usage

perl get_string_position_from_genome.pl <genome.fa> <string> >output

=head1 Example

perl get_string_position_from_genome.pl hg38.fa CG >hg38.CpG.pos.txt

=head1 Version

v1.1 Ming-an Sun, 2020-03-23

Last modified: 2020-03-24

=cut

die `pod2text $0` unless @ARGV==2;

my $infile = shift;
my $pattern = shift;

$pattern = uc($pattern);

warn "Processing $infile ...\n";
open(IN,$infile)||die"Cannot read $infile\n";
$/ = ">";
while(<IN>){
	s/>//;
	if(/^(\S+)[^\n]*\n(.*)$/s){
		my ($name, $seq) = ($1, $2);
        warn "$name\n";
		$seq =~ s/\W//g;
		$seq = uc $seq;
	    foreach my $pos (&get_str_pos($seq, $pattern)){
		    print "$name\t$pos\n";
        }
	}
}
$/ = "\n";
close IN;

warn "Done.\n";

##### subroutines #####
#
# read fasta file
sub read_fasta{
	my $faFile = shift;
	my @info;
	open(IN,$faFile)||die"Cannot open $faFile\n";
	$/ = ">";
	while(<IN>){
		s/>//;
		if(/^(\S+)[^\n]*\n(.*)$/s){
			my ($name, $seq) = ($1, $2);
            warn "$name\n";
			$seq =~ s/\W//g;
			$seq = uc $seq;
			push(@info, \{'name'=>$name, 'seq'=>$seq});
		}
	}
	$/ = "\n";
	close IN;
	return @info;
}

# get string position from sequence
sub get_str_pos{
	my ($seq, $str) = @_;
	my $strLen = length($str);
	my @pos;
	my $s = 0;
	my $pos;
	while( ($pos = index($seq, $str, $s) ) >= 0){
		my $pos2 = $pos+1;
		push(@pos, $pos2);
		$s = $pos2;
	}
	return @pos;
}

