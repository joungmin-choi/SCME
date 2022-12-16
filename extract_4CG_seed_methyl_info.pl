#!/usr/bin/perl -w
# extract fragments (of 4 CpG dyads) with at least one unmethylated (C_C_C_C) and one methylated (mC_mC_mC_mC) continuous CpG dyads
# Ming-an Sun, May 21, 2013
# Ming-an Sun, Oct 02, 2013: add a column to show the number of full reads (4CG with no missing data)
# Ming-an Sun, Oct 09, 2013: read one line at each time to reduce memory requirement
# Ming-an Sun, Apr 25, 2014: remove useless columns
use warnings;
use strict;

=head1 Description

This script is used to extract 4CG seed information from merged CpG results. The input should be 
sorted by chrs and coordinates already. The output of merge_CpG_methyl_info.pl can be used
directly).

=head1 Usage

perl extract_4CG_seed_methyl_info.pl <CpG coord file> <CpG methyl info file> > output

=head1 Example

perl extract_4CG_seed_methyl_info.pl hg38.CPG brain.CpG.inf.txt >brain.CpG.4CG.txt

=head1 Version

Ver 1.0, Ming-an Sun, May 19, 2014

Ver 1.1, Ming-an Sun, Mar 3, 2015 (support pair-end results now)

Ver 1.2, Ming-an Sun, 2020-03-24

Last modified: 2020-03-24

=cut

die `pod2text $0` unless @ARGV==2;

my $coord_file = shift;
my $siteInfoFile = shift;

my $frag_C_num = 4;
my @info;

warn "Reading CpG coordinate from $coord_file ...\n";
my %coord_arr;
my %coord_idx;
open(COORD, $coord_file)||die"Cannot open $coord_file\n";
while(<COORD>){
    warn $. . "\n" if $. % 1000000 == 0;
    if(/^(\S+)\s+(\d+)/){
        my ($chr, $pos) = ($1, $2);
        push(@{$coord_arr{$chr}}, $pos);
    }
}
close COORD;

warn "Sorting CpG coordinates for each chromosome ...\n";
foreach my $chr (sort keys %coord_arr){
    warn "  $chr\n";
    @{$coord_arr{$chr}} = sort {$a <=> $b} @{$coord_arr{$chr}};
    for(my $i=0; $i<@{$coord_arr{$chr}}; $i++){
        $coord_idx{$chr}->{$coord_arr{$chr}->[$i]} = $i;
    }

}

warn "Parse $siteInfoFile ...\n";
open(SITE, $siteInfoFile)||die"Cannot open $siteInfoFile\n";
my $l = <SITE>;
push(@info, $l);

# print header
print "#chr\tposition\tcount_all\tcount_full\tcount_1111\tcount_0000\tmeth_pattern\n";
while(my $ln = <SITE>){
    warn $. . "\n" if $. % 100000==0;

    if(&getChr($ln) eq &getChr($info[-1])){
        push(@info, $ln);
        
        my $chr = &getChr($info[0]);
        # array not full (size=4), continue
        unless(scalar(@info) == 4){
            next;
        }
        
        # the 4CGs in the array are not continuous, shift the 1st item, then continue
        my $pos0    = &getPos($info[0], $coord_idx{$chr});
        my $pos1    = &getPos($info[1], $coord_idx{$chr});
        my $pos2    = &getPos($info[2], $coord_idx{$chr});
        my $pos3    = &getPos($info[3], $coord_idx{$chr});
        my $pos1exp = $coord_arr{$chr}->[$coord_idx{$chr}->{$pos0}+1];
        my $pos2exp = $coord_arr{$chr}->[$coord_idx{$chr}->{$pos0}+2];
        my $pos3exp = $coord_arr{$chr}->[$coord_idx{$chr}->{$pos0}+3];
        
        unless($pos1 == $pos1exp && $pos2 == $pos2exp && $pos3 == $pos3exp){
            shift @info;
            next;
        }
        
        # get 4CG seed information
        my %pattern;
        my $cgPos = '';
        for(my $j=1; $j<=$frag_C_num; $j++){
            my $idx = $j-1;
            my ($thisChr, $pos, $posStr, $negStr) = split(/\s+/, $info[$idx]);
            $chr = $thisChr;
            $cgPos .= "$pos;";
            my @posArr = split(';', $posStr);
            my @negArr = split(';', $negStr);
            foreach my $readName (@posArr){
                next if $readName eq '.';
                $pattern{$readName}  = defined($pattern{$readName}) ? &fillPattern($pattern{$readName}, $j-1) : &fillPattern("", $j-1);
                $pattern{$readName} .= 1;
            }
            foreach my $readName (@negArr){
                next if $readName eq '.';
                $pattern{$readName}  = defined($pattern{$readName}) ? &fillPattern($pattern{$readName}, $j-1) : &fillPattern("", $j-1);
                $pattern{$readName} .= 0;
            }
        }
        shift(@info);
    
        my $all_num      = 0;
        my $all_full_num = 0;
        my $all_umC_num  = 0;
        my $all_mC_num   = 0;
        foreach my $readName (keys %pattern){
            $pattern{$readName} = &fillPattern($pattern{$readName}, $frag_C_num);
            $all_num++;
            $all_full_num++ if $pattern{$readName} =~ /^[01]+$/;
            $all_umC_num++  if $pattern{$readName} =~ /^0+$/;
            $all_mC_num++   if $pattern{$readName} =~ /^1+$/;
        }
    
        next if $all_full_num == 0;
        print "$chr\t$cgPos\t$all_num\t$all_full_num\t$all_mC_num\t$all_umC_num\t";
        foreach my $readName (sort { $pattern{$a} cmp $pattern{$b} } keys %pattern){
            print "$readName\:$pattern{$readName};";
        }
	#ADD
	print "\n"
    }
    else{
        @info = ($ln);
    }
        
}

close SITE;


warn "Done.\n";

######################################################
################ subroutines #########################
# fill missing data with ?
sub fillPattern{
    my ($str, $len) = @_;
    if(length($str) < $len){
        $str .= '?' x ($len-length($str));
    }
    return $str;
}

sub getNum{
    my $str = shift;
    my $num = 0;
    while($str =~ /[01]/g){
        $num++;
    }
    return $num;
}

sub getChr{
    my $str = shift;
    my $chr = '';
    if($str =~ /^(\S+)/){
        $chr = $1;
    }
    return $chr;
}

sub getPos{
    my ($str, $hashP) = @_;
    my $newPos;
    if($str =~ /^(\S+)\s+(\d+)/){
        my ($chr, $pos) = ($1, $2);
        if(defined $hashP->{$pos}){
            $newPos = $pos;
        }
        elsif(defined $hashP->{$pos-1}){
            $newPos = $pos-1;
        }
        else{
            $newPos = '';
        }
    }
    return $newPos;
}
__END__
chr22    16050703    9    4    0.444444444444444    SRR478986.26503380;SRR478988.11084554;SRR478988.32330423;SRR478989.33851660    SRR478983.166545;SRR478983.166762;SRR478986.72053982;SRR478987.12565932;SRR478989.52275411;
chr22    16051245    9    7    0.777777777777778    SRR478987.27708289;SRR478988.22934565;SRR478988.34714815;SRR478988.38921509;SRR478988.66032570;SRR478989.1780925;SRR478989.12319698    SRR478987.18231494;SRR478989.62561892
chr22    16051936    2    2    1    SRR478987.7180754;SRR478987.52363258    .
chr22    16051963    1    1    1    SRR478987.52363258    .
chr22    16052079    7    2    0.285714285714286    SRR478988.58477654;SRR478988.63843962    SRR478986.78820689;SRR478988.178318;SRR478989.11906916;SRR478989.57457115;SRR478989.67383192
