#!/usr/bin/perl -w
# convert seed file to the format for Prof. Wu's script for CSM inference
# Ming-an Sun, Jan 23, 2014
use strict;

=head1 Description

This script is used to convert 4CG seed file format. The input is the file generated using
extract_4CG_seed_methyl_info.pl, which contains rich information such as read names and 
missing info (such as ??11). The output only keeps patterns with full information,
and also includes the calculated ML and ME for each 4CG seed.

=head1 Usage

perl convert_4CG_seed_format.pl <detailed seed file>  > new_seed_file

=head1 Example

perl convert_4CG_seed_format.pl brain.CpG.4CG.txt >brain.CpG.4CG.cnv.txt

=head1 Version

v1.0, Ming-an Sun, Aug 12, 2014

v1.1, Ming-an Sun, 2020-03-24: add header line

Last modified: 2020-03-24

=cut

die `pod2text $0` unless @ARGV;


my $seed_file = shift;

my $entropyBit = log(10)/log(2)/4; # so the value of ME is between 0 and 1. But for weighted_WE, it can never be so high because it is weighted.

my %penalty_4CG = (
    '0000'=>0, '1111'=>0,
    '0111'=>1/3, '1000'=>1/3, '0001'=>1/3, '1110'=>1/3, '0011'=>1/3, '1100'=>1/3,
    '0100'=>2/3, '0010'=>2/3, '1001'=>2/3, '0110'=>2/3, '1011'=>2/3, '1101'=>2/3,
    '0101'=>1, '1010'=>1
);
                    
#IN: 
#chr10	3107671;3107723;3107731;3107754;	5	5	0	5	SRR921832.69500918:1111;SRR921841.28066836:1111;SRR921835.57540705:1111;SRR921836.96218669:1111;SRR921841.5413196:1111;
#OUT:
#chr19	89762;89781;89784;89790;	11	5	0	0	0.821917808219178	0.492936814282045	0011:1;1011:1;1101:3;1110:1;1111:5;?011:1;?100:2;?101:1;?111:5;??11:1;
warn "Processing $seed_file ...\n";
print "Chrom\tPosition\tCount_full\tCount_1111\tCount_0000\tBipolar_status\tML\tME\tPattern\n";
open(SEED, $seed_file)||die"Cannot open $seed_file\n";
while(<SEED>){
    warn "$.\n" if $. % 1000000 == 0;
    next if /^\#/;
    my @a = split;
    my $chr = $a[0];
    my $pos = $a[1];
    my $count_full = 0;
    my $count_0000 = 0;
    my $count_1111 = 0;
    my %pattern;
    foreach my $x (split(';', $a[-1])){
        if($x =~ /\:([01]{4})$/){
            my $p = $1;
            $pattern{$p}++;
            $count_full++;
            if($p eq '0000'){
                $count_0000 ++;
            }elsif($p eq '1111'){
                $count_1111 ++;
            }else{
                ;
            }
        }
    }
    my $bipolar_status = $count_0000>0 && $count_1111>0 ? 1 : 0;
    my $ml = &getML(\%pattern);
    my $me = &getME(\%pattern);
    #my $weighted_me = &get_trans_penalty_ME(\%pattern);
    print "$chr\t$pos\t$count_full\t$count_1111\t$count_0000\t$bipolar_status\t$ml\t$me\t";
    foreach my $p (sort keys %pattern){
        print "$p\:$pattern{$p};";
    }
    print "\n";
    
 }
 close SEED;

warn "Done.\n";

############ subroutines  #####################

## get ML
sub getML{
	my $hashP = shift;
    
	my $m_C_count = 0;
    my $all_C_count = 0;
    foreach my $key (keys %{$hashP}){
        while($key =~ /1/g){
            $m_C_count += ${$hashP}{$key};
        }
        $all_C_count += length($key) * ${$hashP}{$key};
    }
    
    my $ml = $all_C_count > 0 ? $m_C_count/$all_C_count : 'NA';
   
    return $ml;
}


## get ME
sub getME{
	my $hashP = shift;
    
	my $me = 0;
    my $depth = 0;
    foreach my $patt (keys %{$hashP}){
	    $depth += ${$hashP}{$patt};
    }
   
    foreach my $patt (keys %{$hashP}){
	    $me -= (${$hashP}{$patt}/$depth)*log10(${$hashP}{$patt}/$depth);
	}
	$me *= $entropyBit;
    return $me;
}

## get weighted ME based on transition numbers
sub get_trans_penalty_ME{
	my $hashP = shift;
    
	my $me = 0;
    my $depth = 0;
    foreach my $patt (keys %{$hashP}){
	    $depth += ${$hashP}{$patt};
    }
   
    foreach my $patt (keys %{$hashP}){
	    $me -= (${$hashP}{$patt}/$depth)*log10(${$hashP}{$patt}/$depth)*$penalty_4CG{$patt};
	}
	$me *= $entropyBit;
    return $me;
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

__END__


==> Input: ../01.seed_merging/6wk.seed <==
chr10	4	3107671;3107723;3107731;3107754;	5	5	0	5	SRR921832.69500918:1111;SRR921841.28066836:1111;SRR921835.57540705:1111;SRR921836.96218669:1111;SRR921841.5413196:1111;
chr10	4	3112573;3112593;3112597;3112652;	3	3	0	2	SRR921833.45341529:0010;SRR921841.2543280:1111;SRR921841.2095162:1111;
chr10	4	3112929;3113001;3113005;3113020;	7	3	0	2	SRR921840.79143805:1110;SRR921833.32009299:1111;SRR921832.85688663:1111;SRR921835.59726952:?111;SRR921839.98256580:?111;SRR921833.32544226:?111;SRR921838.2179321:?111;
chr10	4	3113001;3113005;3113020;3113081;	7	4	0	4	SRR921840.79143805:110?;SRR921835.59726952:1111;SRR921839.98256580:1111;SRR921833.32544226:1111;SRR921838.2179321:1111;SRR921833.32009299:111?;SRR921832.85688663:111?;
chr10	4	3118973;3118993;3119040;3119042;	3	3	0	1	SRR921835.51821600:0011;SRR921835.91851117:0111;SRR921832.16924448:1111;
chr10	4	3139694;3139740;3139770;3139776;	12	2	0	2	SRR921841.11551426:1111;SRR921836.33553237:1111;SRR921833.36697159:?011;SRR921832.18030911:?011;SRR921839.53168596:?111;SRR921836.1464399:?111;SRR921835.56719469:?111;SRR921836.91747149:?111;SRR921832.98356137:?111;SRR921840.24940430:??11;SRR921832.94681478:??11;SRR921840.39546009:??11;
chr10	4	3139740;3139770;3139776;3139782;	12	8	0	6	SRR921833.36697159:0111;SRR921832.18030911:0111;SRR921839.53168596:1111;SRR921836.1464399:1111;SRR921835.56719469:1111;SRR921836.33553237:1111;SRR921836.91747149:1111;SRR921832.98356137:1111;SRR921841.11551426:111?;SRR921840.24940430:?111;SRR921832.94681478:?111;SRR921840.39546009:?111;
chr10	4	3139770;3139776;3139782;3139854;	12	3	0	0	SRR921840.24940430:1110;SRR921832.94681478:1110;SRR921840.39546009:1110;SRR921839.53168596:111?;SRR921836.1464399:111?;SRR921835.56719469:111?;SRR921836.33553237:111?;SRR921833.36697159:111?;SRR921836.91747149:111?;SRR921832.98356137:111?;SRR921832.18030911:111?;SRR921841.11551426:11??;
chr10	4	3146188;3146228;3146236;3146250;	4	4	0	1	SRR921840.66394888:0110;SRR921839.17642848:0111;SRR921832.16212256:0111;SRR921835.13340666:1111;
chr10	4	3146441;3146453;3146501;3146525;	1	1	0	0	SRR921838.1971600:1101;

==> Output: infile.txt <==
chr19	89762;89781;89784;89790;	11	5	0	0	0.821917808219178	0.492936814282045	0011:1;1011:1;1101:3;1110:1;1111:5;?011:1;?100:2;?101:1;?111:5;??11:1;
chr19	93344;93350;93372;93377;	10	8	0	0	0.91044776119403	0.230482023721841	1100:1;1110:1;1111:8;?111:1;??00:1;??10:1;??11:9;???1:2;
chr19	93350;93372;93377;93388;	11	5	0	0	0.841463414634146	0.419184257513033	1001:1;1101:1;1110:4;1111:5;?000:1;?100:1;?111:9;??11:2;???0:1;
chr19	93372;93377;93388;93395;	22	14	1	1	0.846938775510204	0.418233257968342	0000:1;0011:1;1000:1;1011:1;1101:4;1111:14;?111:2;??01:1;???1:2;
chr19	93377;93388;93395;93402;	24	15	2	1	0.849056603773585	0.410788969616458	0000:2;0111:2;1011:4;1110:1;1111:15;?011:1;??11:2;???1:3;
chr19	93388;93395;93402;93411;	21	11	0	0	0.871559633027523	0.459229495225082	0001:2;0111:5;1101:1;1110:2;1111:11;111?:4;?111:2;??11:3;???1:1;
chr19	93395;93402;93411;93419;	19	12	0	0	0.91588785046729	0.416987799299661	0011:2;1011:1;1101:2;1110:2;1111:12;111?:4;11??:4;?111:3;??11:1;
chr19	93402;93411;93419;93427;	22	14	0	0	0.902912621359223	0.409654990898302	0111:3;1011:2;1100:2;1110:1;1111:14;11??:4;1???:4;?111:1;
chr19	93411;93419;93427;93433;	22	17	0	0	0.91578947368421	0.302507128508691	0111:2;1000:1;1001:1;1101:1;1111:17;111?:1;1???:4;
chr19	93419;93427;93433;93449;	18	11	0	0	0.852272727272727	0.352711968701591	0001:1;001?:1;1010:1;1110:5;1111:11;111?:3;11??:1;???0:1;???1:1;
