#! /usr/bin/perl

use strict;
use warnings;

# die with stack trace
use Carp;
$SIG{__DIE__} = sub { confess $_[0] };

use Data::Dump qw(dump);

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $FA_IN = $ARGV[0];

main();

sub main {

    if(!defined($FA_IN) or $FA_IN eq ''){
        die("Input read fasta is required!");
    }

    print "reading fasta\n";
    my $seqs = read_fasta($FA_IN);

    print "seq nbr = " . scalar(@$seqs) . "\n";
}

sub read_fasta{
    my($FA_IN) = @_;

    open(my $FH, "<$FA_IN") or die("Unable to open $FA_IN");
    my @lines = <$FH>;
    chomp(@lines);
    my $c=0;
    my $s="";
    my $h="";
    my @seqs;
    foreach my $li (@lines){
        if($li =~/^\>/){
            if($seq ne ""){
                push(@seqs,{
                    id => $h,
                    seq => $s
                });
            }
            $h = $li;
            $s = "";
        }
        else{
            $s = $s . $li;
        }
    }

    return \@seqs;
}
