#! /usr/bin/perl

use strict;
use warnings;

# die with stack trace
use Carp;
$SIG{__DIE__} = sub { confess $_[0] };

use Data::Dump qw(dump);
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $BED_IN = $ARGV[0];
main();

sub main {

    if(!defined($BED_IN) or $BED_IN eq ''){
        die("Input Bed annotation file is required!");
    }

    open(my $FH, "<$BED_IN");

    my @lines = <$FH>;
    chomp(@lines);

    my %annot;
    my @all_subs;
    foreach my $l (@lines) {
        my($chr,$start,$end,$name,@rest) = split("\t",$l);
        if(!exists($annot{$chr})){
            $annot{$chr} = {
                name_str => '',
                coord_str => ''
            };
        }

        if(!exists($annot{$chr}->{$name})){
            $annot{$chr}->{$name} = 0;
            push(@all_subs,$name);
        }

        $annot{$chr}->{$name}++;
        $annot{$chr}->{name_str} .= "-" . $name;
        $annot{$chr}->{coord_str} .= "-(" . $start . ":" . $end . ")";
    }

    my @unique_subs = uniq @all_subs;

    # use Data::Dump qw(dump);
    # print STDERR dump(%annot) . "\n";
    # print STDERR dump(sort {$a <=> $b} keys(%annot)) . "\n";
    my $header_print = 0;
    foreach my $chr_name (sort {$a cmp $b} keys(%annot)) {
        my $out_l = $chr_name . "\t";

        my $header_str = "chr\t";
        foreach my $sub (sort {$a cmp $b} @unique_subs){
            if($sub eq "name_str" or $sub eq "coord_str"){
                next;
            }

            if(!$header_print){
                $header_str .= $sub . "_cnt\t";
            }

            my $cnt = 0;
            if(exists($annot{$chr_name}->{$sub})){
                $cnt = $annot{$chr_name}->{$sub};
            }
            $out_l .= $cnt . "\t";
        }

        if(!$header_print){
            $header_str .= "organisation_str\torganisation_coords";
        }
        $out_l .= substr($annot{$chr_name}->{name_str},1) . "\t" . substr($annot{$chr_name}->{coord_str},1);

        if(!$header_print){
            print $header_str . "\n";
            $header_print = 1;
        }
        print $out_l . "\n";
    }

}
