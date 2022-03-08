#! /usr/bin/perl

use strict;
use warnings;

# die with stack trace
use Carp;
$SIG{__DIE__} = sub { confess $_[0] };

use Data::Dump qw(dump);
use List::Util qw( min max );

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $NHMMER_IN = $ARGV[0];
main();

sub main {

    if(!defined($NHMMER_IN) or $NHMMER_IN eq ''){
        die("Input NHMMER report tblout file is required!");
    }

    open(my $FH, "<$NHMMER_IN") or die("Unable to open $NHMMER_IN");
    my @lines = <$FH>;
    chomp(@lines);

    my $curr_target = "";
    my $curr_query = "";
    my @curr_hmmfrom = ();
    my @curr_hmmto = ();
    my $curr_alifrom;
    my $curr_alito;
    my $curr_seqlen;
    my $curr_strand;
    my @curr_evalue = ();
    my @curr_score = ();
    my @curr_bias = ();
    foreach my $li (@lines){
        #56b84e33-119a-40d3-863f-2e765a6c27c2	25838	25847	CACACACCA	0	+
        my($target,$taccession,$query,$qaccession,$hmmfrom,$hmmto,$alifrom,$alito,$envfrom,$envto,$seqlen,$strand,$evalue,$score,$bias,$desc) = split("\t",$li);
        # print "new_id=$new_id\n";
        if($target ne $curr_target){
            # new target, print last stretch
            print
                $curr_target . "\t-\t"
                $curr_query . "\t-\t"
                min(@curr_hmmfrom) . "\t"
                max(@curr_hmmto) . "\t"
                $curr_alifrom . "\t"
                $curr_alito . "\t"
                ($curr_alifrom + 1). "\t"
                ($curr_alito + 1) . "\t"
                $curr_seqlen . "\t"
                $curr_strand . "\t"
                max(@curr_evalue) . "\t"
                max(@curr_score) . "\t"
                max(@curr_bias) . "\t-\n";

            $curr_target = $curr_target;
            $curr_query = $curr_query;
            @curr_hmmfrom = ($hmmfrom);
            @curr_hmmto = ($hmmto);
            $curr_alifrom = $alifrom;
            $curr_alito = $alito;
            $curr_seqlen = $seqlen;
            $curr_strand = $strand;
            @curr_evalue = ($evalue);
            @curr_score = ($score);
            @curr_bias = ($bias);

        }
        elsif($strand ne $curr_strand){
            # same target but change strand
            # not the same stretch, print last stretch
            print
                $curr_target . "\t-\t"
                $curr_query . "\t-\t"
                min(@curr_hmmfrom) . "\t"
                max(@curr_hmmto) . "\t"
                $curr_alifrom . "\t"
                $curr_alito . "\t"
                ($curr_alifrom + 1). "\t"
                ($curr_alito + 1) . "\t"
                $curr_seqlen . "\t"
                $curr_strand . "\t"
                max(@curr_evalue) . "\t"
                max(@curr_score) . "\t"
                max(@curr_bias) . "\t-\n";

            $curr_target = $curr_target;
            $curr_query = $curr_query;
            @curr_hmmfrom = ($hmmfrom);
            @curr_hmmto = ($hmmto);
            $curr_alifrom = $alifrom;
            $curr_alito = $alito;
            $curr_seqlen = $seqlen;
            $curr_strand = $strand;
            @curr_evalue = ($evalue);
            @curr_score = ($score);
            @curr_bias = ($bias);
        }
        elsif($strand eq '+' and $alifrom < $curr_alifrom){
            # overlapping stretch
            push(@curr_hmmfrom,$hmmfrom);
            push(@curr_hmmto,$hmmto);
            $curr_alito = $alito
            push(@curr_evalue,$evalue);
            push(@curr_score,$score);
            push(@curr_bias,$bias);
        }
        elsif($strand eq '-' and $alito < $curr_alifrom){
            # overlapping stretch
            push(@curr_hmmfrom,$hmmfrom);
            push(@curr_hmmto,$hmmto);
            $curr_alito = $alito
            push(@curr_evalue,$evalue);
            push(@curr_score,$score);
            push(@curr_bias,$bias);
        }
    }


}

sub flatten_telo_match{
    my($all_telo_match) = @_;

    my %struct;

    foreach my $id (keys(%$all_telo_match)){
        my @sorted = @{$all_telo_match->{$id}};
        # my @sorted = reverse(sort { $a <=> $b } @coords);

        my $is_ongoing_stretch = 0;
        my $stretch_start = 0;
        my $stretch_end = 0;
        $struct{$id} = [];

        for(my $i=0; $i < scalar(@sorted); $i++){
            my $curr_start = $sorted[$i];
            my $curr_end = $curr_start + $SEED_LEN + 1;

            if($i == 0){
                # $is_ongoing_stretch = 1;
                $stretch_start = $curr_start;
                $stretch_end = $curr_end;
            }

            # check if this is last pos in array
            if($i == scalar(@sorted) - 1){
                # last item
                # if($is_ongoing_stretch){
                    push(@{$struct{$id}}, [$stretch_start,$curr_end]);
                # }
                last;
            }

            my $next_start = $sorted[$i+1];
            my $next_end = $sorted[$i+1] + $SEED_LEN + 1;

            if($next_start > $curr_end + $SEED_LEN + 1){
                # end of stretch
                # $is_ongoing_stretch = 0;
                push(@{$struct{$id}}, [$stretch_start,$curr_end]);
                $stretch_start = $next_start;
                $next_end = $curr_end;
            }
            else{
                # steetch keep going
                # $is_ongoing_stretch = 1;
                $stretch_end = $curr_end;
            }
        }

    }

    return \%struct
}

sub read_seqkit{
    my($SEQKIT_IN) = @_;

    open(my $FH, "<$SEQKIT_IN") or die("Unable to open $SEQKIT_IN");
    my @lines = <$FH>;
    chomp(@lines);

    my $curr_id = "";
    my $curr_start;
    my $curr_end;
    my %struct;
    foreach my $li (@lines){
        #56b84e33-119a-40d3-863f-2e765a6c27c2	25838	25847	CACACACCA	0	+
        my($new_id,$new_start,$new_end,$new_motif,$new_x,$new_strand) = split("\t",$li);
        # print "new_id=$new_id\n";
        if($new_id ne $curr_id){
            # init entry
            $curr_id = $new_id;
            $struct{$curr_id} = [];
        }

        push(@{$struct{$curr_id}},$new_start);
        $curr_start = $new_start;
        $curr_end = $new_end;
    }

    # sort coords
    foreach my $k (keys(%struct)){
        my @vals = @{$struct{$k}};
        my @s = sort { $a <=> $b } @vals;
        $struct{$k} = \@s;
    }
    return \%struct;
}
