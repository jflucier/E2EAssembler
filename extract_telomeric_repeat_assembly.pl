#! /usr/bin/perl

use strict;
use warnings;

# die with stack trace
use Carp;
$SIG{__DIE__} = sub { confess $_[0] };

use Data::Dump qw(dump);

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $SEQKIT_IN = $ARGV[0];
my $SEED_LEN = $ARGV[1];

my $MIN_TELO_LEN = 10;

if(!defined($SEED_LEN)){
    print STDERR "Telomeric seed length not defined, will set to default 6 nt\n";
    $SEED_LEN = 6;
}

main();

sub main {

    if(!defined($SEQKIT_IN) or $SEQKIT_IN eq ''){
        die("Input SEQKIT locate report BED file is required!");
    }

    my $all_telo_match = read_seqkit($SEQKIT_IN);

    # print dump($all_telo_match) . "\n";
    # die();

    my $telo_region = flatten_telo_match($all_telo_match);

    # print dump($telo_region->{'tig00000001'}) . "\n";
    # die();


    foreach my $id (keys(%$telo_region)){
        my @sorted = @{$telo_region->{$id}};
        # print dump( @sorted) . "\n";
        foreach my $coords (@sorted){
            # print dump($coords) . "\n";
            my $curr_start = $coords->[0];
            my $curr_end = $coords->[1];

            if(($curr_end-$curr_start) > $SEED_LEN + 1){
                print $id . "\t"
                    . $curr_start . "\t" . $curr_end . "\t+\n";
            }
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
