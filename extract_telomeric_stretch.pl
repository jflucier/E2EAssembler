#! /usr/bin/perl

use strict;
use warnings;

# die with stack trace
use Carp;
$SIG{__DIE__} = sub { confess $_[0] };

use Data::Dump qw(dump);

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $EXTR = $ARGV[0];
my $SEQKIT_IN = $ARGV[1];
my $SEED_LEN = $ARGV[2];
# using reads with TELOTAG found within 200 nt. Telomeric sequence must be found within 200 nt of the telotag
my $STRETCH_START_WINDOW = 200;

if(!defined($SEED_LEN)){
    print STDERR "Telomeric seed length not defined, will set to default 6 nt\n";
    $SEED_LEN = 6;
}

if(!defined($EXTR)){
    print STDERR "You must define a sequence extremity to analyse: 5 or 3\n";
    print STDERR "Example call: extract_telomeric_stretch.pl <<extremity>> <<seqkit bed>> <<seed lenght>>\n";
    die();
}


main();

sub main {

    if(!defined($SEQKIT_IN) or $SEQKIT_IN eq ''){
        die("Input SEQKIT locate report BED file is required!");
    }

    print "reading BED\n";
    my $all_telo_match = read_seqkit($SEQKIT_IN);

    my $telo_region = flatten_telo_match($all_telo_match);

    my $longest_telo = filter_longest_telostretch($telo_region);

    foreach my $id (keys(%$longest_telo)){
        my($nid) = $id =~ /^(.*)\_\d+$/;
        if(!defined($nid)){
            die("Problem extracting id without read sequence length. Problematic id is $id");
        }
        print $EXTR . "\t" . $nid . "\t"
            . $longest_telo->{$id}->{start} . "\t"
            . $longest_telo->{$id}->{end} . "\t"
            . $longest_telo->{$id}->{len} . "\n";
    }
}

sub filter_longest_telostretch {
    my($telo_region) = @_;

    my %struct;

    foreach my $id (keys(%$telo_region)){
        my $stretches = $telo_region->{$id};

        foreach my $st (@$stretches){
            my $len = $st->[1] - $st->[0];

            if(!exists($struct{$id})){
                $struct{$id} = {
                    start => $st->[0],
                    end => $st->[1],
                    len => $st->[1] - $st->[0]
                };
            }
            elsif($struct{$id}->{len} < $len){
                $struct{$id} = {
                    start => $st->[0],
                    end => $st->[1],
                    len => $st->[1] - $st->[0]
                };
            }
        }
    }
    return \%struct;
}

sub flatten_telo_match{
    my($all_telo_match) = @_;

    my $telo_region;
    if($EXTR == 3){
        $telo_region = flatten_telo_match3($all_telo_match);
    }
    else{
        $telo_region = flatten_telo_match5($all_telo_match);
    }
    return $telo_region;
}

sub flatten_telo_match3{
    my($all_telo_match) = @_;

    my %struct;

    foreach my $id (keys(%$all_telo_match)){
        my @coords = @{$all_telo_match->{$id}};
        my @sorted = reverse(sort { $a <=> $b } @coords);

        my $is_ongoing_stretch = 0;
        my $stretch_start = 0;
        my $stretch_end = 0;
        $struct{$id} = [];
        my($len_seq) = $id =~ /^.*\_(\d+)$/;
        if(!defined($len_seq)){
            die("Problem extracting read sequence length from read id. Problematic id is $id");
        }

        for(my $i=0; $i < scalar(@sorted); $i++){
            my $curr_start = $sorted[$i];
            my $curr_end = $curr_start + $SEED_LEN + 1;

            if($i == 0 and $curr_end > ($len_seq - $STRETCH_START_WINDOW)){
                $is_ongoing_stretch = 1;
                $stretch_start = $curr_start;
                $stretch_end = $curr_end;
            }

            # check if this is last pos in array
            if($i == scalar(@sorted) - 1){
                # last item
                if($is_ongoing_stretch){
                    push(@{$struct{$id}}, [$stretch_start,$stretch_end])
                }
                last;
            }

            my $next_end = $sorted[$i+1] + $SEED_LEN + 1;

            if($is_ongoing_stretch){
                # check if stetch ends
                # keep ongoing strech if 2 pos are not overlapping by less then 2x seed length
                if($next_end < $curr_start - (2 * ($SEED_LEN + 1)) ){
                    $is_ongoing_stretch = 0;
                    # record stretch
                    push(@{$struct{$id}}, [$curr_start,$stretch_end])
                }
                else{
                    # stretch keeps on going!
                    $stretch_start = $curr_start;
                }
            }
            else{
                if($curr_end > ($len_seq - $STRETCH_START_WINDOW)){
                    $is_ongoing_stretch = 1;
                    $stretch_start = $curr_start;
                    $stretch_end = $curr_end;
                }
            }
        }

        if(scalar(@{$struct{$id}}) == 0){
            delete($struct{$id});
        }
        # old [[40679, 41040]]
        # new [40634, 41040]]
        # if($id eq "b750794c-84f4-48ee-9e2e-8a4d6ce38157_41065"){
        #     print "#### id $id ######\n";
        #     print "init coords=" . dump(@sorted) . "\n";
        #     print "final stretch=" . dump($struct{$id}) . "\n";
        # }

    }

    return \%struct
}

sub flatten_telo_match5{
    my($all_telo_match) = @_;

    my %struct;

    foreach my $id (keys(%$all_telo_match)){
        my @coords = @{$all_telo_match->{$id}};
        my @sorted = sort { $a <=> $b } @coords;
        my $is_ongoing_stretch = 0;
        my $stretch_start = 0;
        my $stretch_end = 0;
        $struct{$id} = [];

        for(my $i=0; $i < scalar(@sorted); $i++){
            my $curr_start = $sorted[$i];
            my $curr_end = $curr_start + $SEED_LEN + 1;

            if($i == 0 and $curr_start < $STRETCH_START_WINDOW){
                $is_ongoing_stretch = 1;
                $stretch_start = $curr_start;
                $stretch_end = $curr_end;
            }

            # check if this is last pos in array
            if($i == scalar(@sorted) - 1){
                # last item
                if($is_ongoing_stretch){
                    push(@{$struct{$id}}, [$stretch_start,$stretch_end])
                }
                last;
            }

            my $next_start = $sorted[$i+1];

            if($is_ongoing_stretch){
                # check if stetch ends
                # keep ongoing strech if 2 pos are not overlapping by less then 2x seed length
                if($next_start > $curr_end + (2 * ($SEED_LEN + 1)) ){
                #if($next_start > $curr_end ){
                    $is_ongoing_stretch = 0;
                    # record stretch
                    push(@{$struct{$id}}, [$stretch_start,$curr_end])
                }
                else{
                    # stretch keeps on going!
                    $stretch_end = $curr_end;
                }
            }
            else{
                if($curr_start < $STRETCH_START_WINDOW){
                    $is_ongoing_stretch = 1;
                    $stretch_start = $curr_start;
                    $stretch_end = $stretch_start + $SEED_LEN + 1;
                }
            }
        }

        if(scalar(@{$struct{$id}}) == 0){
            delete($struct{$id});
        }

        #old final stretch=[[4, 203]]
        #new final stretch=[[4, 373]]
        # if($id eq "17a4f741-215e-41fa-a324-c8c8e145f983_61589"){
        #     print "#### id $id ######\n";
        #     print "init coords=" . dump(@sorted) . "\n";
        #     print "final stretch=" . dump($struct{$id}) . "\n";
        # }
    }

    return \%struct
}


sub read_seqkit{
    my($SEQKIT_IN) = @_;

    print "opening $SEQKIT_IN\n";
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

        if($new_id ne $curr_id){
            # init entry
            $curr_id = $new_id;
            $struct{$curr_id} = [];
        }

        push(@{$struct{$curr_id}},$new_start);
        $curr_start = $new_start;
        $curr_end = $new_end;
    }

    return \%struct;
}
