#!/usr/bin/perl
# abodrug 31.08.2016
# This script detects overlaps at contig edges in the summary output of quickmerge and suggests trimming solutions

use strict;
use warnings;

my ($cordfile) = @ARGV ;

open(my $table, '<:encoding(UTF-8)', $cordfile)
or die "Could not open file '$cordfile' !";

my @arrovlps; # array of arrays
my @subovlps; # array of arrays subset
my @bv3_id; # queries ids
my @bv4_id; # reference ids

sub get_overlap_length {
  my @sorted_coordinates = sort { $a <=> $b } @_;
  my $overlap_length = $sorted_coordinates[2] - $sorted_coordinates[1]; # Here, you sorted all the four coordinates on a overlap on ref                                                                                                                         # the middle two are the ones which should be overlapping
  return $overlap_length
}

sub get_trim {
  my ($id01,$id02,$len01,$len02,$c101,$c201,$d102, $d202,$trim_dist) = @_ ;
  my @sorted_ccoords =  sort { $a <=> $b } ($c101, $c201); # this is for splitting queries
  my @sorted_dcoords =  sort { $a <=> $b } ($d102, $d202);
  if ( $len01 >= $len02 ){
    return ($id01, 1, $sorted_ccoords[1]-$trim_dist);
  }else{
    return ($id02, $sorted_dcoords[0]+$trim_dist, $len02);
  }
}

sub overlap2trim {
  # There are four cases of overlapping coordinates with a different overlap length calculus each time. Get ready for the ride. 
  # First case: queries are oriented both to the right: a1<=a2 && b1<=b2. There is an overlap if and only if a1<=b2 && b1<=a2. Overlap lenght = |b1-a2|
  # Second case: queries align face to face -><- on the reference. There is an overlap if and only if a1>=b2 && a2<=b1
  # Third case: queries align backing each other on the reference <-->. There is an overlap if and only if a1>=b1 && a2<=b2
  # Fourth case: both queries are oriented to the left on the reference contig. There is an oberlap if and only if a1<=b1 && a2>=b2
  # Coordinate constrains is different if it is queries or references which overlap
  my($ovlp01, $ovlp02, $type) = @_;
  my($ref_id01, $qwr_id01, $ref_c101, $ref_c201, $qwr_a101, $qwr_a201, $ori01, $innie01, $olen01, $oprop01, $novpend01, $overhang01) = @$ovlp01;
  my($ref_id02, $qwr_id02, $ref_d102, $ref_d202, $qwr_b102, $qwr_b202, $ori02, $innie02, $olen02, $oprop02, $novpend02, $overhang02) = @$ovlp02;
  my @ref_idarr01 = split(/_len=/, $ref_id01, 2); my $ref_len01 = $ref_idarr01[1];
  my @ref_idarr02 = split(/_len=/, $ref_id02, 2); my $ref_len02 = $ref_idarr02[1];
  my @qwr_idarr01 = split(/_len=/, $qwr_id01, 2); my $qwr_len01 = $qwr_idarr01[1];
  my @qwr_idarr02 = split(/_len=/, $qwr_id02, 2); my $qwr_len02 = $qwr_idarr02[1];
  my @trim = ("id", 0, 0);
  my $overlap_length = -1;
  if ($type eq "query"){
    if ( ($qwr_a101 < $qwr_a201 and $qwr_b102 < $qwr_b202) and ($qwr_a101 <= $qwr_b202 and $qwr_b102 <= $qwr_a201) ){
      $overlap_length = get_overlap_length($qwr_a101, $qwr_a201, $qwr_b102, $qwr_b202); # OVERLAP LENGTH
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($ref_id01,$ref_id02,$ref_len01,$ref_len02,$ref_c101, $ref_c201,$ref_d102, $ref_d202,$trim_dist);
      @trim = ($tid, $tx, $ty); # TRIM SUGGESTION
    }elsif ( ($qwr_a101 > $qwr_a201 and $qwr_b102 > $qwr_b202) and ($qwr_a101 >= $qwr_b202 and $qwr_b102 >= $qwr_a201) ){
      $overlap_length = get_overlap_length($qwr_a101, $qwr_a201, $qwr_b102, $qwr_b202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($ref_id01,$ref_id02,$ref_len01,$ref_len02,$ref_c101, $ref_c201,$ref_d102, $ref_d202,$trim_dist);
      @trim = ($tid, $tx, $ty);
    }elsif ( ($qwr_a101 > $qwr_a201 and $qwr_b102 < $qwr_b202) and ($qwr_a101 >= $qwr_b102 and $qwr_a201 <= $qwr_b202) ) {
      $overlap_length = get_overlap_length($qwr_a101, $qwr_a201, $qwr_b102, $qwr_b202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($ref_id01,$ref_id02,$ref_len01,$ref_len02,$ref_c101, $ref_c201,$ref_d102, $ref_d202,$trim_dist);
      @trim = ($tid, $tx, $ty);
    }elsif ( ($qwr_a101 > $qwr_a201 and $qwr_b102 > $qwr_b202) and ($qwr_a101 >= $qwr_b202 and $qwr_a201 <= $qwr_b102) ){
      $overlap_length = get_overlap_length($qwr_a101, $qwr_a201, $qwr_b102, $qwr_b202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($ref_id01,$ref_id02,$ref_len01,$ref_len02,$ref_c101, $ref_c201,$ref_d102, $ref_d202,$trim_dist);
      @trim = ($tid, $tx, $ty);
    }
    return ($qwr_id01, $overlap_length, @trim);
  }elsif ($type eq "reference"){
    if ( ($ref_c101 < $ref_c201 and $ref_d102 < $ref_d202) and ($ref_c101 <= $ref_d202 and $ref_d102 <= $ref_c201) ){
      $overlap_length = get_overlap_length($ref_c101, $ref_c201, $ref_d102, $ref_d202); # OVERLAP LENGTH
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($qwr_id01,$qwr_id02,$qwr_len01,$qwr_len02,$qwr_a101, $qwr_a201,$qwr_b102, $qwr_b202,$trim_dist);
      @trim = ($tid, $tx, $ty); # TRIM SUGGESTION
    }elsif ( ($ref_c101 > $ref_c201 and $ref_d102 > $ref_d202) and ($ref_c101 >= $ref_d202 and $ref_d102 >= $ref_c201) ){
      $overlap_length = get_overlap_length($ref_c101, $ref_c201, $ref_d102, $ref_d202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($qwr_id01,$qwr_id02,$qwr_len01,$qwr_len02,$qwr_a101, $qwr_a201,$qwr_b102, $qwr_b202,$trim_dist);
      @trim = ($tid, $tx, $ty);
    }elsif ( ($ref_c101 > $ref_c201 and $ref_d102 < $ref_d202) and ($ref_c101 >= $ref_d102 and $ref_c201 <= $ref_d202) ) {
      $overlap_length = get_overlap_length($ref_c101, $ref_c201, $ref_d102, $ref_d202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($qwr_id01,$qwr_id02,$qwr_len01,$qwr_len02,$qwr_a101, $qwr_a201,$qwr_b102, $qwr_b202,$trim_dist);
    }elsif ( ($ref_c101 > $ref_c201 and $ref_d102 > $ref_d202) and ($ref_c101 >= $ref_d202 and $ref_c201 <= $ref_d102) ){
      $overlap_length = get_overlap_length($ref_c101, $ref_c201, $ref_d102, $ref_d202);
      my $trim_dist = $overlap_length + 2;
      my ($tid, $tx, $ty) = get_trim($qwr_id01,$qwr_id02,$qwr_len01,$qwr_len02,$qwr_a101, $qwr_a201,$qwr_b102,$qwr_b202,$trim_dist);
    }
    return ($ref_id01, $overlap_length, @trim);
  }
}

sub is_at_edge {
  # This subroutine checks if an overlap is situated at the edges of contigs. If both reference and query overlap at edges, returns 1.
  my $bool = 0;
  my($ref_id, $qwr_id, $ref_a1, $ref_a2, $qwr_b1, $qwr_b2, $ori01, $innie01, $olen01, $oprop01, $novpend01, $overhang01) = @_;
  my @ref_idarr = split(/_len=/, $ref_id, 2); my $ref_len = $ref_idarr[1];
  my @qwr_idarr = split(/_len=/, $qwr_id, 2); my $qwr_len = $qwr_idarr[1];
  my $maxoverhang = 500;
  my @wanted_edges_ref = (1+$maxoverhang, $ref_len-$maxoverhang);
  my @wanted_edges_qwr = (1+$maxoverhang, $qwr_len-$maxoverhang);
  my @edges_ref = ($ref_a1, $ref_a2);
  my @edges_qwr = sort { $a <=> $b } ($qwr_b1, $qwr_b2);
  #print $ref_a1, "\n";
  if ( ( $edges_ref[0] < $wanted_edges_ref[0] or $edges_ref[1] > $wanted_edges_ref[1]) and ($edges_qwr[0] < $wanted_edges_qwr[0] or $edges_qwr[1] > $wanted_edges_qwr[1])  ){
    #print $ref_a1, " HEYHEYHEYYES?\n";
    $bool = 1;
  }

  return $bool; 
}

while (my $line = <$table>) {
  chomp $line;
  my @overlap = (split /	/, $line);
  push @arrovlps, \@overlap;
  push @bv3_id, $overlap[1] unless $overlap[1] ~~ @bv3_id ;
  push @bv4_id, $overlap[0] unless $overlap[0] ~~ @bv4_id ;
}

close $table; # closed the file

my @onearr = (1,1); # doesn't work directly so I put it in a variable

# Looping over ids to group 
my @treated_suboverlap_ids;
foreach my $id3 (@bv3_id){
  foreach my $ovlp (@arrovlps){
    my $ovlp_id3 = @$ovlp[1];
    if ($id3 eq $ovlp_id3){
      push @subovlps, \@$ovlp;
    }
  }
  if (scalar(@subovlps) > 1) { # For an id, there should be more than one alignment to be able to compare... two alignments. 
    foreach my $ovlp01 (@subovlps){
      foreach my $ovlp02 (@subovlps){
        # I do not want to compare same alignments or alignment to identique query
        my @edgy = (is_at_edge(@$ovlp01), is_at_edge(@$ovlp02));
        unless (@$ovlp01 ~~ @$ovlp02 or @$ovlp01[0] eq @$ovlp02[0] or @$ovlp02[0] ~~ @treated_suboverlap_ids){
          if ( @edgy ~~ @onearr ){
            my ($bv3_contig, $overlap_length, @trim) = overlap2trim(\@$ovlp01, \@$ovlp02, "query") ;
            if ($overlap_length ne -1) {
              print "####\nOn the QUERY contig: ", $bv3_contig, "\nTwo reference contigs overlap:\n", "\t\t@$ovlp01\n", "\t\t@$ovlp02\n","Overlap length is:", $overlap_length ;
              print "\nThe solution provided is the follwing trim: ", "@trim\n", "\n @edgy\n";
            }
          }
        }
      }
      push @treated_suboverlap_ids, @$ovlp01[0]; # Overlap 01 was treated, it will not be looked at a second time as Overlap 02 to avoid duplicated output.  
    }
  }
  undef(@subovlps);
}

foreach my $id4 (@bv4_id){
  foreach my $ovlp (@arrovlps){
    my $ovlp_id4 = @$ovlp[0];
    if ($id4 eq $ovlp_id4){
      push @subovlps, \@$ovlp;
    }
  }
  if (scalar(@subovlps) > 1) {
    foreach my $ovlp01 (@subovlps){
      foreach my $ovlp02 (@subovlps){
        my @edgy = (is_at_edge(@$ovlp01), is_at_edge(@$ovlp02));
        unless (@$ovlp01 ~~ @$ovlp02 or @$ovlp01[1] eq @$ovlp02[1] or @$ovlp02[1] ~~ @treated_suboverlap_ids){
          if ( @edgy ~~ @onearr ){
            my ($bv4_contig, $overlap_length, @trim) = overlap2trim(\@$ovlp01, \@$ovlp02, "reference") ;
            if ($overlap_length ne -1) {
              print "####\nOn the REFERENCE contig: ", $bv4_contig, "\nTwo query contigs overlap:\n", "\t\t@$ovlp01\n", "\t\t@$ovlp02\n","Overlap length is:", $overlap_length ;
              print "\nThe solution provided is the follwing trim: ", "@trim\n", "\n @edgy\n";
            }
          }
        }
      }
      push @treated_suboverlap_ids, @$ovlp01[1];
    }
  }
  undef(@subovlps);
}
