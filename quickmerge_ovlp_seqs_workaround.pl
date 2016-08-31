#!/usr/bin/perl
#
#

use strict;
use warnings;

my ($cordfile) = @ARGV ;

open(my $table, '<:encoding(UTF-8)', $cordfile)
or die "Could not open file '$cordfile' !";

my @arrovlps; # array of arrays
my @subovlps; # array of arrays subset
my @bv3_id;
my @bv4_id;

sub splitid {
  my $id = @_;
  my @idlen = split(/_len=/, $id, 2);
  return $idlen[1];
}

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
#print "@trim\nhey";
#return @trim;
}

sub overlap2trim_query {
  my($ovlp01, $ovlp02) = @_;
  my($bv4_id01, $bv3_id01, $bv4_c101, $bv4_c201, $bv3_a101, $bv3_a201, $ori01, $innie01, $olen01, $oprop01, $novpend01, $overhang01) = @$ovlp01;
  my($bv4_id02, $bv3_id02, $bv4_d102, $bv4_d202, $bv3_b102, $bv3_b202, $ori02, $innie02, $olen02, $oprop02, $novpend02, $overhang02) = @$ovlp02;
  my @bv4_idarr01 = split(/_len=/, $bv4_id01, 2); my $bv4_len01 = $bv4_idarr01[1];
  my @bv4_idarr02 = split(/_len=/, $bv4_id02, 2); my $bv4_len02 = $bv4_idarr02[1];  
  
  my @trim = ("id", 0, 0);
  my $overlap_length = -1;
  # There are four cases of overlapping coordinates with a different overlap length calculus each time. Get ready for the ride. 
  # First case: queries are oriented both to the right: a1<=a2 && b1<=b2. There is an overlap if and only if a1<=b2 && b1<=a2. Overlap lenght = |b1-a2|
  if ( ($bv3_a101 < $bv3_a201 and $bv3_b102 < $bv3_b202) and ($bv3_a101 <= $bv3_b202 and $bv3_b102 <= $bv3_a201) ){
    $overlap_length = get_overlap_length($bv3_a101, $bv3_a201, $bv3_b102, $bv3_b202); # OVERLAP LENGTH
    my $trim_dist = $overlap_length + 2; 
    # I want to cut the longest contig. Cutting at head coordinate if first contigs is larger, cutting at tail coordinate if second contig larger. 
    # The coordinates for query sequence are always first_coord < second_coord.
    my ($tid, $tx, $ty) = get_trim($bv4_id01,$bv4_id02,$bv4_len01,$bv4_len02,$bv4_c101, $bv4_c201,$bv4_d102, $bv4_d202,$trim_dist);
    @trim = ($tid, $tx, $ty); # TRIM SUGGESTION
  # Second case: queries align face to face -><- on the reference. There is an overlap if and only if a1>=b2 && a2<=b1
  }elsif ( ($bv3_a101 > $bv3_a201 and $bv3_b102 > $bv3_b202) and ($bv3_a101 >= $bv3_b202 and $bv3_b102 >= $bv3_a201) ){
    $overlap_length = get_overlap_length($bv3_a101, $bv3_a201, $bv3_b102, $bv3_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv4_id01,$bv4_id02,$bv4_len01,$bv4_len02,$bv4_c101, $bv4_c201,$bv4_d102, $bv4_d202,$trim_dist);
    @trim = ($tid, $tx, $ty);
  # Third case: queries align backing each other on the reference <-->. There is an overlap if and only if a1>=b1 && a2<=b2
  }elsif ( ($bv3_a101 > $bv3_a201 and $bv3_b102 < $bv3_b202) and ($bv3_a101 >= $bv3_b102 and $bv3_a201 <= $bv3_b202) ) {
    $overlap_length = get_overlap_length($bv3_a101, $bv3_a201, $bv3_b102, $bv3_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv4_id01,$bv4_id02,$bv4_len01,$bv4_len02,$bv4_c101, $bv4_c201,$bv4_d102, $bv4_d202,$trim_dist); 
    @trim = ($tid, $tx, $ty);  
  # Fourth case: both queries are oriented to the left on the reference contig. There is an oberlap if and only if a1<=b1 && a2>=b2
  }elsif ( ($bv3_a101 > $bv3_a201 and $bv3_b102 > $bv3_b202) and ($bv3_a101 >= $bv3_b202 and $bv3_a201 <= $bv3_b102) ){
    $overlap_length = get_overlap_length($bv3_a101, $bv3_a201, $bv3_b102, $bv3_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv4_id01,$bv4_id02,$bv4_len01,$bv4_len02,$bv4_c101, $bv4_c201,$bv4_d102, $bv4_d202,$trim_dist);
    @trim = ($tid, $tx, $ty);
  }

  return ($bv3_id01, $overlap_length, @trim);
}

sub overlap2trim_reference {
  my($ovlp01, $ovlp02) = @_;
  my($bv4_id01, $bv3_id01, $bv4_a101, $bv4_a201, $bv3_c101, $bv3_c201, $ori01, $innie01, $olen01, $oprop01, $novpend01, $overhang01) = @$ovlp01;
  my($bv4_id02, $bv3_id02, $bv4_b102, $bv4_b202, $bv3_d102, $bv3_d202, $ori02, $innie02, $olen02, $oprop02, $novpend02, $overhang02) = @$ovlp02;
  my @bv3_idarr01 = split(/_len=/, $bv3_id01, 2); my $bv3_len01 = $bv3_idarr01[1];
  my @bv3_idarr02 = split(/_len=/, $bv3_id02, 2); my $bv3_len02 = $bv3_idarr02[1];

  my @trim = ("id", 0, 0);
  my $overlap_length = -1;

  if ( ($bv4_a101 < $bv4_a201 and $bv4_b102 < $bv4_b202) and ($bv4_a101 <= $bv4_b202 and $bv4_b102 <= $bv4_a201) ){
    $overlap_length = get_overlap_length($bv4_a101, $bv4_a201, $bv4_b102, $bv4_b202); # OVERLAP LENGTH
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv3_id01,$bv3_id02,$bv3_len01,$bv3_len02,$bv3_c101, $bv3_c201,$bv3_d102, $bv3_d202,$trim_dist);
    @trim = ($tid, $tx, $ty); # TRIM SUGGESTION
  }elsif ( ($bv4_a101 > $bv4_a201 and $bv4_b102 > $bv4_b202) and ($bv4_a101 >= $bv4_b202 and $bv4_b102 >= $bv4_a201) ){
    $overlap_length = get_overlap_length($bv4_a101, $bv4_a201, $bv4_b102, $bv4_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv3_id01,$bv3_id02,$bv3_len01,$bv3_len02,$bv3_c101, $bv3_c201,$bv3_d102, $bv3_d202,$trim_dist);
    @trim = ($tid, $tx, $ty);
  }elsif ( ($bv4_a101 > $bv4_a201 and $bv4_b102 < $bv4_b202) and ($bv4_a101 >= $bv4_b102 and $bv4_a201 <= $bv4_b202) ) {
    $overlap_length = get_overlap_length($bv4_a101, $bv4_a201, $bv4_b102, $bv4_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv3_id01,$bv3_id02,$bv3_len01,$bv3_len02,$bv3_c101, $bv3_c201,$bv3_d102, $bv3_d202,$trim_dist);
  }elsif ( ($bv4_a101 > $bv4_a201 and $bv4_b102 > $bv4_b202) and ($bv4_a101 >= $bv4_b202 and $bv4_a201 <= $bv4_b102) ){
    $overlap_length = get_overlap_length($bv4_a101, $bv4_a201, $bv4_b102, $bv4_b202);
    my $trim_dist = $overlap_length + 2;
    my ($tid, $tx, $ty) = get_trim($bv3_id01,$bv3_id02,$bv3_len01,$bv3_len02,$bv3_c101, $bv3_c201,$bv3_d102,$bv3_d202,$trim_dist);
  }
  
  return ($bv4_id01, $overlap_length, @trim);
}

sub is_at_edge {
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

my @onearr = (1,1);

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
            my ($bv3_contig, $overlap_length, @trim) = overlap2trim_query(\@$ovlp01, \@$ovlp02) ;
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
            my ($bv4_contig, $overlap_length, @trim) = overlap2trim_reference(\@$ovlp01, \@$ovlp02) ;
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



=pod
foreach my $overlap (@arrovlps) {
  my $bv3_id = @$overlap[1]; 
  foreach my $id3 (@bv3_id) {
    if ($id3 eq $bv3_id){
      push @subsetovl, \@$overlap; 
    }
  }
  # Now you have a subset of @arrovlps for which all id3 are the same
  # You need to look for overlaps within this subset
  foreach my $sub01_overlap (@subsetovl){   
    my($sub01_bv4_id, $sub01_bv3_id, $sub01_bv4_x, $sub01_bv4_y, $sub01_bv3_x, $sub01_bv3_y, $sub01_orientation, $sub01_innie, $sub01_olen, $sub01_oprop, $sub01_novpend, $sub01_overhang) = @$sub01_overlap;
    my $sub01_bv4_len = splitid($sub01_bv4_id);
    foreach my $sub02_overlap (@subsetovl){
      my($sub02_bv4_id, $sub02_bv3_id, $sub02_bv4_x, $sub02_bv4_y, $sub02_bv3_x, $sub02_bv3_y, $sub02_orientation, $sub02_innie, $sub02_olen, $sub02_oprop, $sub02_novpend, $sub02_overhang) = @$sub02_overlap;
      my $sub02_bv4_len = splitid($sub02_bv4_id);
      print $sub01_bv4_len,"\t", $sub02_bv4_len, "\n";
    }
  }
}


    my @idlen = split(/_len=/, $subbv4_id, 2); 
    my $subbv4_len = $idlen[1];
    if ($subbv4_id ne $bv4_id and ($subbv3_x ~~ [$bv3_x..$bv3_y] xor $subbv3_y ~~ [$bv3_x..$bv3_y]) ){
      my $overlapsize = abs($bv3_x - $subbv3_x);
      my @to_trim;
      if ($bv4_len > $subbv4_len){
        my $trimlim = $bv4_y-$overlapsize-2;
        push(@to_trim, $bv4_id, $bv4_x, $trimlim); 
      }else{
        my $trimlim = $subbv4_y-$overlapsize-2;
        push(@to_trim,$subbv4_id, $subbv4_x, $trimlim);
      }
      #print "@to_trim\n";
      push @to_trim_array, \@to_trim unless \@to_trim ~~ @to_trim_array;
      print "@$overlap\n", "@$suboverlap\n","\n" unless \@to_trim ~~ @to_trim_array;
    } elsif ($subbv4_id ne $bv4_id and ($subbv3_x ~~ [$bv3_y..$bv3_x] xor $subbv3_y ~~ [$bv3_y..$bv3_x]) ){
      my $overlapsize = abs($bv3_y - $subbv3_y);
      my @to_trim;
      if ($bv4_len > $subbv4_len){
        my $trimlim = $bv4_x-$overlapsize;
        push(@to_trim, $bv4_id, $trimlim, $bv4_y);
      }else{
        my $trimlim = $subbv4_y-$overlapsize;
        push(@to_trim,$subbv4_id, $subbv4_x, $trimlim);
      }
      push @to_trim_array, \@to_trim unless \@to_trim ~~ @to_trim_array;
      print "@$overlap\n", "@$suboverlap\n","\n" unless \@to_trim ~~ @to_trim_array;    
    }
  }

}

=cut

close $table;
