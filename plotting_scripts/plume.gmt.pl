use Math::Trig;
use warnings;
use strict;
use Carp ();
 local $SIG{__WARN__} = \&Carp::cluck;
 
BEGIN { $| = 1 }

#my $args = @ARGV;
#if ($args < 5) {
#  print STDERR "USAGE: 
#  perl plume.pl <MAX column height (m)> <vent elevation (m)> <column steps> <alpha> <beta>\n\n";
#}
open DAT, "<", $ARGV[0] or die "cannot open < $ARGV[0]: $!";
my @Data = <DAT>;
(our $max_plume_elevation, our $vent_height, our $col_steps, our $alpha, our $beta) = split " ", $Data[0];

# GMT 
my $out = $ARGV[1];

#our $max_plume_elevation = $Data[0];
#our $vent_height = $Data[1];
#our $col_steps = $Data[2];
#our $alpha = $Data[3];
#our $beta = $Data[4];

print STDERR "MAX column height = $max_plume_elevation\nvent elevation = $vent_height\nColumn steps = $col_steps\nalpha = $alpha\nbeta = $beta\n";

# GMT
my $west = 0;
my $east = 1;
my $south = 0;
my $north = $max_plume_elevation+1000; 
my $max_prob = 0;

my $ht_section_width = $max_plume_elevation - $vent_height; 
my $ht_step_width = $ht_section_width / $col_steps;
my $step_norm = $ht_step_width / $ht_section_width;
my $sum_prob = 0.0;
my $x_norm = 0.0;
my $cum_prob = 0.0;
my $x1;
my $p_out = "column_probabilities";
open(PROB, " >$p_out") || die("cannot open $p_out: $!");
printf(STDERR "Vent elevation = %.0f  Max Plume Ht = %.0f  step = %g\n", $vent_height, $max_plume_elevation, $step_norm); 

#First calculate the normalization constant, so that the total probability integrates to 1
my $prob = plume_pdf2($x_norm, $step_norm, $alpha, $beta, $sum_prob);
$sum_prob = $prob;
#printf(STDERR "\nSum of prob: %g\n\n", $sum_prob);
system "echo '' > points";
#Now calculate the probability
my $x = $vent_height;
$x_norm = 0.0;
print(PROB "# x | probability |  Normalized-x | cummulative-probability\n");
for (my $i = 0; $i < $col_steps; $i++) {
     $x += $ht_step_width;
     $x_norm += $step_norm;
     if ($x_norm eq 0) {
       $x1 = $x_norm + 0.001;
       $prob = plume_pdf2(($x1), $step_norm, $alpha, $beta, $sum_prob);
     } elsif ($x_norm eq 1) {
       $x1 = $x_norm - 0.001;
       $prob = plume_pdf2($x1, $step_norm, $alpha, $beta, $sum_prob);
     } else {
       $x1 = $x_norm;
       $prob = plume_pdf2($x1, $step_norm, $alpha, $beta, $sum_prob);
     } 
       $prob /= $sum_prob;
     
       # This is the total probability and must be equal to 1
       $cum_prob += $prob;
       if ($max_prob < $prob) { $max_prob = $prob; }
       printf(PROB "%0.1f\t\t%0.3g\t\t%0.3g\t\t%0.3g\n", $x, $prob, $x1, $cum_prob);
       # GMT - output plot points
       $prob *= 100;
       system "echo $prob $x >> points";
       
       
       #system "psxy -Sc.2c -R -Gred -JX -W.25p,0 -N -O -K -V >> $out << eof
#$prob $x 
#eof
#";
}
$east = $max_prob * 100 + 1;
$east = 3;
# GMT - draw basemap
`gmt psbasemap  -JX2i/3i -X1i -Y1i -R$west/$east/$south/$north --MAP_FRAME_TYPE=plain --MAP_FRAME_PEN=0.5 --MAP_TICK_PEN=0.5 --FONT_LABEL=10 --FONT_ANNOT_PRIMARY=8 --FORMAT_FLOAT_OUT=%g -Ba1:'probability (%)':/a2000:'Column  (masl)':/:.'':WS -P -K > $out`;

`gmt psxy points -R -JX -W1p,0 -O -K >> $out`;
`gmt psxy points -Sc.1c -R -G150 -JX -W.25p,0 -O >> $out`;
#system "psbasemap --BASEMAP_TYPE=plain --FRAME_PEN=0.5 --LABEL_FONT_SIZE=10 --ANNOT_FONT_SIZE=8 --D_FORMAT=%g -JX -R -Ba0.1:'Probability':/a2000:'Column  (masl)':/:.'':WS -V -O >> $out";

`gmt ps2raster $out -A -Tg`;

sub plume_pdf2 {
	my $x_norm = $_[0];
	my $step = $_[1];
	my $alpha = $_[2];
	my $beta = $_[3];
	my $sum_prob = $_[4];
	
	my $probability = 0.0;
  my $x, my $x1;
  my $prob;
  my $i;
  
  $x = $x_norm;
  if ($sum_prob eq 0) {
    for ($i = int 0; $i < $col_steps; $i++) {
      #step is the small slice of the column as a fraction of the whole
      $x += $step;
      if ($x eq 0) {
        $x1 = $x + 0.001;
        #prob = pow(x1, (alpha - 1.0)) * pow((1.0 - x1), (beta - 1.0));
        $prob = ($x1**($alpha - 1.0)) * ((1.0 - $x1)**($beta - 1.0));
      } elsif ($x eq 1) {
        $x1 = $x - 0.001;
        #prob = pow(x1, (alpha - 1.0)) * pow((1.0 - x1), (beta - 1.0));
        $prob = ($x1**($alpha - 1.0)) * ((1.0 - $x1)**($beta - 1.0));
      } else {
         $x1 = $x;
        #prob = pow(x, (alpha - 1.0)) * pow((1.0 - x), (beta - 1.0));
        $prob = ($x1**($alpha - 1.0)) * ((1.0 - $x1)**($beta - 1.0));
        #printf(STDOUT "[%d]  x=%0.3g\tprob=%0.3g\talpha=%g\tbeta=%g\n", $i, $x, $prob, $alpha, $beta);
      }
      #printf(STDERR "[%d]  x=%g\tstep = %g\tprob=%g\n", $i, $x, $step,$prob); 
      $probability += $prob;
	  }
  } else {
    #Just calculate the probability for one column step
    #probability = pow(x, (alpha - 1.0) ) * pow((1.0 - x), (beta - 1.0));
    $probability = ($x**($alpha - 1.0)) * ((1.0 - $x)**($beta - 1.0));
    #printf(STDOUT "x=%g prob=%g\n", $x, $probability); 
  }
  return $probability;
}



