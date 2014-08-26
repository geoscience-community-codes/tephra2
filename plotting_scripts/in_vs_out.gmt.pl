$points_in = $ARGV[0];
$points_out = "tephra.out";
#$points_out = $ARGV[1];
$points = "plot.data";
$pointsxy = "plotxy.data";
$diff_less = "diff.less";
$diff_more = "diff.more";
#$out = "in_vs_out.eps";
$out = $ARGV[1];
$min=1e12;
$max=0;
$rmse = 0;
$ct = 0;
$max_diff = 0;
$min_diff = 1e6;


@data_in;
@data_out;
open(DATA_IN, "<$points_in") || die("cannot open $points_in: $!");
open(DATA_OUT, "<$points_out") || die("cannot open $points_out: $!");
open(DATA, ">$points") || die("cannot open $points: $!");
open(DATAXY, ">$pointsxy") || die ("cannot open $pointsxy: $!");
open (DM, ">$diff_more") || die ("cannot open $diff_more: $!");
open (DL, ">$diff_less") || die ("cannot open $diff_less: $!");
$line2 = <DATA_OUT>;
while (<DATA_IN>) {
  $line1 = $_;
  $line2 = <DATA_OUT>;
  ($e1, $n1, $el1, $data_in[$ct]) = split(" ", $line1);
  ($e2, $n2, $el2, $data_out[$ct]) = split(" ", $line2);
  if ($data_in[$ct] < $min) {$min = $data_in[$ct];}
  if ($data_in[$ct] > $max) {$max = $data_in[$ct];}
  if ($data_out[$ct] < $min) {$min = $data_out[$ct];}
  if ($data_out[$ct] > $max) {$max = $data_out[$ct];}
  $s_data_in = sqrt($data_in[$ct]);
  $s_data_out = sqrt($data_out[$ct]);
  printf DATA "%.3f %.3f 0\n", $s_data_in, $s_data_out;

#$data_in /= 10;
#$data_out /= 10;
#if ($data_out < $min) {$min = $data_out;}
# if ($data_out > $max) {$max = $data_out;}
# $diff = sqrt(($data_in-$data_out)*($data_in-$data_out));
#if ($data_in > 0) {$diff /= $data_in;}

#system "psxy -JX -R -N -Sc.15i -G0 -V -O -K <<eof>> $out
#$s_data_in $s_data_out
#eof
#";

  

  # DOes the model over or under estimate the observed data

  $diff = $data_out[$ct] - $data_in[$ct];
  $diff2 = abs($diff);
  $rmse += ($diff*$diff);
#if ($data_in > $max_diff) {$max_diff = $data_in;}
  if ($diff2 > $max_diff) {$max_diff = $diff2;}
#if ($data_in < $min_diff) {$min_diff = $data_in;}
  if ($diff2 < $min_diff) { $min_diff = $diff2;}

# IF model over-estimates
  if ($diff > 0) {
  	if ($data_in[$ct] != 0) {$diff /= $data_in[$ct];}
  	printf DM "%.0f %.0f %.2f %g\n", $e1, $n1, $data_in[$ct], $diff;
  # ELSE IF model under-estimates
  } elsif ($diff < 0) {
    if ($data_in[$ct] != 0) {$diff /= $data_in[$ct];}
    printf DL "%.0f %.0f %.2f %g\n", $e1, $n1, $data_in[$ct], abs($diff);
  } else { printf "$e1 $e2 $data_in[$ct]";}
  $ct++;
}

$pmin = sqrt($min);
$pmax = sqrt($max);
printf stderr "Min=$pmin  MAX=$pmax\n";
system "psbasemap --HEADER_FONT_SIZE=14 -JX6i -R$pmin/$pmax/$pmin/$pmax -Ba2g2:'sqrt(Observed [kg m\@+-2\@+])':/a2g2:'sqrt(Calculated [kg m\@+-2\@+)]':/:.'':WS -K -P > $out ";
system "psxy -JX -R -W2p/255/0/0 -O -K <<eof>> $out
$pmin $pmin
$pmax $pmax
eof
";

foreach my $dat (@data_in) {
    $log_diff = 0;
    #$out_cm = $data_out[$ct]/10.43;
    #$in_cm = $data_in[$ct]/10.43;
    $diff = abs($data_out[$ct] - $data_in[$ct]);
    if ($diff != 0) { $log_diff = log($diff)/log(10); }
    #$diff2 = ($diff/$max_diff) * $diff;
    
    printf DATAXY "%.0f %.0f %.2f %g\n",$e1, $n1, $data_in[$ct], $log_diff;
}

#The normalized root mean squared error (NRMSE) is the RMSE divided by the range of observed values
system "minmax $points";
$rmse /= $ct;
$rmse = sqrt($rmse);
$rmse /= ($max - $min);
$title = sprintf "NRMSE = %.2f", $rmse;
printf "$title\n";
system "psbasemap --HEADER_FONT_SIZE=14 -JX -R -B:'':/:'':/:.'$title':wnes -K -O >> $out ";
system "psxy $points -JX -N -Sc.18i -G0 -R -O >> $out";

system "ps2raster $out -A -Tg";
system "echo $rmse > nrmse";
#system "perl write_conf.pl parameters.README tephra2.conf 400 400 $rmse";

