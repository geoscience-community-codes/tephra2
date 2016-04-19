#$in = "psrose.in";
#$out = "psrose0.eps";
$in = $ARGV[1];
$out = $ARGV[2];

open (IN, ">$in") || die ("can't open : $!");
$min_ht = 10e9;
$max_ht = 0;
$min_speed = 10e9;
$max_speed = 0;
$max_dir = 0;
$min_dir = 361;
open (DATA, "<$ARGV[0]") || die ("can't open : $!");
$total = 0;
$ct = 0;
while (<DATA>) {
  ($ht, $sp, $dir) = split " ", $_;

  if ($ht =~ m/^#/) {}
  else {
  	
  	$ht /= 1000;
  	
  	#if ($dir >= 180) {$dir -=360;}
    if ($ht >= $max_ht) { $max_ht = int($ht) +1; }
  	if ($ht <= $min_ht) { $min_ht = int($ht); }
  	if ($sp >= $max_speed) { $max_speed = int($sp); }
  	if ($sp <= $min_speed) { $min_speed = int($sp); }
  	if ($dir >= $max_dir) { $max_dir = int($dir); }
  	if ($dir <= $min_dir) { $min_dir = int($dir); }
  	
  	$sp *= .16;  
  	$total += $sp;	
    $height[$ct] = $ht;
    $speed[$ct] = $sp;
    $degree[$ct] = $dir;
    $ct++;
  }
}

$i = 0;
for (0..$ct-1) {
	$sca = $speed[$i]/$total *15;
	printf IN "%.0f 0 %f %f %f %f\n", $degree[$i], $height[$i], $sca, $degree[$i], $speed[$i];
	$i++;
}



printf STDERR "Height\tSpeed\tDirection\n";
printf STDERR "%d\t%.1f\t%d\n", $min_ht, $min_speed, $min_dir;
printf STDERR "%d\t%.1f\t%d\n", $max_ht, $max_speed, $max_dir;
#$min_speed *= .5;
#$max_speed *= .5;
$max_dir += 5;
$min_dir -= 5;
#system "psrose $in -Y2i -R0/30/0/360/0/30 -G0 -M0.04i/0.12i/0.1i/255/0/0 -W1p/50 -Ba1g1:'Wind Speed (m/sec)':/g25:'':/:.'Wind Direction (wind blowing toward)': -S3 -C -V > $out";
`gmt makecpt -Cpanoply -T0/$max_ht/5 -I > cpt`;
`gmt psxy --FONT_ANNOT_PRIMARY=8p,2,0 --FONT_LABEL=8p,2,0 $in -R$min_dir/$max_dir/0/$max_speed -Sc0.1c -SV+e -N -Jpa0.16/0.16 -Bxa30g30 -Bya10f1g5+u"ms^-1" -W0.5p,150 -Ccpt -K -P > $out`;
`gmt psscale --FONT_LABEL=8p,2,0 --FONT_ANNOT_PRIMARY=8p,1,0 -D-1.0/1.5/4/.2c -Ccpt -E -Baf+l"km" -O >> $out`;
`gmt ps2raster $out -A -Tg`;
