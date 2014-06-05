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

while (<DATA>) {
#foreach $data  (@Data) {
  ($ht, $sp, $dir) = split " ", $_;
  #print "$data\n";
  if ($ht =~ m/^#/) {}
  else {
  	
  	$ht /= 1000;
  	if ($dir >= 180) {$dir -=360;}
    if ($ht >= $max_ht) { $max_ht = int($ht) +1; }
  	if ($ht <= $min_ht) { $min_ht = int($ht); }
  	if ($sp >= $max_speed) { $max_speed = int($sp); }
  	if ($sp <= $min_speed) { $min_speed = int($sp); }
  	if ($dir >= $max_dir) { $max_dir = int($dir); }
  	if ($dir <= $min_dir) { $min_dir = int($dir); }
  	$sp *= .2;
    printf IN "%.0f 0 %f %f %f\n", $dir, $ht, $dir, $sp;
  }
}
printf STDERR "%d %.1f %d\n", $min_ht, $min_speed, $min_dir;
printf STDERR "%d %.1f %d\n", $max_ht, $max_speed, $max_dir;
#$min_speed *= .2;
#$max_speed *= .2;
$max_dir += 5;
$min_dir -= 5;
#system "psrose $in -Y2i -R0/30/0/360/0/30 -G0 -M0.04i/0.12i/0.1i/255/0/0 -W1p/50 -Ba1g1:'Wind Speed (m/sec)':/g25:'':/:.'Wind Direction (wind blowing toward)': -S3 -C -V > $out";
system "makecpt -Cpanoply -T0/$max_ht/5 -I > cpt";
system "psxy $in -Y2i -R$min_dir/$max_dir/0/$max_speed -Sc.1c -SV -N -Jpa.2/.2 -Ba30g30/f1g5 -Gred -Ccpt  -V -K -P > $out";
system "psscale -D12/6/5/.2c -Ccpt -Li -B:'':/:'Wind Levels (kmasl)': -V -O >> $out";
system "ps2raster $out -A -Tg -V";
