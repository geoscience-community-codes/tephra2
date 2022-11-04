#!/usr/bin/perl

# devise an input grid that scales with distance from the
# volcano - this is crucial for producing smooth contours
# on the map - but otherwise not important

$conf = "create_grid.conf";
#$out = "grid.in";

open(DATA, "<$conf") || die("can't open $conf: $!");
#open(OUT, ">$out") || die("can't open $out: $!"); 

while (<DATA>) {
  ($key, $value) = split(" ", $_);
  $grid{$key} = $value;
  print STDERR "$key = $grid{$key}\n";
}


$DEG2RAD = 0.017453293;
$z = $grid{ELEVATION};

for (my $deg = 0; $deg<360; $deg+=5) {
  $theta = $deg * $DEG2RAD;
  $r = 0;
  for (my $i = 5; $i < $grid{NUM_RINGS}; $i++){
    $r += $grid{GRID_STEP} * $i;
    $x = $r * cos($theta) + $grid{VOLCANO_EAST};
    $y = $r * sin($theta) + $grid{VOLCANO_NORTH}; 
    printf "%.0f %.0f %.0f\n",$x,$y, $z;    
  }
}

for (my $x = $grid{MIN_EAST}; $x <= $grid{MAX_EAST}; $x += $grid{SPACING}) {
  for (my $y = $grid{MIN_NORTH}; $y <= $grid{MAX_NORTH}; $y += $grid{SPACING}) {
    printf "%.0f %.0f %.0f\n", $x, $y, $z;
  }
}     

