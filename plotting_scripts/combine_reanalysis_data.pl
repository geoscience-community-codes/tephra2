use Math::Trig;
use Scalar::Util qw(looks_like_number);

###############################
# How to run the 'ncdump' code:
###############################
#print "Running: ncdump -v vwnd 'vwnd.2006.nc' > vwnd.2006.nc.out\n";
#system "ncdump -v vwnd 'vwndxxxxx.nc' > vwnd.nc.out";

#print "Running: ncdump -v uwnd 'vwnd.2006.nc' > uwnd.2006.nc.out\n";
#system "ncdump -v uwnd 'uwndxxxxx.nc' > uwnd.nc.out";

#print "Running: ncdump -v hgt 'hgt.2006.nc' > hgt.2006.nc.out\n";
#system "ncdump -v hgt 'hgtxxxx.nc' > hgt.nc.out";

my $conf = $ARGV[0];
print STDERR "Opening configuration file: $conf\n";
open CONF, "< $conf" or die "Can't open $conf : $!";

my %Param;
while(<>) {
	if (/^$/ or /^#/) { next; }
       ($key, $value) = split "=",$_;
        chomp($value);
		$Param{$key} = $value;
        print STDERR "$key=$Param{$key}\n";
} 
$level_file = $Param{LEVEL_FILE};
$uwind_file = $Param{UWIND_FILE};
$vwind_file = $Param{VWIND_FILE};

open  (LEVEL, "<$level_file") || die ("$!");
open (UWIND, "<$uwind_file") || die ("$!");
open (VWIND, "<$vwind_file") || die ("$!");

@level = <LEVEL>;
@uwind = <UWIND>;
@vwind = <VWIND>;

my $data_ct = -1;
my $scale = 1;
my $offset = 0;
my $num_levels = 0;
Loop1: foreach $lev (@level) {
	($name, $eq, $value) = split " ", $lev;
	#print "$name $eq $value\n";
	if ($data_ct >= 0) {
		($val) = split ",", $name;
		#chop($name);
		#print "$name\n";
		if ($val =~ m/;/ or $val =~ m/}/ or $val =~ m/,/ or $val =~ m/^$/) {
			#print "hgt=$val\n";
			next Loop1;
		}
		$level[$data_ct] = $val * $scale + $offset;
		#print " $name:$level[$data_ct] ";
		$data_ct++;
  }
	if ($name eq "hgt") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
		if ($name eq "hgt:scale_factor") {
		#print STDERR "$name $eq $value\n";
		chomp($value);
		chop($value);
		print STDERR "$name $eq $value\n";
		$scale = $value;
    }
    if ($name eq "hgt:add_offset") {
      #print STDERR "$name $eq $value\n";
      chomp($value);
      chop($value);
      print STDERR "$name $eq $value\n";
      $offset = $value;
    }
    if ($name eq "level") {
      chomp($value);
      print STDERR "$name $eq $value\n";
      $num_levels = $value;
    }
}

$data_ct = -1;
$scale = 1;
$offset = 0;
Loop2: foreach $uw (@uwind) {
	($name, $eq, $value) = split " ", $uw;
	if ($data_ct >=0) {
		($val) = split ",", $name;
		#print "$val\n";
		#chop($name);
		if ($val =~ m/;/ or $val =~ m/}/ or $val =~ m/,/ or $val =~ m/^$/) {
			print "uwind=$val\n";
			next Loop2;
		}
		$uwind[$data_ct] = $val * $scale + $offset;
		#print "$name : $uwind[$data_ct]\n";
		$data_ct++;
  }
	if ($name eq "uwnd") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
	if ($name eq "uwnd:scale_factor") {
	    #print STDERR "$name $eq $value\n";
	    chomp($value);
		chop($value);
		print STDERR "$name $eq $value\n";
		$scale = $value;
    }
    if ($name eq "uwnd:add_offset") {
      #print STDERR "$name $eq $value\n";
      chomp($value);
      chop($value);
      print STDERR "$name $eq $value\n";
      $offset = $value;
    } 
}

$data_ct = -1;
$scale = 1;
$offset = 0;
Loop3: foreach $vw (@vwind) {
	($name, $eq, $value) = split " ", $vw;
	if ($data_ct >=0) {
		($val) = split ",", $name;
		#chop($name);
		#print "$val\n";
		if ($val =~ m/;/ or $val =~ m/}/ or $val =~ m/,/ or $val =~ m/^$/) {
			print "vwind=$val\n";
			next Loop3;
		}
		$vwind[$data_ct] = $val * $scale + $offset;
		$data_ct++;
     }
	if ($name eq "vwnd") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
		if ($name eq "vwnd:scale_factor") {
		 #print STDERR "$name $eq $value\n";
		  chomp($value);
		  chop($value);
		  print STDERR "$name $eq $value\n";
		  $scale = $value;
    }
    if ($name eq "vwnd:add_offset") {
     # print STDERR "$name $eq $value\n";
        chomp($value);
        chop($value);
        print STDERR "$name $eq $value\n";
        $offset = $value;
    } 
}
# $data_ct--;
mkdir "wind_db", 0777 unless -d "wind_db";
my $file_no = 0;
my $j = $num_levels;
print "$num_levels $data_ct\n";
for ($i = 0; $i < $data_ct; $i++) {
	if ($j == $num_levels) {
		close OUT;
		$j = 0;
		my $file = "wind_db/wind$file_no";
		$file_no++;
		open OUT, "> $file" or die "Can't open $file: $!";
	}
	if ($uwind[$i] == 0) {
		#print "[uwind=0] ";
		if ($vwind[$i] < 0) {
			$wind_direction = 180;
		}
		elsif ($vwind[$i] > 0) {
			$wind_direction =0;
		}
		#print "$wind_direction\n";
	}
	else {
		$wind_direction = -180/pi *atan($vwind[$i]/$uwind[$i]);
			if ($uwind[$i] > 0) {
				$wind_direction += 90;
				#print "$wind_direction\n";
			}
			else {
				$wind_direction += 270;
				#print "$wind_direction\n";
			}
	}
	$wind_speed = sqrt($uwind[$i]*$uwind[$i]+$vwind[$i]*$vwind[$i]);
	if ($wind_speed > 200) {
		printf STDERR "%d %d %f %f\n", $i, $level[$i], $uwind[$i], $vwind[$i];
	}
  printf OUT "%.0f %.1f %.1f\n", $level[$i], $wind_speed, $wind_direction;
  $j++;
}
$file_no--;
print STDERR "Wrote out $file_no wind files.\n";
print STDERR "Done\n";  
