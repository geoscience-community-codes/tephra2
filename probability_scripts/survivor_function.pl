use Math::Random qw(random_uniform random_uniform_integer);
use POSIX;
use Geo::Proj4;
use FileHandle;
use Carp;    
    # warn user (from perspective of caller)
    #carp "string trimmed to 80 chars";
    # die of errors (from perspective of caller)
    #croak "We're outta here!";use Geo::Proj4;
    # die of errors with stack backtrace
    #confess "not implemented";
    # cluck not exported by default
    #use Carp qw(cluck);
    #cluck "This is how we got here!";
use strict;

my $args = @ARGV;
if ($args < 7) {
  print STDERR "USAGE: perl survivor_function.pl tephra2_exe Volcano_file SITE_file wind_directory NUM_simulations utm_zone vei3 [vei4] [vei5] [vei6]\n\n";
  exit;
}
else {
	print STDERR "Running .....\n";
}


our %H;

$H{TEPHRA2_EXE} = $ARGV[0];
$H{VOLCANO} = $ARGV[1];
$H{SITE} = $ARGV[2];
$H{WIND_DB} = $ARGV[3];
$H{END} = $ARGV[4];
$H{ZONE} = $ARGV[5];
print STDERR "$H{VOLCANO} : $H{SITE} : ";
my $arg_num = @ARGV;
my @VEI;
for (my $i = 6; $i < $arg_num; $i++) {
	push(@VEI, $ARGV[$i]);
	print STDERR "$ARGV[$i] ";
} 
print STDERR "\n";
#$H{VEI} = [@VEI];

# Convert long / lat to utm 
my $proj = Geo::Proj4->new(proj => "utm", ellps => "WGS84", datum => "WGS84", zone => $H{ZONE} );

# Open the volcano file for retrieving easting, northing, elevation
my $volc = open_or_die("<", $H{VOLCANO});
while (<$volc>) {
	unless ($_ =~ /^#/ || $_ =~ /^\n/) {
		my ($lon, $lat, $el) = split " ";
		($H{VENT_EASTING}, $H{VENT_NORTHING}) = $proj->forward($lat, $lon);	
		$H{VENT_EASTING} = floor($H{VENT_EASTING});
		$H{VENT_NORTHING} = floor($H{VENT_NORTHING});
		$H{VENT_ELEVATION} = $el;
	}
}

# Open site file and convert location to UTM
my $site = open_or_die("<", "$H{SITE}");
my $site_out = open_or_die(">", "site");
while (<$site>) {
	unless ($_ =~ /^#/ || $_ =~ /^\n/) {
		my ($lon, $lat, $el) = split " ";
		my ($e, $n) = $proj->forward($lat, $lon);
		$n = floor($n);
		$e = floor($e);
		print $site_out "$e $n $el"; 
	}
}

# Open the wind database for volcano chosen
our $DIR = "$H{WIND_DB}";
opendir(DIR, $DIR) || croak "Not happening: [$DIR] $!";
our @CONTENTS = readdir(DIR);

# number of simulation to RUN
our $END = $H{END};

our $min_ms;
our $max_ms;

our $min_col_ht;
our $max_col_ht;

######################################################################
# Open the VEI configuration file and read in each key value pair
# This is also stored in the %H hash array
my $conf = open_or_die("<", "vei.conf");
my $key;
my $value;
while (<$conf>) {
  unless ($_ =~ /^#/ || $_ =~ /^\n/) {
  	my ($key, $value) = split "=",$_;
  	chomp($value);
	  $H{$key} = $value;
  }
}

my $tephra2_output;
my $survivor_output;
my $plot_inputs;
my $data;
my $fixed = 0;
foreach my $v (@VEI) {
		if ($fixed > 0) {
	 		if ($fixed ==1) {
	 			# Fixed mass and column height 
				$min_ms = $v;
				$max_ms= $v;
				$fixed++;
				next;
			}
			elsif ($fixed == 2) {
				$min_col_ht = log($v)/log(10);
				$max_col_ht = log($v)/log(10);
			  $fixed = 0;
			}
		}
	 elsif ($v =~m/fixed/) {
	 	$H{VEI} = $v;
			$fixed++;
			$tephra2_output = 	"output/fixed_tephra2.out";	
			$survivor_output = "output/fixed_survivor.out";
			#$H{BETA} = 1.0;
			next;
		}
		elsif ($v =~ /vei2/) {
		
			# VEI-2 data specified in the config file
			$tephra2_output = 	"output/vei2_tephra2.out";	
			$survivor_output = "output/vei2_survivor.out";
			$H{VEI} = $v;
			
			# mass in kg/m^2
			$min_ms = log(1.0e9)/log(10);
			$max_ms = log(1.0e10)/log(10);
			# column height in km
			$min_col_ht = log(1.0)/log(10);
			$max_col_ht = log(5.0)/log(10);
			#$H{BETA} = 1.0;
		}
		elsif ($v =~ /vei3/) {
		
			# VEI-3 data specified in the config file
			$tephra2_output = 	"output/vei3_tephra2.out";	
			$survivor_output = "output/vei3_survivor.out";
			$H{VEI} = $v;
			
		# mass in kg/m^2
			$min_ms = log(1.0e10)/log(10);
			$max_ms = log(1.0e11)/log(10);
			# column height in km
			$min_col_ht = log(3.0)/log(10);
			$max_col_ht = log(15.0)/log(10);
			#$H{BETA} = 1.0;
		}
		elsif ($v =~ /vei4/) {
	
			# VEI-4 data specified in the config file
			$tephra2_output = 	"output/vei4_tephra2.out";	
			$survivor_output = "output/vei4_survivor.out";
			$H{VEI} = $v;
			
			# mass in kg/m^2
			$min_ms = log(1.0e11)/log(10);
			$max_ms = log(1.0e12)/log(10);
			# column height in km
			$min_col_ht = log(10.0)/log(10);
			$max_col_ht = log(25.0)/log(10);
		
			#$H{BETA} = 1.0;
		}
		elsif ($v =~ /vei5/) {
		
			# VEI-5 data specified in the config file
			$tephra2_output = 	"output/vei5_tephra2.out";	
			$survivor_output = "output/vei5_survivor.out";
			$H{VEI} = $v;
			
			# mass in kg/m^2
			$min_ms = log(1.0e12)/log(10);
			$max_ms = log(1.0e13)/log(10);
			# column height in km
			$min_col_ht = log(20.0)/log(10);
			$max_col_ht = log(35.0)/log(10);
		
			#$H{BETA} = 1.0;
		}
		elsif ($v =~ /vei6/) {
	
			# VEI-6 data specified in the config file
			$tephra2_output = 	"output/vei6_tephra2.out";	
			$survivor_output = "output/vei6_survivor.out";

			# mass in kg/m^2
			$min_ms = log(1.0e13)/log(10);
			$max_ms = log(1.0e14)/log(10);
			# column height in km
			$min_col_ht = log(25.0)/log(10);
			$max_col_ht = log(40.0)/log(10);
		
			#$H{BETA} = 1.0;
		}
		elsif ($v =~ /vei7/) {
			# VEI-7 data specified in the config file
			$tephra2_output = 	"output/vei7_tephra2.out";	
			$survivor_output = "output/vei7_survivor.out";

			# mass in kg/m^2
			$min_ms = log(1.0e14)/log(10);
			$max_ms = log(1.0e15)/log(10);
			# column height in km
			$min_col_ht = log(30.0)/log(10);
			$max_col_ht = log(50.0)/log(10);
		
			#$H{BETA} = 1.0;
		}

	`> $survivor_output`;
	$data = open_or_die(">>", $survivor_output);
	#while (my ($key,$value) = each %H) {
	#	print "$key = $value\n";
	#}
	run_tephra2($tephra2_output, $data);
	$plot_inputs .= $survivor_output." ";
} # END foreach
	
print STDERR "perl survivor_plot.gmt.pl $plot_inputs\n";
`perl survivor_plot.gmt.pl $plot_inputs`;

##############################################################
sub run_tephra2 {

my $tephra2_output = $_[0];
my $data_output = $_[1];

for (my $i = 0; $i < $END; $i++) {
  my $col_ht = random_uniform(1, $min_col_ht, $max_col_ht);
  $col_ht = 10.0**$col_ht;
  my $erupt_mass = random_uniform(1, $min_ms, $max_ms);
  $erupt_mass = 10.0**$erupt_mass;
  $H{PLUME_HEIGHT} = $col_ht * 1000.0;
  $H{ERUPTION_MASS} = $erupt_mass;

  # Mass(kg) = duration(s) x 1000(kg/m^3) x [col_ht(km)/1.67]^4 
  #$H{ERUPTION_MASS} =  $duration * 3.6e6 * ($col_ht/1.67)**4.0;

	my $conf = open_or_die(">", "tephra2.conf");
	while ( my ($key, $value) = each(%H)){
 	print $conf "$key $value\n";
 }

 my $num_contents = @CONTENTS;
 my $wind = random_uniform(1, 2, $num_contents);
 $wind = floor($wind);
 print STDERR "wind file used = $DIR$CONTENTS[$wind]\n";

 # printf DATA "%d %d %s %.0f %.1e %\n", $H{VENT_EASTING}, $H{VENT_NORTHING}, $contents[$wind], $H{PLUME_HEIGHT}, $H{ERUPTION_MASS};

 # the command line for running tephra should go between quotes such as:
 #system "tephra2 tephra2.conf npp.in wind.in >> tephra2.out";
 my $run_command = "$H{TEPHRA2_EXE} tephra2.conf site $DIR$CONTENTS[$wind] > $tephra2_output";
 print STDERR "Now Running: $run_command\n\n";
 system "$run_command";
 my $tephra = open_or_die("<", "$tephra2_output");
 my @Line = <$tephra>;
 foreach my $line (@Line) {
 	my ($east, $north, $elev, $mass) = split " ", $line;
 	if ($east =~ m/^#/ ) {next;}
	#print STDERR "$east, $north, $elev, $mass\n"; 
 	printf $data_output "%d %d %s %.0f %0.1e %s\n", $H{VENT_EASTING}, $H{VENT_NORTHING}, $CONTENTS[$wind], $H{PLUME_HEIGHT}, $H{ERUPTION_MASS}, $mass;
	}
} # END for (my $i = 0; $i < $END; $i++)
} # END sub


sub open_or_die {
	my ($mode, $filename) = @_;
	open my $h, $mode, $filename
		or croak "Could not open '$filename': $!";
return $h;}
