use Carp;  

# Extract 5th column; Sort mass in ascending order; 
# output: mass vs. N-i/N, where N = total number of mass values, i = sort ranking (i.e. index, starts at 0)

my $arg_num = @ARGV;
my @in;
my @output;
my @sorted;
my $max_out = 100;
my $total;

my $i = 0;
foreach my $arg (@ARGV) {
	# Open output file on command line for reading
	$in[$i] = open_or_die("<", $arg);
	my @Mass;
	# Read each line from file and save mass values
	my $j = 0;
	my $f = $in[$i];
	while (<$f>) {
	($east, $north, $wind, $colHT, $tMass, $Mass[$j++]) = split " ";
	}
	
	# Sort mass values from smallest to largest
	my @Sorted_data = sort {$a <=> $b} @Mass;
	$total = scalar(@Sorted_data);
	print STDERR "Num of data values = $total\n";
	
	# Create a file for sorted output and open for writing
	$output[$i] = "$arg.sorted";
	$sorted[$i] = open_or_die(">", $output[$i]);	
	print "Opening $output[$i] for plotting\n";
	for (my $j = 0; $j < $total; $j++) {
		$out = ($total-$j) / $total * 100;
		print {$sorted[$i]} "$Sorted_data[$j] $out\n";
		if ($Sorted_data[$j] > $max_out) {
			$max_out = $Sorted_data[$j];
		}
	}
	
	$i++;
}	

print STDERR "$total $max_out\n";
#$total = 1/$total;

$eps = "survivor.eps";
my @colors = qw(lightblue green2 yellowgreen yellow2 orange1 red3 black);
#system "psbasemap -R.1/$max_out/$total/1 -Jx2l/4l -Ba-1pg3p:'Mass Loading (kg/m\@+2\@+)':/a1pg3p:'Exceedence Probability':WS -P -K -V > $eps";
system "psbasemap -R.1/$max_out/0/100 -Jx2l/.06 -Ba-1pg3p:'Mass Loading (kg/m\@+2\@+)':/a10g10f10:'Exceedence Probability (%)':WS -P -K -V > $eps";
for ($i = 0; $i < $arg_num; $i++) { 
	system "psxy $output[$i] -R -Jx -W2p,$colors[$i] -O -K -V >> $eps";
}
system "psbasemap -R -Jx -B:.'': -O -V >> $eps";

`ps2raster $eps -A -P -V -Tg`;

sub open_or_die {
	my ($mode, $filename) = @_;
	open my $h, $mode, $filename
		or croak "Could not open '$filename': $!";
return $h;
}
