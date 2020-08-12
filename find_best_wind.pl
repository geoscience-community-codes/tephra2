# Some configuration variables
my $EXEI="tephra2-inversion_2020";
my $INPUT="ararat_7samples.in";
my $WIND;
my $ICONF="tephra2-inversion.conf";
my $VENT="ararat.wgs84z38.utm";
my $MACHINE="machines2";
my $COLUMN=200;
my $PARTICLE=100;
my $NODES=2;
my $RUN;

# The results file
my $out = "results.txt";

# Open the wind database for list of wind files
my $DIR = "wind_db";
opendir(DIR, $DIR) || die "Not happening: [$DIR] $!";
my @CONTENTS = readdir(DIR);
foreach (@CONTENTS) {
	unless ($_ =~ /\./) {
		#$Wind[$w] = $_;
		$WIND = "$DIR/$_";
		# Run the inversion code for each wind in the wind database
		$RUN="mpirun -np $NODES -hostfile $MACHINE $EXEI $ICONF $INPUT $WIND";
		system "$RUN";
		# Record the results for each inversion, one line for each run
		$RUN="perl scripts/write_stats.pl parameters.README $_";
		system "$RUN";		
	}
}
