my $c_steps = $ARGV[2];
my $p_steps = $ARGV[3];

my $Col_Ht;
my $Alpha;
my $Beta;
my $Diff_coef;
my $Fall_tt;
my $Max_phi;
my $Min_phi;
my $Median_phi;
my $Std_phi;
my $Total_mass;

open(FIT, "<$ARGV[4]") || die ("Cannot open $ARGV[4]:$!");
($nrmse) = split("\n"); 

open(README, "<$ARGV[0]") || die("Cannot open $ARGV[0]:$!");
    
LINE:	while(<README>) {

      # Record the results from this run
      ($a, $b) = split(":", $_);
      #print "$_\n";
      
      if ($a =~ m/Parameter Ranges/) {
      	last LINE;
      }
      
      if ($a =~ m/Max Column Height/) { 
      	($c, $d) = split(" ", $b); 
      	$Col_Ht = $c; 
      	#print "Max Column Height = $Col_Ht\n";
      }
      elsif ($a =~ m/Alpha Param/) { 
      	($c, $d) = split(" ", $b);
      	$Alpha = $c; 
      	#print "Alpha Param = $Alpha";
      }
      elsif ($a =~ m/Beta Param/) { 
      	($c, $d) = split(" ", $b);
      	$Beta = $c; 
      	#print "Beta Param = $Beta";
      }
      elsif ($a =~ m/Diffusion Coefficient/) { 
      	($c, $d) = split(" ", $b); 
      	$Diff_coef = $c; 
      	#print "Diffusion Coefficient = $Diff_coef\n";
      }
      elsif ($a =~ m/Fall Time Threshold/) { 
        ($c, $d) = split(" ", $b); 
        $Fall_tt = $c;
        #print "Fall Time Threshold = $Fall_tt\n";
      }
      elsif ($a =~ m/Total Mass Ejected/) { 
      	($c, $d) = split(" ", $b);
      	$Total_mass = $c; 
      	#print "Total Mass Ejected = $Total_mass\n";
      }
      elsif ($a =~ m/Max Particle Size/) { 
      	($c, $d) = split(" ", $b);
      	$Max_phi = $c; 
      	#print "Max Particle Size = $Max_phi\n";
      }
      elsif ($a =~ m/Min Particle Size/) { 
      	($c, $d) = split(" ", $b);
      	$Min_phi = $c; 
      	#print "Min Particle Size = $Min_phi\n";
      }
      elsif ($a =~ m/Median Size/) { 
      	($c, $d) = split(" ", $b);
      	$Median_phi = $c; 
      	#print "Median Size = $Median_phi\n";
      }
      elsif ($a =~ m/Std. Dev./) { 
      	($c, $d) = split(" ", $b);
      	$Std_phi = $c; 
      	#print "Std. Dev. = $Std_phi\n";
      }
    }
close(README);

open(RESULTS, ">>results") || die("Cannot open results:$!");
printf RESULTS "%.2f %.0f %.1f %.1f %.0f %.0f %.2f %.2f %.2g\n", $nrmse, $Col_Ht,  $Alpha,  $Beta,  $Diff_coef,  $Fall_tt, $Median_phi,  $Std_phi,  $Total_mass;
close RESULTS;

# Read old configuration file and update, writing out a new file
$ret = open(CONF, "<$ARGV[1]");
if ($ret) {
	@lines = <CONF>;
	close(CONF);
}
else {
	print stderr "No tephra2 configuration file to update!\n";
	exit;
}

open(CONF, ">$ARGV[1]") || die("Cannot open $ARGV[1]:$!");
foreach $line (@lines) {
	chomp $line;
	($key, $val) = split" ", $line;
  if ($key =~ m/PLUME_HEIGHT/) { print CONF "$key $Col_Ht\n"; }
  elsif ($key =~ m/ALPHA/) { print CONF "$key $Alpha\n"; }
  elsif ($key =~ m/BETA/) { print CONF "$key $Beta\n"; }
  elsif ($key =~ m/ERUPTION_MASS/) { print CONF "$key $Total_mass\n"; }
  elsif ($key =~ m/MAX_GRAINSIZE/) { print CONF "$key $Max_phi\n"; }
  elsif ($key =~ m/MIN_GRAINSIZE/) { print CONF "$key $Min_phi\n"; }
  elsif ($key =~ m/MEDIAN_GRAINSIZE/) { print CONF "$key $Median_phi\n"; }
  elsif ($key =~ m/STD_GRAINSIZE/) { print CONF "$key $Std_phi\n"; }
  elsif ($key =~ m/COL_STEPS/) { print CONF "$key $c_steps\n"; }
  elsif ($key =~ m/PART_STEPS/) { print CONF "$key $p_steps\n"; }
  elsif ($key =~ m/DIFFUSION_COEFFICIENT/) { print CONF "$key $Diff_coef\n"; }
  elsif ($key =~ m/FALL_TIME_THRESHOLD/) { print CONF "$key $Fall_tt\n"; }
  else { printf CONF "$line\n"; }
}
close(CONF);
