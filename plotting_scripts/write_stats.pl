my $wind = $ARGV[1];

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

#open(FIT, "<$ARGV[2]") || die ("Cannot open $ARGV[2]:$!");
#$nrmse = <FIT>;
#close FIT;
#chomp $nrmse;
#($nrmse) = split("\n"); 

$ct = 0;
open(README, "<$ARGV[0]") || die("Cannot open $ARGV[0]:$!");

LINE:	while(<README>) {
   ($none, $fit_temp) = split "=", $_; 
if ($none =~ m/FIT/) {
	$fit = $fit_temp;
	chomp($fit_temp);
}
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
printf RESULTS "%.0f %s %.0f %.1f %.1f %.0f %.0f %.2f %.2f %.2g\n", $fit, $wind, $Col_Ht,  $Alpha,  $Beta,  $Diff_coef,  $Fall_tt, $Median_phi,  $Std_phi,  $Total_mass;
close RESULTS;


