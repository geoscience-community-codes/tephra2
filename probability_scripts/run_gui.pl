#!/usr/bin/perl
# Edit the pathname for the variable: our $Tephra2_exe

use strict;
use Tk;
use Tk::PNG;
use Image::Resize;
#use Image::Imlib2;
use Tk::MsgBox;
use Tk::Pane;
use Tk::ProgressBar;
use Tk::PathEntry;
#use Tk::BrowseEntry;
use Tk::LogScale;
use Data::Dumper;
use Cwd;

###########################################################################################
# Make sure that the location of the tephra2 executable is correctly specified below!
############################################################################################
our $Tephra2_exe = "../tephra2-2012";
our $run_line = "perl survivor_function.pl $Tephra2_exe ";
#our $num_runs = 1000;
#our $wind_dir = "wind_db";
#our $zone = 13;

our @colors = (  0, '#ff002a',  1, '#ff0014',  2, '#ff000a',  3, '#ff0500',  4, '#ff1000',
	     5, '#ff1b00',  6, '#ff3000',  7, '#ff3b00',  8, '#ff4600',  9, '#ff5100',
	    10, '#ff6100', 11, '#ff7600', 12, '#ff8100', 13, '#ff8c00', 14, '#ff9700',
	    15, '#ffa100', 16, '#ffbc00', 17, '#ffc700', 18, '#ffd200', 19, '#ffdd00',
	    20, '#ffe700', 21, '#fffd00', 22, '#f0ff00', 23, '#e5ff00', 24, '#dbff00',
	    25, '#d0ff00', 26, '#baff00', 27, '#afff00', 28, '#9fff00', 29, '#95ff00',
	    30, '#8aff00', 31, '#74ff00', 32, '#6aff00', 33, '#5fff00', 34, '#54ff00',
	    35, '#44ff00', 36, '#2eff00', 37, '#24ff00', 38, '#19ff00', 39, '#0eff00',
	    40, '#03ff00', 41, '#00ff17', 42, '#00ff21', 43, '#00ff2c', 44, '#00ff37',
	    45, '#00ff42', 46, '#00ff57', 47, '#00ff67', 48, '#00ff72', 49, '#00ff7d',
	    50, '#00ff87', 51, '#00ff9d', 52, '#00ffa8', 53, '#00ffb8', 54, '#00ffc3',
	    55, '#00ffcd', 56, '#00ffe3', 57, '#00ffee', 58, '#00fff8', 59, '#00faff',
	    60, '#00eaff', 61, '#00d4ff', 62, '#00c9ff', 63, '#00bfff', 64, '#00b4ff',
	    65, '#00a9ff', 66, '#008eff', 67, '#0083ff', 68, '#0079ff', 69, '#006eff',
	    70, '#0063ff', 71, '#004eff', 72, '#003eff', 73, '#0033ff', 74, '#0028ff',
	    75, '#001dff', 76, '#0008ff', 77, '#0200ff', 78, '#1200ff', 79, '#1d00ff',
	    80, '#2800ff', 81, '#3d00ff', 82, '#4800ff', 83, '#5300ff', 84, '#5d00ff',
	    85, '#6e00ff', 86, '#8300ff', 87, '#8e00ff', 88, '#9900ff', 89, '#a300ff',
	    90, '#ae00ff', 91, '#c900ff', 92, '#d400ff', 93, '#df00ff', 94, '#e900ff',
	    95, '#f400ff', 96, '#ff00f3', 97, '#ff00e3', 98, '#ff00d9', 99, '#ff00ce' );
	    
# Creation of Widgets ---------------------------------------------------------	    
	    
our $mw = MainWindow->new(-title=>"Survivor Function Inputs"); # Main Window
our $f = $mw->Scrolled("Frame", -scrollbars=>'ose', -width  => 500, -height => 750) ;
#$f->gridPropagate(0);
$mw->resizable( 1, 1 );
$f->form();

my @Step;
my @Labl;
my $s=0;
my $i=0;

my $wlabel = $f->Label(-text=>"Input Parameters") ->form(-top=>'%0', -right=>'%55');
$Step[$s] = $f->Label(-text=>"Step 1:", -fg=>'red')->form(-top =>$wlabel);
$Labl[$s] = $f->Label(-text=>"Choose a volcano:")->form(-top=>$Step[$s]);
# Open the volcano database for list of volcanos
my $DIR = "data/volcanos/";
opendir(DIR, $DIR) || die "Not happening: [$DIR] $!";
my @CONTENTS = readdir(DIR);
our @Volcano;
#our @Volcano = qw(Cerro_Negro Colima Concepcion Masaya Momotombo Nejapa San_Cristobal Telica);
my $v = 0;
foreach (@CONTENTS) {
	unless ($_ =~ /\./) {
		$Volcano[$v++] = $_;
	}
}
our $wvolcano = $f->Scrolled("Listbox", -scrollbars=>'oe', -selectmode=> "single" );
$wvolcano->bind('<ButtonRelease>' => \&selectVolcano);
#opendir DIR, "data/volcanos/";
#$wvolcano->insert('end', grep { -f $_ && -r $_ } readdir DIR);
#close DIR;
$wvolcano->form(-top => $Labl[$s]);
foreach my $item (@Volcano) {
	$wvolcano->insert('end', $item);
}


$s++;
$Step[$s] = $f->Label(-text=>"Step 2:", -fg=>'red')->form(-top=>$wvolcano);
$Labl[$s] = $f->Label(-text=>"Complete pathname to wind database:")->form(-top=>$Step[$s]);
our $path = cwd();
our $wpath = $f->Scrolled("PathEntry", -scrollbars=>'os', -width=>50,  -textvariable=>\$path, -highlightcolor=>'black', -background=>'white' );
$wpath->configure(-state => 'disabled');
$wpath->form(-top=>$Labl[$s], -left=>'%0', -right=>'%100');

$s++;
our @Cities;
$Step[$s] = $f->Label(-text=>"Step 3:", -fg=>'red')->form(-top=>$wpath);
$Labl[$s] = $f->Label(-text=>"Choose a city:")->form(-top=>$Step[$s]);
our $wcities = $f->Scrolled("Listbox", -scrollbars=>'oe', -selectmode=> "single");
$wcities->bind('<ButtonRelease>' => \&whichCity);
$wcities->form(-top=>$Labl[$s]);

$s++;
$Step[$s] = $f->Label(-text=>"Step 4:", -fg=>'red')->form(-top=>$wcities);
$Labl[$s] = $f->Label(-text=>"Choose one or more eruption styles:")->form(-top=>$Step[$s]);

our @cb;
our @vei = qw(vei2 vei3 vei4 vei5 vei6 vei7 fixed); # check all
$i=0;
our @ck_vei;
foreach (@vei) {
	my $tex = "";
	if ($i eq 0) { $tex .= "VEI-2    Mass Range: 1e9 -- 1e10 kg       Column Range:  1 -- 5 km"; }
	elsif ($i eq 1) { $tex .= "VEI-3    Mass Range: 1e10 -- 1e11 kg     Column Range:  3 -- 15 km"; }
	elsif ($i eq 2) { $tex .= "VEI-4    Mass Range: 1e11 -- 1e12 kg     Column Range:  10 -- 25 km"; }
	elsif ($i eq 3) { $tex .= "VEI-5    Mass Range: 1e12 -- 1e13 kg     Column Range:  20 -- 35 km"; }
	elsif ($i eq 4) { $tex .= "VEI-6    Mass Range: 1e13 -- 1e14 kg     Column Range:  25 -- 40 km"; }
	elsif ($i eq 5) { $tex .= "VEI-7    Mass Range: 1e14 -- 1e15 kg     Column Range:  30 -- 50 km"; }
	elsif ($i eq 6) {$tex .= "Fixed Mass and Column Height";}
	$cb[$i] = $f->Checkbutton(-text => $tex, -variable =>\$ck_vei[$i] );
	if ($i == 6) { $cb[6]->bind('<ButtonRelease>' => \&fixed); }
	($i > 0) ? $cb[$i]->form(-top=>$cb[$i-1]) : $cb[$i]->form(-top=>$Labl[$s]);
		$i++;
}

my $Labl_col = $f->Label(-text=>"Eruption Mass as log(Kg):")->form(-top=>$cb[$i-1], -left=>'%10');
my $Labl_ht = $f->Label(-text=>"Column Height (km):")->form(-top=>$cb[$i-1], -left=>'%60');
our $emass = $f->Scale(-orient=>'horizontal', -from=>'9', -to=>'15', -resolution=>1);
$emass->set(1);
$emass->configure(-state => 'disabled');
$emass->form(-top=>$Labl_col, -left=>'%10');
our $cht = $f->Scale(-orient=>'horizontal', -from=>'1', -to=>'50', -resolution=>1);
$cht->set(1);
$cht->configure(-state => 'disabled');
$cht->form(-top=>$Labl_col, -left=>'%60');

$s++;
$Step[$s] = $f->Label(-text=>"Step 5:", -fg=>'red')->form(-top=>$emass);
$Labl[$s] = $f->Label(-text=>"Number of iterations:")->form(-top=>$Step[$s]);
our $wruns = $f->Scale(-orient=>'horizontal', -from=>'100', -to=>'10000', -resolution=>100);
$wruns->set(100);
$wruns->form(-top=>$Labl[$s], -left=>'%0', -right=>'%100');


our $wrun_line = $f->Text(-width=>50, -height=>5, -exportselection=>1, -highlightcolor=>'black', -background=>'white', -wrap=>'char');
$wrun_line->form(-top=>$wruns, -left=>'%0', -right=>'%100');


our $but_run = $f->Button(-text => "Run", -command =>\&push_button);  
$but_run->form(-top=>$wrun_line);
our $but_image = $f->Button(-text=>'Show image', -command=>\&showImage);
$but_image->form(-top=>$wrun_line, -left=>'%38');
$but_image->configure(-state => 'disabled');
our $but_quit = $f->Button(-text => "Quit", -command => \&exitProgram);
$but_quit->form(-top=>$wrun_line, -right=>'%100');
#my $but_quit = ->MsgBox(-title => "Title", -type => "okcancel");
#$but_run->grid($but_quit);   

our $progress = $f->ProgressBar(-borderwidth=>2, -blocks=>100, -gap=>0,
			  -troughcolor=>'white',-colors=>\@colors,
			  -length=>106 );
			  
$progress->form(-top=>$but_quit, -left=>'%0',  -right=>'%100');

my $stuff = $wrun_line->Contents();
print STDERR "$stuff\n";
$stuff = $wpath->get();
print STDERR "$stuff\n"; 

MainLoop;

#This subroutine is executed when a volcano is selected
sub selectVolcano {
	
	my $text = "";
	(my $volc) = $wvolcano->curselection;
	$wvolcano->activate($volc);
	$text = "Volcano selected: $Volcano[$volc]\n";
	print STDERR $text;
	$wrun_line -> insert('end',$text);
	$wpath->configure(-state => 'normal');
	$path = cwd();
	$path .= "/data/volcanos/$Volcano[$volc]/wind_db/";
	print STDERR "$path\n";
	
	$wcities->delete(0, 'end');
	# Open the volcano database for list of volcanos
my $DIR = "data/volcanos/$Volcano[$volc]/cities";
opendir(DIR, $DIR) || die "Not happening: [$DIR] $!";
my @CONTENTS = readdir(DIR);
our @Cities = [] ;
my $v = 0;
foreach (@CONTENTS) {
	unless ($_ =~ /\./) {
		$Cities[$v++] = $_;
	}
}
#	 if ($Volcano[$volc] =~ m/Colima/) {
#	   @Cities = qw(Aguascalientes Ameca Atotonilco Colima Ciudad_Guizman Cuquio Encarnacion_de_Diaz Fresnillo Guadalajara Guanajuato La_Barca Lagos_de_Moreno Manzanillo Mazamitla Poncitlan Saltillo San_Luis_Potosi Sayula Teocaltiche Tizapan_El_Alto Tonila Yahualica Yurecuaro Zacatecas Zacoalco);
#	 }
#	 else {
#	 @Cities = qw(Blue_Fields Chichigalpa Chinandega Corinto Esteli Granada Leon Managua Masaya Matagalpa Nagarote Rivas);
#	 
#	 }
 	my $num = @Cities;
 	
	foreach my $item (@Cities) {
		$wcities->insert('end', $item);
	}
	
}

sub whichCity {

	my $text;
	(my $city) = $wcities->curselection;
	$wcities->activate($city);
	$text = "City chosen: $Cities[$city]\n";
	print STDERR $text;
	$wrun_line -> insert('end',$text);
}

sub fixed {
	print STDERR "Fixed wind and Column Height\n";
	$emass->configure(-state => 'normal');
	$cht->configure(-state => 'normal');
}

sub showImage {
	print STDERR "Open image: survivor.png\n";
#	`xv survivor.png`;
		my $iw = new MainWindow(-title=>"Survivor Function Output");
		my $image = Image::Resize->new('survivor.png');		
		my $ht = $image->height();		
		my $wd = $image->width();
		$ht *= .5;
		$wd *= .5;
		my $resize = $image->resize($wd, $ht);
		open (FILE, ">surviv.png");
		print FILE $resize->png();
		close FILE;
      my $canvar = $iw ->Canvas(-width  => $wd, -height => $ht);
      $canvar->pack;
      my $img =
       $canvar->Photo( 'IMG',
                       -file =>"surviv.png" );

      $canvar->create( 'image',0,0,
                       '-anchor' => 'nw',
                       '-image'  => $img );

	$but_run->configure(-state => 'active');
#	$but_image->configure(-state => 'disabled');
}

#This is executed when the Quit button is pressed ------------------------------------------
sub exitProgram {

    my $d = $mw->MsgBox(-title => "Quit?", -type => "okcancel", -message=>"Ok to close?" );
    my $ans = $d->Show();
    print STDERR "$ans\n";
    if ($ans =~ /ok/) { exit; }
}

#This is executed when the Run button is pressed ------------------------------------------
sub push_button {
	
	my $zone = 0;
	my $text = "";
	
	(my $volc) = $wvolcano->index("active");
	$text = "Chosen volcano: $volc $Volcano[$volc]\n";
	print STDERR $text;
	$wrun_line -> insert('end',$text);
	my $file = "data/volcanos/$Volcano[$volc]/location.txt";
	open FILE, "<$file" or die "Unable to open $file: $!";
	my $line = <FILE>;
	(my $lon, my $lat, my $el, $zone) = split " ", $line; 
	$run_line .= "$file ";
	close FILE;
	 
	(my $city) = $wcities->index("active");
	print STDERR "Chosen city: $city $Cities[$city]\n";
	$run_line .= "data/volcanos/$Volcano[$volc]/cities/$Cities[$city] ";
	
	my $wind_dir = $wpath->get;
	print STDERR "Wind database: $wind_dir\n";
	$run_line .= "$wind_dir ";
	
	my $n = $wruns->get;
	print STDERR "Number of runs: $n\n";
	$run_line .= "$n ";
	
	#my $zone = $wzone->get;
	$run_line .= "$zone ";
	print STDERR "UTM zone: $zone\n";
	
	for (my $i = 0; $i < (@cb); $i++) {
		if ($ck_vei[$i]) {
			print STDERR "Chosen eruption style: $vei[$i] \n";
			$run_line .= "$vei[$i] ";
			if ($vei[$i] =~ m/fixed/) {
				my $mass = $emass->get;
				$run_line .= "$mass ";
				my $col = $cht->get;
				$run_line .= "$col ";
			}
		}
	}
	$run_line .= "2> survivor.out";

	$wrun_line -> insert('end',"Running: $run_line\n");
	print STDERR "Running: $run_line\n";
	`$run_line`;
	
	$but_run->configure(-state => 'disabled');
	$but_image->configure(-state => 'active');
	$run_line = "";
	
}

sub reset {

	
}
