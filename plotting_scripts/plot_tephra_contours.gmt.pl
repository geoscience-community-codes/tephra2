#$out = "contours.eps";
#$in = "tephra2.out.xyz";
$out = "$ARGV[2]";
$in = "$ARGV[1]";
$samples = $ARGV[0];
$vent = $ARGV[3];

$cellsize = 100;
$contour = "contours2";
$west = 540000;
$east = 760000;
$south = 2120000;
$north = 2340000;

$west = 620000;
$east = 680000;
$south = 2150000;
$north = 2200000;

@contours = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 10, 20, 30, 40, 50);

open(CONT, ">", "$contour") or die ("cannot open $contour: $!");
foreach my $c (@contours) {
	print CONT "$c a\n";
}
close CONT;
#system "nearneighbor --D_FORMAT=%f $in -Gtephra.grd -N8 -S2000 -I$cellsize -R$west/$east/$south/$north ";
system "surface --D_FORMAT=%g $in -Gtephra.grd -I$cellsize -R$west/$east/$south/$north";
system "psbasemap --BASEMAP_TYPE=plain --ANNOT_FONT_PRIMARY=2  --ANNOT_FONT_SIZE=12 --D_FORMAT=%.0f -Y2i -X1i -Jx1:350000 -R -Ba10000:'':/a10000:'':/:.'':WSen -P -K > $out";

system "psxy $samples -Jx -R -Sc8p  -Gred -Wthinnest,0 -O -K >> $out";
system "grdcontour tephra.grd -C$contour -A+f0+s8+kblack -Wthinner,0 -Gn2 -Jx -O -K >> $out";
system "psxy $vent -Jx -R -St2p  -Gyellow -Wthinnest,0 -O >> $out";
system "ps2raster $out -A -Tg";

#$out = "contours_more.eps";
$out = "$out.more.eps";
$west = 600000;
$east = 700000;
$south = 2140000;
$north = 2280000;
system "surface $in -Gtephra.grd -I$cellsize -R$west/$east/$south/$north -V";
system "psbasemap --BASEMAP_TYPE=plain --ANNOT_FONT_PRIMARY=2 --ANNOT_FONT_SIZE=10 --D_FORMAT=%.0f -Y2i -X1i -Jx1:600000 -R -Ba20000:'':/a20000:'':/:.'':WSen -P -K > $out";

system "psxy $samples -Jx -R -Sc8p  -Gred -Wthinnest,0 -O -K >> $out";
#system "grdcontour tephra.grd -C$contour -A10+f0+s8+kblack -Wthinnest,0 -Gn2 -Jx -R -O -K >> $out";
system "grdcontour tephra.grd -C0.25 -A2+f6+r6+s6 -Q25 -Wa0.75p,0 -Wc.25p,150 -Gn1 -Jx -R -O -K -V >> $out";
system "psxy $vent -Jx -R -St20p  -Gyellow -Wthinnest,0 -O >> $out";
system "ps2raster $out -A -Tg";
system "rm $contour";
