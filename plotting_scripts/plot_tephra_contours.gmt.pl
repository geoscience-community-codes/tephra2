#$out = "contours.eps";
#$in = "tephra2.out.xyz";
$samples = $ARGV[0];
$in = $ARGV[1];
$out = $ARGV[2];
$vent = $ARGV[3];
$cellsize = 1000;
$contour = "contours2";
@contours = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 10, 20, 30, 40, 50);
open(CONT, ">", "$contour") or die ("cannot open $contour: $!");
foreach my $c (@contours) {
print CONT "$c a\n";
}
close CONT;
$west= 660000;
$east = 740000;
$south = 1577000;
$north = 1627000;
system "surface $in -R$west/$east/$south/$north -Gtephra.grd -I$cellsize -V";
system "psbasemap --BASEMAP_TYPE=plain --ANNOT_FONT_PRIMARY=2 --ANNOT_FONT_SIZE=8 --D_FORMAT=%.0f -Y2i -X1i -Jx1:500000 -R -Ba10000:'':/a10000:'':/:.'':WSen -P -K -V > $out";
#system "psxy $samples -Jx -R -Sp -Gred -Wgray -O -K -V >> $out";
system "psxy $samples -Jx -R -Sc5p -Gred -Wthinnest,0 -O -K >> $out";
#system "grdcontour tephra.grd -C$contour -A10+f0+s8+kblack -Wthinnest,0 -Gn2 -Jx -R -O -K >> $out";
system "grdcontour tephra.grd -C0.25 -A2+f6+r6+s6 -Q25 -Wa0.75p,0 -Wc.25p,150 -Gn1 -Jx -R -O -K -V >> $out";
system "psxy $vent -Jx -R -St20p -Gyellow -Wthinnest,0 -O -V >> $out";
system "ps2raster $out -A -Tg -V";
system "rm $contour";

    Status
    API
    Training
    Shop
    Blog
    About

    © 2015 GitHub, Inc.
    Terms
    Privacy
    Security
    Contact

