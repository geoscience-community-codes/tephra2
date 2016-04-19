#$out = "contours.eps";
#$in = "tephra2.out.xyz";
$samples = $ARGV[0];
$in = $ARGV[1];
$out = $ARGV[2];
$vent = $ARGV[3];
$cellsize = 1000;
$contour = "contours2";
@contours = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80);
open(CONT, ">", "$contour") or die ("cannot open $contour: $!");

foreach my $c (@contours) {
if ($c%10 == 0) {print CONT "$c a\n";}
else {print CONT "$c c\n";}
}
close CONT;
$west= 400000;
$east = 600000;
$south = 4350000;
$north = 4500000;
`gmt surface $in \`gmtinfo -I- $in \` -Gtephra.grd -I$cellsize -V`;
`gmt psbasemap --MAP_FRAME_TYPE=plain --FONT_ANNOT_PRIMARY=8p,2,0 --FONT_LABEL=10p,2,0 --FORMAT_FLOAT_OUT=%.0f -Y2i -X1i -Jx1:1200000 -R$west/$east/$south/$north -Bxa50000+l'Easting (m)' -Bya50000+l'Northing (m)' -BWS -P -K -V > $out`;
#system "psxy $samples -Jx -R -Sp -Gred -Wgray -O -K -V >> $out";
#system "grdcontour tephra.grd -C$contour -A10+f0+s8+kblack -Wthinnest,0 -Gn2 -Jx -R -O -K >> $out";
`gmt grdcontour tephra.grd -C$contour -Q25 -Wa0.75p,0 -Wc.25p,150 -Gn1 -Jx -R -O -K -V >> $out`;
`gmt psxy $samples -Jx -R -Sc5p -Gred -Wthinnest,0 -O -K >> $out`;
`gmt psxy $vent -Jx -R -St20p -Gyellow -Wthinnest,0 -O -V >> $out`;
`gmt ps2raster $out -A -Tg -V`;
`rm $contour`;

