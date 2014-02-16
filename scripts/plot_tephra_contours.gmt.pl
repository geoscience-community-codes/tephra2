#$out = "contours.eps";
#$in = "tephra2.out.xyz";
$out = "$ARGV[2]";
$in = "$ARGV[1]";
$samples = $ARGV[0];
$vent = "../colima_vent.xy";

$cellsize = 100;
$contours = "../scripts/contours2";
$west = 540000;
$east = 760000;
$south = 2120000;
$north = 2340000;

$west = 620000;
$east = 680000;
$south = 2150000;
$north = 2200000;

#system "nearneighbor --D_FORMAT=%f $in -Gtephra.grd -N8 -S2000 -I$cellsize -R$west/$east/$south/$north ";
system "surface --D_FORMAT=%f $in -Gtephra.grd -I$cellsize -R$west/$east/$south/$north -T.05 -V";
system "psbasemap --BASEMAP_TYPE=plain --ANNOT_FONT_PRIMARY=2  --ANNOT_FONT_SIZE=12 --D_FORMAT=%.0f -Y2i -X1i -Jx1:16000 -R -Ba10000:'':/a10000:'':/:.'':WSen -P -V -K > $out";

system "psxy $samples -Jx -R -Sc8p  -Gred -Wthinnest,0 -V -O -K >> $out";
system "grdcontour tephra.grd -C$contours -A10+f0+s8+kblack -Wthinner,0 -Gn2 -Jx -R -O -K >> $out";
system "psxy $vent -Jx -R -St20p  -Gyellow -Wthinnest,0 -V -O >> $out";
system "ps2raster $out -A -Tg";

#$out = "contours_more.eps";
$out = "$out.more.eps";
$west = 600000;
$east = 700000;
$south = 2140000;
$north = 2280000;
system "surface --D_FORMAT=%f $in -Gtephra.grd -I$cellsize -R$west/$east/$south/$north -T.05 -V";
system "psbasemap --BASEMAP_TYPE=plain --ANNOT_FONT_PRIMARY=2 --ANNOT_FONT_SIZE=12 --D_FORMAT=%.0f -Y2i -X1i -Jx1:26000 -R -Ba20000:'':/a20000:'':/:.'':WSen -P -V -K > $out";

system "psxy $samples -Jx -R -Sc8p  -Gred -Wthinnest,0 -V -O -K >> $out";
system "grdcontour tephra.grd -C$contours -A10+f0+s8+kblack -Wthinner,0 -Gn2 -Jx -R -O -K >> $out";
system "psxy $vent -Jx -R -St20p  -Gyellow -Wthinnest,0 -V -O >> $out";
system "ps2raster $out -A -Tg";
