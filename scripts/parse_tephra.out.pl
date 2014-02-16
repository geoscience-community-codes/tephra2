while (<>) {
  ($east, $north, $elev, $mass) = split " ", $_;
  		$thickness = $mass/10.43;
                #$thickness = $mass;
     printf "$east $north $thickness\n";
}
