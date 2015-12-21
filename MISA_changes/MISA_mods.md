#Modifications to MISA scripts.

The following modifications to the [scripts provided by the MISA developers](http://pgrc.ipk-gatersleben.de/misa/primer3.html) were made for this study.

## p3_in.pl

These changes were for compatability with primer3 versions 2 and newer.

- Line 28, change to:
>  print OUT "SEQUENCE_ID=$id"."_$ssr_nr\nSEQUENCE_TEMPLATE=$seq\n";

- Line 30, change to:
> print OUT "SEQUENCE_TARGET=",$start-3,",",$size+6,"\n";

## p3_out.pl
Several changes were made to p3_out.pl and the new file is provided as ps_out.MAG.pl.
  1. Modified tags for compatability with primer3 versions 2 and newer 
  2. Removed printing of failr prirmer desing attempts
  3. Add output a file of names of scaffolds with successful primer developmet to be used in other scripts developed here.
