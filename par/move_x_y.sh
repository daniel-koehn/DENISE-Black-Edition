#!/bin/csh

set x=1
while ( $x < 98)

mv su/DENISE_MARMOUSI_x.su.shot$x.it1 su/MARMOUSI_spike/DENISE_MARMOUSI_x.su.shot$x
mv su/DENISE_MARMOUSI_y.su.shot$x.it1 su/MARMOUSI_spike/DENISE_MARMOUSI_y.su.shot$x

set x = `expr $x + 1`

end
