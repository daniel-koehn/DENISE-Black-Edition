#!/bin/csh

set nshots=14
set x=1
while ( $x < $nshots)

supswigp label1="Time [s]"  label2="Trace No." title="Shot no. $x" < DENISE_mod_y.su.shot$x.it1 > test_$x.ps

set x = `expr $x + 1`

end

# merge all shots in one PS-file
cat test_1.ps > allshots.ps 

set x=2
while ( $x < $nshots)

cat test_$x.ps >> allshots.ps

set x = `expr $x + 1`

end

# convert PS->PDF
ps2pdf allshots.ps

# tidy up the place
rm *.ps
