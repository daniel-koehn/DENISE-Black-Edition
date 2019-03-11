#!/bin/bash

# output of gx-header values for first shot
sugethw < su/DENISE_wagrien_y.su.shot1.it1 output=geom key=gx > su/hdrfile

x=1
while [ $x -le 400 ]
do 

	y=$(awk 'NR=='$x'' su/hdrfile)

	# extract traces at shot positions
	suwind < su/DENISE_wagrien_y.su.shot$x.it1 key=gx min=$y max=$y > su/ZO/tmp_$x.su
        (( x++ ))

done

# merge individual traces
cat su/ZO/tmp_* > su/ZO/tmp_cat.su
susort < su/ZO/tmp_cat.su > su/ZO/zero_offset.su gx

rm su/ZO/tmp*
rm su/hdrfile

