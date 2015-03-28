#!/bin/csh

set x=1
while ( $x < 101)

#sudiff su/MARMOUSI_10Hz_time_1/DENISE_MARMOUSI_x.su.shot$x su/MARMOUSI_10Hz_time_0/DENISE_MARMOUSI_x.su.shot$x > su/MARMOUSI_10Hz_data_diff/DENISE_MARMOUSI_x.su.shot$x
sudiff su/CAES_spike_time_1/DENISE_CAES_y.su.shot$x su/CAES_spike_time_0/DENISE_CAES_y.su.shot$x > su/CAES_spike_data_diff/DENISE_CAES_y.su.shot$x

set x = `expr $x + 1`

end

