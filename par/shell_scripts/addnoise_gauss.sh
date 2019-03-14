#!/bin/csh

set x=1
while ( $x < 101)

suaddnoise sn=25 seed=34520 < su/CAES_spike_time_1_sg_80/DENISE_CAES_y.su.shot$x f=30,60 amps=1,0 > su/CAES_spike_time_1/DENISE_CAES_y.su.shot$x
suaddnoise sn=25 seed=87414 < su/CAES_spike_time_0_no_noise/DENISE_CAES_y.su.shot$x f=30,60 amps=1,0 > su/CAES_spike_time_0/DENISE_CAES_y.su.shot$x

set x = `expr $x + 1`

end

