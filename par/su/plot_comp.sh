#!/bin/csh

set nshots=2
set x=1
while ( $x < $nshots)

# set tracl header
sushw < DENISE_tastrup_y.su.shot$x.it1 > tmp2.su key=tracl a=1 b=1

# extract data 
suwind < tmp2.su > tmp.su key=tracl j=10
suwind < tastrup/DENISE_tastrup_y.su.shot$x > tmp1.su key=tracl j=10

# plot field and model data
sugain qbal=1 < tmp.su | supswigp label1="Time [s]" label2="Trace No." title="Shot no. $x" fill=0 tracecolor=blue > true_$x.ps
sugain qbal=1 < tmp1.su | supswigp label1="Time [s]" label2="Trace No." title="Shot no. $x" fill=0 tracecolor=red > mod_$x.ps

psmerge in=true_$x.ps in=mod_$x.ps > test_$x.ps

#psmerge in=true_$x.ps in=mod_$x.ps > true_mod_$x.ps

# calculate and plot data residuals
#sudiff cachtice/profil_1/cachtice_data_FWI_Prof1_y.su.shot$x DENISE_mod_y.su.shot$x.it1 > diff.su

#supswigp label1="Time [s]"  label2="Trace No." title="Data Residuals shot no. $x" < diff.su > res_$x.ps

# merge field/model data comparison and residuals into one PS-file 
#merge2 true_mod_$x.ps res_$x.ps > test_$x.ps

set x = `expr $x + 1`

end

# merge all shots in one PS-file
cat test_1.ps > comp_shots.ps 

set x=2
while ( $x < $nshots)

cat test_$x.ps >> comp_shots.ps

set x = `expr $x + 1`

end

# convert PS->PDF
ps2pdf comp_shots.ps

# tidy up the place
rm *.ps
rm tmp.su
rm tmp1.su
rm tmp2.su
#rm diff.su
