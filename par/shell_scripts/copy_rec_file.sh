#!/bin/bash

x=1
while [ $x -le 100 ]
do

    cp ../receiver/receiver_OBC.dat ../receiver/receiver_OBC_shot_$x.dat
    (( x++ ))

done
