#!/bin/bash

suaddhead < jacobian/jacobian_Test_P_image ns=174 > jacobian/jacobian_Test_P_image.su
suaddhead < jacobian/jacobian_Test_P_image_shot_$1 ns=174 > jacobian/jacobian_Test_P_image_shot_.su

sudiff jacobian/jacobian_Test_P_image.su jacobian/jacobian_Test_P_image_shot_.su| suximage


