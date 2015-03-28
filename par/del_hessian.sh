#!/bin/csh
set NT=124

set x=1
while ( $x < 41 )
  
  rm jacobian/jacobian_Test_hessian_$x.*
  rm jacobian/jacobian_Test_hessian_u_$x.*
  rm jacobian/jacobian_Test_hessian_rho_$x.*
  
  echo "delete distributed hessian for frequency ... $x"
  set x = `expr $x + 1`
end
