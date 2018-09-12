#Armijo linesearch
include("ForwardBackward2.jl")


function  search( funct, x0, d, alpha_guess )
  # find step that satisfy Armijo condition
  # if ok = false search failed
  alpha0,alpha1,nn,mm,ierr = ForwardBackward(funct, x0, d, alpha_guess )
  ok = true
  alpha_ott = false
  if ierr == 0 # only alpha0 satisfy Armijo
    alpha_ott = alpha0
  elseif ierr== 1 || ierr == 2 # only alpha0 and alpha1 satisfy Armijo, take average
    alpha_ott = (alpha0+alpha1)/2
  else
    ok = false
  end
  return alpha_ott,ok
end
