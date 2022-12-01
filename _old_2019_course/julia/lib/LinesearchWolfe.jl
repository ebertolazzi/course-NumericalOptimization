include("ForwardBackward2.jl")

function strongWolfe_on( )
  strongWolfe = true
end

function strongWolfe_off( )
  strongWolfe = false
end

function  search( funct, x0, d, alpha_guess )
  # find step that satisfy Armijo condition
  # of ok = false search failed
  aLO,aHI,fLO,fHI,ierr = ForwardBackward(funct, x0, d, alpha_guess )
  if ierr == 0 # only alpha0 satisfy Armijo, check if minimum is on [0,aLO]
    if eval_D( funct, x0, d, aLO) > 0
      alpha_ott, ok = Zoom( funct, x0, d, 0, eval_1Dcut( funct, x0, d, 0.0 ), aLO, fLO, strongWolfe )
    else
      alpha_ott, ok = Zoom( funct, x0, d, aLO, fLO, aHI, fHI, strongWolfe )
    end
  elseif ierr == 1 #
     alpha_ott, ok = Zoom( funct, x0, d, aLO, fLO, aHI, fHI, strongWolfe )
  elseif ierr == 2 # LO and HI exchanged
     alpha_ott, ok = Zoom( funct, x0, d, aHI, fHI, aLO, fLO, strongWolfe )
  else
    ok        = false
    alpha_ott = 0
  end
end
