# this file contains the standard methods
#for every 1D functions

# return the value of the functions computed in x
function y(funct, x)
    return funct(x)
end

# return the value of the first derivative computed in x
function Dy(funct, x)
    h = (max(1,abs(x))*(eps()^(1/3)))[1]
    fp = funct(x+h)
    fm = funct(x-h)
    return (fp - fm)/2*h
end

# return the value of the first derivative computed in x
function DDy(funct, x)
    h = max(1 , abs( x )) * (eps()^(1/3))
    fp = funct( x+h )
    fm = funct( x-h )
    fc = funct( x )
    return ( fp + fm -2*fc ) / h^2
end
