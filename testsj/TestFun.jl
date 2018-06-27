
#import Function1D
include("../libj/Function1D.jl")
function TestFun()
    return x -> exp(-3*x)+x^2*sin(3*x)
end

#y(fun,1.57)
