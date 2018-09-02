
#import Function1D
include("../libj/Function1D.jl")
function SinFun()
    return x -> sin(x)
end

#y(fun,1.57)

#Dy(fun,1.57)

#DDy(fun,1.57)
@assert y(SinFun(),pi/2) == 1
@assert y(SinFun(),0) == 0
@assert Dy(SinFun(),pi/2) == 0
@assert DDy(SinFun(),0) == 0

#print(y(fun, [0:0.2:3]))
