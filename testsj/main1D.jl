workspace()

#include("SinFun.jl")
#include("CosFun.jl")
include("FunMia.jl")
#include("TestFun.jl")
include("../libj/LinesearchGoldenSection.jl")

#f = SinFun()
#f = CosFun()
f = FunMia()
#f = TestFun()

#minfunct = f

#a,b = minimize(minfunct, [π,2*π])
a,b = minimize(f, [2,6])
fa = f(a) ;
print("function evaluated in ",a , " = ", fa, "\n")
fb = f(b) ;
print("function evaluated in ",b , " = ", fb, "\n")
