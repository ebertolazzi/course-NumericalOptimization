clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

r = Barrier1();
disp(r.arity());

search_method   = LinesearchGoldenSection();
%%search_method = LinesearchArmijo();
%search_method = LinesearchWolfe();

%search_method.debug_on();

dir_method = MinimizationGradientMethod( r, search_method );
%dir_method = MinimizationBFGS(r,search_method);

dir_method.setMaxIteration( int32(10) );
dir_method.setTolerance(1e-6);
dir_method.debug_on();

x0 = r.guess(int32(1));
[xs,converged] = dir_method.minimize( x0 );
r.contour([-1.5 1.5],[-1.5 1.5],@(z) log(1+z), 80)
axis equal;
dir_method.plotIter();

xs
