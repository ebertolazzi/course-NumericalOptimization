clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%r = Rosenbrock();
r = Rastrigin();

disp(r.arity());

r.contour([-1,1],[-1,1],@(z) log(1+z), 80)
axis equal ;

%search_method   = GoldenSearch();
search_method   = LinesearchArmijo();
gradient_method = MinimizationGradientMethod(search_method,1e-6,10000,true);

search_method.setDebug();
gradient_method.no_FD_D();
gradient_method.setFunction( r ) ;


[xs,converged] = gradient_method.minimize( [0.5;0.8] ) ;

xs


