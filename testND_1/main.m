clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

r = Rosenbrock();

disp(r.arity());

x     = linspace(-1,3,1000);
y     = linspace(-1,2,1000);
[X,Y] = meshgrid(x,y);
%
% X and Y are matrices, to evaluate in [X(i,j);Y(i,j)]
% vector
%
XY        = zeros(2,size(X,1),size(X,2)) ;
XY(1,:,:) = X ;
XY(2,:,:) = Y ;
Z         = r.eval(XY);
contour(X,Y,log(1+Z),100);
axis equal ;

%search_method   = GoldenSearch(1e-5,100);
search_method   = LinesearchArmijo();
gradient_method = MinimizationGradientMethod(search_method,1e-6,10000,true);

gradient_method.setFunction( r ) ;

[xs,converged] = gradient_method.minimize( [1;2] ) ;

xs

