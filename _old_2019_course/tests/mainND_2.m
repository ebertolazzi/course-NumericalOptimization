clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%r = Rosenbrock();
r = Rastrigin();

disp(r.arity());

r.contour([-1,1],[-1,1], 80);
axis equal;

%ls = 'Armijo';
ls = 'MT';
%ls = 'GS';
switch ls
case 'Armijo'
  search_method = LinesearchArmijo();
case 'MT'
  search_method = LinesearchMoreThuente();
case 'GS'
  search_method = LinesearchGoldenSection();
  search_method.setMaxIteration( int32(10) );
  search_method.setTolerance(1e-8);
end

minimization_method = MinimizationConjugateGradient( r, search_method );
minimization_method.selectByName('PRP');
minimization_method.save_iterate_on();
minimization_method.setTolerance(1e-12);

[xs,converged] = minimization_method.minimize( [0.5;0.8] );
minimization_method.plotIter();
xs
