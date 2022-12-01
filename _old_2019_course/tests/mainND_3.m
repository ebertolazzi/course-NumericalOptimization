clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');
addpath('../maps');

fun_name = 'Rosen';

switch fun_name
case 'Quad'
  r = Quadratic2D();
  disp(r.arity());
  r.contour([-1.5 1.5],[-1.5 1.5], 80)
  axis equal;
  x0 = [-1; 1 ];
case 'Rosen'
  r = Rosenbrock();
  disp(r.arity());
  r.contour([-1.5 1.5],[-0.4 1.5], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'Bro'
  r = Brown_bsf();
  disp(r.arity());
  r.contour( [10^6 -2000,10^6 + 2000],[-1000,1000], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'SchafferF6'
  r = SchafferF6();
  disp(r.arity());
  r.contour([-20,20],[-20,20], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'Boha1'
  r = BohachevskyN1();
  disp(r.arity());
  r.contour([-4 4],[-4 4], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'Boha2'
  r = BohachevskyN2();
  disp(r.arity());
  r.contour([-4 4],[-4 4], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'Boha3'
  r = BohachevskyN3();
  disp(r.arity());
  r.contour([-4 4],[-4 4], 80)
  axis equal;
  x0 = r.guess(int32(1));
case 'Hilb'
  r = Hilbert();
  disp(r.arity());
  r.contour([-4 4],[-4 4], 80)
  axis equal;
  x0 = r.guess(int32(1));
end


search_method   = LinesearchGoldenSection();
%search_method = LinesearchArmijo();
%search_method = LinesearchWolfe();

%dir_method = MinimizationGradientMethod(r,search_method);
%dir_method = MinimizationBFGS(r,search_method);
dir_method = MinimizationQuasiNewton(r,search_method);
%dir_method = MinimizationCG(r,search_method);

dir_method.setMaxIteration( int32(100) );
dir_method.setTolerance(1e-6);
dir_method.save_iterate_on();
dir_method.selectByName('BFGS');
%dir_method.selectByName('DFP');

[xs,converged] = dir_method.minimize( x0 );
dir_method.plotIter();

fprintf('\nThe solution found by the algorithm is:\n');
fprintf('       %d\n',xs)
exact = r.exact_solutions;
fprintf('\nThe analytical solution is:\n');
fprintf('       %d\n',exact)
