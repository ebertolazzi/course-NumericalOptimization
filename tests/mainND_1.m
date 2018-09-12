clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');
addpath('../maps');

%r = Rosenbrock();
%r = SchafferF7();
%r = BoxThreeDimensionalFunction(3);
%r = KowalikAndOsborne();
%r = Leon();
%r = Box3();
%r = PowellSingular();
%r = Han();
%r = Helix(); %%%
%r = SixHumpCamelBack(); %%%
%r = PenaltyN1();
%r = PenaltyN2();
%r = Colville();
r = ShekelSQRN5();

disp(r.arity());

%linesearch_method = LinesearchArmijo();
linesearch_method = LinesearchMoreThuente();
%linesearch_method = LinesearchGoldenSection();

if false
  minimization_method = MinimizationConjugateGradient( r, linesearch_method );
  minimization_method.selectByName('FR');
else
  minimization_method = MinimizationQuasiNewton( r, linesearch_method );
end

minimization_method.setMaxIteration( int32(200) );
minimization_method.setTolerance(1e-9);
minimization_method.verbose_on();
minimization_method.save_iterate_on();
minimization_method.no_FD_D();

fprintf('method = %s\n',minimization_method.activeMethod());

guess = r.guess(int32(1));
[ xs, converged ] = minimization_method.minimize( guess );

subplot(2,1,1);
[xmin,ymin,xmax,ymax] = minimization_method.iterRange();
r.contour( [xmin,xmax],[ymin,ymax], 40, @(x) sqrt(sqrt(sqrt(x))) );
axis equal;
minimization_method.plotIter();

subplot(2,1,2);
minimization_method.plotResidual();
xs

fprintf('method = %s\n',minimization_method.activeMethod());
