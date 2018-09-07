clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%r = Rosenbrock();
%r = SchafferF7();
r = SixHumpCamelBack();

disp(r.arity());

linesearch_method = LinesearchArmijo();
%linesearch_method = LinesearchWolfe();
%linesearch_method = LinesearchGoldenSection();

minimization_method = MinimizationConjugateGradient( r, linesearch_method );
minimization_method.selectByName('FR');

minimization_method.setMaxIteration( int32(1000) );
minimization_method.setTolerance(1e-9);
minimization_method.verbose_on();
minimization_method.save_iterate_on();
minimization_method.no_FD_D();

fprintf('method = %s\n',minimization_method.activeMethod());

guess = r.guess(int32(1));
[xs,converged] = minimization_method.minimize( guess );

subplot(2,1,1);
[xmin,ymin,xmax,ymax] = minimization_method.iterRange();
r.contour( [xmin,xmax],[ymin,ymax], 40 );
axis equal;
minimization_method.plotIter();

subplot(2,1,2);
minimization_method.plotResidual();
xs

fprintf('method = %s\n',minimization_method.activeMethod());
