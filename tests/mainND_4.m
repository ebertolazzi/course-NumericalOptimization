clear all;
close all;
clc;

warning('on');

addpath('../lib');
addpath('../functions');

r = Barrier1();

%search_method   = LinesearchGoldenSection();
search_method   = LinesearchMoreThuente();
%search_method = LinesearchArmijo();
%search_method = LinesearchWolfe();

%search_method.debug_on();
dir_method = MinimizationConjugateGradient( r, search_method );
%dir_method = MinimizationQuasiNewton(r, search_method);


dir_method.setMaxIteration( int32(400) );
dir_method.setTolerance(1e-6);
dir_method.save_iterate_on();

t = 0:2*pi/1000:2*pi;

subplot(1,2,1);
r.contour([-1.5 1.5],[-1.5 1.5], 80, true);
hold on;
plot(cos(t),sin(t),'-b','Linewidth',2);
axis equal;
title('CG x0 = [0,0.9999]');

subplot(1,2,2);
r.contour([-1.5 1.5],[-1.5 1.5], 80, true);
hold on;
plot(cos(t),sin(t),'-b','Linewidth',2);
axis equal;
title('CG x0 = [0,0.8]');

if true
  for kkk=24:24
<<<<<<< HEAD
    fprintf('\n\n\n\n\n\n\nmethod N.%d\n\n\n\n\n\n\n',kkk);
=======
    fprintf('\n\n\n\n\n\n\n%d\n\n\n\n\n\n\n',kkk);
>>>>>>> 777dc57db9fd3e478fd721447f3381155df7b8c5
    dir_method.selectByNumber(kkk);
    %x0 = r.guess(int32(1));
    x0 = [0; 0.999993944545];
    [xs,converged] = dir_method.minimize( x0 );
    xs
    subplot(1,2,1);
    dir_method.plotIter();
    hold on;
    x0 = [0; 0.8];
    [xs,converged] = dir_method.minimize( x0 );
    xs
    subplot(1,2,2);
    dir_method.plotIter();
    hold on;
  end
else
  %x0 = r.guess(int32(1));
  x0 = [0; 0.999993944545];
  [xs,converged] = dir_method.minimize( x0 );
  subplot(1,2,1);
  dir_method.plotIter();
  xs
  x0 = [0; 0.8];
  [xs,converged] = dir_method.minimize( x0 );
  subplot(1,2,2);
  dir_method.plotIter();
  xs
end
