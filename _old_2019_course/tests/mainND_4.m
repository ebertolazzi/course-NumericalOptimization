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
minimization_method = MinimizationConjugateGradient( r, search_method );
minimization_method.setMaxIteration( int32(400) );
minimization_method.setTolerance(1e-8);
minimization_method.save_iterate_on();

t = 0:2*pi/1000:2*pi;

subplot(1,2,1);
r.contour([-1.5 1.5],[-1.5 1.5], 80, @(x) x);
hold on;
plot(cos(t),sin(t),'-b','Linewidth',2);
axis equal;
title('CG x0 = [0,0.9999]');

subplot(1,2,2);
r.contour([-1.5 1.5],[-1.5 1.5], 80, @(x) x);
hold on;
plot(cos(t),sin(t),'-b','Linewidth',2);
axis equal;
title('CG x0 = [0,0.8]');

INFO1 = {};
INFO2 = {};
nok = 0;
nno = 0;

for kkk=1:24
  fprintf('\n\n\n\n\n\n\nmethod N.%d\n\n\n\n\n\n\n',kkk);
  minimization_method.selectByNumber(kkk);
  
  INFO1{kkk}.name = minimization_method.activeMethod();
  INFO2{kkk}.name = minimization_method.activeMethod();

  %x0 = r.guess(int32(1));
  x0 = [0; 0.999993944545];
  [xs,converged] = minimization_method.minimize( x0 );
  xs
  INFO1{kkk}.iter = minimization_method.getIteration();
  if converged
    INFO1{kkk}.converged = 'YES';
    nok = nok + 1;
  else
    INFO1{kkk}.converged = 'NO ';
    nno = nno + 1;
  end

  subplot(1,2,1);
  minimization_method.plotIter();
  hold on;
  x0 = [0; 0.8];
  [xs,converged] = minimization_method.minimize( x0 );
  xs
  INFO2{kkk}.iter = minimization_method.getIteration();
  if converged
    INFO2{kkk}.converged = 'YES';
    nok = nok + 1;
  else
    INFO2{kkk}.converged = 'NO ';
    nno = nno + 1;
  end
  subplot(1,2,2);
  minimization_method.plotIter();
  hold on;
end


for kkk=1:length(INFO1)
  fprintf('%4d  converged %s iteration %3d converged %s iteration %3d : %s\n', ...
          kkk, INFO1{kkk}.converged, INFO1{kkk}.iter, ...
          INFO2{kkk}.converged, INFO2{kkk}.iter, INFO1{kkk}.name );
end

fprintf('converged = %d, NOT converged = %d\n', nok, nno);

