clear all;
close all;
clc;

warning('on');

addpath('../lib');
addpath('../functions');

r = Barrier1();
disp(r.arity());

search_method   = LinesearchGoldenSection();
%search_method   = LinesearchMoreThuente();
%search_method = LinesearchArmijo();
%search_method = LinesearchWolfe();

%search_method.debug_on();
%dir_method = MinimizationCG( r, search_method );
%dir_method = MinimizationBFGS(r, search_method);
dir_method = MinimizationQuasiNewton(r, search_method);
%dir_method = MinimizationGradientMethod( r, search_method );


dir_method.setMaxIteration( int32(400) );
dir_method.setTolerance(1e-6);
dir_method.save_iterate_on();

r.contour([-1.5 1.5],[-1.5 1.5], 80);

t = 0:2*pi/1000:2*pi;
hold on;
plot(cos(t),sin(t),'-b','Linewidth',2);
axis equal;

if false
  for kkk=2:24
    fprintf('\n\n\n\n\n\n\n%d\n\n\n\n\n\n\n',kkk);
    dir_method.selectByNumber(kkk);
    x0 = r.guess(int32(1));
    x0 = [0.999993944545; 0];
    [xs,converged] = dir_method.minimize( x0 );
    xs
    dir_method.plotIter(true);
    hold on;
  end
else
  x0 = r.guess(int32(1));
  x0 = [0.999993944545; 0];
  %x0 = [-0.447147077055711;-0.894293265304308];
  [xs,converged] = dir_method.minimize( x0 );
  dir_method.plotIter();
  xs
end
