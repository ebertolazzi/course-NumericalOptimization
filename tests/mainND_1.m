clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');
addpath('../maps');

function_list
map_list

INFO = {};
for kkk=1:(length(FUNCTION_LIST)+length(MAP_LIST))
  if kkk <= length(MAP_LIST)
    N = MAP_LIST{kkk};
  else
    N = FUNCTION_LIST{kkk-length(MAP_LIST)};
  end
  INFO{kkk}.name = N;
  fprintf('\n\n\nfunction %s\n\n\n',N);
  r = feval(N);
  disp(r.arity());

  %linesearch_method = LinesearchArmijo();
  %linesearch_method = LinesearchMoreThuente();
  linesearch_method = LinesearchGoldenSection();

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
  INFO{kkk}.iter = minimization_method.getIteration();
  if converged
    INFO{kkk}.converged = 'YES';
  else
    INFO{kkk}.converged = 'NO ';
  end

  if r.arity() == 2
    subplot(2,1,1);
    [xmin,ymin,xmax,ymax] = minimization_method.iterRange();
    r.contour( [xmin,xmax],[ymin,ymax], 40, @(x) sqrt(sqrt(sqrt(x))) );
    axis equal;
    minimization_method.plotIter();
    subplot(2,1,2);
    minimization_method.plotResidual();
  else    
    minimization_method.plotResidual();
  end
  %xs

  %%fprintf('method = %s\n',minimization_method.activeMethod());
end

for kkk=1:length(INFO)
  fprintf('converged %s iteration %3d : %s\n', ...
          INFO{kkk}.converged, INFO{kkk}.iter, INFO{kkk}.name );
end