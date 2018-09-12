clear all;
% 
%    1  converged YES iteration  13 : Beale
%    2  converged YES iteration   8 : BiggsEXP6
%    3  converged YES iteration  16 : Box3
%    4  converged YES iteration  26 : BoxThreeDimensionalFunction
%    5  converged NO  iteration 200 : BrownBadlyScaled
%    6  converged YES iteration   4 : Gauss
%    7  converged YES iteration  20 : GaussianFunction
%    8  converged YES iteration  11 : Han
%    9  converged YES iteration  27 : Helix
%   10  converged YES iteration   8 : Himmelblau
%   11  converged YES iteration   2 : JennrichAndSampson
%   12  converged YES iteration  26 : KowalikAndOsborne
%   13  converged YES iteration  41 : Leon
%   14  converged YES iteration   1 : PowellSingular
%   15  converged NO  iteration 200 : PowellBadlyScaled
%   16  converged YES iteration  35 : Rosenbrock
%   17  converged YES iteration   2 : TrivialQuadratic
%   18  converged YES iteration   6 : Barrier1
%   19  converged YES iteration  10 : BohachevskyN1
%   20  converged YES iteration   8 : BohachevskyN2
%   21  converged YES iteration   8 : BohachevskyN3
%   22  converged YES iteration   9 : BraninRCOS
%   23  converged YES iteration  35 : Colville
%   24  converged YES iteration   4 : DeJongF1
%   25  converged YES iteration  37 : DeJongF2
%   26  converged YES iteration   2 : DeJongF3
%   27  converged YES iteration  18 : DeJongF4
%   28  converged YES iteration   1 : DeJongF5
%   29  converged YES iteration  10 : Easy_function_3D
%   30  converged YES iteration   8 : GoldsteinPrice
%   31  converged YES iteration   3 : GregoryAndKarney
%   32  converged YES iteration   3 : Hilbert
%   33  converged YES iteration 166 : PenaltyN1
%   34  converged YES iteration  11 : PenaltyN2
%   35  converged YES iteration  10 : Powell3D
%   36  converged YES iteration   5 : Quadratic2D
%   37  converged NO  iteration   3 : Rastrigin
%   38  converged YES iteration   3 : SchafferF6
%   39  converged YES iteration   4 : SchafferF7
%   40  converged NO  iteration 200 : ShekelSQRN10
%   41  converged YES iteration  11 : ShekelSQRN5
%   42  converged YES iteration  10 : ShekelSQRN7
%   43  converged NO  iteration 200 : Shubert
%   44  converged YES iteration   7 : SixHumpCamelBack
%   45  converged YES iteration  19 : mckinnon
%
% converged = 40, NOT converged = 5
%
close all;
clc;

addpath('../lib');
addpath('../functions');
addpath('../maps');

function_list
map_list

INFO = {};
nok = 0;
nno = 0;

for kkk=1:(length(FUNCTION_LIST)+length(MAP_LIST))
%for kkk=43
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
  xs
  INFO{kkk}.iter = minimization_method.getIteration();
  if converged
    INFO{kkk}.converged = 'YES';
    nok = nok + 1;
  else
    INFO{kkk}.converged = 'NO ';
    nno = nno + 1;
  end

  if r.arity() == 2
    subplot(2,1,1);
    [xmin,ymin,xmax,ymax] = minimization_method.iterRange();
    if xmax-xmin > ymax-ymin
      m    = (ymax+ymin)/2;
      ymax = m + (xmax-xmin)/2;
      ymin = m - (xmax-xmin)/2;
    else
      m    = (xmax+xmin)/2;
      xmax = m + (ymax-ymin)/2;
      xmin = m - (ymax-ymin)/2;        
    end
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
  fprintf('%4d  converged %s iteration %3d : %s\n', ...
          kkk, INFO{kkk}.converged, INFO{kkk}.iter, INFO{kkk}.name );
end

fprintf('converged = %d, NOT converged = %d\n', nok, nno);
