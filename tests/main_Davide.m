
clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%% Select test options %%

M = 20;
fun_num   = '16';    % function number to be minimized  ( 10 16 da sistemare!!!)
ls        = 'Armijo';   % linesearch algorithm
minimizer = 'BFGS'; % minimization method

%% Assing function %%

function_assign

%% Assing linesearch %%

switch ls
    case 'GS'
        linesearch_method = LinesearchGoldenSection();
    case 'Armijo'
        linesearch_method = LinesearchArmijo();
    case 'Wolfe'
        linesearch_method = LinesearchWolfe();
end

if ls == "GS"
    linesearch_method.setTolerance( 1e-8 );
    linesearch_method.setMaxIteration( int32(200) );
end
%linesearch_method.debug_on() ;

%% Assing minimization method %%

switch minimizer
    case 'Grad'
        min_method = MinimizationGradientMethod(r,linesearch_method);
    case 'BFGS'
        min_method = MinimizationBFGS( r, linesearch_method );
    case 'CG'
        min_method = MinimizationCG( r, linesearch_method );
end

min_method.setMaxIteration( int32(1000) );
min_method.setTolerance(1e-8);
min_method.debug_on();
min_method.no_FD_D();

%% Minimization algorithm execution

[xs,converged] = min_method.minimize( x0 ) ;

fprintf('\nFound solution:\n')
fprintf('%5.5f\n',xs)
fprintf('\nTest actual solution:\n')
fprintf('%5.5f\n',goal)
fprintf('\nComputed function minima:\n')
fprintf('%5.5f\n',r.eval(xs))
fprintf('\nAnalitycal or approximated minima:\n')
if isnan(r.exact_solutions(1))
    fprintf('%5.5f\n',approx_minima)
else
    fprintf('%5.5f\n',r.eval(goal))
end

fprintf('\n\nAll done folks!\n\n') ;

return

% first = true ;
% %for kk=2:14
% for kk=2
%
%     min_method.selectByNumber( int32(kk) );
%
%     fprintf('\n\nMethod %s\n\n',min_method.activeMethod() );
%
%
%     if first
%
%         subplot(2,1,1) ;
%         r.contour(RX,RY,fplot,80);
%         min_method.plotIter();
%         axis equal ;
%
%         subplot(2,1,2) ;
%         min_method.plotResidual();
%         first = false ;
%     else
%         subplot(2,1,2) ;
%         hold on ;
%         min_method.plotResidual();
%     end
% end




