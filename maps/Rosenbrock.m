classdef Rosenbrock < FunctionMap
  %
  % The Rosenbrock Function
  %
  % An Automatic Method for Finding the Greatest or Least Value of a Function
  % H. H. Rosenbrock
  % The Computer Journal, Volume 3, Issue 3, 1 January 1960, Pages 175–184
  % DOI: 10.1093/comjnl/3.3.175
  %
  % see also in reference test N.1
  %
  % @article{More:1981,
  %   author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},
  %   title   = {Testing Unconstrained Optimization Software},
  %   journal = {ACM Trans. Math. Softw.},
  %   volume  = {7},
  %   number  = {1},
  %   year    = {1981},
  %   pages   = {17--41},
  %   doi     = {10.1145/355934.355936}
  % }
  %
  % Author: Enrico Bertolazzi
  %
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Rosenbrock()
      self@FunctionMap(int32(2),int32(2));
      self.exact_solutions = [1;1];   % one known solution
      self.guesses         = [-1.2;1]; % one guess
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate Rosenbrock (2D) function.
      X = x(1);
      Y = x(2);
      F = [ sqrt(100)*(Y-X^2); 1-X ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1);
      J = [ -2*sqrt(100)*X, sqrt(100); -1, 0 ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T        = zeros(2,2,2);
      T(1,:,:) = [ -2*sqrt(100), 0; 0, 0 ];
      T(2,:,:) = [ 0, 0; 0, 0 ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
