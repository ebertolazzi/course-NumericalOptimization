 classdef FreudensteinAndRoth < FunctionMap
  %
  % Freudenstein and Roth function
  %
  % Freudstein F., and Roth B. "Numerical solutions of systems of nonlinear
  % equations". J ACM 10, Oct. 1963.
  % 
  %
  % see also in reference test N.2
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
  % Author: Eleonora Isotta - University of Trento

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = FreudensteinAndRoth()
      self@FunctionMap(int32(2),int32(2)); 
      self.exact_solutions        = [5;4]; % one exact solution  
      self.approximated_solutions = [11.41;-0.8968]; % one approximated solution 
      self.guesses                 = [0.5;-2];           
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      X1   = x(1);
      X2   = x(2);
      F    = zeros(2,1);
      F(1) = -13 + X1 + ((5 - X2)*X2 - 2)*X2;
      F(2) = -29 + X1 + ((X2 + 1)*X2 - 14)*X2;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X2 = x(2);
      J = [ 1, +10*X2 - 3*X2^2 -2; 1, 3*X2^2 +2*X2 -14 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X2 = x(2);
      T = zeros(2,2,2);
      % D J / D X1
      T(1,:,:) = [ 0, 0; 0, 0 ];
      % D J / D X2
      T(2,:,:) = [ 0, 10 -6*X2; 0, 6*X2 + 2 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
