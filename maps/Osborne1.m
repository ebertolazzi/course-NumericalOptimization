 classdef Osborne1 < FunctionMap
  %
  % Osborne 1 function
  %
  % Osborne, M. R. "Some aspects of nonlinear least squares calculations".
  % Numerical Methods for Nonlinear Optimization, F. A. Lootsma (Ed),
  % Academic Press, New York, 1972, pp 171-189.
  % 
  % see also in reference test N.17
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
  
  properties (SetAccess = private, Hidden = true)
   y
   end
  
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Osborne1() 
      self@FunctionMap(int32(5),int32(33)); 
      % min f = 5.46489... 10^(-5)
      self.guesses         = [0.5; 1.5; -1.0; 0.01; 0.02]; 
      self.y = [0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818,...
                0.784, 0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558,...
                0.538, 0.522, 0.506, 0.490, 0.478, 0.467, 0.457, 0.448, 0.438,...
                0.431, 0.424, 0.420, 0.414, 0.411, 0.406];
      end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      X4   = x(4);
      X5   = x(5);
      F    = zeros(1,33);
      t    = zeros(1,33);
      for i=1:33
          t(i) = 10*(i-1);
          F(i) = self.y(i) - (X1 + X2*exp(-t(i)*X4) + X3*exp(-t(i)*X5));
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X2   = x(2);
      X3   = x(3);
      X4   = x(4);
      X5   = x(5);
      J    = zeros(33,5);
      t    = zeros(1,33);
      for i=1:33
            t(i) = 10*(i-1);
            % DF / DX1
            J(i,1) = -1;
            % DF / DX2
            J(i,2) = -exp(-t(i)*X4);
            % DF / DX3
            J(i,3) = -exp(-t(i)*X5);
            % DF / DX4
            J(i,4) = X2*t(i)*exp(-t(i)*X4);
            % DF / DX5
            J(i,5) = X3*t(i)*exp(-t(i)*X5);
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X2   = x(2);
      X3   = x(3);
      X4   = x(4);
      X5   = x(5);
      t    = zeros(1,33);
      T    = zeros(33,5,5);
      for i=1:33
          t(i) = 10*(i-1);
          T(i,1,:) = [0, 0, 0, 0, 0];
          T(i,2,:) = [0, 0, 0, t(i)*exp(-t(i)*X4), 0];
          T(i,3,:) = [0, 0, 0, 0, t(i)*exp(-t(i)*X5)];
          T(i,4,:) = [0, t(i)*exp(-t(i)*X4), 0, -X2*t(i)^2*exp(-t(i)*X4), 0];
          T(i,5,:) = [0, 0, t(i)*exp(-t(i)*X5), 0, -X3*t(i)^2*exp(-t(i)*X5)];
      end       
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
