 classdef Wood < FunctionMap
  %
  % Wood function
  %
  % Colville, A. R. "A comparative study of nonlinear programming codes".
  % Rep. 320-2949, IBM New York Scientific Center, 1968.
  % 
  % see also in reference test N.14
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
    function self = Wood() 
      self@FunctionMap(int32(4),int32(6)); 
      self.exact_solutions = [1.0; 1.0; 1.0; 1.0]; % one exact solution   
      self.guesses         = [-3.0; -1.0; -3.0; -1.0]; 
      end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      X4   = x(4);
      F    = [10*(X2-X1^2);...
              1-X1;...
              90^(0.5)*(X4-X3^2);...
              1-X3;...
              10^(0.5)*(X2+X4-2);...
              10^(-0.5)*(X2-X4)];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1   = x(1);
      X3   = x(3);
      J    = zeros(6,4);
      % DF / DX1
      J(:,1) = [-20*X1, -1, 0, 0, 0, 0];
      % DF / DX2
      J(:,2) = [10, 0, 0, 0, 10^(0.5), 10^(-0.5)];
      % DF / DX3
      J(:,3) = [0, 0, -180^(0.5)*X3, -1, 0, 0];
      % DF / DX4
      J(:,4) = [0, 0, 90^(0.5), 0, 10^(0.5), -10^(-0.5)];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(6,4,4);
      T(1,:,:) = [-20, 0, 0, 0;...
                    0, 0, 0, 0;...
                    0, 0, 0, 0;...
                    0, 0, 0, 0];
      T(2,:,:) = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0];
      T(3,:,:) = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, -180^(0.5), 0;...
                  0, 0, 0, 0]; 
      T(4,:,:) = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0];  
      T(5,:,:) = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0];
      T(6,:,:) = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 0, 0];        
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
