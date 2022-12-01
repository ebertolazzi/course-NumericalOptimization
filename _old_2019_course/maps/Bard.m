 classdef Bard < FunctionMap
  %
  % Bard function
  %
  % Bard, Y. "Comparison of gradient methods for the solution of nonlinear
  % parameter estimation priblems". SIAM J. Numer. Anal. 7, 1970, 157-186.
  % 
  %
  % see also in reference test N.8
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
  c
  end
  
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Bard()
      self@FunctionMap(int32(3),int32(15)); 
      % self.exact_solutions        = [] % no exact solution  
      self.approximated_solutions = [0.8406;-inf;-inf];  
      self.guesses                 = [1.0;1.0;1.0]; 
      self.c = [0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, ...
          0.58, 0.73, 0.96, 1.34, 2.10, 4.39];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      F    = zeros(15,1);
      u    = zeros(15,1);    
      v    = zeros(15,1);
      w    = zeros(15,1);
      for i = 1:15
         u(i) = i;
         v(i) = 16-i;
         w(i) = min(u(i),v(i));
         F(i) = self.c(i) - (X1+u(i)/(v(i)*X2+w(i)*X3));
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X2   = x(2);
      X3   = x(3);
      J    = zeros(15,3);
      u    = zeros(15,1);    
      v    = zeros(15,1);
      w    = zeros(15,1);
      for i = 1:15
         u(i) = i;
         v(i) = 16-i;
         w(i) = min(u(i),v(i));
         J(i,:) = [-1 ,... % dF(i)/dX1
             u(i)*v(i)/(v(i)*X2+w(i)*X3)^2,... % dF(i)/dX2
             u(i)*w(i)/(v(i)*X2+w(i)*X3)^2]; % dF(i)/dX3
      end  
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X2 = x(2);
      X3   = x(3);
      u    = zeros(15,1);    
      v    = zeros(15,1);
      w    = zeros(15,1);
      T    = zeros(15,3,3);
      for i = 1:15
            % D J / D X1
            T(i,1,:) = [ 0, 0, 0 ];
            % D J / D X2
            T(i,2,:) = [ 0, -2*u(i)*v(i)^2/(v(i)*X2+w(i)*X3)^3, -2*u(i)*v(i)*w(i)/(v(i)*X2+w(i)*X3)^3];
            % D J / D X3
            T(i,3,:) = [ 0, -2*u(i)*v(i)*w(i)/(v(i)*X2+w(i)*X3)^3, -2*u(i)*w(i)^2/(v(i)*X2+w(i)*X3)^3];
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
