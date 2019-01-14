 classdef Watson < FunctionMap
  %
  % Watson function
  %
  % Kowalik, J. S., and Osborne, M. R. "Methods for Unconstrained Optimization Problems". 
  % Elsevier North-Holland, New York, 1968.
  % 
  %
  % see also in reference test N.20
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
    function self = Watson( varargin )
        if nargin == 1 
            if varargin{1} < 2 || varargin{1} > 31 || ~isinteger(varargin{1})
                error('Watson: argument n must be an integer between 2 and 31, found %f', varargin{1});
            else    
                N = varargin{1};
            end    
        else
            N = 6;
        end    
      self@FunctionMap(int32(N),int32(31)); 
      % self.exact_solutions = []; % no exact solution 
      % self.approximated_solutions = []; f = 2.28767...10^(-3) if n = 6,
      %                                   f = 1.39976...10^(-6) if n = 9, 
      %                                   f = 4.72238...10^(-10) if n = 12 
      self.guesses         = zeros(N,1); 
      end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      self.check_x( x );
      F    = zeros(31,1);
      t    = zeros(29,1);    
      for i = 1:29
         t(i) = i/29;
         F(i) = -x(1)^2-1;
         for j = 2:N
             F(i) = F(i) + (j-1)*x(j)*t(i)^(j-2) - (x(j)*t(i)^(j-1))^2;
         end    
      end
      F(30) = x(1);
      F(31) = x(2)-x(1)^2-1;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      J    = zeros(31,N);
      t    = zeros(29,1);    
      for i = 1:29
         t(i) = i/29;
         J(i,1) = -2*x(1);
         for j=2:N
         J(i,j) = (j-1)*t(i)^(j-2) -2*(t(i)^(j-1))^2*x(j);%DFi/Dx(j)
         end
      end   
     J(30,1) = 1;
     J(31,1) = -2*x(1);     J(31,2) = 1;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      t    = zeros(29,1);    
      T    = zeros(31,N,N);
      for i = 1:29
         t(i) = i/29;
         for j =1:N
            T(i,j,j) = -2*(t(i)^(j-1))^2;
         end    
      end
      T(31,1,1) = -2;
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
