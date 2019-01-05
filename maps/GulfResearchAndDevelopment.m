 classdef GulfResearchAndDevelopment < FunctionMap
  %
  % Gulf research and development function
  %
  % Cox, R. A. "Comparison of the performance of seven optimization algorithms
  % on twelve unconstrained optimization problems". Gulf Research and 
  % Development Company, Pittsburg, Jan. 1969.
  % 
  %
  % see also in reference test N.11
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
    function self = GulfResearchAndDevelopment( varargin )
        if nargin == 1 
            if varargin{1} < 3 || varargin{1} > 100 || ~isinteger(varargin{1})
                error('GulfResearchAndDevelopmentD: argument must be an integer between 3 and 100, found %f', varargin{1});
            else    
                M = varargin{1};
            end    
        else
            M = 3;
        end    
      self@FunctionMap(int32(3),int32(M)); 
      self.exact_solutions = [50, 25, 1.5] % one exact solution   
      self.guesses         = [5.0;2.5;0.15]; 
      end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, x )
      % evaluate the entries (not squared).
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      F    = zeros(M,1);
      t    = zeros(M,1);    
      y    = zeros(M,1);
      for i = 1:M
         t(i) = i/100;
         y(i) = 25+(-50*log(t(i)))^(2/3);
         F(i) = exp(-abs(y(i)*M*i*X2)^X3/X1)-t(i);
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      J    = zeros(M,3);
      t    = zeros(M,1);    
      y    = zeros(M,1);
      for i = 1:M
         t(i) = i/100;
         y(i) = 25+(-50*log(t(i)))^(2/3);
         J(i,:) = [exp(-abs(y(i)*M*i*X2)^X3/X1)*(abs(y(i)*M*i*X2)^X3/X1^2),...%DFi/DX1
                  exp(-abs(y(i)*M*i*X2)^X3/X1)*(-X3/X1*abs(y(i)*M*i*X2)^X3/X2),...%DFi/DX2
                  exp(-abs(y(i)*M*i*X2)^X3/X1)*(-abs(y(i)*M*i*X2)^X3/X1*log(abs(y(i)*M*i*X2)))];%DFi/DX3
     end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X1   = x(1);
      X2   = x(2);
      X3   = x(3);
      t    = zeros(M,1);    
      y    = zeros(M,1);
      T    = zeros(M,3,3);
      for i = 1:M
            t(i) = i/100;
            y(i) = 25+(-50*log(t(i)))^(2/3);
            % D2F / DX1DX1
            T(i,1,1) = exp(-abs(y(i)*M*i*X2)^X3/X1)*(-2*abs(y(i)*M*i*X2)^X3/X1^3)+...
                exp(-abs(y(i)*M*i*X2)^X3/X1)*(abs(y(i)*M*i*X2)^X3/X1^2)^2;
            % D2F / DX1DX2
            T(i,1,2) = exp(-abs(y(i)*M*i*X2)^X3/X1)*((-X3/X1/X2*abs(y(i)*M*i*X2)^X3)*...
                (1/X1^2*abs(y(i)*M*i*X2)^X3)+(X3/X1^2/X2*abs(y(i)*M*i*X2)^X3));
            % D2F / DX1DX3
            T(i,1,3) = exp(-abs(y(i)*M*i*X2)^X3/X1)*...
                ((abs(y(i)*M*i*X2)^X3/X1*log(abs(y(i)*M*i*X2)))*(1-abs(y(i)*M*i*X2)^X3/X1^2));
            % D2F / DX2DX1
            T(i,2,1) = T(i,1,2);
            % D2F / DX2DX2
            T(i,2,2) = exp(-abs(y(i)*M*i*X2)^X3/X1)*((-X3/X1/X2*abs(y(i)*M*i*X2)^X3)^2+...
                X3/X1/X2^2*abs(y(i)*M*i*X2)^X3*(1-X3));
            % D2F / DX2DX3
            T(i,2,3) = exp(-abs(y(i)*M*i*X2)^X3/X1)*...
                ((-abs(y(i)*M*i*X2)^X3/X1*log(abs(y(i)*M*i*X2)))*...
                (-X3/X1/X2*abs(y(i)*M*i*X2)^X3)+...
                (-abs(y(i)*M*i*X2)^X3/X1/X2-X3/X1/X2*abs(y(i)*M*i*X2)^X3*log(abs(y(i)*M*i*X2))));
            % D2F / DX3DX1
            T(i,3,1) = T(i,1,3);
            % D2F / DX3DX2
            T(i,3,2) = T(i,2,3);
            % D2F / DX3DX3
            T(i,3,3) = exp(-abs(y(i)*M*i*X2)^X3/X1)*...
                ((-abs(y(i)*M*i*X2)^X3/X1*log(abs(y(i)*M*i*X2)))^2+...
                (-(log(abs(y(i)*M*i*X2)))^2/X1*abs(y(i)*M*i*X2)^X3));
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
