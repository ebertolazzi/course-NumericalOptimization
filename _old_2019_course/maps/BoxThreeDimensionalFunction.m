

classdef BoxThreeDimensionalFunction < FunctionMap
  %
  % Gaussian Function N = 3 , M = 15;
  %
  % Box, M.J. A comparison of several current optimization methods,
  % and the use of transformations in constrained problems.
  % Comput. J 9 (1966), 67-77.
  %
  %
  % see also in reference test N.12
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
  % Author: Giammarco Valenti - University of Trento

  properties( SetAccess = private ,  Hidden = true)
    exact_linear_set % Provides a subspace (linear) of the variable space which is a common global minimum
    ii
    ti
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = BoxThreeDimensionalFunction( varargin )
      if nargin == 1
        M = varargin{1};
      else
        M = 3;
      end
      % BoxThreeDimensionalFunction()...............N = 3; M >= N;
      self@FunctionMap(int32(3),int32(M));              % call superclass constructor (initialize M)
      self.exact_solutions  = [ 1, 10, 1; 10, 1,-1 ].'; % more than one (f = 0 those cases)
      self.exact_linear_set = [ 1; 1; 0 ];              % minimum when x1=x2 and x3=0. so it is span(exact_linear_set)
      self.guesses          = [ 0; 10; 20 ];
      self.ti               = 0.1*(1:double(self.M)).';         % t for the function
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = evalMap( self, xx )
      % evaluate the entries (not squared) of the function.
      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      t  = self.ti;   % t for the function
      f  = exp( -t*x1 ) - exp( -t*x2 ) + x3 * (exp( -10*t ) - exp( -t ));
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, xx )
      % use analytic jacobian
      self.check_x( xx );
      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      t  = self.ti;   % t for the function

      % Matlab fashion assignment of Jacobian
      J = [ -t.*exp( -t*x1 ), t.*exp( -t*x2 ),  exp( -10*t ) - exp( -t ) ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, xx )
      % use analytic tensor

      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      % Create the n-matrices of T

      t = self.ti;   % t for the function

      zz = zeros(self.M,1);

      T1 = [ t.^2.*exp( t*x1 ), zz,                 zz ];
      T2 = [ zz,                -t.^2.*exp( t*x2 ), zz ];
      T3 = [ zz,                zz,                 zz ];

      T  = cat(3,T1,T2,T3);

    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end

%  #######################################################
%  #  _______   __  _____ ______ _______ _ ____   _____  #
%  #         \ |  |   ___|   ___|__   __| |    \ |       #
%  #          \|  |  __| |  __|    | |  | |     \|       #
%  #       |\     | |____| |       | |  | |   |\         #
%  #  _____| \____| _____|_|       |_|  |_|___| \______  #
%  #                                                     #
%  #######################################################
