% DA FARE


classdef KowalikAndOsborne < FunctionMap
  %
  % Kowalink and Osborne function
  %
  % - ref [1]
  %
  %  @article{kowalik1968methods,
  %    title={Methods for unconstrained optimization problems},
  %    author={Kowalik, Janusz S and Osborne, Michael Robert},
  %    year={1968},
  %    publisher={North-Holland}
  % }
  %
  % see also in ref [2] test N.15
  %
  % ref[2]
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
    yi
    ui
  end

  methods

    function self = KowalikAndOsborne()
      % KowalikAndOsborne( )...............N = 4; M = 11;
      self@FunctionMap(int32(4),int32(11));       % call superclass constructor
      %exact_solutions        = [];                % no exacts solution provided
      self.approximated_solutions      = [ +inf , -14.07 , -inf , -inf ].';                 % approximated solution provided
      self.guesses                     = [ 0.25 ,  0.39  , 0.415, 0.39 ].';                 % one guess

      self.yi = [ 0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246 ].';
      self.ui = [ 4      2      1      0.5    0.25   0.1670 0.1250 0.1000 0.0833 0.0714 0.0625 ].';
    end

    function F = evalMap(self,x)
      % evaluate the entries (not squared) of the function.
      X1 = x(1);
      X2 = x(2);
      X3 = x(3);
      X4 = x(4);

      y = self.yi;
      u = self.ui;

      F  =  y + ( X1*(u.^2 + u.*X2)) ./ (u.^2 + u.*X3 + X4 ); % vector of [ f_1(x) ... f_n(x) ] values.
    end

    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1 = x(1);
      X2 = x(2);
      X3 = x(3);
      X4 = x(4);

      y = self.yi;
      u = self.ui;

      % Define common denominator and other useful quantities:
      UXX = (u.^2 + u.*X3 + X4 );
      UXX2 = UXX.^2;
      u2   = u.^2;

      J  = [ (u2 + u.*X2)./UXX, (u.*X1./UXX), - u.*X1.*( u2+u.*X2 )./UXX2, -X1.*( u2+u.*X2 )./UXX2 ];

    end

    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X1 = x(1);
      X2 = x(2);
      X3 = x(3);
      X4 = x(4);

      y = self.yi;
      u = self.ui;

      % Define common denominator and other useful quantities:
      UXX  = (u.^2 + u.*X3 + X4 );
      UXX2 = UXX.^2;
      UXX3 = UXX.^3;
      u2   = u.^2;

      % Create the n-matrices of T

      % D J / D X1
      T1 = [ 0        , u./UXX   , -u.*(u2+u.*X2)./UXX2          , -(u2+u.*X2)./(UXX2)        ];

      % D J / D X2 % start using simmetry
      T2 = [ T1(:,2)  , 0        , - u2.*X1./(UXX2)              , -u.*X1./(UXX2)             ];

      % D J / D X3 %
      T3 = [ T1(:,3)  , T2(:,3)  , 2.*u2.*X1.*(u2+X2.*u)./UXX3   , 2.*u.*X1.*(u2+X2.*u)./UXX3 ];

      % D J / D X4
      T4 = [ T1(:,4)  , T2(:,4)  , T3(:,4)                       , 2.*X1.*(u2+X2.*u)./UXX3    ];

      % Concatenate the n-matrices of T
      % Dimensions = MxNxN
      T  = cat(3,T1,T2,T3,T4);

    end

    % For tensor and jacobian a wolfram file is available: ask to the author if needed

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
