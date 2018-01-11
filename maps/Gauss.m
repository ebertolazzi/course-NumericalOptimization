classdef Gauss < FunctionMap
  %
  % The Gaussian function.
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Zbigniew Michalewicz,
  %   Genetic Algorithms + Data Structures = Evolution Programs,
  %   Third Edition,
  %   Springer Verlag, 1996,
  %   ISBN: 3-540-60676-9,
  %   LC: QA76.618.M53.
  %
  %   Richard Brent,
  %   Algorithms for Minimization with Derivatives,
  %   Dover, 2002,
  %   ISBN: 0-486-41998-3,
  %   LC: QA402.5.B74.
  %
  % Author:
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %

  properties (SetAccess = protected, Hidden = true)
    y
  end

  methods

    function self = Gauss()
      self@FunctionMap(int32(3),int32(15)) ;
      self.exact_solutions = [ 0 ; 0 ; 0 ];     % one known solution 
      self.guesses         = [ 0.4; 1.0; 0.0 ];

      self.y = [ 0.0009, 0.0044, 0.0175, 0.0540, ...
                 0.1295, 0.2420, 0.3521, 0.3989, ...
                 0.3521, 0.2420, 0.1295, 0.0540, ...
                 0.0175, 0.0044, 0.0009 ] ;
    end

    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      Z = x(3) ;
      F = zeros(15,1) ;
      for i = 1 : 15
        F(i) = X * exp ( - 0.5 * Y * ( 3.5 - 0.5 * (i - 1) - Z )^2 ) - self.y(i);
      end
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      Z = x(3) ;
      J = zeros(15,3);
      for i = 1 : 15
        tmp1 = 3.5 - 0.5 * (i - 1) - Z ;
        tmp2 = exp( - 0.5 * Y * tmp1^2 ) ;
        J(i,:) = tmp2 * [ 1, -0.5 * tmp1^2 * X, - X * Y * tmp1 ] ;
      end
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T = zeros(15,3,3);
      X = x(1) ;
      Y = x(2) ;
      Z = x(3) ;
      for i = 1 : 15
        tmp1 = 2*Z-8+i ;
        tmp2 = tmp1^2/8 ;
        tmp3 = exp(-tmp2) ;
        T31  = -(Y/2) * tmp1/2 ;
        T23  = (X/16)*tmp1*(Y*tmp1^2-8) ;
        T(i,:,:) = tmp3 * [   0,     -tmp2, T31 ; ...
                          -tmp2, -X*tmp2^2, T23 ; ...
                            T31,       T23, (X*Y/4)*(Y*tmp1^2-4) ] ;
      end
    end
  
  end
end
