classdef Han < FunctionMap
  %
  % Han-Sun-Han-SAMPAJO
  %
  % Reference:
  %
  %
  % Author:
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %

  methods

    function self = Han()
      self@FunctionMap(int32(2),int32(4)) ;
      self.exact_solutions = zeros(2,1); 
      self.guesses         = [ 5 ; 3 ];
    end

    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      X = x(1);
      Y = x(2);
      F = [ X^4 ; X ; X*Y ; exp(Y^2/2) ] ;
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1);
      Y = x(2);
      J = [ 4*X^3, 0 ; ...
            1, 0 ; ...
            Y, X ; ...
            0, Y*exp(Y^2/2) ] ;
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T = zeros(3,2,2);
      X = x(1);
      Y = x(2);
      T(1,:,:) = [ 12*X^2, 0 ; 0, 0 ] ;
      T(2,:,:) = [ 0, 0 ; 0, 0 ] ;
      T(3,:,:) = [ 0, 1 ; 1, 0 ] ;
      T(4,:,:) = [ 0, 0 ; 0, (1+Y^2)*exp(Y^2/2) ] ;
    end

  end
end
