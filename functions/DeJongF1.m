classdef DeJongF1 < FunctionND
  %
  % The De Jong Function F1
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
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

  methods

    function self = DeJongF1()
      self@FunctionND(int32(3)) ;
      self.exact_solutions = [ 0 ; 0 ; 0 ];     % one known solution 
      self.guesses         = [ -5.12; 0; 5.12 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = dot ( x, x ) ;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = 2*x' ;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = 2*eye(3);
    end
  end
end
