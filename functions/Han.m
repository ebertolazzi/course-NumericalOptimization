classdef Han < FunctionND
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
      self@FunctionND(int32(2)) ;
      self.exact_solutions = zeros(2,1); 
      self.guesses         = [ 5 ; 3 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = x(1)^8 + x(1)^2 + x(1)^2 * x(2)^ 2 + exp(x(2)^2);
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g   = zeros ( 1, 2 );
      t1  = (x(1) ^ 2);
      t3  = (t1 ^ 2);
      t7  = x(2) ^ 2;
      t12 = exp(t7);
      g(1) = (8 * t3 * t1 * x(1)) + (2 * x(1)) + 0.2e1 * x(1) * t7;
      g(2) = (2 * t1 * x(2)) + 0.2e1 * x(2) * t12;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h   = zeros( 2, 2 ) ;
      t1  = x(1)^2 ;
      t2  = t1^2 ;
      t5  = x(2)^2;
      t9  = 4 * x(1) * x(2);
      t11 = exp(t5) ;
      h(1,1) = (56 * t2 * t1) + 2 + 2 * t5;
      h(1,2) = t9;
      h(2,1) = t9;
      h(2,2) = (2 * t1) + 2 * t11 + 4 * t5 * t11;
    end
  end
end
