classdef Osborne_2_function < FunctionMap
    properties(SetAccess = private , Hidden = true)
        yi;
        ti;
        approximated_minima;
    end
    methods
        function self = Osborne_2_function()
            arity  = 11;
            M      = 65;
            self@FunctionMap(int32(arity) , int32(M)) ;
            self.exact_solutions = NaN*ones(arity,1);
            self.guesses         = [1.3; 0.65; 0.65; 0.7; 0.6; 3; 5; 7; 2; 4.5; 5.5];   % one known solution
            self.approximated_minima = 4.01377e-2;
            self.yi = [1.366 1.191 1.112 1.013 0.991 0.885 0.831 0.847 0.786 0.725...
                       0.746 0.679 0.608 0.655 0.616 0.606 0.602 0.626 0.651 0.724...
                       0.649 0.649 0.694 0.644 0.624 0.661 0.612 0.558 0.533 0.495...
                       0.500 0.423 0.395 0.375 0.372 0.391 0.396 0.405 0.428 0.429...
                       0.523 0.562 0.607 0.653 0.672 0.708 0.633 0.668 0.645 0.632...
                       0.591 0.559 0.597 0.625 0.739 0.710 0.729 0.720 0.636 0.581...
                       0.428 0.292 0.162 0.098 0.054 ];
            self.ti = ((1:M)'-1 )/10;
        end
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
           
            F = self.yi' - ( x(1).*exp( -self.ti*x(5)          ) + x(2).*exp( -(self.ti-x(9) ).^2*x(6) ) +...
                             x(3).*exp( -(self.ti-x(10)).^2*x(7)) + x(4).*exp( -(self.ti-x(11)).^2*x(8) )) ;      
        end
        
        % Use finite difference for grad and hessian
        function g = grad( self, x )
            g = self.FD_grad( x );
        end
        
        function h = hessian( self, x )
            h = self.FD_hessian( x );
        end
        
    end
end