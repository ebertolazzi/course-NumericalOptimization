classdef Meyer_function < FunctionMap
    properties(SetAccess = private , Hidden = true)
        yi;
    end
    methods
        function self = Meyer_function()
            arity  = 3;
            M = 16;
            self@FunctionMap(int32(arity) , int32(M)) ;
            self.exact_solutions = [0.02, 4000, 250];   % one known solution
            self.guesses         = [-1,0,0] ; % one guess
            self.approximated_solutions = 87.9458;
            self.yi = [34780 28610 23650 19630 16370 13720 11540 9744 8261 7030 6005 5147 4427 3820 3307 2872];
        end
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
            F = (x(1).*exp(x(2)./((45+5*(1:16)) + x(3)) )- self.yi)';
        end
        
        % Use finite difference for grad and hessian
        function g = grad( self, x )
            g = self.FD_grad( self, x );
        end
        
        function h = hessian( self, x )
            h = self.FD_hessian( self, x );
        end
        
    end
end