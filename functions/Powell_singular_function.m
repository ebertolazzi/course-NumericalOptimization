classdef Powell_singular_function < FunctionMap
    
    methods
        function self = Powell_singular_function()
            arity = 4;
            M     = 4;
            self@FunctionMap(int32(arity) , int32(M)) ;
            self.exact_solutions = [0;0;0;0];    % one known solution
            self.guesses         = [3;-1;0;1] ;  % one guess
        end
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
            f1 = x(1)+10*x(2);
            f2 = sqrt(5)*(x(3)-x(4));
            f3 = (x(2)-2*x(3))^2;
            f4 = sqrt(10)*(x(1)-2*x(4))^2;
            
            F = [ f1 ; f2 ; f3 ; f4 ];
            
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