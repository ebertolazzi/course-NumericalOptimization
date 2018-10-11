classdef Wood_function < FunctionMap
    methods
        function self = Wood_function()
            arity = 4;
            M     = 6;
            self@FunctionMap(int32(arity) , int32(M)) ;
            self.exact_solutions = [1,1,1,1];       % one known solution
            self.guesses         = [-3,-1,-3,-1] ;  % one guess
        end
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
            f1 = 10*(x(2)-x(1)^2);
            f2 = 1-x(1);
            f3 = sqrt(90)*(x(4) - x(3)^2);
            f4 = 1- x(3);
            f5 = sqrt(10)*(x(2) + x(4) - 2);
            f6 = 10^(-1/2)*(x(2) - x(4));
            
            F = [ f1 ; f2 ; f3 ; f4 ; f5 ; f6 ];
            
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