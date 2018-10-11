classdef Helical_valley < FunctionMap
    methods
        function self = Helical_valley()
            arity = 3;
            M     = 3;
            self@FunctionMap(int32(arity) , int32(M)) ;
            self.exact_solutions = [+1,0,0];   % one known solution
            self.guesses         = [-1,0,0] ; % one guess
        end
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
            f1 = 10*(x(3)-10*theta(x(1) , x(2)) );
            f2 = 10*(sqrt(x(1)^2 + x(2)^2) - 1  );
            f3 = x(3);
            F = [ f1 ; f2 ; f3 ];
            
            function a = theta(X1,X2)
                b = 0;
                if X1 < 0 % NB: Case X1 = 0 is not defined in the paper, here I defined it the same as X1 >0
                    b = 0.5;
                end
                a =  1/(2*pi)*atan(X2/X1) + b;
            end
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