classdef Brown_Dennis_function < FunctionMap
    
%     properties(SetAccess = private , Hidden = true)
%         approximated_solutions
%     end
    
    methods
        function self = Brown_Dennis_function( M_input)
            fixed_arity = 4;
            self@FunctionMap(int32(fixed_arity) , int32(M_input)) ; % arity = 4;
            self.guesses         = [25;5;-5;-1] ;  % one guess
            self.M = M_input;
            self.N = fixed_arity;
            if M_input<fixed_arity
                error('M must be greater than 3')
            end
            
            if M_input == 20
               self.approximated_solutions = 85822.2;
            end
        end 
        
        function F = evalMap(self , x) % return a column vector with the not squared entries f1 f2 f3...
            t = 1/5*(1:self.M)';
            F = (x(1) + t.*x(2) - exp(t) ).^2 + (x(3) + x(4)*sin(t) - cos(t)).^2;
            
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