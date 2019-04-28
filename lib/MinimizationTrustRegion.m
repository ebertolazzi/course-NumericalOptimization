classdef MinimizationTrustRegion < MinimizationND
 %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate function using the Trust Region method.
  % The algorithm, described in the references, mainly consists in approximating the 
  % function with a quadratic model and in finding the minimum of the model within a
  % region where we trust that the model well approximates the function. 
  % Indeed the step is limited (max is called trust radius). If the step is
  % successful in some sense it is accepted, otherwise both the direction
  % and the step length are recalculated, and usually the trust radius is
  % changed. Three different Trust Region algorithms are available: Exact
  % Trust Region, Dogleg and DoubleDogleg.
  %
  % References
  % ----------
  %
  % @article{2006,
  %     author   = {},
  %     journal  = {},
  %     keywords = {},
  %     number   = {},
  %     pages    = {},
  %     title    = {},
  %     volume   = {},
  %     year     = {}
  % }
  % @article{2006,
  %     author  = {},
  %     title   = {},
  %     journal = {},
  %     volume  = {},
  %     number  = {},
  %     year    = {},
  %     pages   = {},
  %     doi     = {},
  % }
  %
  % Authors: Enrico Bertolazzi & Eleonora Isotta
  %
  properties (SetAccess = private, Hidden = true)
    delta % trust radius, first guess
    beta1 % coefficient to compare the actual reduction achieved by the step vs the 
    % reduction predicted by the model. If rho (actual reduction/predicted reduction) is
    % < beta1 the step is not accepted and the trust radius recalculated.
    beta2 % coefficient to compare the actual reduction achieved by the step vs the 
    % reduction predicted by the model. If rho (actual reduction/predicted reduction) is
    % > beta1 the step is accepted and the trust radius can be recalculated.
    epsilon % small number, coefficient for exact trust region algorithm
    algorithm_names
    algorithm % type of Trust Region Algorithm that the user wishes to use, 
    % choose among Exact Trust Region, Dogleg and DoubleDogleg
     %
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    function self = MinimizationTrustRegion( fun )
      self@MinimizationND( fun, [] );
	  self.delta = 0.2;
      self.beta1 = 0.1; 
      self.beta2 = 0.5;
      self.epsilon = 0.01;
      self.algorithm_names = { 'Exact Trust Region', 'Dogleg', 'DoubleDogleg' };
      self.algorithm = self.algorithm_names{2};
    end
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    function setDelta( self , d )
      % Change radius of the Trust Region algorithm:
      if d <= 0
        error( 'MinimizationTrustRegion:setDelta(%2.6g) argument must be positive\n', d );
      end
      self.delta = d;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    function setBeta1( self , beta1 )
      % Change the beta1 in Trust Region algorithm, evaluation coefficient
      % for reduction      
      if beta1 <= 0
        error( 'MinimizationTrustRegion:setBeta1(%2.6g) argument must be positive\n', beta1 );
      end
      if beta1 >= self.beta2
        error( 'MinimizationTrustRegion:setBeta1(%2.6g) beta1 must be smaller than beta2, found: %2.6g >= %2.6g\n', beta1, beta1, self.beta2 );
      end
      self.beta1 = beta1;
    end
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
    function setBeta2( self , beta2 )
      % Change the beta2 in Trust Region algorithm, evaluation coefficient
      % for reduction
      if beta2 <= 0
        error( 'MinimizationTrustRegion:setBeta2(%2.6g) argument must be positive\n', beta2 );
      end
      if beta2 <= self.beta1
        error( 'MinimizationTrustRegion:setBeta2(%2.6g) beta2 must be greater than beta1, found: %2.6g <= %2.6g\n', beta2, beta2, self.beta1 );
      end
      self.beta2 = beta2;
    end
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     
    function setEpsilon( self , epsilon )
      % Change the epsilon in Trust Region algorithm: small number, coefficient for exact trust region algorithm
      if epsilon <= 0
        error( 'MinimizationTrustRegion:setEpsilon(%2.6g) argument must be positive\n', epsilon );
      end
      self.epsilon = epsilon;
    end
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    function n = numOfAlgorithms( self )
      % return the number of Trust Region Algorithms available
      n = length(self.algorithm_names);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function name = activeAlgorithm( self )
      % return active Trust Region Algorithm used in minimization
      name = self.Algorithm;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByNumber( self, k )
      % select the Trust Region Algorithm by number
      if k < 1 || k > length(self.algorithm_names)
        error('MinimizationTrustRegion, selectByNumber, k=%d out of range',k);
      end
      self.algorithm = self.algorithm_names{k};
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByName( self, name )
      % select the Trust Region Algorithm by its name
      if ischar(name)
        for k=1:length(self.algorithm_names)
          if strcmp(name,self.algorithm_names{k})
            self.algorithm = name;
            return;
          end
        end
        error('MinimizationTrustRegion, selectByName, name=%s not found',name);
      else
        error('MinimizationTrustRegion, selectByName, expected string as arument, found %s',class(name));
      end
    end
   
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, converged ] = minimize( self, x0 )
      x           = x0(:);
      converged   = false;

      if self.save_iterate
        self.x_history = x;
      end
      I = eye(length(x));
      rho = 0; 
      
      for iter=1:self.max_iter
        self.iter = iter;
        g1 = self.funND.grad( x ).';
        H = self.funND.hessian ( x );     
      
        % check if converged
        norm_inf_g1 = norm(g1,inf);
        converged   = norm_inf_g1 < self.tol;
        if converged
          if self.verbose
            fprintf( 'solution found, ||grad f||_inf = %g < %g\n', ...
                     norm_inf_g1, self.tol );
          end
          break;
        end
        %
        if self.verbose
          fprintf( 'Iter = %5d, rho = %5d ||grad f||_inf = %12.6g Delta = %12.6g\n', ...
                   self.iter, rho, norm_inf_g1, self.delta );
        end     
        
        % find search direction and step length ------------------------------------------
        switch ( self.algorithm ) % apply the desired Trust Region Algorithm
            
            case 'Exact Trust Region'; %EXACT TRUST REGION ALGORITHM------------------------------------------
                mu = 0;
                s = -inv(H)*g1;
                while (norm(s)-self.delta) > self.epsilon && mu >= 0
                    % compute the model
                    sp = -inv(H+mu*I)*s;
                    phi = norm(s)-self.delta;
                    phip = (s.'*sp)/norm(s);
                    % update mu and s
                    mu = mu - phi/phip*norm(s)/self.delta;
                    s = -inv(H+mu*I)*g1;
                end    
                if mu < 0
                    s = - inv(H)*g1;
                end  
                
            case 'Dogleg'; % DOGLEG ALGORITHM------------------------------------------
                sg = -g1*(norm(g1))^2/(g1.'*H*g1);
                sn = -inv(H)*g1;
                if self.delta <= norm(sg)
                    s = self.delta*sg/norm(sg);
                    fprintf('\nGrad\n');
                elseif self.delta>= norm(sn)
                    s = sn;
                    fprintf('\nNewton\n');
                else
                    a = (norm(sg))^2;
                    b = (norm(sn))^2;
                    c = (norm(sg-sn))^2;
                    d = (a+b-c)/2;
                    alpha = (b-self.delta^2)/(b-d+sqrt(d^2-a*b+self.delta^2*c));
                    s = alpha*sg + (1-alpha)*sn; 
                    fprintf('\nHybrid\n');
                end    
                
            case 'DoubleDogleg'; % DOUBLE DOGLEG ALGORITHM------------------------------------------
                sg = -g1*(norm(g1))^2/(g1.'*H*g1);
                sn = -inv(H)*g1;
                gamma = (norm(sg))^2/(sg.'*sn);
                if self.delta <= norm(sg)
                    s = self.delta*sg/norm(sg);
                elseif self.delta <= gamma*norm(sn)
                    A = gamma^2*(norm(sn))^2-(norm(sg))^2;
                    B = self.delta^2 - (norm(sg))^2;
                    alpha = (A-B)/(A+sqrt(A*B));
                    s = alpha*sg+(1-alpha)*sn;
                elseif self.delta <= norm(sn)
                    s = self.delta*sn/norm(sn);
                else
                    s = sn;
                end    
                
            otherwise;
            error( 'MinimizationTrustRegion, algorithm `%s` not supported', self.algorithm );    
        end
        
        % minimize along search direction and check the reduction------------------------------------------
        x_new = x + s; % new point
        p_red = -(g1.'*s + 1/2*s.'*H*s); % predicted reduction, %lascia così oppure scrivi brutalmente la funzione
        a_red = self.funND.eval(x) - self.funND.eval(x_new); % actual reduction
        rho = a_red/p_red; % ratio between actual and predicted reduction: according to its value we can assess the goodness of the step

        if rho >= self.beta2
           x = x_new;
           self.delta = max(2*norm(s), self.delta);
        elseif rho > self.beta1
           x = x_new;
        else
           self.delta = self.delta/2;
        end    
  
        
      end
      %end iteration
    end
    %end minimize
  end
end
