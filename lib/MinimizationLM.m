classdef MinimizationLM < MinimizationND
 %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate functions using Levenberg-Marquardt method.
  % The algorithm is described in the references.
  %
  % The function to minimize must be of class "functionMAP << FunctionND"
  % And it must contain the instance method "jacobian(x)"
  %
  % Reference algorithm is in ref[1]: algorithm 3.16;
  %
  % References
  % ----------
  %
  % ref [1]
  % @article{madsen1999methods,
  %       title={Methods for non-linear least squares problems},
  %       author={Madsen, Kaj and Nielsen, Hans Bruun and Tingleff, Ole},
  %       year={1999}
  % }
  %
  %
  % Author: Enrico Bertolazzi & Giammarco Valenti
  %
  properties (SetAccess = private, Hidden = true)
    % direction_short %
    % angle_too_small %
    tol2 % the second tolerance
    tau  % parameter to choose mi

  end

  methods

    function self = MinimizationLM( fun )

      self@MinimizationND( fun, [] ); % Linesearch is empty -> no linesearch method in LM

      self.tol2        = 1e-6;
      self.tol         = 1e-3;
      self.max_iter    = 100;
      self.debug_state = false;
      self.FD_D        = true;
      self.tau         = 10^-6;

    end

    function setTau( self , tau )
      % Change tau of the LM algorithm:
      % Note the higher the tau the higher the confidence to start close to the minimum:
      % see ref[1] pag 27

      self.tau = tau;

      fprintf( 1 , 'set tau to %2.6g\n' , self.tau );

    end

    function setEpsilon2( self , tolerance )
      % Change tolerance epsilon 2 (see ref[1]):

      self.tol2 = tolerance;

      fprintf(1,'set epsilon2 to %2.10g\n', self.tol2 );
    end



    function [x,converged] = minimize( self , x0 )
      % Launch the minimization algorithm

      % Initial values -> see ref[1] algorithm 3.16

      x         = x0;                     % initial point
      converged = false;                  % initialization
      tau       = self.tau;               % see ref[1]
      J         = self.funND.jacobian(x); % initial jacobian
      f         = self.funND.evalMap(x);  % initial values
      g         = J.' * f;                % (1/2)gradient'
      A         = J.' * J;                % first term of Hessian
      mi        = tau * max(diag(A));     % First mi (ref[1])
      ni        = 2;                      % ni -> multiplicative factor of mi
      eps2      = self.tol2;              % tolerance on the error

      % check convercence and assign the boolean converged.
      nrm_g     = norm(g,inf);      % infinite norm (max)
      converged = nrm_g < self.tol; % Is it already converged?

      if converged
        if self.debug_state % DEBUG ONLY
          fprintf(1,'starting point is already a minimum, ||grad f||_inf = %12.6g ... \n', iter, nrm_g );
        end
        return;
      end

      if self.debug_state % DEBUG ONLY
        self.x_history = reshape( x, length(x), 1 ); % create the vector for history
      end

      % == START ITERATION ======

      for iter = 1:self.max_iter



        if self.debug_state % DEBUG ONLY

          fprintf(1,'iter = %5d ||grad f||_inf = %12.6g ', iter, nrm_g );
          fprintf(1,'ni = %5d , mi = %12.6g ... \n', ni, mi);

          self.x_history(:,iter) = x; % it change size every iteration, but it is not the bottleneck

        end

        % Solve( (J'J - uI)h = - grad )
        h_lm = -  ( A + mi * eye(length(x)) ) \ g;


        % Condition on the h_lm
        nrm_h_lm      = norm(h_lm);
        eps2_x_eps2   = eps2*( norm(x) + eps2  );

        if ( nrm_h_lm <= eps2_x_eps2 )

          converged = true; % update the boolean converged


          if self.debug_state && converged % DEBUG ONLY

              fprintf(1,'solution found (epsilon 2 reached), || n_lm || = %g < %g \n', nrm_h_lm , eps2_x_eps2 );

          end

          break; % out of the loop

        else

        % Update x
        x_new = x + h_lm;

        % compute L0 - L(h_lm)
        L0_Lh = (1/2) * h_lm.' * ( mi * h_lm - g ); % eq(3.7b) version at page 25 of ref [1]

        % Gain ratio (rho)
        rho   = ( self.funND.eval(x) - self.funND.eval(x_new) ) / ( L0_Lh );
        % In the future this step can be optimized evaluating the function only when x changes

          if rho > 0 % Acceptable step

            x = x_new;                           % x update
            J = self.funND.jacobian(x);          % Jacobian (functionMap)
            A = J.'*J;                           % First term hessian
            f = self.funND.evalMap(x);           % function values
            g = J.' * f;                         % (1/2)gradient'

            nrm_g     = norm(g,inf);       % infinite norm (max)
            converged = nrm_g <= self.tol; % check convergence

            if converged % Converged: stop the algorithm

              if self.debug_state  % DEBUG ONLY

                  fprintf(1,'solution found, ||grad f||_inf = %g < %g \n', nrm_g, self.tol );

              end

            break; % out of the loop

          end


            mi = mi * max( 1/3 , ( 1 - ( 2*rho - 1 ).^3 ) ); % Change damping
            ni = 2;                                          % reset ni

          else

            mi = mi*ni; % Increase damping
            ni = 2*ni;  % increase ni

          end

        end

      end
      % == END ITERATION ========
    end

  end

end


%  #######################################################
%  #  _______   __  _____ ______ _______ _ ____   _____  #
%  #         \ |  |   ___|   ___|__   __| |    \ |       #
%  #          \|  |  __| |  __|    | |  | |     \|       #
%  #       |\     | |____| |       | |  | |   |\         #
%  #  _____| \____| _____|_|       |_|  |_|___| \______  #
%  #                                                     #
%  #######################################################
