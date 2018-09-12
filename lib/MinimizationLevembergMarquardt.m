classdef MinimizationLevembergMarquardt < MinimizationND
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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationLevembergMarquardt( fun )
      self@MinimizationND( fun, [] ); % Linesearch is empty -> no linesearch method in LM
      self.tol2 = 1e-6;
      self.tau  = 1e-6;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setTau( self , tau )
      % Change tau of the LM algorithm:
      % Note the higher the tau the higher the confidence to start close to the minimum:
      % see ref[1] pag 27
      if tau <= 0
        error( 'MinimizationLM::setTau(%2.6g) argument musy be positive\n', tau );
      end
      self.tau = tau;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setEpsilon2( self , tolerance )
      % Change tolerance epsilon 2 (see ref[1]):
      if tolerance <= 0
        error( 'MinimizationLM::setEpsilon2(%2.6g) argument musy be positive\n', tolerance );
      end
      self.tol2 = tolerance;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,A] = eval_FGA( self, x )
      use_map = self.funND.is_a_map();
      f = self.funND.eval(x);
      if use_map
        J = self.funND.jacobian(x);
        g = J.' * self.funND.evalMap(x);
        A = J.' * J;                % first term of Hessian
      else
        g = self.funND.grad(x).';
        A = self.funND.hessian(x);
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, converged ] = minimize( self , x0 )
      % Launch the minimization algorithm
      % Initial values -> see ref[1] algorithm 3.16
      x           = x0;    % initial point
      converged   = false; % initialization
      [ f, g, A ] = self.eval_FGA(x);
      mi          = self.tau * max(abs(diag(A))); % First mi (ref[1])
      ni          = 2;                       % ni -> multiplicative factor of mi
      eps2        = self.tol2;               % tolerance on the error

      % check convercence and assign the boolean converged.
      norm_inf_g = norm(g,inf);      % infinite norm (max)
      converged  = norm_inf_g < self.tol; % Is it already converged?

      if converged
        if self.verbose
          fprintf( 'starting point is already a minimum, ||grad f||_inf = %12.6g\n', ...
                    norm_inf_g );
        end
        return;
      end

      if self.save_iterate
        self.x_history = x(:); % create the vector for history
      end

      % == START ITERATION ======
      for self.iter = 1:self.max_iter

        if self.verbose
          fprintf('[LM] iter = %5d ||grad f||_inf = %12.6g ni = %5d , mi = %12.6g\n', ...
                  self.iter, norm_inf_g, ni, mi );
        end

        if self.save_iterate
          self.x_history = [self.x_history x(:)];
        end

        % Solve( (J'J - uI)h = - grad )
        h = -( A + mi * eye(length(x)) ) \ g;

        % Condition on the h
        norm_2_h = norm(h);
        eps2_x_eps2 = eps2*( norm(x) + eps2  );

        if ( norm_2_h <= eps2_x_eps2 )
          if self.verbose
            fprintf( 2, 'step too small (epsilon 2 reached), ||h|| = %g < %g = epsi*||x||\n', ...
                     norm_2_h, eps2_x_eps2 );
          end
          break; % out of the loop
        end

        % Update x
        x_new = x + h;

        % compute L0 - L(h)
        L0_Lh = -dot(h, 0.5*(A*h) + g ); % eq(3.7b) version at page 25 of ref [1]

        % Gain ratio (rho)
        rho = (f - self.funND.eval(x_new)) / L0_Lh;
        if rho > 0 && isfinite(self.funND.eval(x_new)) % Acceptable step
          x = x_new; % x update
          [ f, g, A ] = self.eval_FGA(x);

          norm_inf_g = norm(g,inf);       % infinite norm (max)
          converged  = norm_inf_g <= self.tol; % check convergence
          if converged % Converged: stop the algorithm
            if self.verbose
              fprintf( 'solution found, ||grad f||_inf = %g < %g\n', ...
                        norm_inf_g, self.tol );
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
      % == END ITERATION ========
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
