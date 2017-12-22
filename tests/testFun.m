classdef testFun < Function1D
  %
  % The class DOES implement eval_D and eval_DD
  % so that it overrides the inherited eval_D and eval_DD.
  %
  methods

    function self = testFun()
    end

    function y = eval( self, x )
      y = exp(-3*x)+x.^2.*sin(3*x);
    end

  end
end
