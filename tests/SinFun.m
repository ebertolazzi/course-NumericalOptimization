classdef SinFun < Function1D
  %
  % The class DOES NOT implement eval_D and eval_DD
  % so that it uses the inherited eval_D and eval_DD
  % which compute derivative by finite difference
  %
  methods

    function self = SinFun()
    end

    function y = eval( self, x )
      y = sin(x);
    end

  end
end
