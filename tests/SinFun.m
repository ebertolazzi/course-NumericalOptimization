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

    function Dy = eval_D( self, x )
      disp('Use FD first derivative for SinFun');
      Dy = self.FD_eval_D(x);
    end

    function DDy = eval_DD( self, x )
      disp('Use FD second derivative for SinFun');
      DDy = self.FD_eval_DD(x);
    end

  end
end
