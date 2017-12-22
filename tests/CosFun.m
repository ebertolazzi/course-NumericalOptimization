classdef CosFun < Function1D
  %
  % The class DOES implement eval_D and eval_DD
  % so that it overrides the inherited eval_D and eval_DD.
  %
  methods

    function self = CosFun()
    end

    function y = eval( self, x )
      y = cos(x);
    end

    function y = eval_D( self, x )
      disp('Use analitic first derivative for CosFun');
      y = -sin(x);
    end

    function y = eval_DD( self, x )
      disp('Use analitic second derivative for CosFun');
      y = -cos(x);
    end

  end
end
