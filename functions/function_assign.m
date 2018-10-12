
%fplot    = @(z) log(1+z);

switch fun_num
    
    case '1'
        r  = Rosenbrock();
    case '4'
        r  = Brown_badly_scaled_function();
    case '7'
        r  = Helical_valley();
    case '10'
        r  = Meyer_function();
        approx_minima = r.approximated_minima;
    case '13'
        r  = Powell_singular_function();
    case '16'
        r  = Brown_Dennis_function(M);
        approx_minima = r.approximated_minima;
    case '19'
        r  = Osborne_2_function();
        approx_minima = r.approximated_minima;
        
end


x0 = r.guesses;
goal = r.exact_solutions;