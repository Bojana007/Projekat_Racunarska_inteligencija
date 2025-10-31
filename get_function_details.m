function [lb, ub, dim, fobj] = get_function_details(func_name)
switch func_name
    case 'Rosenbrock'
        fobj = @(x) (x(1)-1).^2 + 100*(x(2)-x(1).^2).^2;
        lb = [-5, -5];
        ub = [5, 5];
        dim = 2;

    case 'Rastrigin'
        fobj = @(x) 10*length(x) + sum(x.^2 - 10*cos(2*pi*x));
        lb = [-5.12, -5.12];
        ub = [5.12, 5.12];
        dim = 2;

    case 'Schwefel'
        fobj = @(x) 418.9829*length(x) - sum(x .* sin(sqrt(abs(x))));
        lb = [-500, -500];
        ub = [500, 500];
        dim = 2;

    otherwise
        error('Nepoznata funkcija: %s', func_name);
end
end
