clear all; clc; close all;

func_name = 'Schwefel';  

[lb, ub, n, fobj] = get_function_details(func_name);

T = 350;
Tc = 0.993;
Tmin = 0.1;

curr = lb + (ub - lb) .* rand(1, n);

vec = [];

for i = 1:10

    T = 350;
    Tc = 0.993;
    Tmin = 0.1;

    [curr, output, vektor] = SA_4(curr, T, Tc, Tmin, n, lb, ub, fobj);

    vec = [vec; output];
end

figure(1)
plot(vec)
title(['Simulated Annealing - ', func_name])
yline(min(vec), 'r--')
xlabel('Iteracija')
ylabel('Najbolja fitnes vrednost')
hold on;

figure(2)
boxplot(vec)
title(['Distribucija rezultata - ', func_name])

minimalno = min(vec)
maximalno = max(vec)
srednja_vr = mean(vec)
sd = std(vec)

fprintf('\n*** Statistika za funkciju %s ***\n', func_name);
fprintf('Minimalna vrednost: %.6f\n', minimalno);
fprintf('Maksimalna vrednost: %.6f\n', maximalno);
fprintf('Srednja vrednost:   %.6f\n', srednja_vr);
fprintf('Standardna devijacija: %.6f\n', sd);

function [curr, output, vektor] = SA_4(curr, T, Tc, Tmin, n, lb, ub, fobj)

    vektor = [];

    while (T > Tmin)

        curr = max(curr, lb);
        curr = min(curr, ub);

        temp = curr + randn(1, n);
        
        f_temp = fobj(temp);
        f_curr = fobj(curr);

        d_E = f_temp - f_curr;

        if d_E < 0
            curr = temp;
        else
            r = rand;
            if r < exp(-d_E / T)
                curr = temp;
            end
        end

        T = Tc * T;
        output = fobj(curr);
        vektor = [vektor, output];
    end
end
