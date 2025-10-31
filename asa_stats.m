function [best_x, best_f, history_all] = asa_stats(func_name, max_iter, seed, num_runs)
if nargin < 1, func_name = 'Rastrigin'; end 
if nargin < 2, max_iter = 10000; end
if nargin < 3, seed = 1; end
if nargin < 4, num_runs = 10; end

[LB, UB, D, fobj] = get_function_details(func_name);

if isscalar(LB), LB = LB * ones(1,D); end
if isscalar(UB), UB = UB * ones(1,D); end

best_f = nan(num_runs,1);
best_x = nan(num_runs, D);
history_all = cell(num_runs,1);

fprintf('Pokrećem ASA eksperimente za funkciju %s: %d run-ova, max_iter=%d\n', func_name, num_runs, max_iter);

for run = 1:num_runs
    rng(seed + run - 1); 

    T0_param = ones(1,D);         
    c_param  = 1.0 * ones(1,D);   
    T0_cost  = 1.0;               
    c_cost   = 1.0 * ones(1,D);   

    Qi = 1;                       
    reanneal_every = 200;         
    accepted = 0; generated = 0;  

    x = LB + (UB - LB) .* rand(1,D);
    fx = fobj(x);
    best_x_run = x; best_f_run = fx;

    maxStore = max_iter + 1;
    hist_best_f = nan(maxStore,1);
    hist_best_x = nan(maxStore, D);
    idx_hist = 1;
    hist_best_f(idx_hist) = best_f_run;
    hist_best_x(idx_hist,:) = best_x_run;

    k = 1; 

    while generated < max_iter
        Tk_param = T0_param .* exp(-c_param .* (k).^(1/(D*Qi))); 

        T_cost_vec = T0_cost .* exp(-c_cost .* ( (accepted+1).^(1/(D*Qi)) ) );
        T_cost = mean(T_cost_vec);  

        u = rand(1,D);
        sgn = sign(u - 0.5);
        exponent = abs(2*u - 1);
        yi = sgn .* Tk_param .* ((1 + 1./Tk_param) .^ exponent - 1);

        xnew = x + yi .* (UB - LB);
        xnew = min(UB, max(LB, xnew)); 
        fx_new = fobj(xnew);
        delta = fx_new - fx;

        if delta <= 0
            accept = true;
        else
            p = exp(-delta / max(T_cost, 1e-12));
            accept = (rand() < p);
        end

        generated = generated + 1;

        if accept
            x = xnew; fx = fx_new; accepted = accepted + 1;
            if fx < best_f_run
                best_f_run = fx; best_x_run = x;
            end

            if mod(accepted, reanneal_every) == 0
                eps_fd = 1e-4 * (UB - LB);
                grads = zeros(1,D);
                for ii = 1:D
                    xp = best_x_run; xm = best_x_run;
                    xp(ii) = min(best_x_run(ii) + eps_fd(ii), UB(ii));
                    xm(ii) = max(best_x_run(ii) - eps_fd(ii), LB(ii));
                    denom = (xp(ii) - xm(ii));
                    if abs(denom) < eps, denom = eps; end
                    grads(ii) = (fobj(xp) - fobj(xm)) / denom;
                end
                s = abs(grads); smax = max(s) + eps;
                scale = smax ./ (s + eps);
                c_param = c_param ./ scale;
                c_cost  = c_cost  ./ scale;
            end
        end

        idx_hist = idx_hist + 1;
        hist_best_f(idx_hist) = best_f_run;
        hist_best_x(idx_hist,:) = best_x_run;

        k = k + 1;
    end

    hist_best_f = hist_best_f(1:idx_hist);
    hist_best_x = hist_best_x(1:idx_hist,:);

    best_f(run) = best_f_run;
    best_x(run,:) = best_x_run;
    history_all{run}.best_f = hist_best_f;
    history_all{run}.best_x = hist_best_x;
    history_all{run}.accepted = accepted;
    history_all{run}.generated = generated;

    fprintf('Run %2d završen: best_f = %.6e, accepted=%d, generated=%d\n', run, best_f_run, accepted, generated);
end

min_val = min(best_f);
max_val = max(best_f);
mean_val = mean(best_f);
std_val  = std(best_f);

fprintf('\n=== REZULTATI ASA (po run-ovima) ===\n');
fprintf('Funkcija: %s\n', func_name);
fprintf('Broj dimenzija: %d\n', D);
fprintf('Broj eksperimenata: %d\n', num_runs);
fprintf('Broj generisanih rešenja po run-u (max_iter): %d\n', max_iter);
fprintf('-------------------------------------\n');
fprintf('Minimum funkcije (po svim run-ovima): %.6e\n', min_val);
fprintf('Maksimum funkcije: %.6e\n', max_val);
fprintf('Srednja vrednost: %.6e\n', mean_val);
fprintf('Standardna devijacija: %.6e\n', std_val);
fprintf('=====================================\n\n');

figure('Name',['ASA - Results per experiment (', func_name, ')'],'NumberTitle','off','Color','w');
subplot(1,2,1);
plot(1:num_runs, best_f, '-o', 'LineWidth', 1.5, 'MarkerSize',8);
xlabel('Eksperiment');
ylabel('Najbolja vrednost funkcije (final)');
title(['Rezultati po eksperimetu (', func_name, ')']);
grid on;

subplot(1,2,2);
boxplot(best_f);
ylabel('Najbolja vrednost funkcije');
title(['Boxplot rezultata (', func_name, ')']);
grid on;

if ~isempty(history_all) && ~isempty(history_all{1}.best_f)
    figure('Name',['Konvergencija (prvi run) - ', func_name],'NumberTitle','off','Color','w');
    plot(history_all{1}.best_f, 'LineWidth', 1.5);
    xlabel('Iteracija');
    ylabel('Najbolja vrednost funkcije');
    title(['Konvergencija (prvi eksperimenat) - ', func_name]);
    grid on;
end

end
