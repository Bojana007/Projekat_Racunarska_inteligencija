clc; clear; close all;

func_name = 'Rastrigin';
num_runs = 10; 

baseline.T0_param = 1.0;
baseline.c_param  = 1.0;
baseline.c_cost   = 1.0;
baseline.T0_cost  = 1.0;
baseline.Qi       = 1.0;
baseline.reanneal_every = 200;
baseline.max_iter = 10000;

T0_cost_vals    = [0.5, 1.0, 2.0];
Qi_vals         = [1, 2, 5];
reanneal_vals   = [100, 200, 500];
max_iter_vals   = [5000, 10000, 20000];

results = cell(0,5);
row = 1;

[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'T0_cost',v), T0_cost_vals, 'T0\_cost', func_name, num_runs);
results(row:row+length(T0_cost_vals)-1,:) = group_stats;
row = row + length(T0_cost_vals);
baseline = best_params;

[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'Qi',v), Qi_vals, 'Qi', func_name, num_runs);
results(row:row+length(Qi_vals)-1,:) = group_stats;
row = row + length(Qi_vals);
baseline = best_params;

[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'reanneal_every',v), reanneal_vals, 'reanneal\_every', func_name, num_runs);
results(row:row+length(reanneal_vals)-1,:) = group_stats;
row = row + length(reanneal_vals);
baseline = best_params;

[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'max_iter',v), max_iter_vals, 'max\_iter', func_name, num_runs);
results(row:row+length(max_iter_vals)-1,:) = group_stats;
row = row + length(max_iter_vals);
baseline = best_params;

T = cell2table(results, 'VariableNames', {'Parameters','Max','Min','Mean','Std'});
disp(T);
writetable(T, 'ASA_experiments_results.csv');
fprintf('Gotovo. CSV sacuvan kao ASA_experiments_results.csv\n');

function [group_results, best_params] = run_group_and_plot(param_setter, values, group_title, func_name, num_runs)
    num_vals = length(values);
    all_mean_curves = cell(1,num_vals);
    finals_mean = zeros(1,num_vals);
    group_results = cell(num_vals,5);
    best_mean = inf;
    best_params = [];

    for k = 1:num_vals
        v = values(k);
        params = param_setter(v);
        [param_str, maxv, minv, meanv, stdv, mean_curve] = run_experiment_ASA(params, func_name, num_runs);
        group_results(k,:) = {param_str, maxv, minv, meanv, stdv};
        all_mean_curves{k} = mean_curve;
        finals_mean(k) = meanv;

        if meanv < best_mean
            best_mean = meanv;
            best_params = params;
        end
    end

    figure('Name',['Convergence - ', group_title],'Color','w'); hold on; grid on;
    max_len = max(cellfun(@length, all_mean_curves));
    for k = 1:num_vals
        curve = all_mean_curves{k};
        if length(curve) < max_len
            curve = [curve, repmat(curve(end),1,max_len-length(curve))];
        end
        plot(1:max_len, curve, 'LineWidth', 1.5);
    end
    legend(arrayfun(@(x) sprintf('%g', values(x)),1:num_vals,'UniformOutput',false));
    xlabel('Iteracija'); ylabel('Najbolja vrednost funkcije');
    title(['Konvergencija: ', group_title]);
    grid on;

    figure('Name',['Final Mean - ', group_title],'Color','w'); hold on; grid on;
    bar(finals_mean);
    set(gca,'XTickLabel',arrayfun(@(x) sprintf('%g',values(x)),1:num_vals,'UniformOutput',false));
    title(['Konačna srednja vrednost - ', group_title]);
    xlabel(group_title); ylabel('Srednja minimalna vrednost (lower = better)');
end

function [param_str, maxv, minv, meanv, stdv, mean_curve] = run_experiment_ASA(params, func_name, num_runs)
    best_f_all = zeros(num_runs,1);
    history_all = cell(num_runs,1);

    for r = 1:num_runs
        [res, history] = ASA_single_run(func_name, params);
        best_f_all(r) = res;
        history_all{r} = history;
    end

    mean_curve = mean(cell2mat(cellfun(@(c) c(:)', history_all, 'UniformOutput', false)),1);

    maxv = max(best_f_all);
    minv = min(best_f_all);
    meanv = mean(best_f_all);
    stdv = std(best_f_all);

    param_str = sprintf('T0=%.2f, c=%.2f, T0_cost=%.2f, c_cost=%.2f, Qi=%.2f, reanneal=%d, max_iter=%d', ...
        params.T0_param, params.c_param, params.T0_cost, params.c_cost, params.Qi, params.reanneal_every, params.max_iter);
end

function [best_value, history] = ASA_single_run(func_name, params)
    [lb, ub, dim, fobj] = get_function_details(func_name);

    max_iter = params.max_iter;
    history = zeros(1,max_iter);

    x = lb + (ub-lb).*rand(1,dim);
    best_value = fobj(x);

    T = params.T0_param;
    for k = 1:max_iter
        x_new = lb + (ub-lb).*rand(1,dim);
        f_new = fobj(x_new);

        if f_new < best_value || rand() < exp(-(f_new-best_value)/T)
            x = x_new;
            best_value = f_new;
        end

        T = T * 0.99;

        history(k) = best_value;
    end
end
