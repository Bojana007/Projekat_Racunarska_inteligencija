clear all; clc; close all;

func_name = 'Schwefel';
[lb, ub, n, fobj] = get_function_details(func_name);

num_runs = 10;

baseline.T0 = 350;
baseline.Tc = 0.993;
baseline.Tmin = 0.1;

T0_list    = [50, 100, 200, 350, 700];
Tmin_list  = [0.001, 0.01, 0.05, 0.1, 0.5];
Tc_list    = [0.90, 0.95, 0.98, 0.993, 0.997];

results = cell(0,4);
row = 1;

[param_cells, best_params] = run_group(@(v,b)set_T0(v,b), T0_list, 'Inicijalna temperatura T0', baseline, num_runs, lb, ub, n, fobj);
results(row:row+4,:) = param_cells;
row = row + 5;
baseline = best_params;

[param_cells, best_params] = run_group(@(v,b)set_Tmin(v,b), Tmin_list, 'Minimalna temperatura Tmin', baseline, num_runs, lb, ub, n, fobj);
results(row:row+4,:) = param_cells;
row = row + 5;
baseline = best_params;

[param_cells, ~] = run_group(@(v,b)set_Tc(v,b), Tc_list, 'Faktor hlađenja Tc', baseline, num_runs, lb, ub, n, fobj);
results(row:row+4,:) = param_cells;

T = cell2table(results, 'VariableNames', {'Parametri','Najbolja','Srednja','Std'});
disp(T);
writetable(T, 'SA_results_simple.csv');

disp('Gotovo. CSV fajl sacuvan kao SA_results_simple.csv');

function [group_results, best_params] = run_group(setter, values, title_str, baseline, num_runs, lb, ub, n, fobj)
    m = numel(values);
    colors = lines(m);
    figure('Name', title_str, 'Color','w');
    hold on; grid on;
    legend_entries = strings(1,m);
    group_results = cell(m,4);
    best_params = baseline;
    best_mean = inf;

    for i = 1:m
        v = values(i);
        params = setter(v, baseline);
        [param_str, best_val, mean_val, std_val, mean_curve] = run_experiment(params, num_runs, lb, ub, n, fobj);
        group_results(i,:) = {param_str, best_val, mean_val, std_val};
        legend_entries(i) = string(v);

        plot(mean_curve, 'LineWidth', 2.2, 'Color', colors(i,:));

        if mean_val < best_mean
            best_mean = mean_val;
            best_params = params;
        end
    end

    xlabel('Iteracija');
    ylabel('Najbolja fitnes vrednost');
    title(['Konvergencija - ', title_str]);
    legend(legend_entries, 'Location','best');
    hold off;
end

function [param_str, best_val, mean_val, std_val, mean_curve] = run_experiment(params, num_runs, lb, ub, n, fobj)
    max_len = 0;
    all_curves = cell(1,num_runs);
    best_values = zeros(num_runs,1);

    for i = 1:num_runs
        curr = lb + (ub - lb).*rand(1,n);
        curve = SA_run(curr, params.T0, params.Tc, params.Tmin, n, lb, ub, fobj);
        all_curves{i} = curve;
        best_values(i) = min(curve);
        max_len = max(max_len, numel(curve));
    end

    M = zeros(num_runs, max_len);
    for i = 1:num_runs
        c = all_curves{i};
        M(i,1:numel(c)) = c;
        if numel(c) < max_len
            M(i,numel(c)+1:end) = c(end);
        end
    end
    mean_curve = mean(M,1);
    mean_val = mean(best_values);
    std_val = std(best_values);
    best_val = min(best_values);
    param_str = sprintf('T0=%.3g, Tmin=%.3g, Tc=%.3g', params.T0, params.Tmin, params.Tc);
end

function curve = SA_run(curr, T, Tc, Tmin, n, lb, ub, fobj)
    curve = [];
    while T > Tmin
        curr = max(curr, lb);
        curr = min(curr, ub);
        temp = curr + randn(1,n);
        dE = fobj(temp) - fobj(curr);
        if dE < 0 || rand < exp(-dE/T)
            curr = temp;
        end
        T = T * Tc;
        curve(end+1) = fobj(curr);
    end
end

function s = set_T0(v,b), s=b; s.T0=v; end
function s = set_Tmin(v,b), s=b; s.Tmin=v; end
function s = set_Tc(v,b), s=b; s.Tc=v; end
