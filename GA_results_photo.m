clc; clear; close all;

func_name = 'Schwefel';
[lb, ub, dim, fobj] = get_function_details(func_name);

num_runs = 10; 

baseline.n = 100;
baseline.max_gen = 100;
baseline.crossover_pct = 90;  
baseline.mutation_rate = 0.02; 
baseline.n_elite = 5;

pop_sizes = [50, 100, 200, 400, 800];
max_gens  = [50, 100, 200, 400, 800];
crossover_pcts = [50, 70, 80, 90, 100];
mutation_rates  = [0.01, 0.03, 0.05, 0.1, 0.2];
n_elite_list    = [1, 2, 5, 10, 20];

results = cell(0,5);
row = 1;

group_title = 'Veličina populacije';
[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'n',v), pop_sizes, group_title);
results(row:row+4,:) = group_stats;
row = row + 5;
baseline = best_params; 

group_title = 'Maksimalan broj generacija';
[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'max_gen',v), max_gens, group_title);
results(row:row+4,:) = group_stats;
row = row + 5;
baseline = best_params;

group_title = 'Procenat ukrštanja (%)';
[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'crossover_pct',v), crossover_pcts, group_title);
results(row:row+4,:) = group_stats;
row = row + 5;
baseline = best_params;

group_title = 'Verovatnoća mutacije (rate)';
[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'mutation_rate',v), mutation_rates, group_title);
results(row:row+4,:) = group_stats;
row = row + 5;
baseline = best_params;

group_title = 'Broj elitnih hromozoma';
[group_stats, best_params] = run_group_and_plot(@(v) setfield(baseline,'n_elite',v), n_elite_list, group_title);
results(row:row+4,:) = group_stats;
row = row + 5;
baseline = best_params;

T = cell2table(results, 'VariableNames', {'Parameters','Max','Min','Mean','Std'});
disp(T);
writetable(T, 'GA_experiments_results.csv');
fprintf('Gotovo. CSV sacuvan kao GA_experiments_results.csv\n');

function [group_results, best_params] = run_group_and_plot(param_setter, values, group_title)
    num_vals = length(values);
    all_mean_curves = cell(1,num_vals);
    all_std_curves  = cell(1,num_vals);
    finals_mean = zeros(1,num_vals);
    group_results = cell(num_vals,5);
    best_mean = inf;
    best_params = [];
    for k = 1:num_vals
        v = values(k);
        params = param_setter(v);
        [pstr, mx, mn, me, sd, mean_curve, std_curve] = run_experiment(params, evalin('base','num_runs'), evalin('base','lb'), evalin('base','ub'), evalin('base','fobj'));
        group_results(k,:) = {pstr, mx, mn, me, sd};
        all_mean_curves{k} = mean_curve;
        all_std_curves{k}  = std_curve;
        finals_mean(k) = me;
        if me < best_mean
            best_mean = me;
            best_params = params;
        end
    end
    figure('Name', ['Convergence - ', group_title], 'NumberTitle', 'off');
    hold on;
    max_len = 0;
    for k = 1:num_vals
        lenk = length(all_mean_curves{k});
        if lenk > max_len, max_len = lenk; end
    end
    for k = 1:num_vals
        curve = all_mean_curves{k};
        if length(curve) < max_len
            curve = [curve, repmat(curve(end), 1, max_len-length(curve))];
        end
        plot(1:max_len, curve, 'LineWidth', 1.5);
    end
    legend_strings = arrayfun(@(x) sprintf('%g', values(x)), 1:num_vals, 'UniformOutput', false);
    legend(legend_strings, 'Location', 'best');
    title(['Konvergencija: ', group_title]);
    xlabel('Generacija');
    ylabel('Minimalna vrednost (srednja preko ponavljanja)');
    grid on;
    hold off;

    % Bar plot final mean values
    figure('Name', ['Final Mean - ', group_title], 'NumberTitle', 'off');
    bar(finals_mean);
    set(gca, 'XTickLabel', legend_strings);
    title(['Konačna srednja vrednost (posle svih generacija) - ', group_title]);
    xlabel(group_title);
    ylabel('Srednja minimalna vrednost (lower = better)');
    grid on;
end

function [param_str, maxv, minv, meanv, stdv, mean_curve, std_curve] = run_experiment(params, num_runs, lb, ub, fobj)
    curves = [];
    for r = 1:num_runs
        population = coding(lb, ub, params.n);
        curve = GA_function(population, params.n, params.n_elite, params.crossover_pct, params.mutation_rate, params.max_gen, fobj, lb, ub);
        curves(r,1:length(curve)) = curve; %#ok<AGROW>
    end
    max_len = size(curves,2);
    for r = 1:size(curves,1)
        if any(isnan(curves(r,:)))
            last = find(~isnan(curves(r,:)),1,'last');
            curves(r,last+1:end) = curves(r,last);
        end
    end
    mean_curve = mean(curves,1);
    std_curve  = std(curves,0,1);
    final_vals = min(curves,[],2); 
    maxv = max(final_vals);
    minv = min(final_vals);
    meanv = mean(final_vals);
    stdv  = std(final_vals);
    param_str = sprintf('n=%d,gens=%d,cross=%d%%,mut=%.3f,elite=%d', ...
        params.n, params.max_gen, params.crossover_pct, params.mutation_rate, params.n_elite);
end

function [population] = coding(lb, ub, n)
    dim = length(lb);
    population = zeros(n, dim);
    for i = 1:n
        for j = 1:dim
            population(i,j) = lb(j) + (ub(j) - lb(j)) * rand();
        end
    end
end

function [selected_individuals, elite_individuals] = selection(population, n, number_of_elite, fobj)
    dim = size(population,2);
    fitness_values = zeros(1,n);
    for i = 1:n
        fitness_values(i) = fobj(population(i,:));
    end
    worst = max(fitness_values);
    scores = worst - fitness_values + eps;
    probs = scores / sum(scores);
    [~, idx_sorted] = sort(scores, 'descend');
    elite_individuals = zeros(number_of_elite, dim);
    for k = 1:number_of_elite
        elite_individuals(k,:) = population(idx_sorted(k), :);
    end
    to_select = n - number_of_elite;
    cumulative = cumsum(probs);
    selected_individuals = zeros(to_select, dim);
    for k = 1:to_select
        r = rand();
        ind = find(cumulative >= r, 1);
        if isempty(ind), ind = n; end
        selected_individuals(k,:) = population(ind,:);
    end
end

function out = crossover(selected_individuals, new_n, n_crossover_pct, lb, ub)
    dim = size(selected_individuals,2);
    num_of_crossovered = round(new_n * (n_crossover_pct / 100));
    if num_of_crossovered > new_n, num_of_crossovered = new_n; end
    if mod(num_of_crossovered,2) == 1, num_of_crossovered = max(0, num_of_crossovered - 1); end
    out = selected_individuals;
    for i = 1:2:num_of_crossovered
        p1 = out(i,:);
        p2 = out(i+1,:);
        if dim > 1
            cp = randi([1, dim-1]);
            c1 = [p1(1:cp), p2(cp+1:end)];
            c2 = [p2(1:cp), p1(cp+1:end)];
        else
            c1 = p2; c2 = p1;
        end
        out(i,:) = c1;
        out(i+1,:) = c2;
    end
    if num_of_crossovered < new_n
        rem = new_n - num_of_crossovered;
        extra = zeros(rem, dim);
        for i = 1:rem
            for j = 1:dim
                extra(i,j) = lb(j) + (ub(j) - lb(j)) * rand();
            end
        end
        startpos = num_of_crossovered + 1;
        out(startpos:new_n, :) = extra;
    end
end

function mutated_individuals = mutation(crossovered_individuals, new_n, mutation_rate)
    dim = size(crossovered_individuals,2);
    mutated_individuals = crossovered_individuals;
    num_to_mutate = max(1, round(new_n * mutation_rate));
    idxs = randperm(new_n, num_to_mutate);
    for k = 1:length(idxs)
        i = idxs(k);
        num_genes = randi([1, dim]);
        genes = randperm(dim, num_genes);
        for g = genes
            mutated_individuals(i,g) = mutated_individuals(i,g) + randn();
        end
    end
end

function minimalne = GA_function(population, n, n_elite, n_crossover_pct, mutation_rate, max_gen, fobj, lb, ub)
    if n_elite >= n, n_elite = max(1, n-1); end
    current_pop = population;
    minimalne = nan(1, max_gen);
    for gen = 1:max_gen
        [selected_individuals, elite_individuals] = selection(current_pop, n, n_elite, fobj);
        new_n = n - n_elite;
        crossovered_individuals = crossover(selected_individuals, new_n, n_crossover_pct, lb, ub);
        mutated_individuals = mutation(crossovered_individuals, new_n, mutation_rate);
        current_pop = vertcat(mutated_individuals, elite_individuals);
        fitness_values = zeros(1,n);
        for i = 1:n
            fitness_values(i) = fobj(current_pop(i,:));
        end
        minimalne(gen) = min(fitness_values);
    end
end

