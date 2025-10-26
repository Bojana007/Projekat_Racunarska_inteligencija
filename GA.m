clc; clear all; close all;

func_name = 'Rosenbrock';
[lb, ub, dim, fobj] = get_function_details(func_name);

n = 100;
n_elite = 5;
n_crossover = 90; 
n_mutation = 2;

[population] = coding(lb, ub, n);

vec = [];

for i = 1:10
    [minimal_value] = GA_function(population, n, n_elite, n_crossover, n_mutation, fobj, lb, ub);
    vec = [vec; minimal_value];
end

minimalno = min(vec);
maksimalno = max(vec);
srednja_vr = mean(vec);
sd = std(vec);

fprintf('\n*** Statistika za funkciju %s ***\n', func_name);
fprintf('Minimalna vrednost: %.6f\n', minimalno);
fprintf('Maksimalna vrednost: %.6f\n', maksimalno);
fprintf('Srednja vrednost:   %.6f\n', srednja_vr);
fprintf('Standardna devijacija: %.6f\n', sd);

disp(vec)

figure(1)
plot(vec)
yline(0, 'r--') 
title(['Optimizacija funkcije: ', func_name])
xlabel('Iteracija')
ylabel('Minimalna vrednost')
hold on;

figure(2)
boxplot(vec)
title(['Distribucija minimalnih vrednosti - ', func_name])

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

   fitness_values = zeros(1,n);
   probability = zeros(1,n); 
   cumulative_sum = zeros(n,2); 
   sum_val = 0;

   for i = 1:n
        fitness_values(i) = fobj(population(i,:));
        fitness_values(i) = 1 / fitness_values(i);   
        sum_val = sum_val + fitness_values(i);
   end

   for i = 1:n
        probability(i) = fitness_values(i)  / sum_val;
   end

   probability = round(probability, 3);

   elite_individuals = zeros(number_of_elite, 2);

   [~, idx] = sort(fitness_values, "descend");
   
   for i = 1:number_of_elite 
       elite_individuals(i,:) = population(idx(i), :); 
   end
 
   new_size = n - number_of_elite;
   cumulative_sum(:,1) = cumsum(probability);
       
   for i = 1:new_size 
      random_num = rand();
      index = find(cumulative_sum(:,1) >= random_num, 1);
       if ~isempty(index)
        cumulative_sum(index, 2) = cumulative_sum(index, 2) + 1;   
       end
   end
    
  current_index = 1; 

  for i = 1:n
    counter = cumulative_sum(i, 2);
    if counter ~= 0
        for j = 1:counter
            selected_individuals(current_index, :) = population(i, :);
            current_index = current_index + 1;
        end   
    end
  end
end


function [selected_individuals] = crossover(selected_individuals, n, n_crossover, lb, ub)
    num_of_crossovered = round(n * (n_crossover / 100));

    for i = 1:2:num_of_crossovered - 2
        crossover_point = randi([1, 2]);
        temp = selected_individuals(i, crossover_point);
        selected_individuals(i, crossover_point) = selected_individuals(i + 1, crossover_point);
        selected_individuals(i + 1, crossover_point) = temp;
    end
    
    [other_individuals] = coding(lb, ub, n - num_of_crossovered);

    j = 1;
    for i = num_of_crossovered+1:n
        selected_individuals(i,:) = other_individuals(j,:);
        j = j + 1;
    end
end


function [mutated_individuals] = mutation(crossovered_individuals, n, M)
    permutations = randperm(n, M);
    randoms = permutations(1:2);

    for i = 1:M 
        crossovered_individuals(randoms(i),1) = crossovered_individuals(randoms(i),1) + randn; 
        crossovered_individuals(randoms(i),2) = crossovered_individuals(randoms(i),2) + randn;     
    end

    mutated_individuals = crossovered_individuals;
end


function [min_value] = GA_function(population, n, n_elite, n_crossover, n_mutation, fobj, lb, ub)
    counter = 1;
    
    while counter <= 100
        [selected_individuals, elite_individuals] = selection(population, n, n_elite, fobj);
        
        new_n = n - n_elite;
        [crossovered_individuals] = crossover(selected_individuals, new_n, n_crossover, lb, ub);
        [mutated_individuals] = mutation(crossovered_individuals, new_n, n_mutation);
           
        population = vertcat(mutated_individuals, elite_individuals);

        for i = 1:n
            fitness_values(i) = fobj(population(i,:));
        end
         
        min_value = min(fitness_values);
        minimalne(counter) = min_value;

        counter = counter + 1;
    end

    min_value = min(minimalne);
end
