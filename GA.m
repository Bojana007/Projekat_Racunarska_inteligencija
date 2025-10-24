clc; clear all; close all;

syms x y n;
syms i_d i_g;

n = 100;
n_elite = 5;
n_crossover = 90; 
n_mutation = 2;

[population] = coding( -5, 5, n );

vec = [];

for i = 1:10
    
    [minimal_value] = Rosenbrock( population, n, n_elite, n_crossover, n_mutation );

    if minimal_value < -5
        greska = 'MB';
    end

    vec = [vec; minimal_value];
end

vec

figure(1)
plot(vec)
yline(0, 'r--') 
hold on;
figure(2)
boxplot( vec )

function z = fun(x, y)
    z = ( x - 1 ).^2 + 100 * ( y - x.^2 ).^2;
end


function[ population ] = coding( i_d, i_g, n ) 
    
    for i = 1:n
        for j = 1:2
            population(i,j) = i_d + ( i_g - i_d ) * rand();
        end
    end
end


function[ selected_individuals, elite_individuals ] = selection( population, n, number_of_elite ) 

   fitness_values = zeros(1,n);
   probability = zeros(1,n); 
   cumulative_sum = zeros(n,2); 
   sum = 0;

   for i = 1:n
        fitness_values(i) = fun( population(i,1), population(i,2));
        fitness_values(i) = 1 / fitness_values(i);   
        sum = sum + fitness_values(i);
   end

   for i = 1:n
        probability(i) = fitness_values(i)  / sum;
   end

   probability = round( probability, 3 );

   elite_individuals = zeros(number_of_elite,2);

   [ ~, idx ] = sort( fitness_values, "descend" );
   
   for i = 1:number_of_elite 
        
       elite_individuals(i,:) = population( idx(i), : ); 
       
   end
 
   new_size = n - number_of_elite;

   cumulative_sum(:,1) = cumsum(probability);
       
   for i = 1:new_size 
      random_num = rand();
      index = find( cumulative_sum(:,1) >= random_num, 1 );
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

function[ selected_individuals ] = crossover( selected_individuals, n, n_crossover, population )

    num_of_crossovered = round( n * ( n_crossover / 100 ) );

    for i = 1:2:num_of_crossovered - 2
    
        crossover_point = randi([1, 2]);
         
        temp = selected_individuals(i, crossover_point);
        selected_individuals(i, crossover_point) = selected_individuals(i + 1, crossover_point);
        selected_individuals(i + 1, crossover_point) = temp;
        
    end
    
    [other_individuals] = coding(-5,5, n - num_of_crossovered);

    j = 1;
    for i = num_of_crossovered+1:n
        selected_individuals( i, : ) = other_individuals( j, : );
        j = j + 1;
    end

end

function[ mutated_individuals  ] = mutation( crossovered_individuals, n, M ) 
    
    permutations = randperm(n, M); 
    randoms = permutations(1:2);

    for i = 1:M 
        crossovered_individuals(randoms(i),1) = crossovered_individuals(randoms(i),1) + randn; 
        crossovered_individuals(randoms(i),2) = crossovered_individuals(randoms(i),2) + randn;     
    end

    mutated_individuals = crossovered_individuals;

end

function[ min_value ] = Rosenbrock( population, n, n_elite, n_crossover, n_mutation )
    
    counter = 1;
    
    while counter <= 100000
            
        [ selected_individuals, elite_individuals ] = selection( population, n, n_elite );
        
        new_n = n - n_elite;
            
        [ crossovered_individuals ] = crossover( selected_individuals, new_n, n_crossover, population );
        [ mutated_individuals ] = mutation( crossovered_individuals, new_n, n_mutation );
           
        population = vertcat(mutated_individuals, elite_individuals);

        for i = 1:n
            fitness_values(i) = fun( population(i,1), population(i,2));
        end
         
        min_value = min( fitness_values );
        minimalne(counter) = min_value;

        counter = counter + 1;
    end

    min_value = min(minimalne);
end