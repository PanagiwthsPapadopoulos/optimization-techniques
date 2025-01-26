
function pop_matrix = initialize_population(pop_size, E, V)
    num_edges = size(E, 1);
    pop_matrix = zeros(pop_size, num_edges);  % Initialize matrix with zeros

    for i = 1:pop_size
        % Assign random values to source edges (no incoming edges)
        source_edges = [1,2,3,4];
        
        % Generate random values for source edges
        random_values = randi([1, 100], 1, length(source_edges));
        
        % Scale values to ensure their sum equals V
        normalized_values = round(random_values / sum(random_values) * V);
        
        % Assign to the population matrix
        pop_matrix(i, source_edges) = normalized_values;
        pop_matrix(i, length(source_edges)) = pop_matrix(i, 4) + (V - sum(pop_matrix(i, source_edges)));

        % Distribute the flow throughout the network
        for j = length(source_edges)+1:num_edges
            
            % Parent edges are edges that a car is coming from
            parent_edges = find(E(:, j) == 1);

            % Sibling edges are edges that a car can choose to go instead
            % of the edge j being checked
            sibling_edges = find(any(E(parent_edges, :) == 1, 1));
                
            % Calculate the flow of the sibling edges
            flow_of_siblings = sum(pop_matrix(i, sibling_edges));
            
            % Calculate the flow of the parent edges
            flow_of_parents = sum(pop_matrix(i, parent_edges));
            
            % If this is the last sibling edge, equal the flow of the
            % parents. Else, assign random value between 0 and the
            % remaining flow of the intersection
            if j == max(sibling_edges)
                pop_matrix(i, j) = flow_of_parents - flow_of_siblings;
            else
                pop_matrix(i, j) = randi([0, flow_of_parents-flow_of_siblings]);
            end
            
        end
    end
end

E = [
    %1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 1
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 2
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0] % 3
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0] % 4
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0] % 5
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0] % 6
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0] % 7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0] % 8
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0] % 9
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1] % 10
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1] % 11
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 12
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0] % 13
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0] % 14
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 15
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 16
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] % 17
    ];

source = [1,2,3,4];
sink = [12,15,16,17];
num_edges = 17;

C = [54.13,21.56,34.08,49.19,33.03,21.84,29.96,24.87,47.24,33.97,26.89,32.76,39.98,37.12,53.83,61.65,59.73];
a = [1.25,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,1.5,1,1,1,1,1,1,1];
V = 100;
t = 1;

generations = 200;  % Number of generations
pop_size = 100;
mutation_rate = 0.1;  % Probability of mutation

% Fitness function: minimize total travel time
fitness_function = @(x) sum(t + (a .* x) ./ (1 - (x ./ C)), 2);

pop_matrix = initialize_population(pop_size, E, V);
disp('Sample population matrix (first 5 solutions):');
disp(pop_matrix(1, :));




for gen = 1:generations
    % Evaluate fitness
    fitness = fitness_function(pop_matrix);

    % Selection (elitism)
    [sorted_fitness, idx] = sort(fitness);
    pop_matrix = pop_matrix(idx(1:pop_size/2), :);

    % Crossover (uniform crossover)
    offspring = pop_matrix;
    for i = 1:size(pop_matrix,1)/2
        alpha = rand();
        offspring(i,:) = alpha * pop_matrix(i,:) + (1-alpha) * pop_matrix(i+1,:);
    end

    % Mutation
    mutation_idx = rand(size(offspring)) < mutation_rate;
    offspring(mutation_idx) = offspring(mutation_idx) + randn * 0.05;

    % Update population with offspring
    pop_matrix = [pop_matrix; offspring];
    pop_matrix = max(min(pop_matrix, C), 0);
end

% Display best solution
best_solution = population(1,:);
best_fitness = fitness_function(best_solution);
disp('Optimal flow distribution:'), disp(best_solution);
disp(['Minimum travel time:', num2str(best_fitness)]);










