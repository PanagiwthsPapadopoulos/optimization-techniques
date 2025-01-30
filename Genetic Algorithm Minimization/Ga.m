
function pop_matrix = initialize_population(pop_size, E, V, C, source_edges)
    
    num_edges = size(E, 1);
    pop_matrix = zeros(pop_size, num_edges);  % Initialize matrix with zeros
    num_source_edges = length(source_edges);

    for i = 1:pop_size
  
        correct_source_edges = false;
        
        while ~correct_source_edges

            % Assign values for all but one source edge
            for j = 1:num_source_edges-1
    
                % Calculate flow for assigned edges
                intake_flow_occupied = sum(pop_matrix(i, source_edges));
    
                % Calculate available flow for the rest of the edges
                available_flow = V - intake_flow_occupied;
             
                % Assign cars to the edge with the maximum being the available
                % flow or maximum allowed
                upper_limit = min(available_flow, C(j));
                pop_matrix(i, j) = rand * upper_limit;
            end
    
            % Assign value to last edge
            pop_matrix(i, length(source_edges)) = V - sum(pop_matrix(i, source_edges) );

            % Check if last edge is within acceptable range, else start
            % again from scratch
            if pop_matrix(i, length(source_edges)) <= C(num_source_edges)
                correct_source_edges = true;
                break;
            else
                pop_matrix(i,:) = zeros(1,num_edges);
            end
            
        end

        correct_edges = false;

        while ~correct_edges
        
            % Distribute the flow throughout the network
            for j = num_source_edges+1:num_edges
    
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
                    
                    % Check if last edge is within acceptable range, else start
                    % again from scratch
                    if pop_matrix(i, j) > C(j) 
                        pop_matrix(i,num_source_edges+1:num_edges) = zeros(1,num_edges-num_source_edges);
                        break;
                    elseif j == num_edges
                        correct_edges = true;  
                        break;
                    end
                    
                else
                    upper_limit = min(flow_of_parents-flow_of_siblings, C(j));
                    pop_matrix(i, j) = rand * upper_limit;
                end

                
            end
        end
    end
end

function fitness = evaluate_fitness(pop_matrix, t, a, C)
    % Evaluate fitness for each solution in the population
    fitness = sum(t + (a .* pop_matrix) ./ (1 - (pop_matrix ./ C)), 2);
end

function offspring = single_point_crossover(parent1, parent2)
    crossover_point = randi([1, length(parent1)-1]);  % Random crossover point
    offspring = [parent1(1:crossover_point), parent2(crossover_point+1:end)];
end

function pop_matrix = apply_mutation(pop_matrix, mutation_rate, C, source_edges, V)
    num_individuals = size(pop_matrix, 1);  % Number of solutions in population
    num_edges = size(pop_matrix, 2);  % Number of edges

    num_mutations = round(num_individuals * mutation_rate);  % Number of individuals to mutate

    for i = 1:num_mutations
        individual_idx = randi(num_individuals);  % Randomly select an individual to mutate
        
        % Select a random edge, ensuring it's NOT a source edge  
        edge_index = randi([1,11]);

        % Mutation value as a fraction of the edge capacity
        mutation_value = (rand - 0.5) * 0.2 * C(edge_index);  % +/- 10% of capacity

        % Apply mutation while respecting capacity limits
        pop_matrix(individual_idx, edge_index) = max(0, ...
            min(pop_matrix(individual_idx, edge_index) + mutation_value, C(edge_index)));
    end
end

function pop_matrix = fix_population(pop_matrix, source_edges, V, C, E)
    num_source_edges = length(source_edges);

    for i = 3:size(pop_matrix, 1)
        % Fix source edge flows
        total_flow = sum(pop_matrix(i, source_edges));
        deficit = V - total_flow;
        correction = round(deficit/num_source_edges);

        if deficit > 0      
            for j=source_edges
                [~, min_index] = min(pop_matrix(i, source_edges));
                if j ~= source_edges(num_source_edges)
                    pop_matrix(i,min_index) = pop_matrix(i,min_index) + correction;
                else
                    total_flow = sum(pop_matrix(i, source_edges));
                    pop_matrix(i,min_index) = pop_matrix(i,min_index) + V-total_flow;
                end
            end
        end

        if deficit < 0
            for j=source_edges
                [~, max_index] = max(pop_matrix(i, source_edges));
                if j ~= source_edges(num_source_edges)
                    pop_matrix(i,max_index) = pop_matrix(i,max_index) + correction;
                else
                    total_flow = sum(pop_matrix(i, source_edges));
                    pop_matrix(i,max_index) = pop_matrix(i,max_index) + V-total_flow;
                end
            end         
        end
        for j=num_source_edges+1:length(C)
            % disp(['index j: ', num2str(j(1,1))])
            parent_edges = find(E(:, j) == 1);
            
            sibling_edges = find(any(E(parent_edges, :) == 1, 1));
            % disp(['Sibling edges for edge ', num2str(j(1,1)), ': ', num2str(sibling_edges)])          
            flow_of_siblings = sum(pop_matrix(i, sibling_edges));

            flow_of_parents = sum(pop_matrix(i, parent_edges));
            
            if abs(flow_of_parents - flow_of_siblings) > 1e-6
                deficit = flow_of_parents - flow_of_siblings;
                num_siblings_edges = length(sibling_edges);

                if deficit > 0
                    for k=sibling_edges
                        [~, min_relative_index] = min(pop_matrix(i, sibling_edges));
                        min_index = sibling_edges(min_relative_index);
                        % disp(['min_index value: ', num2str(min_index)])
                        pop_matrix(i,min_index) = pop_matrix(i,min_index) + deficit/num_siblings_edges;
                    end
                end
        
                if deficit < 0
                    for k=sibling_edges
                        [~, max_relative_index] = max(pop_matrix(i, sibling_edges));
                        max_index = sibling_edges(max_relative_index);
                        % disp(['max_index value: ', num2str(max_index)])
                        pop_matrix(i,max_index) = pop_matrix(i,max_index) + deficit/num_siblings_edges;
                    end
                end
            end
        end      
        % Ensure no edge exceeds capacity
        pop_matrix(i, :) = min(pop_matrix(i, :), C);
    end
end

function result = check_res(pop_matrix, C, V, source_edges, sink_edges, E)
    result = true;
    num_source_edges = length(source_edges);
    num_sink_edges = length(sink_edges);

    for i = 1:size(pop_matrix, 1)
        % disp(['Checking line ' , num2str(i)])
        source_sum = 0;
        for j_idx = 1:length(source_edges)
            j = source_edges(j_idx);
            disp(['Adding source edge with index: ' , num2str(j(1,1)), ' and value: ', num2str(pop_matrix(i,j))])
            source_sum = source_sum + pop_matrix(i,j);
            disp(['source sum after adding element: ', num2str(source_sum)])
        end
        
        if abs(source_sum - V) > 1e-1
            % disp('Source sum is incorrect')
            % disp([num2str(i), ' ', num2str(source_sum)])
            result = false;
            return;
        end
        
        sink_sum = 0;
        % disp('Adding sink edges')
        for j_idx = 1:length(sink_edges)
            j = sink_edges(j_idx);
            % disp(['Adding sink edge with index: ' , num2str(j(1,1)), ' and value: ', num2str(pop_matrix(i,j))])
            sink_sum = sink_sum + pop_matrix(i,j);
            % disp(['sink sum after adding element: ', num2str(sink_sum)])
        end

        if abs(sink_sum - V) > 1e-6
            % disp('Sink sum is incorrect')
            % disp([num2str(i),' ', num2str(sink_sum)])
            result = false;
            return;
        end

        for j=1:length(pop_matrix(i,1))
            if pop_matrix(i,j) < 0 || pop_matrix(i,j) > C(j)
                result = false;
                return;
            end
        end

        num_edges = length(C);
        
        for j = num_source_edges+1:num_edges
    
            parent_edges = find(E(:, j) == 1);
            
            sibling_edges = find(any(E(parent_edges, :) == 1, 1));
            
            flow_of_siblings = sum(pop_matrix(i, sibling_edges));

            flow_of_parents = sum(pop_matrix(i, parent_edges));
            
            if abs(flow_of_siblings - flow_of_parents) > 1e-6
                % disp('Flow sum is incorrect')
                % disp(sibling_edges)
                result = false;
                return;
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
V = 115;
t = 1;

generations = 500;  % Number of generations
pop_size = 100;
mutation_rate = 0.1;  % Probability of mutation

% Fitness function: minimize total travel time
fitness_function = @(x) sum(t + (a .* x) ./ (1 - (x ./ C)), 2);

pop_matrix = initialize_population(pop_size, E, V, C, source);


best_fitness = inf;
best_solution = [];

best_fitness_per_gen = zeros(1, generations / 10);


for generation = 1:generations
    % Evaluate fitness
    fitness = evaluate_fitness(pop_matrix, t, a, C);

    % Selection of the best two solutions
    [~, idx] = sort(fitness);
    parent1 = pop_matrix(idx(1), :);
    parent2 = pop_matrix(idx(2), :);
    
    % Crossover to create new offspring
    new_population = zeros(pop_size, size(pop_matrix, 2));
    new_population(1, :) = parent1;  % Keep the best solution
    new_population(2, :) = parent2;  % Keep the second best solution
    
    for i = 3:pop_size
        new_population(i, :) = single_point_crossover(parent1, parent2);
    end
       
    % Mutation
    new_population = apply_mutation(new_population, mutation_rate, C, find(sum(E, 1)), V);

    % Fix population
    new_population = fix_population(new_population, source, V, C, E);
    
    % Replace old population
    pop_matrix = new_population;

    % Store and print best fitness every 10 generations
    if mod(generation, 10) == 0
        best_fitness = min(fitness);
        best_fitness_per_gen(generation / 10) = best_fitness;
        disp(['Generation ', num2str(generation), ': Best Fitness = ', num2str(best_fitness)]);
    end

end

% Plot best fitness per tenth generation
figure;
plot(10:10:generations, best_fitness_per_gen, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Generation');
ylabel('Best Fitness');
title('Best Fitness Evolution Over Generations');
grid on;

