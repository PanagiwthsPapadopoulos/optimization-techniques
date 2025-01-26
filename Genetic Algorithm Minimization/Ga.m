
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
    mutation_value = rand*5;  % Small random change

    edge_index = randi([1,4]);
    disp(edge_index)
    disp(max(0, min(pop_matrix(3, edge_index) + mutation_value, C(edge_index))))
    pop_matrix(3, edge_index) = max(0, min(pop_matrix(3, edge_index) + mutation_value, C(edge_index)));
    % pop_matrix(3, 2) = pop_matrix(3, 2) + mutation_value;
end

function pop_matrix = fix_population(pop_matrix, source_edges, V, C, E)
    disp(pop_matrix(3,1:13))
    num_source_edges = length(source_edges);

    for i = 3:size(pop_matrix, 1)

        % Fix source edge flows
        total_flow = sum(pop_matrix(i, source_edges));
        deficit = V - total_flow;
        if deficit > 0
              
            for j=1:num_source_edges
                [~, min_index] = min(pop_matrix(i, source_edges));
                pop_matrix(i,min_index) = pop_matrix(i,min_index) + (deficit/num_source_edges);
            end
        end

        if deficit < 0
            
            for j=1:num_source_edges
                [~, max_index] = max(pop_matrix(i, source_edges));
                pop_matrix(i,max_index) = pop_matrix(i,max_index) + (deficit/num_source_edges);
            end
        end

        for j=num_source_edges+1:length(C)

            parent_edges = find(E(:, j) == 1);
            
            sibling_edges = find(any(E(parent_edges, :) == 1, 1));
            
            flow_of_siblings = sum(pop_matrix(i, sibling_edges));

            flow_of_parents = sum(pop_matrix(i, parent_edges));
            % 
            % disp('parentss')
            % disp(parent_edges)
            % disp('siblings')
            % disp(sibling_edges)
            % 
            % 
            % disp('flow of siblings')
            % disp(flow_of_siblings)
            % disp('flow of parents')
            % disp(flow_of_parents)
            
            if ~ (flow_of_parents == flow_of_siblings)
                % disp('Detected discrepancy')
                deficit = flow_of_parents - flow_of_siblings;
                num_siblings_edges = length(sibling_edges);

                % disp('deficit')
                % disp(deficit)

                if deficit > 0
                    for k=1:num_siblings_edges
                        [~, min_index] = min(pop_matrix(i, sibling_edges));
                        pop_matrix(i,min_index) = pop_matrix(i,min_index) + deficit/num_siblings_edges;
                    end
                end
        
                if deficit < 0
                    for k=1:num_siblings_edges
                        [~, max_index] = max(pop_matrix(i, sibling_edges));
                                        disp(max_index)
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
        source_sum = 0;
        for j=1:num_source_edges
            source_sum = source_sum + pop_matrix(i,j);
        end
        
        if ~ (source_sum == V)
            disp('Source sum is incorrect')
            disp([num2str(i),num2str(source_sum)])
            result = false;
            return;
        end
        
        sink_sum = 0;
        for j=length(pop_matrix(i,:)) - num_sink_edges + 1:length(pop_matrix(i,:))
            sink_sum = sink_sum + pop_matrix(i,j);
        end

        if ~ (sink_sum == V)
            disp('Sink sum is incorrect')
            disp([num2str(i),num2str(sink_sum)])
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
        disp(pop_matrix);
        for j = num_source_edges+1:num_edges
    
            parent_edges = find(E(:, j) == 1);
            
            sibling_edges = find(any(E(parent_edges, :) == 1, 1));
            
            flow_of_siblings = sum(pop_matrix(i, sibling_edges));

            flow_of_parents = sum(pop_matrix(i, parent_edges));
            
            if ~ (flow_of_siblings == flow_of_parents)
                disp('Flow sum is incorrect')
                disp(sibling_edges)
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
V = 100;
t = 1;

generations = 200;  % Number of generations
pop_size = 3;
mutation_rate = 0.1;  % Probability of mutation

% Fitness function: minimize total travel time
fitness_function = @(x) sum(t + (a .* x) ./ (1 - (x ./ C)), 2);

pop_matrix = initialize_population(pop_size, E, V, C, source);
% disp('Sample population matrix (first 5 solutions):');
% disp(pop_matrix(1, :));




best_fitness = inf;
best_solution = [];

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
        disp(new_population(:, 1:13));
    % Mutation
    new_population = apply_mutation(new_population, mutation_rate, C, find(sum(E, 1)), V);
    
    
    disp(new_population(:, 1:13));
    % Fix population
    new_population = fix_population(new_population, source, V, C, E);

    disp(new_population);
    disp(check_res(pop_matrix, C, V, source, sink, E))
    break;
    % Replace old population
    pop_matrix = new_population;

    % Store best solution

end


return;


% Display best solution
best_solution = population(1,:);
best_fitness = fitness_function(best_solution);
disp('Optimal flow distribution:'), disp(best_solution);
disp(['Minimum travel time:', num2str(best_fitness)]);










