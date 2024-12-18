f = @(x, y) (1/3)*x.^2 + 3*y.^2;
grad_f = @(x, y) [(2*x)/3; ...
                  6*y];


function [proj] = projection(a, b, x)
    if x<=a
        proj = a; % lower bound 
    elseif x>=b
        proj = b; % upper bound
    else
        proj = x; % no need for projection
    end
end

function [x_vals, y_vals, f_values] = gradient_descent_projection(f, grad_f, x0, max_iter, epsilon, gamma, s)

    x = x0;
    x_vals = x(1, 1);
    y_vals = x(1,2);
    f_values = [f(x(1), x(2))];
    


    for k = 1:max_iter

        % Υπολογισμός της κλίσης
        grad = grad_f(x(1), x(2));  
        
        if norm(grad) < epsilon      
            fprintf("Exit because of norm");
            break;
        end
        
        
        

        % Προβολή στο εφικτό σύνολο
        x_valid = [projection(-10, 5, x(1)-s*grad(1)), projection(-8, 12, x(2)-s*grad(2))];
        
       

        % Ενημέρωση σημείου
        x = x + gamma * (x_valid - x);

        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];
        f_values = [f_values, f(x(1), x(2))];

        
        fprintf('Point [%d %d]\n', x(1), x(2));
        fprintf('Funtion value %d\n', f(x(1), x(2)));
        fprintf('Grad norm %d\n', norm(grad));
        fprintf('\n');
        
    
    end

end



x0 = {[5,-5], [-5, 10], [8, -10]};
disp(x0);
disp(x0(1))
%%
tol = 0.01;
max_iter = 1500;
s = 15;
gamma_values = [0.5, 0.1, 0.2];
s_values = [5, 10/6, 0.1];


for i = 1:length(gamma_values)
    [x_vals, y_vals, f_values] = gradient_descent_projection(f, grad_f, x0{i}, max_iter, tol, gamma_values(i), s_values(i));
    results{i} = struct('x_vals', x_vals, 'y_vals', y_vals, 'f_values', f_values);
end

close all

for i = 1:length(gamma_values)

% Define the range for x and y
x = linspace(-10, 5, 100); % Range for x
y = linspace(-8, 12, 100); % Range for y

% Create a grid for x and y
[X, Y] = meshgrid(x, y);

% Evaluate the function on the grid
Z = f(X, Y);

% Plot the surface
figure;
surf(X, Y, Z);

% Customize the plot
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('3D Plot of f(x, y) = (1/3)x^2 + 3y^2');
shading interp; % Smooth shading
colormap parula;
grid on;
hold on;
%%

% Plot each method's trajectory with the specified colors and DisplayName
h1 = plot3(results{i}.x_vals, results{i}.y_vals, results{i}.f_values, '-o', 'LineWidth', 1.5, 'Color', [1, 0, 0]); 


% Add labels and title
title('Trajectories for gamma ' + string(gamma_values(i)) + ' and s ' + string(s_values(i)) );
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
grid on;
view(3);  % Set 3D view
axis tight;  % Adjust axis limits for better visualization
hold off;

% % Add legend to identify methods (exclude the function plot from the legend)
% legend([h1, h2, h3], 'Location', 'northeast', 'FontSize', 12);


figure;
plot(1:length(results{i}.f_values-1), results{i}.f_values, '-o', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('f(x,y)');
title('f values vs Iteration Number for fixed step');
grid on;

end