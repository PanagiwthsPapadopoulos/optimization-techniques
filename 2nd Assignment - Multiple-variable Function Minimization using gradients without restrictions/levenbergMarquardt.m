f = @(x, y) x.^5 .* exp(-x.^2-y.^2);
grad_f = @(x, y) [5*x.^4.*exp(-x.^2-y.^2) - 2*x.^6.*exp(-x.^2-y.^2); ...
                  -2*y.*x.^5.*exp(-x.^2-y.^2)];
hessian_f = @(x, y) [ 20*x^3*exp(- x^2 - y^2) - 22*x^5*exp(- x^2 - y^2) + ...
                          4*x^7*exp(- x^2 - y^2), 4*x^6*y*exp(- x^2 - y^2) - 10*x^4*y*exp(- x^2 - y^2) ;
                        4*x^6*y*exp(- x^2 - y^2) - 10*x^4*y*exp(- x^2 - y^2),  4*x^5*y^2*exp(- x^2 - y^2) - 2*x^5*exp(- x^2 - y^2) ];


function [alpha, beta, k, cnt] = bisection_method_derivatives(A, B, l, f)
    % alpha is ana array that stores a1,a2,...,ai
    % beta the same for b1,b2,...,bi
    % k is the total number of iterations required to finish the algo
    % cnt is the amount of times fi(x) was called to be calculated
    
    % Calculate the minimum n that satisfies the tolerance inequation
    n = ceil(log((B - A) / l) / log(2));
    % Derivative of f(x)
    df = diff(f);

    alpha = [];
    beta = [];
    k = 1;
    cnt = 0;

    alpha(end+1) = A;
    beta(end+1) = B;
    for k=1:(n) % k = k + 1 is automatic with `for` loop
        xk = (alpha(end)+beta(end))/2;
        dfxk = df(xk);
        cnt = cnt + 1;
        if dfxk == 0
            break;
        elseif dfxk > 0
            % bhma 2
            alpha(end+1) = alpha(end);
            beta(end+1) = xk;
        else % dfxk < 0
            % bhma 3
            alpha(end+1) = xk;
            beta(end+1) = beta(end);
        end
    end
end




function [x_vals, y_vals, f_values] = lm_method_fixed(f, grad_f, hessian_f, x0, max_iter, epsilon, gamma, reg)


    
    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    
    f_values = [];

    for i = 1:max_iter
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        disp(eig(hess))
        f_values = [f_values, f(x(1), x(2))];
        

        eigenvalues = eig(hess);
        
        m = max(abs(eigenvalues)) + reg;
        hess_reg = hess + m*eye(size(hess));

        disp(x);  % Display the current point
        
        
        if norm(grad_f(x(1), x(2))) < epsilon      
            disp("break because norm(grad) < epsilon")
            return;
        end
        
        d = - hess_reg^(-1) * grad;
        x = x + gamma * d.';  

        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];

        
    end

end

function [x_vals, y_vals, f_values] = lm_method_optimal(f, grad_f, hessian_f, x0, max_iter, epsilon, reg)
    syms g;

    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    f_values = [];


    for i = 1:max_iter
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        disp(eig(hess))

        f_values = [f_values, f(x(1), x(2))];

        if norm(grad_f(x(1), x(2))) < epsilon      
            disp("break because norm(grad) < epsilon")
            return;
        end

        eigenvalues = eig(hess);
        m = max(abs(eigenvalues)) + reg;
        hess_reg = hess + m*eye(size(hess));

        
        phi(g) = exp(-(x(1)-g*grad(1))^2 -(x(2)-g*grad(2))^2)*(x(1)-g*grad(1))^5;
        [a,b,cnt,~] = bisection_method_derivatives(-10,10, 0.001, phi);

        gamma = (a(end)+b(end))/2;

        
        
        d = - hess_reg^(-1) * grad;
        x = x + gamma * d.';  

        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];

        
    end
    
    % Final point
    disp('Final point:');
    disp(x);
    disp(['x = ', num2str(x(1))]);
    disp(['y = ', num2str(x(2))]);

end

function [x_vals, y_vals, f_values] = lm_method_armijo(f, grad_f, hessian_f, x0, max_iter, epsilon, alpha, beta, gamma, reg)

    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    f_values = [];


    for i = 1:max_iter
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        disp(eig(hess))

        f_values = [f_values, f(x(1), x(2))];

        if norm(grad_f(x(1), x(2))) < epsilon      
            disp("break because norm(grad) < epsilon")
            return;
        end

        eigenvalues = eig(hess);
        m = max(abs(eigenvalues)) + reg;
        hess_reg = hess + m*eye(size(hess));

       

        d = - hess_reg^(-1) * grad;
        % Line search using Armijo rule
        gamma = 1;
        while f(x(1) + gamma * d(1), x(2) + gamma * d(2)) > f(x(1), x(2)) + alpha * gamma * norm(grad)^2
            gamma = beta * gamma;  
        end
 
        
        
        x = x + gamma * d.';  

        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];

        
    end
    
    % Final point
    disp('Final point:');
    disp(x);
    disp(['x = ', num2str(x(1))]);
    disp(['y = ', num2str(x(2))]);

end



% Parameters
x0 = [1, -1]; % Initial points
max_iter = 100; % Maximum number of iterations
tol = 1e-4;    % Tolerance for gradient norm
alpha = 0.2;  % Armijo condition parameter
beta = 0.01;  % Armijo step size reduction factor
gamma = 0.5;
armijo_initial = 1;
reg = 10;      % Regularization for non positive defined matrices

[x_vals_fixed, y_vals_fixed, f_values_fixed] = lm_method_fixed(f, grad_f, hessian_f, x0, max_iter, tol, gamma, reg);

[x_vals_optimal, y_vals_optimal, f_values_optimal] = lm_method_optimal(f, grad_f, hessian_f, x0, max_iter, tol, reg);

[x_vals_armijo, y_vals_armijo, f_values_armijo] = lm_method_armijo(f, grad_f, hessian_f, x0, max_iter, tol, alpha, beta, armijo_initial, reg);

% Plot the results in a single 3D plot
close all
figure;

% Define a grid for plotting the function surface
[x_grid, y_grid] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
z_grid = f(x_grid, y_grid);  % Evaluate f(x, y) over the grid

% Plot the function surface
surf(x_grid, y_grid, z_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
colormap parula;
hold on;

% Specify colors
color_fixed = [1, 0, 0];  % Red (Fixed Step Size)
color_optimal = [0, 1, 0];  % Green (Optimal Step Size)
color_armijo = [0, 0, 1];  % Blue (Armijo Rule)

% Plot the function surface (without adding it to the legend)
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z = f(X, Y); % Compute function values
surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Plot the function surface
hold on;

% Plot each method's trajectory with the specified colors and DisplayName
h1 = plot3(x_vals_fixed, y_vals_fixed, f(x_vals_fixed, y_vals_fixed), '-o', 'LineWidth', 1.5, 'Color', color_fixed, 'DisplayName', 'Fixed Step Size'); % Fixed step size (Red)
h2 = plot3(x_vals_optimal, y_vals_optimal, f(x_vals_optimal, y_vals_optimal), '-o', 'LineWidth', 1.5, 'Color', color_optimal, 'DisplayName', 'Optimal Step Size'); % Optimal step size (Green)
h3 = plot3(x_vals_armijo, y_vals_armijo, f(x_vals_armijo, y_vals_armijo), '-o', 'LineWidth', 1.5, 'Color', color_armijo, 'DisplayName', 'Armijo Rule'); % Armijo rule (Blue)

% Add labels and title
title('Levenberg-Marquardt Trajectories for Different Methods');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
grid on;
view(3);  % Set 3D view
axis tight;  % Adjust axis limits for better visualization
hold off;

% Add legend to identify methods (exclude the function plot from the legend)
legend([h1, h2, h3], 'Location', 'northeast', 'FontSize', 12);


figure;
plot(1:length(f_values_fixed), f_values_fixed, '-o', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('f(x,y)');
title('f values vs Iteration Number for fixed step');
grid on;
xticks(1:length(f_values_fixed));

figure;
plot(1:length(f_values_optimal), f_values_optimal, '-o', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('f(x,y)');
title('f values vs Iteration Number for optimal step');
grid on;
xticks(1:length(f_values_optimal));


figure;
plot(1:length(f_values_armijo), f_values_armijo, '-o', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('f(x,y)');
title('f values vs Iteration Number for armijo step');
grid on;
xticks(1:length(f_values_armijo));