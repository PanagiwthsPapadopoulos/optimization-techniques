% Define the function whose root we want to find
f = {@(x) (x-2)^2+x*log(x+3),  
     @(x) exp(-2*x)+(x-2)^2,
     @(x) exp(x)*(x^3-1)+(x-1)*sin(x)
     };

epsilon_list = 0.004:0.00005:0.0049;
l = 0.01;
lambda_list = 0.0021:0.001:0.05;
epsilon = 0.001;
lambda_list_static = [0.0021, 0.005, 0.01];

function [k, a, b, fcounter] = dixotomos(l, epsilon, fun)
    a(1) = -1;            
    b(1) = 3;            
    k = 1;    
    fcounter = 0;
    while b(k) - a(k) >= l
    x1(k) = (a(k)+b(k))/2-epsilon;
    x2(k) = (a(k)+b(k))/2+epsilon;
    fx1 = fun(x1(k));
    fx2 = fun(x2(k));
    fcounter = fcounter + 2;   
    if fx1 < fx2
        a(k+1)=a(k);
        b(k+1) = x2(k);
    else
       a(k+1)=x1(k);
       b(k+1) = b(k);
    end
    k = k+1;
    end
end


for graph_index = 1:3
    figure;
    hold on
    for lambda_index = 1:length(lambda_list_static)
        [k, a, b, fcounter] = dixotomos(lambda_list_static(lambda_index), epsilon, f{graph_index});
        plot(1:1:k, a);
        plot(1:1:k, b);
    end
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('Iterations');
    ylabel('a, b');
    legend('a, l=0.0021','b, l=0.0021','a, l=0.005','b, l=0.005','a, l=0.01','b, l=0.01');
end


% Display results for constant l

 for graph_index = 1:3
     for epsilon_index = 1:length(epsilon_list)
        fprintf('Iteration %d\n', epsilon_index);
        [k, a, b, fcounter] = dixotomos(l, epsilon_list(epsilon_index), f{graph_index});
        fcounter_list(graph_index, epsilon_index) = fcounter;
        fprintf('fcounter %d\n', fcounter_list(graph_index, epsilon_index));
        ks_list(graph_index, epsilon_index) = k;
     end
     figure;
     plot(epsilon_list, fcounter_list(graph_index, :));
     title(extractAfter(func2str(f{graph_index}), '@(x)'));
     xlabel('epsilon values');
     ylabel('Calculations of f(x)');
 end


% Display results for constant e

 for graph_index = 1:3
    for lambda_index = 1:length(lambda_list)
        [k, a, b, fcounter] = dixotomos(lambda_list(lambda_index), epsilon, f{graph_index});
        fcounter_list(graph_index, lambda_index) = fcounter;
        ks_list(graph_index, lambda_index) = k;
    end
    figure;
    plot(lambda_list, fcounter_list(1, :));
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('lambda values');
    ylabel('Calculations of f(x)');
 end

