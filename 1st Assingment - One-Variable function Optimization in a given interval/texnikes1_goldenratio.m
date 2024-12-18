f = {@(x) (x-2)^2+x*log(x+3),  
 @(x) exp(-2*x) + (x-2)^2,
 @(x) exp(x)*(x^3-1)+(x-1)*sin(x)};

l = 0.01;
lambda_list = 0.000001:0.0001:0.01;
lambda_list_static = [0.000001, 0.0001, 0.01];

function [k, a, b, fcounter] = golden(l, fun)
    golden_ratio = 0.618;
    a(1) = -1;            
    b(1) = 3;            
        
    x1(1) = a(1) + (1-golden_ratio)*(b(1)-a(1));
    x2(1) = a(1) + golden_ratio*(b(1)-a(1));
    fx1 = fun(x1(1));
    fx2 = fun(x2(1));
    k = 1;
    fcounter = 2;

    while b(k) - a(k) >= l
        if fx1 >= fx2
            a(k+1)=x1(k);
            b(k+1) = b(k);
            x2(k+1) = a(k+1)+golden_ratio*(b(k+1)-a(k+1));
            x1(k+1) = x2(k);
            fx1 = fx2;
            fx2 = fun(x2(k+1));
            fcounter = fcounter + 1;
            k = k+1;
        else
            a(k+1) = a(k);
            b(k+1) = x2(k); 
            x2(k+1)=x1(k);
            x1(k+1) = a(k+1) + (1-golden_ratio)*(b(k+1)-a(k+1));
            fx2 = fx1;
            fx1 = fun(x1(k+1));
            fcounter = fcounter + 1;
            k = k+1;
        end
    end
end


for graph_index = 1:3
    for lambda_index = 1:length(lambda_list)
        [k, a, b, fcounter] = golden(lambda_list(lambda_index), f{graph_index});
        fcounter_list(graph_index, lambda_index) = fcounter;
        ks_list(graph_index, lambda_index) = k;
    end
    figure;
    plot(lambda_list, fcounter_list(graph_index, :));
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('lambda values');
    ylabel('Calculations of f(x)');
end

for graph_index = 1:3
    figure;
    hold on
    for lambda_index = 1:length(lambda_list_static)
        [k, a, b, fcounter] = golden(lambda_list_static(lambda_index), f{graph_index});
        plot(1:1:k, a);
        plot(1:1:k, b);
    end
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('Iterations');
    ylabel('a, b');
    legend('a, l=0.000001','b, l=0.000001','a, l=0.0001','b, l=0.0001','a, l=0.01','b, l=0.01');
end
