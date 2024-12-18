f = {
    @(x) (x-2)^2+x*log(x+3),  
    @(x) exp(-2*x) + (x-2)^2,
    @(x) exp(x)*(x^3-1)+(x-1)*sin(x)
 };

df = {@(x) 2*(x-2) + log(x+3) + x/(x+3),
      @(x) 2*(x-2) - 2*exp(-2*x),
      @(x) sin(x) + cos(x)*(x-1) + exp(x)*(x^3+3*x^2-1)
};

lambda_list = 0.000001:0.0001:0.01;
lambda_list_static = [0.000001, 0.0001, 0.01];

function [k, a, b, fcounter] = dixotomos_paragwgos(l, dfun)
    a(1) = -1;
    b(1) = 3;
    n = 1;
    k = 1;
    fcounter = 0;

    while 0.5^n >= l / (b(1) - a(1))  
        n = n + 1;
    end

    while k < n
    x(k) = (a(k)+b(k))/2;
    df = dfun(x(k));
    fcounter = fcounter + 1;

    if df == 0
        fprintf('df == 0\n');
        fprintf('xk = %d', x(k));
        return;
    else
        if df > 0
            a(k+1) = a(k);
            b(k+1) = x(k);
        else
            a(k+1) = x(k);
            b(k+1) = b(k);
        end
    end
    k = k+1;
     fprintf('xk = %d\n', x(k-1));
    end
end

for graph_index = 1:3
    for lambda_index = 1:length(lambda_list)
        [k, a, b, fcounter] = dixotomos_paragwgos(lambda_list(lambda_index), df{graph_index});
        fcounter_list(graph_index, lambda_index) = fcounter;
        ks_list(graph_index, lambda_index) = k;
    end
    figure;
    plot(lambda_list, fcounter_list(graph_index, :));
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('lambda values');
    ylabel("Calculations of f'(x)");
 end


for graph_index = 1:3
    figure;
    hold on
    for lambda_index = 1:length(lambda_list_static)
        [k, a, b, fcounter] = dixotomos_paragwgos(lambda_list_static(lambda_index), df{graph_index});
        plot(1:1:k, a);
        plot(1:1:k, b);        
    end
    title(extractAfter(func2str(f{graph_index}), '@(x)'));
    xlabel('Iterations');
    ylabel('a, b');
    legend('a, l=0.000001','b, l=0.000001','a, l=0.0001','b, l=0.0001','a, l=0.01','b, l=0.01');
end
