f = @(x, y) x.^5 .* exp(-x.^2-y.^2);

x = linspace(-2, 2, 100); 
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);
Z = f(X, Y);
figure;
mesh(X, Y, Z);
title('f(x,y) = x^5 * exp(-x^2-y^2)');
xlabel('x'); ylabel('y'); zlabel('f(x, y)');

