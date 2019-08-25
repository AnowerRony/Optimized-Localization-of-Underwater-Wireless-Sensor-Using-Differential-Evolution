
function PlotR()

[X,Y] = meshgrid(-5.12:0.03:5.12,-5.12:0.03:5.12);
N = size(X,1);

for i = 1:N
    for j = 1:N
        z = [X(i,j),Y(i,j)];
        Z(i,j) = Rastrigin(z);
    end
end

figure(1);
mesh(X,Y,Z);
title('Objective Function in 2 dimension');
axis([-5.5 5.5 -5.5 5.5 0 80]);

end