function P = DrawTraj3D(rMatrixECI, color)
% Description: create a 3D Plot of the propagated orbit.

if nargin < 2
    color = '#ff7403';
end

X = rMatrixECI(:, 1);
Y = rMatrixECI(:, 2);
Z = rMatrixECI(:, 3);

[x,y,z]=sphere;

I = imread('earth.jpg');

surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct')

hold on
P = plot3(X,Y,Z,'Color',color, 'Linestyle', '-', 'LineWidth', 1.5);
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([120, 30])


end