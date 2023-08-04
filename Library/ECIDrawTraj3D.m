function ECIDrawTraj3D(rMatrixECI)
% Description:
% Create a 3D Plot of the propagated orbit in the ECI reference frame.

global thetaG

X = rMatrixECI(:, 1);
Y = rMatrixECI(:, 2);
Z = rMatrixECI(:, 3);

[x,y,z]=sphere;

I = imread('earth.jpg');

earth = surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);
earth_axis = [0 0 1];
rotate(earth, earth_axis, rad2deg(thetaG))

hold on
plot3(X,Y,Z,'Color','#ff7403', 'Linestyle', '-', 'linewidth', 1.8)
plot3(X(1), Y(1), Z(1), 'Color', '#ff7403', 'LineStyle','none', 'marker', '.', 'markersize', 15)
plot3(X(end), Y(end), Z(end), 'Color', '#ff2e2e', 'LineStyle','none', 'marker', '.', 'markersize', 15)
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([-30, 20])


end