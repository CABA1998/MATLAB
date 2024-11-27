clear;
close all;
clc;

% 1. import a geometry from pdeModeler or
load('MyModel.mat');

% 2. define fegeometry and plot it
model = createpde();
gm = decsg(gd,sf,ns);
geometryFromEdges(model, gm);

figure 
pdegplot(model,"EdgeLabels","on"); 
axis equal

%%

generateMesh(model,Hedge={[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18],0.01}, ...
                           Hvertex={[1 4 6 8 10 15 16 17 18],0.001});
pdemesh(model)
axis equal

%%

% Specify PDE coefficients
applyBoundaryCondition(model, "dirichlet", "Edge", 1:model.Geometry.NumEdges-1, "u", 0);

c = 1; % Wave speed
specifyCoefficients(model, 'm', 1, 'd', 0, 'c', c^2, 'a', 0, 'f', 0);

% Define initial conditions
origin = [0.3, 0.3];
radius = 0.1;
initVal = -1;
u0 = @(location, state) ...
    initVal * ((location.x - origin(1)).^2 + (location.y - origin(2)).^2 <= radius^2);
ut0 = @(location) 0;
setInitialConditions(model, u0, ut0);

% Time span
Tmax = 1;
tlist = linspace(0, Tmax, 50);

% Solve the PDE
result = solvepde(model, tlist);

% Extract solution at the final time step
u = result.NodalSolution(:, end);

% Plot solution at the final time step
figure;
pdeplot(model, 'XYData', u,"ZData", u, 'ColorMap', 'jet', 'Mesh', 'off');
xlabel('x'); ylabel('y'); zlabel('Displacement');
title(['Wave Equation Solution at t = ', num2str(tlist(end))]);
colorbar;
axis equal tight;

%% 
figure;
for i = 1:length(tlist)
    % Plot the solution for the current time step
    pdeplot(model, 'XYData', result.NodalSolution(:, i));
    colormap([linspace(1, 0, 256)' linspace(1, 0, 256)' zeros(256, 1)]);
    %%plot mesh as well    
    hold on;
    pdemesh(model, "EdgeColor", [0, 0, 0]);
    hold off;
    xlabel('x'); ylabel('y');
    title(['Time: t = ', num2str(tlist(i))]);
    colorbar;
    clim([min(u(:)), max(u(:))]); % Dynamic color range
    axis equal tight;
    pause(0.1); % Adjust speed of animation
end