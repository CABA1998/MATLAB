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

%% Poisson equation with Dirichlet BCs

applyBoundaryCondition(model, "dirichlet", "Edge", 1:model.Geometry.NumEdges-1, "u", 1);
applyBoundaryCondition(model, "dirichlet", "Edge", model.Geometry.NumEdges, "u", -1);
applyBoundaryCondition(model, "dirichlet", "Edge", [10, 11], "u", -1);
applyBoundaryCondition(model, "dirichlet", "Edge", [5, 6], "u", 0);
applyBoundaryCondition(model, "dirichlet", "Edge", [1, 9, 14, 16], "u", 0.5);


origin = [0.3, 0.4];
radius = 0.2;
forcingValue = 10;
f = @(location, state) ...
    forcingValue * ((location.x - origin(1)).^2 + (location.y - origin(2)).^2 <= radius^2);
specifyCoefficients(model,"m",0,"d",0,"c",1,"a",0,"f",f); %see help


%%

figure(4);
generateMesh(model,Hedge={[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18],0.01}, ...
                           Hvertex={[1 4 6 8 10 15 16 17 18],0.001});
pdemesh(model)
axis equal

%%

results = solvepde(model);
u = results.NodalSolution;
figure;
pdeplot(model, "XYData", u);
colormap([linspace(1, 0, 256)' linspace(1, 0, 256)' zeros(256, 1)]);
colorbar;
hold on;
pdemesh(model, "EdgeColor", [0.5, 0.5, 0.5]);
title("Numerical Solution");
xlabel("x");
ylabel("y");
hold off;
