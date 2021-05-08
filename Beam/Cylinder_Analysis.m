clear all 
close all

model = createpde('structural','static-solid');
importGeometry(model,'3D_Slab_3.stl');

figure
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Bracket with Face Labels')

figure
pdegplot(model,'FaceLabels','on')
view(-134,-32)
title('Bracket with Face Labels, Rear View')

structuralProperties(model,'Cell',1,'YoungsModulus',200e9, ...
                                    'PoissonsRatio',0.3);
                                
structuralBC(model,'Face',1,'Constraint','fixed');
distributedLoad = 1e1; % Applied load in Pascals
structuralBoundaryLoad (model,'Face',6,'SurfaceTraction',[0;0;distributedLoad]);

bracketThickness = 5e-1;
generateMesh(model,'Hmax',bracketThickness);
figure
pdeplot3D(model)
title('Mesh with Quadratic Tetrahedral Elements');

result = solve(model);

minUz = min(result.Displacement.uz);
fprintf('Maximal deflection in the z-direction is %g meters.', minUz)

figure
pdeplot3D(model,'ColorMapData',result.Displacement.ux)
title('x-displacement')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',result.Displacement.uy)
title('y-displacement')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',result.Displacement.uz)
title('z-displacement')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',result.VonMisesStress)
title('von Mises stress')
colormap('jet')