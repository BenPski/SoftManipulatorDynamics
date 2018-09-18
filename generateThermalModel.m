function [thermalmodel,msh,thermalResults] = generateThermalModel()
    %generates the thermal model for the tcas
    
    thermalmodel = createpde('thermal','transient');

    % Create a geometry that consists of a square with an embedded diamond.

    C1 = [1;0;0;0.008/2];
    C2 = [1;0;0.002;0.0014/2];
    C3 = [1;0;0.002;0.0002];
    C4 = [1;-0.001732;-0.001;0.0007];
    C5 = [1;-0.001732;-0.001;0.0002];
    C6 = [1;0.001732;-0.001;0.0007];
    C7 = [1;0.001732;-0.001;0.0002];

    gd = [C1,C2,C3, C4, C5, C6, C7];
    sf = 'C1-C3-C5-C7+C2-C3+C4-C5+C6-C7';
    ns = char('C1','C2','C3','C4','C5','C6','C7');
    ns = ns';
    dl = decsg(gd,sf,ns);

    geometryFromEdges(thermalmodel,dl);
    %IEEEfigure();
    %figure(1)
    %pdegplot(thermalmodel,'EdgeLabels','on','FaceLabels','on');
    %xlim([-0.0002,0.006])
    %ylim([-0.005,0.005])
    %axis equal

    %% Set thermal
    thermalProperties(thermalmodel,'Face',1,'ThermalConductivity',0.19, ...
                                   'MassDensity',1040, ...
                                   'SpecificHeat',1050);

    thermalProperties(thermalmodel,'Face',2,'ThermalConductivity',0.46, ...
                                   'MassDensity',1300, ...
                                   'SpecificHeat',1267.7); 
    thermalProperties(thermalmodel,'Face',3,'ThermalConductivity',0.46, ...
                                   'MassDensity',1300, ...
                                   'SpecificHeat',1267.7);
    thermalProperties(thermalmodel,'Face',4,'ThermalConductivity',0.46, ...
                                   'MassDensity',1300, ...
                                   'SpecificHeat',1267.7);

    %Specify Convection on the Boundary

    thermalBC(thermalmodel,'Edge',1:4, ...
                           'ConvectionCoefficient',50, ...
                           'AmbientTemperature',25);

    % Mesh the geometry.
    hmax = .0008; % element size
    msh = generateMesh(thermalmodel,'Hmax',hmax);
    % figure(2);
    % pdeplot(thermalmodel,'NodeLabels','on');
    % hold on;
    % axis equal
    
    thermalIC(thermalmodel,25);
    thermalResults = 25;
    
end