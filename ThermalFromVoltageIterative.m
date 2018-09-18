function [tcaTemps,centerTemps,thermalResults] = ThermalFromVoltageIterative(thermalmodel,msh,U,dt,thermalResults)
%since I can't pass around some of the data in pythong need to take care of globals
% global thermalResults
% global thermalmodel
% global msh

%the goal here is to give one set of voltages to run over a certain time
%step, so that there is a bit more control rather than running over long
%periods of time

%so need to provide the voltages, the time step, and the initial
%temperatures

%also want to provide the model so it doesn't have to keep being
%regenerated and the mesh

%intitialize temperatures
    %for the first step fine to just be =25
    % the rest of steps pass in the returned value
%run the simulation from 0->t (use only 2 points?)

%extract temperatures and determine some of the more useful ones


%f: time -> [V,indices]
%also faces are 2,3,4
%% Specify Internal Heat Generation on Entire Geometry
% Create a transient thermal model.

%plot(0:10,Ufun(0:10))

%%
% thermalmodel = createpde('thermal','transient');
% 
% % Create a geometry that consists of a square with an embedded diamond.
% 
% C1 = [1;0;0;0.008/2];
% C2 = [1;0;0.002;0.0014/2];
% C3 = [1;0;0.002;0.0002];
% C4 = [1;-0.001732;-0.001;0.0007];
% C5 = [1;-0.001732;-0.001;0.0002];
% C6 = [1;0.001732;-0.001;0.0007];
% C7 = [1;0.001732;-0.001;0.0002];
% 
% gd = [C1,C2,C3, C4, C5, C6, C7];
% sf = 'C1-C3-C5-C7+C2-C3+C4-C5+C6-C7';
% ns = char('C1','C2','C3','C4','C5','C6','C7');
% ns = ns';
% dl = decsg(gd,sf,ns);
% 
% geometryFromEdges(thermalmodel,dl);
% %IEEEfigure();
% %figure(1)
% %pdegplot(thermalmodel,'EdgeLabels','on','FaceLabels','on');
% %xlim([-0.0002,0.006])
% %ylim([-0.005,0.005])
% %axis equal
% 
% %% Set thermal
% thermalProperties(thermalmodel,'Face',1,'ThermalConductivity',0.19, ...
%                                'MassDensity',1040, ...
%                                'SpecificHeat',1050);
%                            
% thermalProperties(thermalmodel,'Face',2,'ThermalConductivity',0.46, ...
%                                'MassDensity',1300, ...
%                                'SpecificHeat',1267.7); 
% thermalProperties(thermalmodel,'Face',3,'ThermalConductivity',0.46, ...
%                                'MassDensity',1300, ...
%                                'SpecificHeat',1267.7);
% thermalProperties(thermalmodel,'Face',4,'ThermalConductivity',0.46, ...
%                                'MassDensity',1300, ...
%                                'SpecificHeat',1267.7);
% 
% %Specify Convection on the Boundary
% 
% thermalBC(thermalmodel,'Edge',1:4, ...
%                        'ConvectionCoefficient',50, ...
%                        'AmbientTemperature',25);
% 
% % Mesh the geometry.
% hmax = .0008; % element size
% msh = generateMesh(thermalmodel,'Hmax',hmax);
% % figure(2);
% % pdeplot(thermalmodel,'NodeLabels','on');
% % hold on;
% % axis equal
%%
% Specify that the entire geometry generates heat
tlist=[0,dt];
%clear thermalresults1;
%set ICs
thermalIC(thermalmodel,thermalResults);
%set heat generation
R = 3.5;% resistance of a TCA
Volume = 1e-7; % m^3
q_dot = U.^2/R/Volume;
for j=[2,3,4]
    internalHeatSource(thermalmodel,q_dot(j-1),'Face',j);
end
%do the simulation for the timestep
thermalResults = solve(thermalmodel,tlist);
T = thermalResults.Temperature;

%get useful data points
%the tcas
I1=0; n1=0; I2=0; n2=0; I3=0; n3=0;
%the center
IC = 0; nc = 0;

for j=1:length(msh.Nodes)
    if (msh.Nodes(1,j)^2+(msh.Nodes(2,j)-0.002)^2)<(6.8e-4)^2 % selection nodes at the circle
        I1 = I1+T(j,end); % add temperatures
        n1 = n1+1; % add number of nodes
    elseif ((msh.Nodes(1,j)+0.001732)^2+(msh.Nodes(2,j)+0.001)^2)<(6.8e-4)^2
        I2 = I2+T(j,end);
        n2 = n2+1;
    elseif ((msh.Nodes(1,j)-0.001732)^2+(msh.Nodes(2,j)+0.001)^2)<(6.8e-4)^2
        I3 = I3+T(j,end);
        n3 = n3+1;
    elseif (msh.Nodes(1,j)^2+msh.Nodes(2,j)^2)<1e-3^2 %central temperature
        IC = IC+T(j,end);
        nc = nc+1;
    end
end
Tavg1=I1/n1;
Tavg2=I2/n3;
Tavg3=I3/n3;
tcaTemps = [Tavg1;Tavg2;Tavg3];
centerTemps = IC/nc;