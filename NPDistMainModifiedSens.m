clc; clear; clf; close all; 

%main code for sensitivity analysis- look at change in concentration with
%respect to change in parameters
%use ODE45 to solve

% fill constant matrix
Vp=1.7; %plasma vol (mL)
Vl=1.3; %liver vol (mL)
Vt=1.01; %tumor vol (mL)
Vs=0.1; %volume of spleen (mL)
Qs=1.1378; %flow rate from spleen to liver, prior to modification:1.1348*0.1
Ql=10.868; %flow rate to liver from plasma, prior to modification:9.868*1.3
Qt=0.78; %flow rate to tumor from plasma, prior to modification: 0.78
Rl=10.87; %distribution ratio in liver, prior to modification: 7.87
Rt=7.13; %distribution ratio in other, prior to modification: 10.13
Rlt=14.9; %distribution ratio in liver from tumor, prior to modification: 14.9
Rs=1.17; %distribution ratio in spleen, prior to modification: 1.17
Rst=8.43; %distribution ratio in spleen from tumor, prior to modification: 8.43

%import concentration data to use for conc. values in derivative
%use interp1q 
A = importdata('workspace.mat');
time=A.o.x; %time 
y_liver=A.o.y(1,:); %liver conc
y_tumor=A.o.y(2,:); %tumor conc
y_plasma=A.o.y(3,:); %plasma conc
y_spleen=A.o.y(4,:); %spleen conc 
c=[Vp; Ql; Rl; Qt; Rt; Qs; Rs; Vl; Rlt; Vt;Vs;Rst;y_liver';y_tumor';y_plasma';y_spleen']'; %constant matrix 


%fill y matrix of initial concentrations
y=[0;0;0;0;0;0;0;0;0;0;0]; %mcg/mL
odefun=@(t,y)NPDistSens(time,t,y,c);


%call ode45 
tspan=[0:0.001:2]; %time span 
y0=[0;0;0;0;0;0;0;0;0;0;0];
ode45(odefun,tspan,y0);
%[t,yplot]=
%plot(t,yplot);
%'r','b','y','m','c','g','k','-r','ob','-b','ok');

%modify graph 
legend ('dCp/dVp','dCp/dQl','dCp/dRl','dCp/dQt','dCp/dRt','dCl/dVl','dCl/dQl','dCl/dRl','dCt/dVt','dCt/dQt','dCt/dRt');
title('Sensitivity Analysis; NP Concentration Change with respect to Parameters ');
ylabel('(d/dt)(dCi/dZi)');
xlabel('Time (hours)');
