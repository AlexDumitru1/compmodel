clc; clear;


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

c=[Vp; Ql; Rl; Qt; Rt; Qs; Rs; Vl; Rlt; Vt;Vs;Rst]';

%fill y matrix of initial concentrations

y=[0;0;5;0]; %mcg/mL
odefun=@(t,y)NPDist(t,y,c);

%liver; tumor; plasma; spleen
%call ode45 
tspan=[0:0.0001:2]; %time span 

y0=[0;0;5;0];
o=ode45(odefun,tspan,y0);
legend ('Liver', 'Tumor', 'Plasma','Spleen');
title('Nanoparticle Distribution');
ylabel('NP Concentration (mcg/mL)');
xlabel('Time (hours)');
