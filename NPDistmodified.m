function [dy] = NPDistmodified(t,y,c )
%parameters=[Vp;Vl;Vt;Ql;Qt;Rl;Rt]
% Vp=c(1); %plasma volume
% Vl=c(2); %liver vol
% Vt=c(3); %tumor vol
% Ql=c(4); %liver flow rate
% Qt=c(5); %tumor flow rate
% Rl=c(6); %liver dist. ratio
% Rt=c(7); %tumor dist. ratio
% Cl=c(8); %liver np conc.
% Ct=c(9); %tumor np conc

%first 3 (change) parameters are the ones we are sampling
Vp=c(1); %plasma volume
Rl=c(2); %dist. ratio liver
Vt=c(3); %tumor volume
Ql=c(4); %liver flow rate
Vl=c(5); %liver volume
Rt=c(6); %dist. ratio tumor 
Cl=c(7); %NP Conc. liver
Qt=c(8); %tumor flow rate
Ct=c(9); %conc. NP tumor 

%[liver, tumor, plasma]

dy(1)=(y(1)^2)*(-(Ql)/(Vl*Rl))+(y(2)^2)*(0)+(y(3)^2)*((Ql)/Vl); %liver 
dy(2)=0*(y(1)^2)+(-Qt/(Rt*Vt))*(y(2)^2)+(Qt/Vt)*(y(3)^2); %tumor
dy(3)=(y(1)^2)*(Ql/(Vp*Rl))+(y(2)^2)*(Qt/(Vp*Rt))+(y(3)^2)*(-(Ql+Qt)*(1/Vp));%plasma

dy=dy';
end