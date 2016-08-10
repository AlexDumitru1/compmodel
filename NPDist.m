function [dy] = NPDist(t,y,c )

%function for NP Distribution that considers spleen
% c=[Vp;Vl;Vt;Ql;Qt;Rl;Rt]

Vp=c(1);
Ql=c(2);
Rl=c(3);
Qt=c(4);
Rt=c(5);
Qs=c(6);
Rs=c(7);
Vl=c(8);
Rlt=c(9);
Vt=c(10);
Vs=c(11);
Rst=c(12);

dy=zeros(4,1);

%[liver, tumor, plasma,spleen];

dy(1)=(y(1)^2)*(-(Ql)/(Vl*Rl))+(y(2)^2)*(0)+(y(3)^2)*((Ql-Qs)/Vl)+(y(4)^2)*(-Qs/(Rs*Vl)); %liver 
dy(2)=0*(y(1)^2)+(-Qt/(Rt*Vt))*(y(2)^2)+(Qt/Vt)*(y(3)^2)+0*(y(4)^2); %tumor
dy(3)=(y(1)^2)*(Ql/(Vp*Rl))+(y(2)^2)*(Qt/(Vp*Rt))+(y(3)^2)*(-(Ql+Qt+Qs)*(1/Vp))+(y(4)^2)*(Qs/(Vp*Rs));%plasma
dy(4)=0*(y(1)^2)+(y(2)^2)*(0)+(y(3)^2)*(Qs/Vs)+(y(4)^2)*(-(Qs)/(Vs*Rs)); %spleen

end

