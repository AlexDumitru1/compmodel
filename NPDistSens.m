%NPDistSens
%most up to date code for sens. analysis 

function [dy] = NPDistSens(time,t,y,c )
%the constant matrix is in the form [Vp; Ql; Rl; Qt; Rt; Qs; Rs; Vl; Rlt; Vt;Vs;Rst;y_liver;y_tumor;y_plasma;y_spleen];

%constants 
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


dL_dt_y=c(13:40)'; %liver w time
dT_dt_y=c(41:68)'; %tumor w time
dP_dt_y=c(69:96)'; %plasma w time
dS_dt_y=c(97:124)';  %spleen w time


dy=zeros(11,1); %initialize 


%interpolate the data so that we do not treat the concentrations as
%constant
%time is t=[0:0.0001:2]

dL_dt=interp1q(time,dL_dt_y,t);
dT_dt=interp1q(time,dT_dt_y,t);
dP_dt=interp1q(time,dP_dt_y,t);
dS_dt=interp1q(time,dS_dt_y,t);

%derivatives, change in concentration with respect to parameters 
dy(1)=y(1)*(((Ql+Qt)/(Vp^2))-((Ql+Qt)/Vp)); %dCp/dVp
dy(2)=(1/Vp*Rl)*((dP_dt/Vl)-(dL_dt/(Vl*Rl)))-(2/Vp)*y(2)-((Ql*y(2))/Vp); %dCp/dQl
dy(3)=(-Ql/Vp)*(1/(Rl^2))*((Ql*dL_dt)/(Vl*Rl^2)); %dCp/dRl
dy(4)=((1/Vp)*((dP_dt/Vt)-(dT_dt/(Vt*Rt)))*(1/Rt))-(3/Vp)*y(4)-y(4)*(Qt/Vp); %dCp/dQt
dy(5)=(-Qt/Vt)*(1/Rt^2)*((Qt*dT_dt)/(Vt*Rt^2)); %dCp/dRt
dy(6)=(((Ql*dL_dt)/(Rl*Vl^2))-(Ql/(Vl*Rl)))*y(6); %dCl/dVl
dy(7)=(1/Vl)*((dL_dt/(Vp*Rl))-(dP_dt/Vp))-(2/(Rl*Vl))*y(7)-y(7)*(Ql/(Vl*Rl));%dCl/dQl
dy(8)=((Ql/(Vl*Rl^2))-(Ql/(Vl*Rl)))*y(8); %dCl/dRl
dy(9)=((Qt/(Rt*Vt^2))-(Qt/(Vt*Rt)))*y(9); %dCt/dVt 
dy(10)=(1/Vt)*((dT_dt/(Vp*Rt))-(dP_dt/Vp))-(2/(Rt*Vt))*y(10)-y(10)*(Qt/(Vt*Rt));%dCt/dQt
dy(11)=((Qt/(Vt*Rt^2))-(Qt/(Vt*Rt)))*y(11); %dCt/dRt

end

