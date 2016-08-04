function y_ode=sampl_y(Zeta)

%determine concentration values
%input: parameters
%output: true conc. values, ode45 evaluated at parameter values 

y=[0;0;5]; 
odefun_somthng=@(t,y)NPDistmodified(t,y,Zeta); %NPDistmodified returns concentration for matrix of parameters
tspan=[0:1:12]; %time span (min)
y0=[0;0;5];%inital conditions
[t,y_ode]=ode45(odefun_somthng,tspan,y0); %call ode45 to determine non-linearly conc.

end
