% MCMC posterior sampling: NP bounding parameter estimation 
% **Note: In this simulation, the lung is considered the tumor**

function out=MCMC_BayesParamEst() 
clear; clc; clf; close all; 
%function out= MCMC_BayesParamEst()
%proposal is b+N(a given 0, V (where V is user defined) and b is the current state
%parameters=[Vp;Vl;Vt;Ql;Qt;Rl;Rt;Cl;Ct]

%parameters; consider rat (0.25 kg)
Vp=13.5; %plasma vol (mL), 1.7
Vl=19.6; %liver vol (mL), 1.3
Rl=3.84; %distribution ratio in liver,liver:blood, rat, prior to modification: 7.87, 3.12 
Rt=6.87; %distribution ratio in other,lung:blood,rat prior to modification: 10.13
Cl=13.856; %NP conc in liver 
Vt=2.1; %tumor vol (mL); 0.01
Ql=2.1; %flow rate to liver from plasma (ml/min),rat 1.5 
Qt=3.1; %flow rate to tumor from plasma (ml/min),rat prior to modification: 0.78
Ct=69.591; %NP conc in tumor 

%standard deviations
Vp_sd=0; %plasma vol
Vl_sd=0; %liver vol
Vt_sd=0; %tumor vol
Ql_sd=0; %flow rate, liver
Qt_sd=0; %flow rate, tumor
Rl_sd=0; %dist. ratio liver
Rt_sd=0; %dist. ratio tumor

% Sample and set initial zeta
Zeta_true=[Vp;Rl;Vt;Ql;Vl;Rt;Cl;Qt;Ct];
%Zeta_true=[Vp;Vl;Rl;Rt;Cl];

%User definied:
Noise=1; %noise, might be un-used
M=1000; %number of iterations

%user defined; alter size for parameter(s) of interest
Sigma=0.05*Zeta_true(1:2); %***

%call ode45 to determine conc. at Zeta_true
odefun=@(t,y)NPDistmodified(t,y,Zeta_true); %NPDistmodified returns concentration for matrix of parameters
tspan=[0:1:12]; %time span (min) see if moves on time 
y0=[0;0;5];%inital conditions
[t,y_ode_t]=ode45(odefun,tspan,y0); %call ode45 to determine non-linearly conc.
true=y_ode_t; %true holds conc. values for Zeta_true

Noise_Cov=max(max(abs(true)))/10; %denominator can be altered***
true_pertb=true+normrnd(0,Noise_Cov,13,3); %observation; true perturbation 


%Initial condition, can be modified; alter for number of parameters
%sampling
Zeta=1*ones(2,1); %***

%Zeta=Zeta_true+normrnd(0,1,9,1);
%Zeta=Zeta_true(1:5)+normrnd(0,min(Zeta_true(1:5))/10,5,1);

%draw sample 
for i=1:M 
Zeta_tilde=-1;
while min(Zeta_tilde)<0
Zeta_tilde=mvnrnd(Zeta,diag(Sigma.^2))'; %random vectors norm dist. mu=Zeta, Sigma=user defined, draw zeta_tilde from proposal
end

%alter for number of parameters interested in sampling:
Zeta_old=[ Zeta(1:2);Vt; Ql;Vl;Rt;Cl;Qt;Ct]; %***
Zeta_tilde=[Zeta_tilde(1:2);Vt;Ql;Vl;Rt;Cl;Qt;Ct]; %***

out=Zeta_tilde(1:2) %alter; output ***

alp=min(1,(PostFun(Zeta_tilde,Noise_Cov,true_pertb)./PostFun(Zeta_old,Noise_Cov,true_pertb))); %probability; assume symmetric proposal density 

%to accept with probability a:
U=rand(); %uniform random variable on [0,1] 
if U<alp
Zeta=Zeta_tilde(1:2); %accept new parameter values, else keep old parameter values, change **

%disp('accept');
%show(j)=Zeta_tilde(i);

else
    %disp('reject')
end

end
end


%future work

%figure out ode function
%play with M, NoiseCov and bounds 
%repeat for 100 samples, 500 samples, 1000 samples, etc.
%can change time. tolerance, skip over 3*1
%play with parameters and see what get 

%first derivative;