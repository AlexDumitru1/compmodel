function Posterior=PostFun(Zeta,NoisCov,true_pertb) 
%calculate posterior probability 
%Zeta=matrix of parameters,[Vp;Vl;Vt;Ql;Qt;Rl;Rt;Cl;Ct]
%NoisCov=standard deviation


%make a matrix of upper and lower bound ******

%standard deviation: 
Zeta_dev=[0,0,0,0,0;20*ones(1,5)]';


%call sampl_y function to determine the concentration value (sampy) for the input
%parameter matrix (Zeta)
sampy=sampl_y(Zeta);
    
%compute the likhd; P(C|Zeta)
Likhd=prod(prod(normpdf(true_pertb,sampy,NoisCov))); %likhd is product of normally distributed func evaluated at Cp (given time and parameters), mu=sol, sigma=NoisCov (same length as soln)

%determine the prior using the Prior function, must be altered for number
%of parameters sampling:
Prior=Priorfun(Zeta(1:2),Zeta_dev); %determine the prior prob. at the parameters ***

%compute the posterior probability
Posterior=Likhd*Prior; %ignore constant term; posterior is product of likehd and prior prob.



