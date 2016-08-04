function ProbVal= Priorfun(x,XrangeMat)

%function to determine the prior probability; (P(parameter))
%x is a parameter vector,[Vp;Vl;Vt;Ql;Qt;Rl;Rt;Cl;Ct]
%XrangeMat is a length(x) by 2 each row is lower and upper bounds on x_{i}
%ProbVal is the value of the prior probability

for i=1:length(x)
    if x(i)<XrangeMat(i,1) || x(i)>XrangeMat(i,2) %check to see if x is outside lower or upper bound
    Prob(i) =0; %if outside the range, then prob=0     
    else
    Prob(i)=1/(abs(diff(XrangeMat(i,:)))); %the prob is inverse of the absolute val difference between the upper and lower bound
    end
end

ProbVal=prod(Prob); %the value returned is the product of all the probabilities 
