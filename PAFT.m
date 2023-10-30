function [beta,var_e,var_v]=PAFT(data)
%% Estimation of parameters under normal random effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmf=@test_linear;
m=size(data,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of beta assuming no random effect only with uncensored data
[b,~,~,~,stats]=regress(data(data(:,m)==1,2),[data(data(:,m)==1,m),data(data(:,m)==1,3:(m-1))]);
cons=stats(4)*0.8;
bb=stats(4)/5;
beta=b;
% Estimation of beta with normal random effect using numerical integration
rep=101;
pseudo=norminv((1:rep)/(rep+1),0,1);
pseudo=pseudo/sqrt(mean(pseudo.^2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_cons=myfun_frailty_normal(data,bb,beta,mmf,cons,pseudo);
beta=beta_cons(1:m-2)';
var_e=beta_cons(m-1);
var_v=beta_cons(m);
