function [r_effect,va]=estimated_reffect_PAFT(data,beta_PAFT,PAFT_sigmae,PAFT_sigmav)

rep=101;
pseudo=norminv((1:rep)/(rep+1),0,1);
pseudo=pseudo/sqrt(mean(pseudo.^2)); 

temp=(pseudo+1)./(pseudo+1);
temp=temp/length(temp);
[r_effect,va]=p_reffect(data,pseudo*sqrt(PAFT_sigmav),temp,beta_PAFT,PAFT_sigmae,@test_linear);
