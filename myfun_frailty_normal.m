function [f,fval,exi]=myfun_frailty_normal(data,bb,beta,meanfunction,cons,pseudo)
    
d=size(data,2)-1;
l_beta=length(beta);

function g=loglik(theta)            
    g=0;
    id=unique(data(:,1));
    th_beta=theta(1:l_beta);
    epsilon=theta(l_beta+1);
    th_bb=theta(l_beta+2);
    for h=1:length(id)
        id_index=id(h);
        tmp_data=data(data(:,1)==id_index,:);        
        me=meanfunction(th_beta,tmp_data(:,3:d));    
        
   
        if (prod(tmp_data(:,d+1))==1)
            tmp_mix=mvnpdf(tmp_data(:,2),me,epsilon*eye(length(tmp_data(:,2)))+th_bb);           
        else
            tmp_mix=density_normal_censored(tmp_data(:,2),tmp_data(:,d+1),me,th_bb,epsilon,pseudo); 
        end
        g=g+log(tmp_mix);
    end
    g=-g;    
end

bound=10*max(abs(beta))*beta'./beta';
options = optimoptions(@fmincon,'Display','none','Algorithm','interior-point');
[f,fval,exi]=fmincon(@loglik,[beta;cons;bb],[],[],[],[],[-bound,0,0],[bound,90,90],[],options);
end



