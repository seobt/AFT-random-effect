function [f,fval,exi]=cnm_ms(data,support,weight,beta,meanfunction,sigma)
        
id=unique(data(:,1));
d=size(data,2)-1;
tmp=[];
l_support=length(support);
s_bound=3*abs(support);    
l_beta=length(beta)+1;
function g=loglik(theta)               
    th_beta=theta(1:l_beta-1);
    epsilon=theta(l_beta);       
    supp=theta(l_beta+1:l_beta+l_support);
    wei=theta(l_beta+l_support+1:length(theta));
    tmp=zeros(length(id),1);
    for h=1:length(id)
        g=id(h);
        tmp_data=data(data(:,1)==g,:);
        me=meanfunction(th_beta,tmp_data(:,3:d));
        ydata=tmp_data(:,2);
        dij=tmp_data(:,d+1);
        tmp(h)=log(density_censored_mix(ydata,dij,supp,wei,me,epsilon)); 
    end
    g=-sum(tmp);
end
   
bound=10*max(abs(beta))*beta'./beta';
options = optimoptions(@fmincon,'Display','none','Algorithm','interior-point');
lb=[-bound,0,-s_bound,zeros(1,l_support)];
ub=[bound,10,s_bound,ones(1,l_support)];
A=[zeros(1,l_beta+l_support),ones(1,l_support)];
[f,fval,exi]=fmincon(@loglik,[beta;sigma;support';weight'],A,1,[],[],lb,ub,[],options);
end


