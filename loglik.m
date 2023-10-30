function g=loglik(data,support,weight,beta,meanfunction,sigma)  
id=unique(data(:,1));
d=size(data,2)-1;
tmp=zeros(length(id),1);
    for h=1:length(id)
        g=id(h);
        tmp_data=data(data(:,1)==g,:);
        me=meanfunction(beta,tmp_data(:,3:d));
        ydata=tmp_data(:,2);
        dij=tmp_data(:,d+1);
        tmp(h)=log(density_censored_mix(ydata,dij,support,weight,me,sigma));        
    end
    g=sum(tmp);
end