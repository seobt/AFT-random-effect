function jpdf = density_normal_censored(x,dij,mu,bb,cons,pseudo)

if (size(x,1)==1)
    jpdf=dij*normpdf(x,mu,sqrt(cons+bb))+(1-dij)*(1-normcdf(x,mu,sqrt(cons+bb)));
    return
end

rep=length(pseudo);
pseudo=pseudo*sqrt(bb);
a1=dij(:,1)==1;a2=dij(:,1)==0;
x1=x(a1,:);
mu1=mu(a1,:);
x2=x(a2,:);
mu2=mu(a2,:);
xmat1=repmat(x1,1,rep);    
nuint1=repmat(pseudo,size(x1,1),1);
xmat2=repmat(x2,1,rep);    
nuint2=repmat(pseudo,size(x2,1),1);
mumat1=repmat(mu1,1,rep);
mumat2=repmat(mu2,1,rep);
consmat1=ones(size(x1,1),rep)*sqrt(cons);
consmat2=ones(size(x2,1),rep)*sqrt(cons);
surv=(1-normcdf(xmat2+nuint2,mumat2,consmat2));
den=normpdf(xmat1+nuint1,mumat1,consmat1);
den=max(den,1.0e-20);
surv=max(surv,1.0e-20);
result=[den;surv];
jpdf=mean(prod(result));


    
