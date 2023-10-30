function answer=gradient_support(data,phi,beta,support,weight,sigma,meanfunction)

d=size(data,2)-1;
id=unique(data(:,1));
tmp=zeros(length(id),1);
num=tmp;

for h=1:length(id)
    g=id(h);
    tmp_data=data(data(:,1)==g,:);
    me=meanfunction(beta,tmp_data(:,3:d));
    ydata=tmp_data(:,2);
    dij=tmp_data(:,d+1);
    tmp(h)=density_censored_mix(ydata,dij,support,weight,me,sigma);
    num(h)=density_censored(ydata,dij,phi,me,sigma);
end
answer=sum(num./tmp)-length(id);
end

