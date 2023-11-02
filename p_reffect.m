function [p_reffect,var_reffect]=p_reffect(frailty_data,support,weight,beta,sigmae,mmf2)


%% Predict random effect 
data=frailty_data;
d=size(data,2)-1;
id=unique(data(:,1));
length_id=length(id);
temp={};
meanfunction=mmf2; 
l_support=length(support);

for iid=1:length(id)
    id_index=id(iid);
    temp{id_index}=data(data(:,1)==id_index,:);
end     

tmp=zeros(length_id,l_support);num=tmp;
for k=1:l_support

    for h=1:length_id
        id_index=id(h);
        tmp_data=temp{id_index};
        me=meanfunction(beta',tmp_data(:,3:d));
        ydata=tmp_data(:,2);
        dij=tmp_data(:,d+1);        
        tmp(h,k)=max(density_censored_mix(ydata,dij,support,weight,me,sigmae),1.0e-10);
        num(h,k)=weight(k)*density_censored(ydata,dij,support(k),me,sigmae);
    end
end
pos=num./tmp;

a=0*pos(:,1);

a=0*pos(:,1);b=0*pos(:,1);

for i=1:l_support
    a=a+support(i)*pos(:,i);
    b=b+support(i)^2*pos(:,i);
end
p_reffect=a;
var_reffect=b-a.^2;


