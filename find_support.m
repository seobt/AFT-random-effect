function [newx,ttt]=find_support(data,beta,support,weight,meanfunction,sigma)

d=size(data,2)-1;
id=unique(data(:,1));
length_id=length(id);
options = optimoptions('fmincon','Display','Off','Algorithm','interior-point');
%% Define gradient function
function answer=grad(pointer) 
    tmp=zeros(length_id,1);num=tmp;
    for h=1:length_id
        id_index=id(h);
        tmp_data=data(data(:,1)==id_index,:);  
        me=meanfunction(beta,tmp_data(:,3:d));
        ydata=tmp_data(:,2);
        dij=tmp_data(:,d+1);        
        tmp(h)=density_censored_mix(ydata,dij,support,weight,me,sigma);
        num(h)=density_censored(ydata,dij,pointer,me,sigma);
    end
answer=sum(num./tmp)-length(id);
answer=-answer;      
end
%%%%%%%%%%%%%%%%%%%%%% find new support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intercept=sum(support.*weight);
var_vi=sum(weight.*((support-intercept).^2));
tmp1=data(:,2)-meanfunction(beta,data(:,3:d));
max_range=max(tmp1)+5*max(var_vi,sigma);
min_range=min(tmp1)-5*max(var_vi,sigma);
x=(min_range:(max_range-min_range)/100:max_range);
x=sort([-3*abs(min_range),x,3*abs(max_range),support]);
newx=[];
ran=[];
for ii=1:length(x)        
    br(ii)=(gradient_support(data,x(ii)+1.0e-12,beta,support,weight,sigma,meanfunction)-...
        gradient_support(data,x(ii)-1.0e-12,beta,support,weight,sigma,meanfunction))/2.0e-12;    
    if (ii>1 && br(ii-1)>0 && br(ii)<0 ) 
        newx=[newx,(x(ii-1)+x(ii))/2];
        ran=[ran;x(ii-1),x(ii)];
    end
end


if (br(length(x))>0) 
    newx=[newx,max_range+1];
    ran=[ran;x(size(x,1),2),max(x)];
end

%%%%%%%%%%%%%%%%%%%%%%%%% Find support within a given interval%%%%%%%
final_can=[];
for t=1:length(newx)        
    [answer,fval,]=fmincon(@grad,newx(t),[],[],[],[],ran(t,1),ran(t,2),[],options);
    if (fval<-0.1) 
        final_can=[final_can,answer];
    end
end

if (~isempty(final_can)) 
    newx=final_can;
else
    ttt=0;
    return;
end

if (isempty(newx)) 
    ttt=0;
    return;
end

newx=unique(newx*1.0e+10)/1.0e+10;

newxx=[];
ttt=[];

for jj=1:length(newx)
     ttt=[ttt,gradient_support(data,newx(jj),beta,support,weight,sigma,meanfunction)];
    if (gradient_support(data,newx(jj),beta,support,weight,sigma,meanfunction)>0.0000001) 
        newxx=[newxx,newx(jj)];
    end
end




    
newx=unique(newxx);   
end