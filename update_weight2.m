function [support,weight,exit]=update_weight2(data,support,weight,beta,meanfunction,sigma)

n=length(unique(data(:,1)));
d=size(data,2)-1;
m=length(support);
old_weight=weight;
S=zeros(n,m);
S1=S;S2=S;
lik_old=-1.0e+100;
for itr=1:30     
    id=unique(data(:,1));
    for i=1:length(id)
        id_index=id(i);
        tmp_data=data(data(:,1)==id_index,:);      
        me=meanfunction(beta,tmp_data(:,3:d));
        for j=1:length(weight)
            S1(i,j)=density_censored(tmp_data(:,2),tmp_data(:,d+1),support(j),me,sigma); 
            S2(i,j)=weight(j)*S1(i,j);
        end   
        S(i,:)=S1(i,:)/sum(S2(i,:));
        SS(i,:)=S2(i,:)/sum(S2(i,:));
    end
    gamma=sqrt(n*10^(-8));
    Sc=[gamma*S;ones(1,m)];
    onec=[2*gamma*ones(n,1);1];
    [weight,~,~,exit]=lsqnonneg(Sc,onec);
    weight=weight/sum(weight);
    weight=weight';  
    if (exit==0) 
        weight=mean(SS);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Armijo rule %%%%%%%%%%%%%%%%% 
    eta=weight-old_weight;
    sig=1/2;alpha=0.1;
    old_lik=loglik(data,support,old_weight,beta,meanfunction,sigma);
    for k=0:20        
        test_weight=abs(old_weight+sig^k*eta);
        test_weight=test_weight/sum(test_weight);
        new_lik=loglik(data,support,test_weight,beta,meanfunction,sigma);
        if (new_lik>=old_lik+ones(1,n)*S*eta'*sig^k*alpha) 
            weight=test_weight;
            break
        end
    end   
    lik_new=loglik(data,support,test_weight,beta,meanfunction,sigma);       
    if (lik_new-lik_old<1.0e-14) 
        break;
    end
    old_weight=weight;
    lik_old=lik_new;    
end
support(weight<=eps)=[];
weight(weight<=eps)=[];


