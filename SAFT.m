function [beta,var_e,var_v,r_dist,log_lik]=SAFT(data,ini_beta,ini_sigmae)
p=size(data,2);

%% Proposed method: Model 3
% Set initial value
beta0=ini_beta(2:p-2)';
residual=[data(:,2)-data(:,3:p-1)*beta0,data(:,p)];
support=quantile(residual(:,1),[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]);
weight=support./support/length(support);
mmf2=@test_linear2;
[support,weight]=update_weight2(data,support,weight,beta0,mmf2,ini_sigmae);

for k=1:30  
    %% If the lastgest residual is censored,...       
    residual=[data(:,2)-data(:,3:p-1)*beta0,data(:,p)];
    data(residual(:,1)==max(data(:,2)-data(:,3:p-1)*beta0),p)=1;     

    %% Find new support for the random efffect distribution
    [new_support,grad_value]=find_support(data,beta0,support,weight,mmf2,ini_sigmae);

    %% Update random effect distribution
    support=unique([support,new_support]);
    weight=ones(1,length(support))/length(support);
    [support,weight]=update_weight2(data,support,weight,beta0,mmf2,ini_sigmae);
    
    %% Update parameters
    f=cnm_ms(data,support,weight,beta0,mmf2,ini_sigmae); % CNMMS
    f=f';
    beta0=f(1:p-3)';ini_sigmae=f(p-2);   
    f2=f;f2(1:p-2)=[];
    l_supp=length(f2)/2;
    support=f2(1:l_supp);
    weight=f2(l_supp+1:2*l_supp);
    [support,weight]=update_weight2(data,support,weight,beta0,mmf2,ini_sigmae);
    intercept=sum(support.*weight);
    var_v=sum(weight.*((support-intercept).^2));
    beta=[intercept,beta0'];
    var_e=ini_sigmae;
    r_dist=[support;weight];
    log_lik=loglik(data,support,weight,beta0,mmf2,ini_sigmae);

    if (max(grad_value)<0.1) % if the gradient is less than 0.1, stop
        return
    end
end