load('kindey.mat')
frailty_data=kidney;

load('bladder0.mat')
frailty_data=bladder0;
frailty_data(frailty_data(:,3)==0,:)=[];
frailty_data=[frailty_data,frailty_data(:,4)];
frailty_data(:,[1,4])=[];
frailty_data(:,2)=log(frailty_data(:,2));

load('cgd.mat')
frailty_data=cgd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmf=@test_linear;
data=frailty_data;
[n,m]=size(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation of beta assuming no random effect
[b,bint,r,rint,stats]=regress(frailty_data(frailty_data(:,m)==1,2),[frailty_data(frailty_data(:,m)==1,m),frailty_data(frailty_data(:,m)==1,3:(m-1))]);
cons=stats(4)*0.8;
bb=stats(4)/5;
beta=b;
beta_naive=beta';
data=frailty_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation of beta with normal random effect: Model 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep=101;
pseudo=norminv([1:rep]/(rep+1),0,1);
pseudo=pseudo/sqrt(mean(pseudo.^2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_cons=myfun_frailty_normal(data,bb,beta,mmf,cons,pseudo);
beta=beta_cons(1:3);
cons=beta_cons(4);
bb=beta_cons(5);
weight=1;
beta_normal=beta_cons';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Proposed method: Model 3
% Set initial value
beta0=beta_normal(2:3)';
sigma=beta_normal(4);
residual=[data(:,2)-data(:,3:4)*beta0,data(:,5)];
support=quantile(residual(:,1),[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]);
weight=support./support/length(support);
data=frailty_data;
if (data(residual(:,1)==max(data(:,2)-data(:,3:4)*beta0),5)==0) 
        'Largest observation is censored'
end
data(residual(:,1)==max(data(:,2)-data(:,3:4)*beta0),5)=1;
mmf2=@test_linear2;
[support,weight]=update_weight2(data,support,weight,beta0,mmf2,sigma);

   for k=1:30  
       %% If the lastgest observation is censored,....
       data=frailty_data;
       residual=[data(:,2)-data(:,3:4)*beta0,data(:,5)];
       data(residual(:,1)==max(data(:,2)-data(:,3:4)*beta0),5)=1;
       
        %% Find new support
        [new_support,ttt]=find_support(data,beta0,support,weight,mmf2,sigma);
        max(ttt)
        if (max(ttt)<0.1) 
            intercept=sum(support.*weight);
            var_vi=sum(weight.*((support-intercept).^2));
            beta_model3=[intercept,beta0',sigma,var_vi];
            break
        end
        support=unique([support,new_support]);
        weight=ones(1,length(support))/length(support);
        [support,weight]=update_weight2(data,support,weight,beta0,mmf2,sigma);
        
        [f,fval,exi]=cnm_ms(data,support,weight,beta0,mmf2,sigma); % CNMMS
        f=f';
        beta0=f(1:2)';sigma=f(3);   
        f2=f;f2(1:3)=[];
        l_supp=length(f2)/2;
        support=f2(1:l_supp);
        weight=f2(l_supp+1:2*l_supp);
        [support,weight]=update_weight2(data,support,weight,beta0,mmf2,sigma);
        intercept=sum(support.*weight);
        var_vi=sum(weight.*((support-intercept).^2));
        beta_model3=[intercept,beta0',sigma,var_vi];  
        loglik(data,support,weight,beta0,mmf2,sigma);
        %[beta_normal,beta_model3]
   end
   beta_normal
   beta_model3


  %% predict random effect 
    aaa=(pseudo+1)./(pseudo+1);
    aaa=aaa/length(aaa);
    [r_normal1, r_normal2]=p_reffect2(frailty_data,pseudo*sqrt(beta_normal(5)),aaa,beta_normal,mmf);
    %[r_normal1,r_normal2]=p_reffect(frailty_data,pseudo*sqrt(beta_normal(5))+beta_normal(1),aaa,beta_normal,mmf2);
    %[r_mix1,r_mix2]=p_reffect(data,support,weight,beta_model3,mmf2);
    [r_mix1,r_mix2]=p_reffect2(data,support-beta_model3(1),weight,beta_model3,mmf);
    %writematrix([r_normal1,r_mix1],'M2.csv') 
  %% predict log survival time
    p_mix=predict_y(data,beta_model3(1:3)',mmf,r_mix2);
    p_normal=predict_y(data,beta_normal(1:3)',mmf,r_normal2);
    [sqrt(mean((p_normal-data(:,2)).^2)),sqrt(mean((p_mix-data(:,2)).^2))]
    
  %% plot observed versus predicted  
    uncen=find(data(:,5)==1);
    cen=find(data(:,5)==0);
    figure
    plot(data(uncen,2),p_normal(uncen),'.b',data(cen,2),p_normal(cen),'.r',[min(data(:,2)),max(p_normal)],[min(data(:,2)),max(p_normal)],'-k')
    figure
    plot(data(uncen,2),p_mix(uncen),'.b',data(cen,2),p_mix(cen),'.r',[min(data(:,2)),max(p_normal)],[min(data(:,2)),max(p_normal)],'-k')

  %% CDF plot
    figure
    hold on
    %cdfplot(reffect,ones(length(reffect),1)/length(reffect))
    cdfplott(support'-intercept,weight')
    xx=[-3:0.02:3.5];
    yy=0*xx;
    for i=1:length(xx)
        yy(i)=normcdf(xx(i),0,sqrt(beta_normal(5)));
    end
    plot(xx,yy,'-k')
    legend({'SAFT','PAFT'},'Location','northwest')
    hold off   
%% Gradient plot
    xx=[];
    yy=[];
    for i=1:100
        xx(i)=1+i/10;
        yy(i)=gradient_support(data,xx(i),beta0,support,weight,sigma,@test_linear2);
    end
    plot(xx-intercept,yy)
    xlabel('\xi');ylabel('Gradient function')
    hline = refline([0 0]);hline.Color = 'r';
  
%% plot for estimated random effect
    plot(unique(frailty_data(:,1)),r_normal1,'k--o',unique(frailty_data(:,1)),r_mix1,'r--*')
    xlabel('Observation number');ylabel('Predicted random effect')
    legend({'PAFT','SAFT'},'Location','northwest')
    
%% plot for estimated random effect
    plot(unique(frailty_data(:,1)),r_normal1,'ko',unique(frailty_data(:,1)),r_mix1,'r*',[0,38],[0,0],'-b')
    xlabel('Patient ID');ylabel('Predicted random effect')
    legend({'PAFT','SAFT'},'Location','northwest')    


%% remove 21st data for kidney data
 
frailty_data(frailty_data(:,1)==21,:)=[];

for i=1:(size(frailty_data,1)/2)
    frailty_data(2*i-1,1)=i;
    frailty_data(2*i,1)=i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

