%% This files shows an example to fit SAFT using CGD data

load('cgd.mat')
data=cgd;

[beta_PAFT,PAFT_sigmae,PAFT_sigmav]=PAFT(data)
[beta_SAFT,SAFT_sigmae,SAFT_sigmav,r_dist,log_lik]=SAFT(data,beta_PAFT,PAFT_sigmae)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

   

