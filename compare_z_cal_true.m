function [var_R,fitresult,fitresult_med,fitresult_mean,z_true,z_cal] = compare_z_cal_true(nR,nt,Q,dif,dA,dA_unc,dt,if_plot,L,R_choice)
state_ep=2:22;
%% z_cal
for rch = 1:nR
    for i = 1:length(state_ep)
        z_cal{i}(rch,:) =-(dA(rch,i+1:end)-dA(rch,1:end-i))/(dt*i);
    end
end

%% build H

for i = state_ep
    H{i-1} = zeros(nR,i*nR);
    comb = [ones(nR,1),i*ones(nR,1)];
    for j = 1: size(comb,1)
        if j == 1
            
            b = 1/(L(j)+L(j+1));
            c = b;
            H{i-1}(j,j+nR*(comb(j,1)-1):j+nR*(comb(j,1)-1)+1) = 0.5*[-c, b];
            H{i-1}(j,j+nR*(comb(j,2)-1):j+nR*(comb(j,2)-1)+1) = 0.5*[-c, b];
            
        else if j == nR
                a = 1/(L(j)+L(j-1));
                c = -a;
                H{i-1}(j,j-1+nR*(comb(j,1)-1):j+nR*(comb(j,1)-1)) = 0.5*[-a, -c];
                H{i-1}(j,j-1+nR*(comb(j,2)-1):j+nR*(comb(j,2)-1)) = 0.5*[-a, -c];
            else
                a = 1/(L(j)+L(j-1));
                b = 1/(L(j)+L(j+1));
                c = (L(j-1)-L(j+1))/((L(j)+L(j-1))*(L(j)+L(j+1)));
                H{i-1}(j,j-1+nR*(comb(j,1)-1):j+1+nR*(comb(j,1)-1)) = 0.5*[-a, -c, b];
                H{i-1}(j,j-1+nR*(comb(j,2)-1):j+1+nR*(comb(j,2)-1)) = 0.5*[-a, -c, b];
            end
        end
    end
end

%% build z_true
for i = 1:length(state_ep)
    tmp = [];
    for k = 1:nt-state_ep(i)+1
        x_true = [];
        for j = 1:state_ep(i)
            x_true = [x_true;Q(:,k+j-1)];
        end
        tmp = [tmp,H{i}*x_true];
    end
    z_true{i} = tmp;
end

%% calculate var of R from std of residual
for i = 1:length(state_ep)
    var_R_res{i} = 1*std(z_true{i}-z_cal{i},1,2).^2;
%     var_R_res2{i} = 1*rms(z_cal{i}-z_true{i},2).^2;
end

%% calculate var of R from error propogation of Hx
for i = 1:length(state_ep)
    var_Q_true = [];
    for j = 1:state_ep(i)
        var_Q_true = [var_Q_true;(0.1*mean(Q,2)).^2];
    end
    var_R_err_Hx{i} = diag([H{i}*diag(var_Q_true)*H{i}'])';
end

%% study var_R_res: fit
for i = 1:length(state_ep)
    med_std_res(i) =sqrt(median(var_R_res{i}));
    mean_std_res(i) =sqrt(mean(var_R_res{i}));
end
for j = 1:nR
    tmp1  = [];
    for i = 1:length(state_ep)
        tmp1 = [tmp1,var_R_res{i}(j)];
    end
    var_R_res_comp{j} = sqrt(tmp1);
    [fitresult{j}, gof] = createFit_stdR_res(state_ep-1, var_R_res_comp{j});
end
[fitresult_med, gof] = createFit_stdR_res(state_ep-1, med_std_res);
[fitresult_mean, gof] = createFit_stdR_res(state_ep-1, mean_std_res);

%% study var_R_res: plot
% plot \sigma_{res} vs. fitted one
tmp_name={};
% fac = mean(mean(Q))/mean(L);
if if_plot ==1
% std_res
figure,
for j = 1:nR
    tmp3 = "\sigma_{res} in Reach "+ "{"+num2str(j)+"}";
    tmp_name = [tmp_name;tmp3];
    plot(var_R_res_comp{j},'*','linewidth',1.2,'markersize',10);
    hold on
 
end
% fitted one
for j = 1:nR
    tmp4 = "fitted \sigma_{res} in Reach "+ "{"+num2str(j)+"}";
    tmp_name = [tmp_name;tmp4];
    plot( 0.5:0.1:21.5,fitresult{j}(0.5:0.1:21.5),'linewidth',1.2);
    hold on
end
% median
tmp3 = "median \sigma_{res} of reaches";
tmp4 = "fitted median \sigma_{res} of reaches";
tmp_name = [tmp_name;tmp3;tmp4];
% mean
tmp3 = "mean \sigma_{res} of reaches";
tmp4 = "fitted mean \sigma_{res} of reaches";
tmp_name = [tmp_name;tmp3;tmp4];
plot(med_std_res,'^','linewidth',1.2,'markersize',10);
hold on
 plot( 0.5:0.1:21.5,fitresult_med(0.5:0.1:21.5),'--','linewidth',1.2);

hold on
plot(mean_std_res,'g^','linewidth',1.2,'markersize',10);
 hold on
  plot( 0.5:0.1:21.5,fitresult_mean(0.5:0.1:21.5),'--','linewidth',1.2);
legend(tmp_name)
xlabel('\Deltat')
ylabel('\sigma_{res}')
grid on
xlim([0.5,21.5])
xticks(1:21)
set(gca,'FontSize',15)

%% study var_R_res percentage 
% plot \sigma_{res} vs. fitted one
tmp_name={};

fac_mean = mean(mean(Q,2))./mean(L);
fac = repmat(fac_mean,nR);%mean(Q,2)./L;
% std_res
figure,
for j = 1:nR
    tmp3 = "\sigma_{res}/factor in Reach "+ "{"+num2str(j)+"}";
    tmp_name = [tmp_name;tmp3];
    tmp1  = [];
    for i = 1:length(state_ep)
        tmp1 = [tmp1,var_R_res{i}(j)];
    end
    var_R_res_comp{j} = sqrt(tmp1);
    plot(var_R_res_comp{j}/fac(j)*100,'*','linewidth',1.2,'markersize',10);
    hold on
end
% fitted one
for j = 1:nR
    tmp4 = "fitted \sigma_{res}/factor in Reach "+ "{"+num2str(j)+"}";
    tmp_name = [tmp_name;tmp4];
    plot( 0.5:0.1:21.5,fitresult{j}(0.5:0.1:21.5)/fac(j)*100,'linewidth',1.2);
    hold on
end
% median
tmp3 = "median \sigma_{res}/factor of reaches";
tmp4 = "fitted median \sigma_{res}/factor of reaches";
tmp_name = [tmp_name;tmp3;tmp4];
% mean
tmp3 = "mean \sigma_{res}/factor of reaches";
tmp4 = "fitted mean \sigma_{res}/factor of reaches";
tmp_name = [tmp_name;tmp3;tmp4];
plot(med_std_res/fac_mean*100,'^','linewidth',1.2,'markersize',10);
hold on
 plot( 0.5:0.1:21.5,fitresult_med(0.5:0.1:21.5)/fac_mean*100,'--','linewidth',1.2);

hold on
plot(mean_std_res/fac_mean*100,'g^','linewidth',1.2,'markersize',10);
 hold on
  plot( 0.5:0.1:21.5,fitresult_mean(0.5:0.1:21.5)/fac_mean*100,'--','linewidth',1.2);
legend(tmp_name)
xlabel('\Deltat')
ylabel('\sigma_{res}/factor [%]')
grid on
xlim([0.5,21.5])
xticks(1:21)
set(gca,'FontSize',15)
end

%% calculate var of R from std of residual,using fitted one
for i = state_ep-1
    for j = 1:nR
    var_R_res_fit{i}(j,1) = fitresult{j}(i).^2;
    end
end

%% calculate var of R from std of residual,using fitted median one
for i = state_ep-1
    for j = 1:nR
    var_R_res_fit_med{i}(j,1) = fitresult_med(i).^2;
    end
end

%% calculate var of R from std of residual,using fitted mean one
for i = state_ep-1
    for j = 1:nR
    var_R_res_fit_mean{i}(j,1) = fitresult_mean(i).^2;
    end
end


%% make choice!
if R_choice ==1 % residual
    var_R = var_R_res;
else if R_choice ==2 % error propogation from Hx
        var_R = var_R_err_Hx;
    else if R_choice ==4 % residual but fitted one
            var_R = var_R_res_fit;
        else if R_choice ==5  % residual but fitted median one
                var_R = var_R_res_fit_med;
            else if R_choice ==6 % residual but fitted mean one
                    var_R = var_R_res_fit_mean;
                else
                    
                    var_R ={};
                end
            end
        end
    end
end

for i=1:length(state_ep)
    
    var_R_err_unc{i}=(dA_unc(:,state_ep(i):end).^2+dA_unc(:,1:end-state_ep(i)+1).^2)/(dt*(state_ep(i)-1))^2;
end
%% plot
%% z_true vs. z_obs
nR =3;
if if_plot ==1
    figure('position',[1000,1000,800,4000])
    for rch = 1:nR
        
        set(gca,'FontSize',15)
        title_gauged = [{'\Deltat=1 day','\Deltat=4 days','\Deltat=8 days','\Deltat=12 days','\Deltat=16 days','\Deltat=20 days'}];
        
        for i = 1:length(dif)
            pos = (i-1)*nR+rch;
            subplot(length(dif),nR,pos)
            %                      boundedline(1:length(obs_cal),obs_cal,abs(obs_cal-obs_true));
            plot(z_cal{dif(i)}(rch,:),'b')
            hold on
            plot(z_true{dif(i)}(rch,:),'r')
            mean(z_true{dif(i)}(rch,:))
%             plot(abs(z_true{dif(i)}(rch,:)-z_cal{dif(i)}(rch,:))./abs(z_true{dif(i)}(rch,:))*100,'r.')
%             ylim([0,1000])
            hold on
            title(title_gauged(i))
            
            % plot(time,Q_KF2{i},'b-+', 'Linewidth',1.2,'MarkerSize',3);
            grid on
            
            xlim([1,nt-1])
            if pos==17
                xlabel('day')
            end
            if pos==7
                ylabel('[m^2/s]')
            end
            hold on
            set(gca,'FontSize',15)
        end
        
    end
    
    
    set(gca,'FontSize',15)
    legend({'z_{simu.}','Hx_{true}'},'Position',[0.5 0.05 0.01 0.01])
end
set(gca,'FontSize',15)


%% var_R_res vs. var_R_err_Hx vs. var_R_err_unc
if if_plot ==1
figure
for rch = 1:nR
    set(gca,'FontSize',15)
    title_gauged = [{'\Deltat=1 day','\Deltat=4 days','\Deltat=8 days','\Deltat=12 days','\Deltat=16 days','\Deltat=20 days'}];
    for i = 1:length(dif)
        pos = (i-1)*nR+rch;
        subplot(length(dif),nR,pos)
        %                      boundedline(1:length(obs_cal),obs_cal,abs(obs_cal-obs_true));
        yline(sqrt(var_R_res{dif(i)}(rch)),'b--', 'Linewidth',1.2)
        hold on
        yline(sqrt(var_R_err_Hx{dif(i)}(rch)),'r--', 'Linewidth',1.2)
        hold on
        title(title_gauged(i))
        
        plot(dif(i)+1:nt,sqrt(var_R_err_unc{dif(i)}(rch,:)),'g-s', 'Linewidth',1.2,'MarkerSize',3);
        grid on
        
        xlim([dif(i)+1,nt])
        if pos==17
            xlabel('day')
        end
        if pos==7
            ylabel('[m^2/s]')
        end
        hold on
        ylim auto
        set(gca,'FontSize',15)
    end   
end
legend({"\sigma_{res}","\sigma_{err}^{\Sigmax}","\sigma_{err}^{\sigma{\deltaA}}"})
end
end