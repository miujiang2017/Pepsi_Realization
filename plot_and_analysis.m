close all;
clear all;
clc

addpath('FullUncertainty');
addpath('Realizations');

file_KF1 = dir(['FullUncertainty/*.nc']);
file_KF2 = dir(['Realizations/*.mat']);
for m = 8%2:17
    river_name = file_KF1(m).name;
    grch1 = ncread(river_name,'/River_Info/gdrch');
    H_syn = ncread(river_name,'/Reach_Timeseries/H');
    W_syn = ncread(river_name,'/Reach_Timeseries/W');
    S_syn = ncread(river_name,'/Reach_Timeseries/S');
    nR_g = length(grch1);
    Obs = load(file_KF2(m-1).name);
    Q_read = Obs.Observations.obs.Q;
    S_read = Obs.Observations.obs.S;
    H_read = Obs.Observations.obs.H;
    dA_read = Obs.Observations.obs.dA;
    W_read = Obs.Observations.obs.W;
    
    Htrue = Obs.Observations.true.H;
    Qtrue = Obs.Observations.true.Q;
    Strue = Obs.Observations.true.S;
    dAtrue = Obs.Observations.true.dA;
    Wtrue = Obs.Observations.true.W;
    for j = 1:100
        
        for k = 1:100
            rlz = 100*(j-1)+k;
            Q_obs(:,:,rlz) = Q_read(:,:,j,k);
            S_obs(:,:,rlz) = S_read(:,:,j,k);
            H_obs(:,:,rlz) = H_read(:,:,j,k);
            dA_obs(:,:,rlz) = dA_read(:,:,j,k);
            W_obs(:,:,rlz) = W_read(:,:,j,k);
            
            for i = 1:nR_g
                [idx1] = find(~isnan(Q_obs(i,:,rlz))==1);
                r_Q{m}(i,rlz) = median(abs((Q_obs(i,idx1,rlz)-Qtrue(i,idx1))./Qtrue(i,idx1))*100);
                r_S{m}(i,rlz) = median(abs((S_obs(i,:,rlz)-Strue(i,:))./Strue(i,:))*100);
                r_dA{m}(i,rlz) = median(abs((dA_obs(i,2:end,rlz)-dAtrue(i,2:end))./dAtrue(i,2:end))*100);
                r_W{m}(i,rlz) = median(abs((W_obs(i,:,rlz)-Wtrue(i,:))./Wtrue(i,:))*100);
                r_H{m}(i,rlz) = median(abs((H_obs(i,:,rlz)-Htrue(i,:))./Htrue(i,:))*100);
            end
        end
    end
    
    na_true = Obs.Parameters.Optimal.na;
    b_true = Obs.Parameters.Optimal.b;
    A0_true = Obs.Parameters.Optimal.A0;
    na_obs = Obs.Parameters.Ensemble.na;
    b_obs = Obs.Parameters.Ensemble.b;
    A0_obs = Obs.Parameters.Ensemble.A0;
    
    r_na{m}=abs((na_true- na_obs)./na_true*100);
    r_A0{m}=abs((A0_true- A0_obs)./A0_true*100);
    r_b{m}=abs((b_true- b_obs)./b_true*100);
     
    %% plot real 
    i = 1;%reach
    ep = 8;%epoch
    idx2 = find(isnan(Q_obs(i,ep,:))==0);
    da = 2;
        figure
        
    subplot(1,3,1)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(A0_obs(i,:),10);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        xlabel('A_0 [m^2]')
        h = xline(mean(A0_obs(i,:)),'k-.','linewidth',da);
        k = xline(A0_true(i),'r-','linewidth',da);
        
    subplot(1,3,2)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(na_obs(i,:),10);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        xlabel('n_a [-]')
        h = xline(mean(na_obs(i,:)),'k-.','linewidth',da);
        k = xline(na_true(i),'r-','linewidth',da);
      
        
    subplot(1,3,3)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(b_obs(i,:),5);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        xlabel('b [-]')
        h = xline(mean(b_obs(i,:)),'k-.','linewidth',da);
        k = xline(b_true(i),'r-','linewidth',da);
        legend([k,h],{"true","mean"})
        
          figure
        
    subplot(1,5,1)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(Q_obs(i,ep,idx2),20);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',17)
        xlabel('Q [m^3/s]')
        h = xline(mean(Q_obs(i,ep,idx2)),'k-.','linewidth',da);
        k = xline(Qtrue(i,ep),'r-','linewidth',da);
        
            subplot(1,5,2)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(W_obs(i,ep,:),20);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',17)
        xlabel('W [m]')
          k = xline(Wtrue(i,ep),'r-','linewidth',da);
        h = xline(mean(W_obs(i,ep,:)),'k-.','linewidth',da);
         g = xline(mean(W_syn(i,ep,:)),'g-','linewidth',da);
        
                    subplot(1,5,3)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(H_obs(i,ep,:),20);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',17)
        xlabel('H [m]')
        k = xline(Htrue(i,ep),'r-','linewidth',da);
        h = xline(mean(H_obs(i,ep,:)),'k-.','linewidth',da);       
      g = xline(mean(H_syn(i,ep,:)),'g-','linewidth',da);
                            subplot(1,5,4)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(S_obs(i,ep,:)*1e5,20);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',17)
        xlabel('S [cm/km]')
        k = xline(Strue(i,ep)*1e5,'r-','linewidth',da);
        h = xline(mean(S_obs(i,ep,:))*1e5,'k-.','linewidth',da);
        g = xline(mean(S_syn(i,ep,:))*1e5,'g-','linewidth',da);
          legend([k,h,g],{"true","mean pertubed data","synthetic data"})
                            subplot(1,5,5)
%             [num,bin] = histcounts(Q_obs(i,1,idx2));
%         idx = find(num<1000);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
          histogram(dA_obs(i,ep,:),20);
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',17)
        xlabel('\deltaA [m^2]')
        k = xline(dAtrue(i,ep),'r-','linewidth',da);
        h = xline(mean(dA_obs(i,ep,:)),'k-.','linewidth',da);
        
        
    %%plot relative error
    figure
    title('Q')
    for i = 1:nR_g
        rt = r_Q{m}(i,:);
        subplot(ceil(nR_g/2),2,i)
        [num,bin] = histcounts(r_Q{m}(i,:));
        idx = find(num<100);
        
        for j = 1:length(idx)
            rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
        end
        histogram(rt,'BinLimits',[min(rt),max(rt)]);
        
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        subtitle(['Reach ',int2str(grch1(i))]);
        xlabel('[%]')
        h = xline(median(r_Q{m}(i,:)),'r-.','linewidth',1.2);
    end
    title('Q')
    legend(h,{"median of relative error"})
    
    %plot relative error of dA
    figure
    title('dA')
    for i = 1:nR_g
        rt = r_dA{m}(i,:);
        subplot(ceil(nR_g/2),2,i)
        [num,bin] = histcounts(r_dA{m}(i,:));
        idx = find(num<100);
        
        for j = 1:length(idx)
            rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
        end
        histogram(rt,'BinLimits',[min(rt),max(rt)]);
        
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        subtitle(['Reach ',int2str(grch1(i))]);
        xlabel('[%]')
        h = xline(median(r_dA{m}(i,:)),'r-.','linewidth',1.2);
    end
    title('dA')
    legend(h,{"median of relative error"})
    
    figure
    title('H')
    for i = 1:nR_g
        
        rt = r_H{m}(i,:);
        subplot(ceil(nR_g/2),2,i)
        [num,bin] = histcounts(r_H{m}(i,:));
        idx = find(num<100);
        
        for j = 1:length(idx)
            rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
        end
        histogram(rt,'BinLimits',[min(rt),max(rt)]);
        
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        subtitle(['Reach ',int2str(grch1(i))]);
        xlabel('[%]')
        h = xline(median(r_H{m}(i,:)),'r-.','linewidth',1.2);
    end
    title('H')
    legend(h,{"median of relative error"})
    
    
    figure
    title('S')
    for i = 1:nR_g
        rt = r_S{m}(i,:);
        subplot(ceil(nR_g/2),2,i)
        [num,bin] = histcounts(r_S{m}(i,:));
        idx = find(num<100);
        
        for j = 1:length(idx)
            rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
        end
        histogram(rt,'BinLimits',[min(rt),max(rt)]);
        
        %         xlim([median(r)-200,median(r)+200])
        set(gca,'FontSize',20)
        subtitle(['Reach ',int2str(grch1(i))]);
        xlabel('[%]')
        h = xline(median(r_S{m}(i,:)),'r-.','linewidth',1.2);
    end
    title('S')
    legend(h,{"median of relative error"})
%     figure
%     title('W')
%     for i = 1:nR_g
%         rt = r_W{m}(i,:);
%         subplot(ceil(nR_g/2),2,i)
%         [num,bin] = histcounts(r_W{m}(i,:));
%         idx = find(num<100);
%         
%         for j = 1:length(idx)
%             rt(rt>bin(idx(j)) & rt<bin(idx(j)+1)) = [];
%         end
%         histogram(rt,'BinLimits',[min(rt),max(rt)]);
%         
%         %         xlim([median(r)-200,median(r)+200])
%         set(gca,'FontSize',20)
%         subtitle(['Reach ',int2str(grch1(i))]);
%         xlabel('[%]')
%         h = xline(median(r_W{m}(i,:)),'r-.','linewidth',1.2);
%     end
%     title('W')
%     legend(h,{"median of relative error"})
%     
%     
%     figure
%     title('na')
%     for i = 1:nR_g
%         rt = r_na{m}(i,:);
%         subplot(ceil(nR_g/2),2,i)
%         histogram(rt,20);
%         set(gca,'FontSize',20)
%         subtitle(['Reach ',int2str(grch1(i))]);
%         xlabel('[%]')
%         h = xline(median(r_na{m}(i,:)),'r-.','linewidth',1.2);
%     end
%      title('na')
%     legend(h,{"median of relative error"})
%     
%     figure
%     title('b')
%     for i = 1:nR_g
%         rt = r_b{m}(i,:);
%         subplot(ceil(nR_g/2),2,i)
%         histogram(rt,20);
%         set(gca,'FontSize',20)
%         subtitle(['Reach ',int2str(grch1(i))]);
%         xlabel('[%]')
%         h = xline(median(r_b{m}(i,:)),'r-.','linewidth',1.2);
%     end
%      title('b')
%     legend(h,{"median of relative error"})
%     
%     figure
%    
%     for i = 1:nR_g
%         rt = r_A0{m}(i,:);
%         subplot(ceil(nR_g/2),2,i)
%         histogram(rt,20);
%         set(gca,'FontSize',20)
%         subtitle(['Reach ',int2str(grch1(i))]);
%         xlabel('[%]')
%         h = xline(median(r_A0{m}(i,:)),'r-.','linewidth',1.2);
%     end
%      title('A0')
%     legend(h,{"median of relative error"})
    
    Q_obs = [];
    dA_obs = [];
    S_obs = [];
    W_obs = [];
    H_obs = [];
    Qobs_median = [];
    %     figure
    %     for i = 1:nR
    %         subplot(nR,1,i),
    %         hist( Qobs_median(i,:))
    %         hold on
    %         xline(Qtrue_median(i),'r')
    %     end
end



