function plot_result(Q_tmp,Qest_weight,Qest_arith,Qest_med,nR,nt,Q_true,state_ep,orbit,reach,reserve,Q_withunc)
%% combine and plot Q
% nR=3;
% m=[180 180 180];
% figure('position',[1000,1000,800,4000])
% set(gca,'FontSize',18)
% title_gauged = [{'reach 1','reach 2','reach 3','reach 4','reach 5'}];
% for i = 1:5%nR
%     subplot(5,1,i)
%      set(gca,'FontSize',20)
%     idx  = reach{i};
%     plot(2:nt,Q_true(i,2:end),'k-','Linewidth',3);
%     hold on
%     plot(reserve{i},Q_withunc(i,reserve{i}),'cs','Linewidth',2);
%     hold on
%     for j = 1:state_ep/2%1:7
%         h3=plot(2:nt,Q_tmp{j}(i,:),'-','Linewidth',1);
%         set(h3,'color',m/255)
%         hold on
%     end
%     for j = state_ep/2+1:state_ep%8:14
%         h4=plot(2:nt,Q_tmp{j}(i,:),'-','Linewidth',1);
%         set(h4,'color',m/255)
%         hold on
%     end
%     %     for j = 15:22%state_ep/2+1:state_ep
%     %         plot(2:nt,Q_tmp{j}(i,:),'--');
%     %         hold on
%     %     end
%     
%     plot(2:nt,Qest_weight{1}(i,:),'b-','Linewidth',1.2);
% %     set(h1,'color',m/255)
%     hold on
%     plot(2:nt,Qest_arith{1}(i,:),'g-','Linewidth',1.2);
% %     set(h2,'color',m/255)
%     hold  on
%     plot(2:nt,Qest_med{1}(i,:),'r-','Linewidth',1.2);
% %     set(h2,'color',m/255)
%     hold  on
%     for j = 1:length(idx)
%         xline(orbit{idx(j)},'m-.','Linewidth',1)
%         hold on
%     end
%     xlim([2,nt])
%     ylim([min(Qest_med{1}(i,:))-1000,max(Qest_med{1}(i,:))+1000])
%     set(gca,'FontSize',20)
%     title(title_gauged(i))
%     xlabel('day')
%     if i ==3
%         ylabel('Discharge [m/s^3]')
%     end
%     if i == 5
%         tmp_name={};
%         for j = 1:state_ep
%             tmp = "x_{est}^"+ "{"+num2str(j)+"}";
%             tmp_name = [tmp_name;tmp];
%             
%             
%         end
%         tmp_name = ["x_{true}";"Q_{obs}";...
%             tmp_name;"x^{weight}_{est}";"x^{arith.}_{est}";"x^{median}_{est}";"Meas. day"];
%         %         tmp_name = [tmp_name];
%         legend(tmp_name,'Position',[0.5 0.05 0.01 0.01])
%        
%     end
% end

m=[180 180 180];
figure('position',[1000,1000,800,4000])
set(gca,'FontSize',18)
title_gauged = [{'reach 1','reach 2','reach 3'}];
for i = 1:3
    subplot(3,1,i)
     set(gca,'FontSize',15)
    idx  = reach{i};
    plot(2:nt,Q_true(i,2:end),'k-','Linewidth',3);
    hold on

    
    plot(2:nt,Qest_med{1}(i,:),'r-','Linewidth',3);
        plot(reserve{i},Q_withunc(i,reserve{i}),'bs','Linewidth',3);
    hold on
%     set(h2,'color',m/255)
    hold  on
    for j = 1:length(idx)
        xline(orbit{idx(j)},'-','color',[180/255 180/255 180/255],'Linewidth',1)
        hold on
    end
    xlim([2,nt])
   % ylim([0,1000])
    %ylim([min(Q_true(i,:))-1000,max(Q_true(i,:))+1000])
    set(gca,'FontSize',20)
%    title(title_gauged(i))
    xlabel('day')
    if i ==2
        ylabel('Discharge [m/s^3]')
    end
    if i == 3
        tmp_name = ["Q_{true}";"Q_{est}";"Q_{obs}";"Meas. day"];
        %         tmp_name = [tmp_name];
        legend(tmp_name,'Position',[0.5 0.05 0.01 0.01])
       
    end
end
%%
% best = [6,7,8,12];
% figure('position',[1000,1000,800,4000])
% set(gca,'FontSize',18)
% title_gauged = [{'reach 1','reach 2','reach 3','reach 4','reach 5'}];
% for i = 1:nR
%     subplot(nR,1,i)
%     plot(2:nt,Q_true(i,2:end),'-', 'Linewidth',1.2);
%     hold on
%     for j = 1:length(best)%1:7
%         plot(2:nt,Q_tmp{best(j)}(i,:),'-','Linewidth',1.2);
%         hold on
%     end
%
%     idx  = reach{i};
%     for j = 1:length(idx)
%         xline(orbit{idx(j)},'-.')
%         hold on
%     end
%     xlim([2,nt])
%     set(gca,'FontSize',20)
%     title(title_gauged(i))
%     xlabel('day')
%     if i ==3
%     ylabel('Discharge [m/s^3]')
%     end
%     if i == nR
%         tmp_name={};
%         for j = 1:length(best)
%             tmp = "x_{est}^"+ "{"+num2str(best(j))+"}";
%             tmp_name = [tmp_name;tmp];
%
%
%         end
%         tmp_name = ["x_{true}";tmp_name;"Meas."];
%         %         tmp_name = [tmp_name];
%         legend(tmp_name,'Position',[0.5 0.05 0.01 0.01])
%     end
% end

end

