function plot_validation9(state_ep,med_corr,med_NSE,med_rRMSE,med_corr_unc,med_NSE_unc,med_rRMSE_unc,med_corr_obs,med_NSE_obs,med_rRMSE_obs)
tmp_name={};
grenze = [inf,10,100];
for j = 1:state_ep
    tmp = num2str(j);
    tmp_name = [tmp_name;tmp];  
end
tmp_name = [tmp_name;"weight";"arith.";"median"];
figure,

subplot(3,1,1)
m = 1;
for i = 1
    
    idx = find(abs(med_corr{1,i})<=grenze(1));
    if isempty(idx)
       plot(1:25,100*ones(25,1),'-*', 'Linewidth',1.2,'MarkerSize',10)
       hold on
    else
                max_val(m) = max(med_corr{1,i}(idx));
                min_val(m) = min(med_corr{1,i}(idx));
        m = m+1;
        plot(idx,med_corr{1,i}(idx)','-*', 'Linewidth',1.2,'MarkerSize',10)
        hold on
    end
end
yline(med_corr_obs{1},'g','Linewidth',2)
min_val(m) = med_corr_unc{1};
max_val(m) = med_corr_unc{1};
xticks(1:25)
xticklabels(tmp_name)

% subtitle('Mean correlation')
% ylim([min(min_val)-0.1,max(max_val)+0.1])
ylim([-1,1])
xlabel('x_{est}')
ylabel('[-]')
set(gca,'FontSize',20)
max_val=[];
min_val=[];

subplot(3,1,2)

m = 1;
for i =  1
    idx = find(abs(med_NSE{1,i})<=grenze(2));
        if isempty(idx)
        plot(1:25,100*ones(25,1),'-*', 'Linewidth',1.2,'MarkerSize',10)
        hold on
        else
        max_val(m) = max(med_NSE{1,i}(idx));
        min_val(m) = min(med_NSE{1,i}(idx));
        m = m+1;
    plot(idx,med_NSE{1,i}(idx)','-*', 'Linewidth',1.2,'MarkerSize',10)
    hold on
        end
end

yline(med_NSE_obs{1},'g','Linewidth',2)
min_val(m) = med_NSE_unc{1};
max_val(m) = med_NSE_unc{1};
xticks(1:25)
xticklabels(tmp_name)
subtitle('Median NSE (0<NSE<1)')% 
% subtitle('Mean NSE')
ylim([-10,1])
% ylim([min(min_val)-0.1,max(max_val)+0.1])
xlabel('x_{est}')
ylabel('[-]')
set(gca,'FontSize',20)
max_val=[];
min_val=[];
subplot(3,1,3)
m=1;
for i =  1
    idx = find(abs(med_rRMSE{1,i})<=grenze(3));
    if isempty(idx)
       plot(1:25,100*ones(25,1)','-*', 'Linewidth',1.2,'MarkerSize',10)
       hold on
    else
                max_val(m) = max(med_rRMSE{1,i}(idx));
                min_val(m) = min(med_rRMSE{1,i}(idx));
        m = m+1;
        plot(idx,med_rRMSE{1,i}(idx)','-*', 'Linewidth',1.2,'MarkerSize',10)
        hold on
    end
end
yline(med_rRMSE_obs{1},'g','Linewidth',2)
min_val(m) = med_rRMSE_unc{1};
max_val(m) = med_rRMSE_unc{1};
subtitle('Median rRMSE (0%<rRMSE<10%)')
% subtitle('Mean rRMSE')
xlabel('x_{est}')
ylabel('[%]')
ylim([0,20])
% ylim([min(min_val)-0.1,max(max_val)+0.1])
set(gca,'FontSize',20)
xticks(1:25)
xticklabels(tmp_name)
 legend({'Q_{est} from 10k realizations','Q_{obs} from the median of 10k realizations'})