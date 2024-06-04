function [Phi_st,Q_st] = build_Phi_SWOT_month(nR,L_center_pos,c,D,tau,Q,coeff,state_ep,choice_Q)
R_crs_mat = zeros(state_ep*nR,state_ep*nR);
R_auto_mat = zeros(state_ep*nR,state_ep*nR);
% R_crs_mat1 = zeros(2*nR,2*nR);
% R_auto_mat1 = zeros(2*nR,2*nR);
% R_crs_mat2 = zeros(2*nR,2*nR);
% R_auto_mat2 = zeros(2*nR,2*nR);
% % Q_unc1 = std(Q(:,1:end-2),0,2);
% Q_unc2 = std(Q(:,2:end-1),0,2);
% Q_unc3 = std(Q(:,3:end),0,2);
nt = size(Q,2);
Q_unc1 = std(Q(:,1:end),1,2);

% 30 days a std
% idx_len = ceil(nt/30);
% for i = 1:idx_len
%     if i<idx_len
%         Q_unc2(:,i) = std(Q(:,1+30*(i-1):30*i),1,2);
%     else
%         Q_unc2(:,i) = std(Q(:,1+30*(i-1):end),1,2);
%     end
% end

% length according to month
fir = datenum(2022,1,1);
las = fir+nt-1;
dates = datevec(fir:las);
jie = max(dates(:,2));
for i = 1:jie
    idx_mon = find(i == dates(:,2));
    Q_unc2(:,i) = std(Q(:,idx_mon),1,2);
end
%%
for i = 1:nR
    idx1 = nR-i+1;
    %cross
    x_lags1 = repmat(L_center_pos,state_ep,1)-L_center_pos(i);
    idx2 = find(x_lags1(1:nR)<0);
    idx2 = repmat(idx2,state_ep,1);
    idx3 = find(x_lags1<0);
    tmp=1:state_ep;
    t_lags_cr=[];
    for k = 1:state_ep
        lag = reshape(repmat(tmp,nR,1),state_ep*nR,1);
        t_lags_cr{k} = lag.*86400;
        tmp = tmp-1;
        r_cr = corrModel_AD_eq( [c(i) D(i) tau(i)] , t_lags_cr{k},x_lags1);
        r_tmp_cr = [];
        if ~isempty(idx2)
            for j = 1:length(idx3)
                r_tmp_cr(j,1) = corrModel_AD_eq( [c(idx2(j)) D(idx2(j)) tau(idx2(j))] , (-t_lags_cr{k}(idx3(j))), (-x_lags1(idx3(j))));
            end
        end
        R_crs_mat(:,i+nR*(k-1)) = r_cr;
        R_crs_mat(idx3,i+nR*(k-1)) = r_tmp_cr;
    end
    
    
    %auto
    x_lags1 = repmat(L_center_pos,state_ep,1)-L_center_pos(i);
    idx2 = find(x_lags1(1:nR)<0);
    idx2 = repmat(idx2,state_ep,1);
    idx3 = find(x_lags1<0);
    tmp=0:state_ep-1;
    t_lags_at = [];
    for k = 1:state_ep
        lag = reshape(repmat(tmp,nR,1),state_ep*nR,1);
        t_lags_at{k} = lag.*86400;
        tmp = tmp-1;
        r_at = corrModel_AD_eq( [c(i) D(i) tau(i)] , t_lags_at{k}, x_lags1);
        r_tmp_at = [];
        if ~isempty(idx2)
            for j = 1:length(idx3)
                r_tmp_at(j,1) = corrModel_AD_eq( [c(idx2(j)) D(idx2(j)) tau(idx2(j))] ,(-t_lags_at{k}(idx3(j))),(-x_lags1(idx3(j))));
            end
        end
        R_auto_mat(:,i+nR*(k-1)) = r_at;
        R_auto_mat(idx3,i+nR*(k-1)) = r_tmp_at;
    end
    
end

% figure,imagesc(R_crs_mat)
% xlabel('reach')
% ylabel('reach')
% set(gca,'FontSize',20)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% %
% %
% figure,imagesc(R_auto_mat)
% xlabel('reach')
% ylabel('reach')
% set(gca,'FontSize',20)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])

%     sig = Q_unc;
%     R_crs_mat = R_crs_mat(1:nR,1:nR);
%     R_auto_mat = R_auto_mat(1:nR,1:nR);
% sig1 = abs([Q(:,1);Q(:,2)]*0.1);
% sig2 = abs([Q(:,2);Q(:,3)]*0.1);
for i = 1:jie%idx_len
    sig1 = repmat(Q_unc2(:,i),state_ep,1);
    sig2 = repmat(Q_unc2(:,i),state_ep,1);
    Sigma = diag(sig2)*R_auto_mat*diag(sig2);
    Sigma_delta = diag(sig1)*R_crs_mat*diag(sig2);
    Phi_st{i} = Sigma_delta*inv(Sigma);
    
    % Phi_st= diag(sig2)*R_crs_mat*pinv(R_auto_mat)*pinv(diag(sig2));
    % Q_st = diag([repmat(zeros(nR,1),state_ep-1,1);(mean(Q,2).^2)*(coeff)^2]);
    if choice_Q == 1
        Q_st{i} = Sigma -  Sigma_delta * pinv(Sigma) *Sigma_delta';
    else if choice_Q == 2
            Q_st = diag(repmat(mean(Q,2),state_ep,1).^2)*(coeff)^2;
        else if choice_Q == 3
                Q_st = diag([repmat(zeros(nR,1),state_ep-1,1);(mean(Q,2).^2)*(coeff)^2]);
            else
                Q_st = diag([repmat(zeros(nR,1),state_ep,1)]);
            end
        end
    end
end
% Q_st_plot(Q_st_plot<1e-5&Q_st_plot>-1e-5)=0.00001;

% figure
% imagesc((Q_st_plot))
% c=colorbar
% c.Label.String = 'log_{10}';
% xlabel('reach')
% ylabel('reach')
% set(gca,'FontSize',20)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])

% figure
% imagesc(Q_st)
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% % title('Q_st')
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% %
% figure
% imagesc(Phi_st)
% % xlabel('reach')
% % ylabel('reach')
% set(gca,'FontSize',30)
% axis equal
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% % title('Phi_st')
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% % % clim([0,1])
% figure
% imagesc(R_auto_mat)
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% axis equal
%
% % title('Sigma_delta_transpose')
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% figure
% imagesc(R_crs_mat)
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
%
% figure
% imagesc(inv(R_auto_mat))
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% axis equal
% % title('Sigma_delta_transpose')
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% figure,
% plot(D2,'*-','linewidth',1.2)
%  xlabel('reach')
% ylabel('D [m^2/s]')
% hold on
% plot(D,'o-','linewidth',1.2)
% legend('D_{NaN}')
% legend('D_{NaN}','D_{corr}')
% set(gca,'FontSize',20)
% figure
% imagesc( Sigma_delta * inv(Sigma) *Sigma_delta')
% xlabel('reach')
% ylabel('reach')
% set(gca,'FontSize',20)
% title('Sigma_delta * inv(Sigma) *Sigma_delta')
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% figure
% imagesc(diag(sig2))
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% figure
% imagesc(inv(diag(sig2)))
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'FontSize',30)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% U = [Q(:,1);Q(:,2);Q(:,2);Q(:,3)];
% figure
% imagesc(Sigma)
% xlabel('reach')
% ylabel('reach')
% % title('Sigma')
% set(gca,'FontSize',20)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% figure
% imagesc(inv(Sigma))
% xlabel('reach')
% ylabel('reach')
% title('inv(Sigma)')
% set(gca,'FontSize',20)
% axis equal
% axis([0.5,state_ep*nR+0.5,-inf,inf])
% for i = 1:nR
%     for j = 1:nR
%         mmmm(i,j) = R_auto_mat(i,j)-R_auto_mat(j,i);
%     end
% end
% figure
% for i = 1:5
%     mmm=corrModel_AD_eq( [c(i) D(i) tau(i)] , ones(1,250001)*86400 , 0:250000);
%     plot(mmm)
%     hold on
% %     plot(L_center_pos-L_center_pos(1),mmm(L_center_pos-L_center_pos(1)+1),'*')
%
% end
% legend({'1','2','3','4','5'})
end