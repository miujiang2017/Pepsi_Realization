%% new process model:
%22/12-epoch state vector
%concept in ppt:24_Mar, page 2

close all;
clear all;
clc
addpath("boundedline-pkg-master/Inpaint_nans")
addpath("boundedline-pkg-master/boundedline")
addpath("boundedline-pkg-master/catuneven")
addpath("boundedline-pkg-master/singlepatch")
addpath('FullUncertainty');
addpath('Realizations');
% load discharge
file_KF = dir(['FullUncertainty/*.nc']);

% delete_reach{1} = [];
% delete_reach{2} = [1];
% delete_reach{3} = [3];
% delete_reach{4} = [5];
% delete_reach{5} = [1 3];
% delete_reach{6} = [3 5];
% delete_reach{7} = [1 3 5];
% delete_reach{8} = [1 2 3 5];
% delete_reach{9} = [1 3 4 5];
% delete_reach{10} = [1 2 3 4 5];
% delete_reach{1} = [];
% delete_reach{2} = [4];
% delete_reach{3} = [4 8];
% delete_reach{4} = [1 4 8];
% delete_reach{5} = [1 4 6 8];
% delete_reach{6} = [1 2 4 6 8];
% delete_reach{7} = [1 2 4 5 6 8];
% delete_reach{8} = [1 2 3 4 5 6 8];
% delete_reach{9} = 1:8;
% percent_dl = [0.05, 0.1, 0.2, 0.3, 0.4];
rvr = [4,10,11,13,14,15,16];
use_Q_obs =1;
use_Q_rlz =1 ;
for n = 1%1:length(rvr)
    m =16;%;rvr(n)
    % for m = 1:5%length(file_KF)
    
    
    %% read data
    river_name = file_KF(m).name
    
    %         grch = 1:9;
    Q = ncread(river_name,'/Reach_Timeseries/Q');
    nR = size(Q,1);
    nt = size(Q,2);
    grch = 1:nR;
    grch1 =ncread(river_name,'/River_Info/gdrch');
    Q = Q(grch,:);
    h = ncread(river_name,'/Reach_Timeseries/H');
    h = h(grch,:);
    W = ncread(river_name,'/Reach_Timeseries/W');
    W = W(grch,:);
    A = ncread(river_name,'/Reach_Timeseries/A');
    A = A(grch,:);
    S = ncread(river_name,'/Reach_Timeseries/S');
    S = S(grch,:);
    xs_rch =ncread(river_name,'/XS_Timeseries/xs_rch');
    n_node =ncread(river_name,'/XS_Timeseries/n');
    QWBM = ncread(river_name,'/River_Info/QWBM');
    Strue= ncread(river_name,'/Reach_Timeseries/Strue');
    Strue = Strue(grch,:);
    htrue = ncread(river_name,'/Reach_Timeseries/Htrue');
    htrue = htrue(grch,:);
    Wtrue = ncread(river_name,'/Reach_Timeseries/Wtrue');
    Wtrue = Wtrue(grch,:);
    L = diff(ncread(river_name,'/River_Info/rch_bnd'));
    
    L_center_pos = cumsum(L)-L/2; % reach center position
    L_center_pos = L_center_pos(grch,:);
    L_center = [L_center_pos(1);diff(L_center_pos)]; % reach distance between centers
    
    if length(L)>length(grch)
        m
        'dams'
    end
    % uncertainty per reach
    Q_unc = std(Q,1,2);
    
    
    %% swot simulator
    
    [orbit,reach,delete,reserve] = passplan(m,nR,nt);
    
    
    %% Observations
    file_KF = dir(['Realizations/*.mat']);
    Obs = load(file_KF(m-1).name);
    dA_read = Obs.Observations.obs.dA;
    Q_read = Obs.Observations.obs.Q;
    dAerr_obs = Obs.ObsError.dA;
    Qerr_obs = Obs.Qerror;
    %     na_obs = Obs.Parameters.Ensemble.na;
    %     b_obs = Obs.Parameters.Ensemble.b;
    %     A0_obs = Obs.Parameters.Ensemble.A0;
    S_read = Obs.Observations.obs.S;
    W_read = Obs.Observations.obs.W;
    
    for j = 1:100
        for k = 1:100
            rlz = 100*(j-1)+k;
            dA_obs(:,:,rlz) = dA_read(:,:,j,k);
            Q_obs(:,:,rlz) = Q_read(:,:,j,k);
            W_obs(:,:,rlz) = W_read(:,:,j,k);
            S_obs(:,:,rlz) = S_read(:,:,j,k);
        end
    end
    
    tmp_Q = isnan(Q_obs(:,:,1));
    for j = 2:10000
        tmp_Q = tmp_Q+isnan(Q_obs(:,:,j));
    end
    str_name1 = ['Obs.dAopt',num2str(m),'.csv'];
    str_name2 = ['Obs.dAunopt',num2str(m),'.csv'];
    dA_sg = table2array(readtable(str_name1));
    dA_unc_sg = table2array(readtable(str_name2));
    %     dA_cell_sg = table2cell(readtable(str_name1));
    %     dA_unc_cell_sg = table2cell(readtable(str_name2));
    
    %% calculate median Q, median dA, var_Q, and var_dA
    med_Q_obs = [];
    med_dA_obs = [];
    med_W_obs = [];
    med_S_obs = [];
    for i = 1:length(grch1)
        for j = 1:nt
            [idx1] = find(~isnan(Q_obs(i,j,:))==1);
            med_Q_obs(i,j) =median(Q_obs(i,j,idx1));
            med_dA_obs(i,j) =median(dA_obs(i,j,:));
            med_W_obs(i,j) =median(W_obs(i,j,:));
            med_S_obs(i,j) =median(S_obs(i,j,:));
%             sigma_Q_obs(i,j) =sqrt(sum(Qerr_obs(i,j,idx1).^2)/10000);
             sigma_Q_obs(i,j) = iqr(Q_obs(i,j,idx1))/1.349;
            sigma_dA_obs(i,j) = iqr(dA_obs(i,j,:))/1.349;
%                       sigma_dA_obs(i,j) =sqrt(sum(dAerr_obs(i,j,:).^2)/10000);
        end
    end
    

    %% Q  with 30% as unc
    load(['Q',num2str(m),'.mat'])
    med_dA_obs = dA_sg;
    sigma_dA_obs =dA_unc_sg ;
        
    Q_withunc_1 = Q_withunc; % synthetic one
     % realization contains Nan for 10k rlz
        [idx_r,idx_c]= find(isnan(med_Q_obs)==1);

    med_Q_obs(idx_r,idx_c) = Q_withunc_1(idx_r,idx_c);
    
    
    Q_withunc_2(grch1,:) = med_Q_obs;
    diff_idx = setdiff(grch,grch1);
    Q_withunc_2(diff_idx,:) = Q_withunc(diff_idx,:);
    sigma_Q_withunc = [];
    if use_Q_rlz == 1
        Q_withunc = Q_withunc_2;
        sigma_Q_obs(idx_r,idx_c) = (abs(Q(idx_r,idx_c)-Q_withunc_2(idx_r,idx_c)));
        sigma_Q_withunc(grch1,:) = sigma_Q_obs;
        sigma_Q_withunc(diff_idx,:) = min(sigma_Q_obs,[],'all');
    else
        Q_withunc  = Q_withunc(grch,:);
    end
    Q_res_unc = std(Q(:,2:end)-Q_withunc(:,2:end),1,2);
    save(['Q_rlz',num2str(m),'.mat'],'Q_withunc_2')
    save(['Q_syn',num2str(m),'.mat'],'Q_withunc_1')
    save(['reserve',num2str(m),'.mat'],'reserve')
    
    %% validate Q_obs
    for i = 1:length(L_center)
        Q_true = Q(i,reserve{i});
        Q_comp_unc = Q_withunc_1(i,reserve{i});
        Q_comp_unc2 = Q_withunc_2(i,reserve{i});
        Q_obs_mean1(i,1) = mean(Q_comp_unc);%syn
        Q_obs_mean2(i,1) = mean(Q_comp_unc2);%ptb
        corr_unc(i) = corr(Q_true',Q_comp_unc');
        corr_obs(i) = corr(Q_true',Q_comp_unc2');
        RMSE_unc(i) = sqrt(sum((Q_true-Q_comp_unc).^2));
        RMSE_obs(i) = sqrt(sum((Q_true-Q_comp_unc2).^2));
        rRMSE_unc(i) = RMSE_unc(i)./sqrt(sum(Q_true.^2))*100;
        rRMSE_obs(i) = RMSE_obs(i)./sqrt(sum(Q_true.^2))*100;
        NSE_unc(i) = 1-sum((Q_true-Q_comp_unc).^2)/sum((Q_true-mean(Q_true)).^2);
        NSE_obs(i) = 1-sum((Q_true-Q_comp_unc2).^2)/sum((Q_true-mean(Q_true)).^2);
        
        Q_true = Q(i,:);
        Q_comp_unc = Q_withunc_1(i,:);
        Q_comp_unc2 = Q_withunc_2(i,:);
        corr_unc_daily(i) = corr(Q_true',Q_comp_unc');
        corr_obs_daily(i) = corr(Q_true',Q_comp_unc2');
        RMSE_unc_daily(i) = sqrt(sum((Q_true-Q_comp_unc).^2));
        RMSE_obs_daily(i) = sqrt(sum((Q_true-Q_comp_unc2).^2));
        rRMSE_unc_daily(i) = RMSE_unc_daily(i)./sqrt(sum(Q_true.^2))*100;
        rRMSE_obs_daily(i) = RMSE_obs_daily(i)./sqrt(sum(Q_true.^2))*100;
        NSE_unc_daily(i) = 1-sum((Q_true-Q_comp_unc).^2)/sum((Q_true-mean(Q_true)).^2);
        NSE_obs_daily(i) = 1-sum((Q_true-Q_comp_unc2).^2)/sum((Q_true-mean(Q_true)).^2);
    end
    med_corr_unc{1} = median(corr_unc(grch1))
    med_NSE_unc{1} = median(NSE_unc(grch1))
    med_rRMSE_unc{1} =median(rRMSE_unc(grch1))
     med_NSE_obs_daily{1} = median(NSE_obs_daily(grch1))
     med_rRMSE_obs_daily{1} =median(rRMSE_obs_daily(grch1))
     med_corr_obs{1} = median(corr_obs(grch1))
    med_NSE_obs{1} = median(NSE_obs(grch1))
    med_rRMSE_obs{1} =median(rRMSE_obs(grch1))
%     med_NSE_obs_daily{1} = median(NSE_obs_daily(grch1))
%     med_rRMSE_obs_daily{1} =median(rRMSE_obs_daily(grch1))
    %     reserve = reserve(grch1);
    %     delete = delete(grch1);
    
    W1 = med_W_obs;
    S1 = med_S_obs;
    [c,D,tau] = build_cDtau((Q),A,W1,S1,QWBM,grch1,W,S,nR,Q_withunc);
    dA_cell = cell(nR,nt);
    dA_cell_full = cell(nR,nt);
    dA_unc_cell = cell(nR,nt);
    ns = 0;
    for i = 1:nR
        if ~isnan(find(i==grch1))
            ns = ns+1;
            for j = 1:length(reserve{i})
                dA_cell{i,reserve{i}(j)}=med_dA_obs(ns,reserve{i}(j));
                dA_unc_cell{i,reserve{i}(j)}=sigma_dA_obs(ns,reserve{i}(j));
            end
%                         else
%                             for j = 1:length(reserve{i})
%                                 dA_cell{i,reserve{i}(j)}=dA_sg(i,reserve{i}(j));
% %                                dA_unc_cell{i,reserve{i}(j)}=max(sigma_dA_obs,[],'all');
% %                                 dA_unc_cell{i,reserve{i}(j)}=max(dA_unc_sg,[],'all');
%                                  dA_unc_cell{i,reserve{i}(j)}=dA_unc_sg(i,reserve{i}(j));
%                             end
        end
        for j = 1:length(reserve{i})
            dA_cell_full{i,reserve{i}(j)}=1;
        end
    end
    
    dt = 86400;
    dif_meas = [1,4,8,12,16,20];
    
    
    
    %% choice
    % 1. col -> Process noise: calculated=1,diag(Q)=2,diag(Q_end)=3,(Q=0)=3
    % 2. col -> cov. of obs.: R_res=1,R_err(SIGMA)=2,R_err(dA)=3,R_res_fit = 4
    %                         R_res_fit_med = 5, R_res_fit_mean = 6
    choice = [1,3;1,2;1,1;
        2,3;2,2;2,1;
        3,3;3,2;3,1;
        4,3;4,2;4,1;];
    %         choice = [3,1;3,4;3,5;3,6];
    nR = size(Q,1);
    for cho =1%size(choice,1)
        
        R_choice = choice(cho,2);
        %% compare z_cal and z_true
        if_plot = 2; % yes=1,no=2
        % 2. col of R_choice -> cov. of obs.:
        %  R_res=1,R_err(SIGMA)=2,R_err(dA)=3,R_res_fit = 4
        %  R_res_fit_med = 5, R_res_fit_mean = 6
        if R_choice ==3
            var_R ={};
        end
        
        %[var_R,fit_std,fit_std_med,fit_std_mean,z_true,z_cal] = compare_z_cal_true(nR,nt,Q,dif_meas,dA,dA_unc,dt,if_plot,L,R_choice);
        
        if use_Q_rlz == 1
            Q_obs_mean = Q_obs_mean2;
        else
            Q_obs_mean = Q_obs_mean1;
        end
        %% Qprocess percentage
        %         pct = [10;100;1e4;1e7];
        pct = [10;30;40;50;60;70];
        %% KF
        for k = 1%1:length(pct) % percentage for process noise
            state_ep = 22;
            %% begin values
            xn1n1 = [];
            for i = 1:state_ep
               % xn1n1 = [xn1n1;Q(:,i)-mean(Q_withunc,2)];
                    xn1n1 = [xn1n1;Q_withunc(:,i)-Q_obs_mean2];%Q_obs_mean2];%mean(Q_withunc,2)];;%
            end
              %   xn1n1 = [Q(:,1)];
            if use_Q_rlz == 1
                 Pn1n1 = diag(repmat(std(Q_withunc,1,2).^2,state_ep,1));
                %  Pn1n1 = diag(repmat(std(Q,1,2).^2,state_ep,1));
%                   tmp = [];
%                   for i = 1:state_ep
%                       tmp =  [tmp;sigma_Q_withunc(:,i)];
%                       Pn1n1 = diag(tmp.^2);
%                   end
            else
%                 tmp = [];
%                 for i = 1:state_ep
%                                         tmp =  [tmp;(Q_withunc(:,i)-Q(:,i)).^2];
%                    % tmp =  [tmp;(0.3*Q(:,i)).^2];
% %                     tmp =  [tmp;sigma_Q_withunc(:,i)];
%                     Pn1n1 = diag(tmp);
%                 end
             %                    Pn1n1 = diag(repmat(std(Q,1,2).^2,state_ep,1));
              %            Pn1n1 = diag(repmat((0.3*Q(:,1)).^2,state_ep,1));
               Pn1n1 = diag(repmat((Q_withunc(:,1)-Q(:,1)).^2,state_ep,1));
              % Pn1n1 = diag(repmat(std(Q,1,2).^2,state_ep,1));
            end
            
            %%  1. iteration
            i=1;
            
            % build transition matrix Phi_st and process noise
            choice_Q = choice(cho,1);
            [Phi,Q_st] = build_Phi_SWOT(nR,L_center_pos,c,D,tau,Q,pct(k)*0.01,state_ep,choice_Q);
            
            % Berechnungen vorab:
            xnn1(:,i) =  Phi *  xn1n1;
            Pnn1{i} = (Phi * Pn1n1 * Phi') + Q_st;
            
            [H,zn,R,z_idx,idx] = build_H_obs_SWOT(L(grch),dA_cell,dA_unc_cell,i,Q,state_ep);
            %                 dt_samp = [dt_samp;idx(:,2)-idx(:,1)];
            if ~isempty(zn)
                if R_choice ~=3
                    dif = idx(:,2)-idx(:,1);
                    tmp = [];
                    for j = 1:length(z_idx)
                        tmp = [tmp;var_R{dif(j)}(z_idx(j))];
                    end
                    R = diag(tmp);
                end
                %     R = diag(rms(z_cal-z_true,2).^2);
                if use_Q_obs == 1
                    
                    [H_Q,z_Q,R_Q] = build_H_obs_SWOT_Q(L(grch),dA_cell,i,Q,state_ep,Q_withunc,Q_res_unc,R_choice,use_Q_rlz,sigma_Q_withunc,Q_obs_mean);
                    H= [H;H_Q];
                    zn = [zn;z_Q];
                    size_R = size(R,1);
                    size_RQ = size(R_Q,1);
                    R = [R,zeros(size_R,size_RQ);zeros(size_RQ,size_R),R_Q];
                end
                % Berechnung Kalman Gain:
                Kn = Pnn1{i} * H' * pinv(H * Pnn1{i} * H' + R);
                
                
                % Update durch Beobachtung:
                xnn{i} = xnn1(:,i) + Kn*(zn - H * xnn1(:,i));
                
                % Berechnung der a-posteriori Kovarianz:
                Pnn{i} = (eye(state_ep*nR) - Kn * H) * Pnn1{i};
                
            else
                xnn{i} = xnn1(:,i);
                Pnn{i} = Pnn1{i};
            end
            
            %% other iterations
            for i = 2:nt-state_ep
                %% build transition matrix Phi_st and process noise
                xnn1(:,i) =  Phi *  xnn{i-1};
                Pnn1{i} = (Phi * Pnn{i-1} * Phi') + Q_st;
                
                [H,zn,R,z_idx,idx] = build_H_obs_SWOT(L(grch),dA_cell,dA_unc_cell,i,Q,state_ep);
                %                     dt_samp = [dt_samp;idx(:,2)-idx(:,1)];
                
                if ~isempty(zn)
                    if R_choice ~=3
                        
                        dif = idx(:,2)-idx(:,1);
                        tmp = [];
                        for j = 1:length(z_idx)
                            tmp = [tmp;var_R{dif(j)}(z_idx(j))];
                        end
                        R = diag(tmp);
                    end
                    %         R = diag(rms(z_cal-z_true,2).^2);
                    if use_Q_obs == 1
                        [H_Q,z_Q,R_Q] = build_H_obs_SWOT_Q(L(grch),dA_cell,i,Q,state_ep,Q_withunc,Q_res_unc,R_choice,use_Q_rlz,sigma_Q_withunc,Q_withunc);
                        H= [H;H_Q];
                        zn = [zn;z_Q];
                        size_R = size(R,1);
                        size_RQ = size(R_Q,1);
                        R = [R,zeros(size_R,size_RQ);zeros(size_RQ,size_R),R_Q];
                    end
                    % Berechnung Kalman Gain:
                    Kn = Pnn1{i} * H' * pinv(H * Pnn1{i} * H' + R);
                    %                         figure,plot(sqrt(diag(H * Pnn1{i} * H')))
                    %                         hold on
                    %                         plot(sqrt(diag(R)))
                    %                         legend({'HpH','R'})
                    % Update durch Beobachtung:
                    xnn{i} = xnn1(:,i) + Kn*(zn - H * xnn1(:,i));
                    
                    % Berechnung der a-posteriori Kovarianz:
                    Pnn{i} = (eye(state_ep*nR) - Kn * H) * Pnn1{i};
           %         schurma = Phi-Kn*H*Phi;
                else
                    xnn{i} = xnn1(:,i);
                    Pnn{i} = Pnn1{i};
                end
                
            end
            
            %% combine xnn
            [Q_est,Qest_weight,Qest_arith,Qest_med] = combine_xnn_SWOT(xnn,Pnn,nR,nt,state_ep,Q_obs_mean2);
            
            %% plot
                           plot_result(Q_est,Qest_weight,Qest_arith,Qest_med,nR,nt,Q,state_ep,orbit,reach,reserve,Q_withunc)
            %
            
            [corr_est,~,rRMSE_est,NSE_est,med_corr_est,med_NSE_est,med_rRMSE_est]=validate(Q_est,Q,nR,nt);
            [corr_weight,~,rRMSE_weight,NSE_weight,med_corr_weight,med_NSE_weight,med_rRMSE_weight]=validate(Qest_weight,Q,nR,nt);
            [corr_arith,~,rRMSE_arith,NSE_arith,med_corr_arith,med_NSE_arith,med_rRMSE_arith]=validate(Qest_arith,Q,nR,nt);
            [corr_med,~,rRMSE_med,NSE_med,med_corr_med,med_NSE_med,med_rRMSE_med]=validate(Qest_med,Q,nR,nt);
            xnn={};
            Pnn={};
            Pnn1={};
            xnn1=[];
        end
        %% combine validation results
        %             Q_with_unc{1} = Q_withunc(:,2:end);
        %             [~,~,~,~,med_corr_unc,med_NSE_unc,med_rRMSE_unc]=validate(Q_with_unc,Q,nR,nt)
        med_corr{cho} = [cell2mat(med_corr_est);cell2mat(med_corr_weight);cell2mat(med_corr_arith);cell2mat(med_corr_med)];
        med_NSE{cho} = [cell2mat(med_NSE_est);cell2mat(med_NSE_weight);cell2mat(med_NSE_arith);cell2mat(med_NSE_med)];
        med_rRMSE{cho} = [cell2mat(med_rRMSE_est);cell2mat(med_rRMSE_weight);cell2mat(med_rRMSE_arith);cell2mat(med_rRMSE_med)];
        
        median(NSE_med{1}(grch1))
        median(rRMSE_med{1}(grch1))

        Qest_syn = [Q_est,Qest_weight{1},Qest_arith{1},Qest_med{1}];
        corr_syn = [corr_est,corr_weight{1},corr_arith{1},corr_med{1}];
        NSE_syn = [NSE_est,NSE_weight{1},NSE_arith{1},NSE_med{1}];
        rRMSE_syn = [rRMSE_est,rRMSE_weight{1},rRMSE_arith{1},rRMSE_med{1}];
        %             for n = 1:state_ep-1
        %                 idx_dt = find(n==dt_samp);
        %                 sum_dt(n) = length(idx_dt);
        %
        %             end
        %             figure,plot(cumsum(sum_dt./xx),'-*')
        %             ylim([0,1])
        %             xlim([0,22])
        %             xticks(1:state_ep-1)
        %             grid on
        %             set(gca,'FontSize',15)
        %             xlabel('\Deltat')
        %             ylabel('CDF')
        %             figure,plot((sum_dt./xx),'-s')
        %             xlim([0,22])
        %             xticks(1:state_ep-1)
        %             grid on
        %             set(gca,'FontSize',15)
        %             xlabel('\Deltat')
        %             ylabel('CDF')
    end
    
    
    %         clearvars -except Q Q_withunc percent_dl delete_reach med_corr med_NSE med_rRMSE file_KF state_ep
    %%
    plot_validation9(state_ep,med_corr,med_NSE,med_rRMSE,med_corr_unc,med_NSE_unc,med_rRMSE_unc,med_corr_obs,med_NSE_obs,med_rRMSE_obs)
    
    %%
    
   
             mm=5;
%                 save(['corr_syn_tau',num2str(mm),'.mat'],'corr_syn')
%         save(['NSE_syn_tau',num2str(mm),'.mat'],'NSE_syn')
%         save(['rRMSE_syn_tau',num2str(mm),'.mat'],'rRMSE_syn')
    save(['corr_unc',num2str(m),'.mat'],'corr_unc')
    save(['NSE_unc',num2str(m),'.mat'],'NSE_unc')
    save(['rRMSE_unc',num2str(m),'.mat'],'rRMSE_unc')
    save(['corr_obs',num2str(m),'.mat'],'corr_obs')
    save(['NSE_obs',num2str(m),'.mat'],'NSE_obs')
    save(['rRMSE_obs',num2str(m),'.mat'],'rRMSE_obs')
        save(['corr_unc_daily',num2str(m),'.mat'],'corr_unc_daily')
    save(['NSE_unc_daily',num2str(m),'.mat'],'NSE_unc_daily')
    save(['rRMSE_unc_daily',num2str(m),'.mat'],'rRMSE_unc_daily')
    save(['corr_obs_daily',num2str(m),'.mat'],'corr_obs_daily')
    save(['NSE_obs_daily',num2str(m),'.mat'],'NSE_obs_daily')
    save(['rRMSE_obs_daily',num2str(m),'.mat'],'rRMSE_obs_daily')
     if use_Q_rlz == 1
         save(['Qest_Qrlz',num2str(m),'.mat'],'Qest_syn') 
                 save(['corr_Qrlz',num2str(m),'.mat'],'corr_syn')
        save(['NSE_Qrlz',num2str(m),'.mat'],'NSE_syn')
        save(['rRMSE_Qrlz',num2str(m),'.mat'],'rRMSE_syn')
     else
         save(['Qest_syn',num2str(m),'.mat'],'Qest_syn')  
        save(['corr_syn',num2str(m),'.mat'],'corr_syn')
        save(['NSE_syn',num2str(m),'.mat'],'NSE_syn')
        save(['rRMSE_syn',num2str(m),'.mat'],'rRMSE_syn')

     end
    Qest_rlz=Qest_syn;
end
% %% plot median validation values
% state_ep = 22;
% % a = [];
% % for i = 1:9
% %    [ a(i)] = max(med_corr{i}(7:9));
% %  [b(i)] = max(med_NSE{i}(7:9));
% %  [c(i)] = min(med_rRMSE{i}(7:9));
% % end
% % [h,d]=max(a)
% % [hh,e]=max(b)
% % [hhh,f]=min(c)


