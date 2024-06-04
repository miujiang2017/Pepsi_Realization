function [c,D,tau] = build_cDtau(Q,A,W,S,QWBM,grch1,W_z,S_z,nR,Q_withunc)
%QWBM = median(Q,'all');
%% hydraulic diffusivity: hydraulic diffusivity D0 was calculated as
%% the median of the values computed for in all days using equation:
%% D = Q/(2WS)
%     D = abs(mean((abs(Q)./(2*W.*(S))),2));
%     D = abs(QWBM./(2*median(abs(W),2).*median(abs(S),2)));
%  D = abs(median(abs(Q),2)./(2*median(W,2).*median((S),2)));
%     D = abs(QWBM./(2*mean(abs(W),2).*mean((S),2)));
%     D = abs(mean((Q),2)./(2*mean((W),2).*mean((S),2)));
% QWBM=median(Q,'all');

idx1 = setdiff(1:nR,grch1);
ns = 0;
for i = 1:nR
    if isempty(find(idx1==i))
        ns = ns+1;
        idx = find(S(ns,:)~=0);
          % D(i,1) = abs(mean(Q(ns,idx))./(2*mean(W(ns,idx)).*mean(abs(S(ns,idx)))));
            %  D(i,1) = abs(median(Q(ns,idx))./(2*median(W(ns,idx)).*median(abs(S(ns,idx)))));
        %        D(i,1) = abs(QWBM./(2*mean(W(ns,idx)).*mean((S(ns,idx)))));
        %        D(i,1) = abs(QWBM./(2*(W_z(i)).*abs(S_z(i))));
        D(i,1) = abs(mean(QWBM./(2*W(ns,idx).*(S(ns,idx)))));
      %  D(i,1) = abs(QWBM./(2*mean(W(ns,idx)).*mean(abs(S(ns,idx)))));
       %   D(i,1) = abs(median(QWBM./(2*W(ns,idx).*(S(ns,idx)))));
    else
        idx = find(S_z(i,:)~=0);
       % D(i,1) = abs(QWBM./(2*mean(W_z(i,idx)).*mean(abs(S_z(i,idx)))));
        D(i,1) = abs(mean(QWBM./(2*W_z(i,idx).*(S_z(i,idx)))));
    end
    % D(i,1) = abs(QWBM./(2*median((W(i,idx))).*median(abs(S(i,idx)))));
    % D(i,1) = abs(QWBM./(2*mean(W(i,idx)).*mean((S(i,idx)))));
    %         D(i,1) = abs(median((Q(i,idx))./(2*W(i,idx).*(S(i,idx)))));
end


%             figure,
%         for i = 1:size(S,1)
%             subplot(ceil(size(S,1)/2),2,i)
%             plot((Q(i,:))./(2*W(i,:).*(S(i,:))),'linewidth',1.2)
%             hold on
%             yline(median((Q(i,:))./(2*W(i,:).*(S(i,:)))),'r')
%             hold on
%             yline(mean((Q(i,:))./(2*W(i,:).*(S(i,:)))),'g')
%             %             plot(abs((Strue(i,:)-S(i,:))./Strue(i,:))*100)
%             name = "Reach "+num2str(i);
%             x = mean((Q(i,:))./(2*W(i,:).*(S(i,:))));
%             ylabel('[m/m]')
%             xlabel('day')
%             ylim([x-1e5,x+1e5])
%             title(name)
%             grid on
%             set(gca,'FontSize',15)
%         end

%% diffussive wave celerity and decorrelation length
%% diffusive wave celerity c0 was estimated by fitting a potential
%% equation Q=a(A-A0)^b to A and Q values (a, b, and A0 are parameters)
%% and taking c0 as the mean derivative of this function
c = [];
c_rrmse = [];
tau = [];
for i = 1:size(A,1)
     y =Q(i,:)';
    x = A(i,:)';
    
    %% diffussive wave celerity
    %                 % begin value
    %                 A0_0 = min(x,[],'all')-0.1; %Cannot fit Power functions to data where X has nonpositive values.
    %                 input_x = x-A0_0;
    %                 [fitresult, gof] = createFit_power(input_x, y);
    %                 a0 = fitresult.a;
    %                 b0 = fitresult.b;
    %                 a=[a0 A0_0 b0];
    %
    %
    %                 % fitting
    %                 % f=@(a,x)a(1)*sin(x).*exp(x)-a(2)./log(x);
    %                 f1=@(coef,inp)coef(1).*(inp-coef(2)).^coef(3);
    %                 % lsqcurvefit (resnorm: quadrate sum of residuals)
    %                 [unk,resnorm]=lsqcurvefit(f1,a,x,y,[0,0,0]);
    %                 c_rrmse(i,1) = sqrt(resnorm/sum(y.^2));
    %
    %                 % c = dQ/dA = a*b*(A-A_0)^(b-1)
    %                 c(i,1) = mean(unk(1).*unk(3).*(x-unk(2)).^(unk(3)-1));
    
    [fit_dir, gof] = createFit_direct(x, y, min(x));
    %         c = dQ/dA = a*b*(A-A_0)^(b-1);
    c(i,1) = mean(fit_dir.a.*fit_dir.c.*(x-fit_dir.b).^(fit_dir.c-1));
    
    %         if gof.sse<resnorm
    %             'yes'
    %         end
    
    
    %% decorrelation length: estimated by fitting r=e^(-abs(dt)/tau) to
    %% the empirical autocorrelation function estimated using SWOT
    %% discharge estimates for each reach.
     y =Q(i,:)';
    [cval,lags] = xcorr(y,'normalized');
    corr_y = cval(length(x)-1:-1:1);
    dt = [1:length(x)-1]';
    % %
    %         obs = log(corr_y);
    %         A_mat = -dt;
    %         tau(i,1) = 1./(inv(A_mat'*A_mat)*A_mat'*obs)*86400;
    %         tau(i,1) = -sum(dt.^2)/sum(log(corr_y).*dt)  *86400;
    tau(i,1) = -sum(log(corr_y).*dt)/sum(log(corr_y).^2)  *86400;
    %         m =  tau(i,1)/86400;
    %         a_start = 4*D(i)/c(i)^2/86400;
    %         a_start = tau(i,1)/86400;
    %         [fit_exp, gof] = createFit_exp(dt, corr_y, a_start);
    %         tau(i,1) = fit_exp.a*86400;
    %                         figure
    %                 plot(dt,exp(-dt./(tau(i,1)/86400)),'r.')
    %                 hold on
    %                 plot(dt,exp(-dt./m),'g.')
    %                 hold on
    %                 plot(dt,corr_y,'b')
    
    % obs = 1./log(corr_y);
    % A_mat = -1./dt;
    % x = inv(A_mat'*A_mat)*A_mat'*obs
    % x_dach1 = -sum(1./dt.*1./log(corr_y))./sum((1./dt).^2)
    
    %     x_dach2 = -sum(log(corr_y).*dt)/sum(log(corr_y).^2)
    %     x_dach3 = -sum(dt.^2)/sum(log(corr_y).*dt)
end
% %  c = abs(random('Normal',c,c*2));
%  D = abs(random('Normal',D,D*200));
% tau = abs(random('Normal',tau,tau*5));
end

