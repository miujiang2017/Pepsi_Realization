function [corr_comp_true,RMSE_comp_true,rRMSE_comp_true,NSE_comp_true,med_corr,med_NSE,med_rRMSE]=validate(Q_est,Q,nR,nt)
for i = 1:length(Q_est)
    for j = 1:nR
        Q_comp = Q_est{i}(j,:);
        Q_true = Q(j,2:end);
        corr_comp_true{i}(j,1) = corr(Q_comp',Q_true');
        RMSE_comp_true{i}(j,1) = sqrt(sum((Q_true-Q_comp).^2));
        rRMSE_comp_true{i}(j,1) = RMSE_comp_true{i}(j)./sqrt(sum(Q_true.^2))*100;
        NSE_comp_true{i}(j,1) = 1-sum((Q_true-Q_comp).^2)/sum((Q_true-mean(Q_true)).^2);
        
        %         xxx2(j) = corr(Q_comp',Q_rdm(j,2:end)');
        %         xxx3(j) = corr(Q(j,2:end)',Q_rdm(j,2:end)');
        
    end
    med_corr{i,1} =  median(corr_comp_true{i});
    med_NSE{i,1} = median(NSE_comp_true{i});
    med_RMSE{i,1} = median(RMSE_comp_true{i});
    med_rRMSE{i,1} =  median(rRMSE_comp_true{i});
%     med_corr{i,1} =  mean(corr_comp_true{i});
%     med_NSE{i,1} = mean(NSE_comp_true{i});
%     med_RMSE{i,1} = mean(RMSE_comp_true{i});
%     med_rRMSE{i,1} =  mean(rRMSE_comp_true{i});
end

end