function [H,z,R] = build_H_obs_SWOT_Q(L,dA,ep,Q,state_ep,Q_withunc,Q_res_unc,R_choice,use_Q_rlz,sigma_Q_withunc,Q_obs_mean2)
nR = length(L);
H = zeros(1,state_ep*nR);
R = [];
z_idx = [];
% H = zeros(nR_dA,2*nR);
% R = zeros(nR_dA,nR_dA);
z = [];
emp = cellfun(@isempty,dA(:,ep+1:ep+state_ep));
tmp = reshape(emp,state_ep*nR,1);

idx = find(tmp ==0);
Q_dl = Q_withunc-mean(Q_withunc,2);
Q_tmp = reshape(Q_dl(:,ep+1:ep+state_ep),state_ep*nR,1);
if R_choice == 1
    Q_tmpunc = repmat(Q_res_unc,state_ep,1);
end
%     Q_tmpunc = repmat(Q_res_unc,state_ep,1);
if use_Q_rlz==1
  Q_tmpunc = reshape(sigma_Q_withunc(:,ep+1:ep+state_ep),state_ep*nR,1);  
else
    
Q_tmpunc = reshape(0.3*Q(:,ep+1:ep+state_ep),state_ep*nR,1);

end
for j = 1:length(idx)
        H(j,idx(j)) = 1;
        z = [z;Q_tmp(idx(j))];
        R = [R;Q_tmpunc(idx(j)).^2];
end
R = diag(R);
