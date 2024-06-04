
function [Q_KF,Q_KF_var] = build_weighted_arith(xnn,Pnn,nR,nt,state_ep,mean_func)
% combines Q and its variance of the same month
Qest{1} = xnn{1,1}(1:nR);
Qest_cov{1} = Pnn{1,1}(1:nR,1:nR);
Qest{nt-1} = xnn{1,end}(end-(nR-1):end);
Qest_cov{nt-1} = Pnn{1,end}(end-(nR-1):end,end-(nR-1):end);

% combines Q and its variance of the same month 2-21
for i = 2:state_ep-1
    A = repmat(eye(nR),i,1);
    y = [];
    tmp = zeros(i*nR,i*nR);
    for j=1:i
        y = [y ;xnn{1,j}((i-j)*nR+1:(i-j)*nR+nR)];
        tmp((j-1)*nR+1:(j-1)*nR+nR,(j-1)*nR+1:(j-1)*nR+nR) = Pnn{1,j}((i-j)*nR+1:(i-j)*nR+nR,(i-j)*nR+1:(i-j)*nR+nR);
        
    end
    if mean_func == 1
        P = pinv(tmp);
    else
        P = eye(i*nR,i*nR);
    end
    x_dach = pinv(A'*P*A)*A'*P*y;
    e_dach = y-A*x_dach;
    sigma0 = sqrt((e_dach'*P*e_dach)/(size(A,1)-size(A,2)));
    SIGMA_dach_x_dach = sigma0^2*pinv(A'*P*A);
    
    Qest{i} = x_dach;
    Qest_cov{i} = SIGMA_dach_x_dach;
    
end


% combines Q and its variance of the same month end-21:end-2
for i = 2:state_ep-1
    A = repmat(eye(nR),i,1);
    y = [];
    tmp = zeros(i*nR,i*nR);
    for j=1:i
        y = [y ;xnn{1,end-i+j}((state_ep-j)*nR+1:(state_ep-j)*nR+nR)];
        tmp((j-1)*nR+1:(j-1)*nR+nR,(j-1)*nR+1:(j-1)*nR+nR) = ...
            Pnn{1,end-i+j}((state_ep-j)*nR+1:(state_ep-j)*nR+nR,(state_ep-j)*nR+1:(state_ep-j)*nR+nR);
        
    end
    if mean_func == 1
        P = pinv(tmp);
    else
        P = eye(i*nR,i*nR);
    end
    x_dach = pinv(A'*P*A)*A'*P*y;
    e_dach = y-A*x_dach;
    sigma0 = sqrt((e_dach'*P*e_dach)/(size(A,1)-size(A,2)));
    SIGMA_dach_x_dach = sigma0^2*pinv(A'*P*A);
    
    Qest{nt-i} = x_dach;
    Qest_cov{nt-i} = SIGMA_dach_x_dach;
    
end


% combines Q and its variance of the same month
A = repmat(eye(nR),state_ep,1);

for i = 1:length(xnn)-(state_ep-1)
    y = [];
    tmp = zeros(state_ep*nR,state_ep*nR);
    for j=1:state_ep
        y = [y ;xnn{1,i+j-1}((state_ep-j)*nR+1:(state_ep-j)*nR+nR)];
        tmp((j-1)*nR+1:(j-1)*nR+nR,(j-1)*nR+1:(j-1)*nR+nR) = Pnn{1,i+j-1}((state_ep-j)*...
            nR+1:(state_ep-j)*nR+nR,(state_ep-j)*nR+1:(state_ep-j)*nR+nR);
        
    end
    if mean_func == 1
%          P = eye(state_ep*nR,state_ep*nR);
        P = pinv(tmp);
    else
        P = eye(state_ep*nR,state_ep*nR);
    end
    
    x_dach = pinv(A'*P*A)*A'*P*y;
    e_dach = y-A*x_dach;
    sigma0 = sqrt((e_dach'*P*e_dach)/(size(A,1)-size(A,2)));
    SIGMA_dach_x_dach = sigma0^2*pinv(A'*P*A);
    
    Qest{i+state_ep-1} = x_dach;
    Qest_cov{i+state_ep-1} = SIGMA_dach_x_dach;
    
end

for j = 1:nR
    tmp1 =  [];
    tmp2 =  [];
    for i = 1:nt-1
        tmp1 = [ tmp1,Qest{i}(j)];
    end
    Q_KF{1}(j,:) = tmp1;
    
    for i = 1:nt-1
        tmp2 =[tmp2;Qest_cov{1,i}(j,j)];
    end
    Q_KF_var{1}(j,:)  = tmp2;
end

                