function [ lambda, cnt, SINR ] = Algo_percell_low(SP, H, gamma, initlambda)

b = SP.b;
if b == inf
    alpha = 1;
else
    bTable = [0.3634, 0.1175, 0.03454, 0.009497, 0.002499];
    if b > 5
        alpha = (1 - pi*sqrt(3)/2*2.^(-2*b));
    else
        alpha = (1 - bTable(floor(b)));
    end
end

Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;
th = SP.th;
algoMax = SP.algoMax;
algoMax2 = SP.algoMax2;
lambda_c = initlambda; % initial guess
lambda_p = zeros(Nu,Nc);

cnt = 0;
count = 0;
count2 = 0;
while (sum( abs((lambda_c(:) - lambda_p(:))./lambda_c(:)) > th) ~= 0 && count2 < algoMax2) % all cells  % relative error, was abs(lambda_c - lambda_p)
    
    lambda_p = lambda_c;
    for i = 1:Nc
        
        % Noise estimate
        lambda_out = lambda_p(:);
        lambda_out(Nu*(i-1)+1:Nu*(i-1)+Nu) = 0; % Users in current cell
        noise = diag(alpha^2*H(:,:,i)*diag(lambda_out)*H(:,:,i)') + alpha*eye(Nr); % was alpha^2 * ICI of other cells
        noise_q = alpha*(1-alpha)*diag(H(:,:,i)*diag(lambda_p(:))*H(:,:,i)');
        
        lambda_cell_c = lambda_c(:,i); 
        lambda_cell_p = zeros(Nu,1);
        while (sum( abs((lambda_cell_c - lambda_cell_p)./lambda_cell_c) > th) ~= 0 && count < algoMax )
            lambda_cell_p = lambda_cell_c;
            
            lambda = lambda_p; % the values generated before getting into the inner loop
            lambda(:,i) = lambda_cell_p; % Replace the value of corresponding cell only
            Lambda_percell = diag(lambda(:));
            
            ICI = H(:,:,i)*Lambda_percell*H(:,:,i)';
            H_incell = H(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
            lambda_incell = diag(lambda_cell_p);
            
            W = pinv(alpha^2*H_incell*lambda_incell*H_incell' + diag(noise + noise_q))*H_incell; % incell MMSE
            %             W = pinv(H_incell*lambda_incell*H_incell')*H_incell; %ZF
            %             W = pinv(ICI+eye(Nr))*H_incell; % MMSE
            for u = 1:Nu
                SigPow = abs( alpha*W(:,u)'*H_incell(:,u) )^2;
                NoisePow = W(:,u)'* (alpha*eye(Nr) ...
                    + alpha^2*ICI ...
                    - alpha^2*H_incell(:,u)*lambda_incell(u,u)*H_incell(:,u)'...
                    + alpha*(1-alpha)*diag(diag(H(:,:,i)*Lambda_percell*H(:,:,i)'))) *W(:,u);
                lambda_cell_c(u) = gamma*real(NoisePow/SigPow);
            end
            cnt = cnt+1;

        end
        lambda_c(:,i) = lambda_cell_c;
        count = count + 1;
    end
    count2 = count2 + 1;
end
lambda = lambda_c;

% Compute SINR
SINR = zeros(Nu,Nc);
for i = 1:Nc
    
    % Noise estimate
    lambda_out = lambda(:);
    lambda_out(Nu*(i-1)+1:Nu*(i-1)+Nu) = 0; % Users in current cell
    noise = diag(alpha^2*H(:,:,i)*diag(lambda_out)*H(:,:,i)') + alpha*eye(Nr);
    noise_q = alpha*(1-alpha)*diag(H(:,:,i)*diag(lambda(:))*H(:,:,i)');
    
    % SINR
    H_incell = H(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
    lambda_incell = diag(lambda(:,i));
    W = pinv(alpha^2*H_incell*lambda_incell*H_incell' + diag(noise + noise_q))*H_incell; % incell MMSE
    ICI = H(:,:,i)*diag(lambda(:))*H(:,:,i)';
    for u = 1:Nu
        SigPow = lambda_incell(u,u)*abs( alpha*W(:,u)'*H_incell(:,u) )^2;
        NoisePow = W(:,u)'* (alpha*eye(Nr) ...
            + alpha^2*ICI - alpha^2*H_incell(:,u)*lambda_incell(u,u)*H_incell(:,u)'...
            + alpha*(1-alpha)*diag(diag(H(:,:,i)*diag(lambda(:))*H(:,:,i)'))) *W(:,u);
        SINR(u,i) = SigPow/NoisePow;
    end
end



