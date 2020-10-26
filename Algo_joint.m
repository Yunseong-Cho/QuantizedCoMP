function [ lambda, count , SINR ] = Algo_joint(SP, H, gamma, initlambda, res)
% actual

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
lambda_c = initlambda; % initial guess
lambda_p = zeros(Nu,Nc);


count = 0;
% count2 = 1;
Kz = zeros(Nr,Nr);
while (sum( abs((lambda_c(:) - lambda_p(:))./lambda_c(:)) > th) ~= 0 && count < algoMax)   % relative error, was abs(lambda_c - lambda_p)
                                                                        %  && count < algoMax
    count = count + 1;
    
    lambda_p = lambda_c;
    Lambda = diag(lambda_p(:));
    
    %
    for i = 1:Nc
        

        ICI = H(:,:,i)*Lambda*H(:,:,i)';
        Kz = alpha*ICI ...
             + eye(Nr) ... 
             + (1-alpha)*diag(diag(H(:,:,i)*Lambda*H(:,:,i)')); % ICI + noise + ßQuantization noise
        
        
        for u = 1:Nu
            D = alpha*(1+1/gamma)*H(:,Nu*(i-1)+u,i)'*(Kz\H(:,Nu*(i-1)+u,i));
            lambda_c(u,i) = real(1/D);
        end
    end
    
    abs((lambda_c(:) - lambda_p(:)));
    
end
count;
lambda = lambda_c;

SINR = zeros(Nu,Nc);
for i = 1:Nc
    
    % SINR
    ICI = H(:,:,i)*diag(lambda(:))*H(:,:,i)';
    H_incell = H(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
    lambda_incell = diag(lambda(:,i));
    W = pinv(alpha^2*ICI + alpha*eye(Nr) + alpha*(1-alpha)*diag(diag(ICI)))*H_incell; % incell MMSE
    
    for u = 1:Nu
        SigPow = lambda_incell(u,u)*abs( alpha*W(:,u)'*H_incell(:,u) )^2;
        NoisePow = W(:,u)'* (alpha*eye(Nr) ...
            + alpha^2*ICI - alpha^2*H_incell(:,u)*lambda_incell(u,u)*H_incell(:,u)'...
            + alpha*(1-alpha)*diag(diag(H(:,:,i)*diag(lambda(:))*H(:,:,i)'))) *W(:,u);
        SINR(u,i) = SigPow/NoisePow;
    end
end
10*log10(mean(SINR(:)));
10*log(mean(lambda(:)));
end












