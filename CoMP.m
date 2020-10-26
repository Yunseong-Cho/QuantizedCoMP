function [P_Tc, SINR, count_Tc, lambda, PP] = CoMP( SP, H, n, method, initlambda)


% First, computation of allocated power
% Second, computation of SINR

% - Input
%     SP: System Parameters
%     H: Channel
%     n: index of target SINR array
%     method: specifies how to calculate the powers ('joint','percell','det')
%     res: Resolution ('infinite', 'low')
%
% - Output
%     P_Tc: Average power
%     SINR: SINR of each user in each cell
%     count_Tc: Average number of iterations to converge

gamma = SP.gamma(n); % targe SINR in linear scale

switch method
    case 'joint'
        [lambda, count,SINR] = Algo_joint(SP, H, gamma, initlambda);
    case 'percell'
        [lambda, count,SINR] = Algo_percell_low(SP, H, gamma, initlambda);            
end

    P_it = sum(lambda(:));
    


SINR = abs(SINR);
P_Tc = mean(P_it);
PP = lambda(:);
count_Tc = nanmean(count);

end


