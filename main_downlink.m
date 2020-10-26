clear all; clc;

% 2 & 2 plot
% Equation check -> extension to 4 & 3 snr_db = [-10 10]

SP.Kz_est = 'disable';

Nr_range = [16 32 64];     % # Rx Antennas per BS
SP.Nc = 3;  Nc = SP.Nc;      % # Cells
SP.Nu = 3;  Nu = SP.Nu;      % # Users per cell

SP.pu_dBm = 0;
SP.pu = 10.^(SP.pu_dBm/10); % mW


SP.gamma_dB = linspace(-5, 15, 20);
SP.gamma = 10.^(SP.gamma_dB/10);

SP.th = 1e-5;              % Iteration Threshold
SP.initpower  = 100;
SP.TcMax = 1;
SP.algoMax = 2000;
SP.iterMax = 10;


% Channel Parameters
SP.BW = 10*10^6;
SP.NF = 5;
SP.r_cell = 1000;            % In meters
SP.dmin = 100;
SP.d0 = 100;
SP.rho = 3;
SP.s = 8.7;
SP.LD = 4*10^8/(2.4*10^9);
SP.a = 0.9; % First-order Gauss-Markov model for small scale fading
%%

%%
ADC = [3, inf];

P_low_ul_result = zeros(length(SP.gamma), length(Nr_range), length(ADC));
P_low_dl_result = zeros(length(SP.gamma), length(Nr_range), length(ADC));

%%
rng(14)
f = waitbar(0);
now = 0;
for iter = 1:SP.iterMax 
    iter
    for nb = 1:length(Nr_range)    
        SP.Nr = Nr_range(nb); 
        Nr = SP.Nr;
        
        coeff = hexcell(SP);
        Hcurr = sqrt(1/2)*(randn(Nr, Nu*Nc, Nc) + 1j*randn(Nr, Nu*Nc, Nc));

        H = zeros(Nr, Nu*Nc, Nc);
        
        for c = 1:Nc
            H(:,:,c) = Hcurr(:,:,c)*diag(sqrt(coeff(:,c)));
        end

        for n = 1:length(SP.gamma)   
            
            init = SP.initpower*ones(Nu,Nc);
            

            for i = 1:length(ADC)  
                SP.b = ADC(i);
                [P_jl_ul, ~, ~, lambda_jl] = CoMP(SP, H, n, 'joint', 'low', init);
                [P_jl_dl] = Precoder( SP, H, n, 'low', lambda_jl);

                P_low_ul_result(n,nb,i) = P_low_ul_result(n,nb,i) + (P_jl_ul - P_low_ul_result(n,nb,i))/iter;
                P_low_dl_result(n,nb,i) = P_low_dl_result(n,nb,i) + (P_jl_dl - P_low_dl_result(n,nb,i))/iter;

            end
        end
    end
end

%%
figure
for nn = 1:length(Nr_range)

P = [P_low_ul_result(:,nn,:)];
P = reshape(P, [], length(ADC));
P_d = [P_low_dl_result(:,nn,:)];
P_d = reshape(P_d, [], length(ADC));

hold on
plot(SP.gamma_dB, 10*log10(P_d), 'o', 'MarkerSize', 8)
plot(SP.gamma_dB, 10*log10(P), '-', 'MarkerSize', 8)

end
grid on
legend('Downlink ($b=3$)','Downlink ($b=\infty$)', 'Uplink ($b=3$)','Uplink ($b=\infty$)', 'Interpreter','latex')
xlabel('Target SINR [dB]');
ylabel('Total Transmit Power [dBm]');
axis([-inf inf -inf inf])
title(['(Nr, Nc, Nu)=', num2str([Nc, Nu])])