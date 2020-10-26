clear all; clc;
close all;




SP.Nr = 16;  Nr = SP.Nr;     % # Rx Antennas per BS
SP.Nc = 2;  Nc = SP.Nc;      % # Cells
SP.Nu = 2;  Nu = SP.Nu;      % # Users per cell

SP.pu_dBm = 0;
SP.pu = 10.^(SP.pu_dBm/10); % mW


dB_min = -5;
dB_max = 2;

SP.gamma = linspace(10^(dB_min/10),  10^(dB_max/10), 5);

SP.gamma_dB = 10*log10(SP.gamma);

SP.th = 1e-4;              % Iteration Threshold
SP.initpower  = 100;
initlambda = SP.initpower*ones(Nu,Nc);
SP.algoMax = 400;
SP.algoMax2 = 300;
SP.iterMax = 5;

% Channel Parameters
SP.BW = 10*10^6;
SP.NF = 5;
SP.dmax = 1000;            % In meters
SP.dmin = 100;
SP.d0 = 100;
SP.rho = 3;
SP.s = 8.7;
SP.LD = 4*10^8/(2.4*10^9);
SP.a = 0.8; % First-order Gauss-Markov model for small scale fading
%%

ADC = [2 3 4 inf];

P_joint_result = zeros(length(SP.gamma), length(ADC));
SINR_joint_result = zeros(length(SP.gamma), length(ADC));
Rate_infinite_result = zeros(length(SP.gamma), length(ADC));
count_infinite_result = zeros(length(SP.gamma), length(ADC));

P_percell_result = zeros(length(SP.gamma), length(ADC));
SINR_percell_result = zeros(length(SP.gamma), length(ADC));
Rate_percell_result = zeros(length(SP.gamma), length(ADC));
count_percell_result = zeros(length(SP.gamma), length(ADC));


%%
rng(16)
for iter = 1:SP.iterMax
    coeff = hexcell(SP);
    Hcurr = sqrt(1/2)*(randn(Nr, Nu*Nc, Nc) + 1j*randn(Nr, Nu*Nc, Nc));
    H = zeros(Nr, Nu*Nc, Nc);
    for c = 1:Nc
        H(:,:,c) = Hcurr(:,:,c)*diag(sqrt(coeff(:,c)));
    end
    
    for n = 1:length(SP.gamma)   



        for i = 1:length(ADC)  
            SP.b = ADC(i);

            [P_ji, SINR_ji, count_ji, lambda_ji] = CoMP(SP, H, n, 'joint', initlambda);
            [P_pi, SINR_pi, count_pi, lambda_pi] = CoMP(SP, H, n, 'percell', initlambda);

            % Running average

            P_joint_result(n, i) = P_joint_result(n, i) + (P_ji - P_joint_result(n, i))/iter;
            SINR_joint_result(n, i) = SINR_joint_result(n, i)+ (mean(SINR_ji(:)) - SINR_joint_result(n, i))/iter;
            Rate_infinite_result(n, i) = Rate_infinite_result(n, i) + (mean(log2(1+10.^(SINR_ji(:)/10))) - Rate_infinite_result(n, i))/iter;
            count_infinite_result(n, i) = count_infinite_result(n, i) + (count_ji - count_infinite_result(n, i))/iter;


            P_percell_result(n, i) = P_percell_result(n, i) + (P_pi - P_percell_result(n, i))/iter;
            SINR_percell_result(n, i) = SINR_percell_result(n, i)+ (mean(SINR_pi(:)) - SINR_percell_result(n, i))/iter;
            Rate_percell_result(n, i) = Rate_percell_result(n, i) + (mean(log2(1+10.^(SINR_pi(:)/10))) - Rate_percell_result(n, i))/iter;
            count_percell_result(n, i) = count_percell_result(n, i) + (count_pi - count_percell_result(n, i))/iter;

        end    
    end
end
%%
figure
hold on

plot(SP.gamma_dB, 10*log10(P_joint_result(:,1)), 'bo:')
plot(SP.gamma_dB, 10*log10(P_joint_result(:,2)), 'bo-.')
plot(SP.gamma_dB, 10*log10(P_joint_result(:,3)), 'bo--')
plot(SP.gamma_dB, 10*log10(P_joint_result(:,4)), 'bo-')

plot(SP.gamma_dB, 10*log10(P_percell_result(:,1)), 'ro:')
plot(SP.gamma_dB, 10*log10(P_percell_result(:,2)), 'ro-.')
plot(SP.gamma_dB, 10*log10(P_percell_result(:,3)), 'ro--')
plot(SP.gamma_dB, 10*log10(P_percell_result(:,4)), 'ro-')

xlabel('Target SINR [dB]')
ylabel('Total Transmit Power [dBm]')

legend('Q-iCoMP ($b=2$)','Q-iCoMP ($b=3$) ','Q-iCoMP ($b=4$)','Q-iCoMP ($b=\infty$)',...
    'Q-Percell ($b=2$)','Q-Percell ($b=3$) ','Q-Percell ($b=4$)','Q-Percell ($b=\infty$)',...
    'Interpreter','latex')
grid on


    

