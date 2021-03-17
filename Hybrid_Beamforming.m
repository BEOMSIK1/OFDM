clc, clear, close all;
%% Parameters
Nt = 64;                    % Tx antenna
Nr = 1;                     % Rx antenna per user
N_rf = 8;                   % RF chain
K = 8;                      % User
Modulation_Order = 4;
Symbol_Size = 128;
Data_Size = Symbol_Size*Modulation_Order;
P = -10:5:30;               % Tx power
noise_power = 1;            % 0dB
iteration = 1000;
%% Phase-Zero forcing(PZF)
for power_index=1:length(P)
    for iter=1:iteration
        %% Channel generate
        H = (randn(K,Nt)+j*randn(K,Nt))/sqrt(2);
        %% Data
        X=randi([0 1], [K, Data_Size]);                                    % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, Modulation_Order);                             % modulation
        %% Analog precoding
        F_rf = exp(j.*angle(H'));
        %% Digital precoding
        H_eq = H*F_rf;
        tx_power=1*(10^(P(power_index)/10));
        F_bb_hat = (H_eq'*inv(H_eq*H_eq'));
        gamma = trace(F_rf*F_bb_hat*F_bb_hat'*F_rf')./tx_power;
        F_bb = (H_eq'*inv(H_eq*H_eq'))./sqrt(gamma);
        %% Tx-Rx
        precoded_X = F_rf*F_bb*X_mod;                                                      % 사용자 신호 전처리
        Y_received = awgn_noise(H*precoded_X,0);                    % H*전처리된 각 사용자 신호/sqrt(gamma) + noise
        Y = Y_received.*sqrt(gamma);                                                       % 정규화 factor를 곱해주어 사용자 신호만 추출
        Y_demod_ZF = base_demod(Y,Modulation_Order);
        %% Error
        num_error(iter,power_index)=biterr(X,Y_demod_ZF);
        %% Sum rate
        Sum_rate(iter,power_index)=log2(det(eye(K)+(H*F_rf*F_bb*F_bb'*F_rf'*H')./noise_power));
    end
end
error_rate=(sum(num_error,1)/(Data_Size*K))/iteration;
sum_rate=(sum(abs(Sum_rate),1))/iteration;
%% graph
semilogy(P,error_rate,'-o')
title('BER Performance'), xlabel('Transmit Power(dB)'),ylabel('BER'),legend('PZF')
grid on

figure
plot(P,sum_rate,'-o')
title('Sum Rate'), xlabel('Transmit Power(dB)'),ylabel('Sum Rate(bps/Hz)'),legend('PZF')
grid on
