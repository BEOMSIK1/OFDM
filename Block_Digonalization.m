clc, clear, close all;
%% Parameters
Nt = 4;                                    % 송신안테나 개수
Nr = 2;                                    % 수신안테나 개수 / 유저
K = Nt/Nr;                                 % 유저의 수
modulation_order = 4;                      % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
noise_power = 1;
P=0:3:30;                                  % transmit power
SNR=0;
iteration=10000;
%% Block Digonalization
for power_index=1:length(P)
    for iter=1:iteration
        %% Channel genrate
        H=(randn(Nr*K,Nt)+j*randn(Nr*K,Nt))/sqrt(2);
        %% Data
        X=randi([0 1], [Nt, data_size]);                                    % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, modulation_order);                              % modulation
        %% BD
        [Y_demod,S]=MU_BD(P,power_index,X_mod,H,SNR,modulation_order,Nt,Nr,noise_power);
        num_error(iter,power_index)=biterr(X,Y_demod);                    % x,y의 총 error 발생 횟수
        %% Sum rate
        SR(iter,power_index)=S;    
    end
end
error_rate=(sum(num_error,1)/(data_size*Nt))/iteration;
sum_rate=(sum(abs(SR),1))/iteration;
%% graph
semilogy(P,error_rate,'-o')
title('BER Performance'), xlabel('Transmit Power(dB)'),ylabel('BER'),legend('BD')
grid on
figure
plot(P,sum_rate,'-o')
title('Sum Rate'), xlabel('Transmit Power(dB)'),ylabel('Sum Rate(bps/Hz)'),legend('BD')
grid on