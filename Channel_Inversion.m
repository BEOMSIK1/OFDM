clc, clear, close all;
%% Parameters
Nt = 16;                                     % 송신안테나 개수
Nr = 1;                                     % 수신안테나 개수 / 유저
K = 8;                                  % 유저의 수
modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
P=0:3:30;                                   % power
SNR=0;
iteration=1000;
%% Channel inversion
for power_index=1:length(P)
    for iter=1:iteration
        %% Channel genrate
        H=(randn(Nr*K,Nt)+j*randn(Nr*K,Nt))/sqrt(2);
        %% Data
        X=randi([0 1], [K, data_size]);                                    % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, modulation_order);                              % modulation
        %% ZF
        Y_demod_ZF=MU_ZF(P,power_index,X_mod,H,SNR,modulation_order,Nt);
        %% MMSE
        Y_demod_MMSE=MU_MMSE(P,power_index,X_mod,H,SNR,modulation_order,K,Nt);
        %% Error
        num_error_ZF(iter,power_index)=biterr(X,Y_demod_ZF);                      % x,y의 총 error 발생 횟수
        num_error_MMSE(iter,power_index)=biterr(X,Y_demod_MMSE); 
    end
end
error_rate_ZF=(sum(num_error_ZF,1)/(data_size*Nt))/iteration;
error_rate_MMSE=(sum(num_error_MMSE,1)/(data_size*Nt))/iteration;
%% graph
semilogy(P,error_rate_ZF,'-o')
hold on
semilogy(P,error_rate_MMSE,'-*')
title('BER Performance'), xlabel('Transmit Power(dB)'),ylabel('BER'),legend('ZF','MMSE')
grid on