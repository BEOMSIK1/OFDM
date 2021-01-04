clc, clear, close all;
%% Parameters
N_t = 4;     % 송신안테나 개수
N_r = 1;     % 수신안테나 개수 / 유저
K = N_t/N_r; % 유저의 수

modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;

SNR=0:3:30;
Iteration=5000;

%% channel generate
for i=1:K
    if N_r==1                                                               % 수신 안테나 1개인 경우
        H(i,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    else                                                                    % 수신 안테나 2개 이상인 경우
        idx = N_r*i;                                                        % 인덱싱을 위한 상수
        H(idx-(N_r-1):idx,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    end
end

%% channel inversion
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        X = randi([0 1], [N_t, data_size]);              % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, modulation_order);           % modulation
        X_power = norm(X_mod,'fro');                     % rho 값 (modulation된 신호의 power, 프로베니우스 놈)
        gamma = trace(inv(H*H'))./X_power;               % power 정규화 factor
        precoded_matrix = (H'*inv(H*H'))./sqrt(gamma);   % 전처리 행렬
        precoded_X = precoded_matrix*X_mod;              % 사용자 신호 전처리
        Y_received = awgn_noise(H*precoded_X,SNR_index); % H*전처리된 각 사용자 신호/sqrt(gamma) + noise
        Y = Y_received.*sqrt(gamma);                     % 정규화 factor를 곱해주어 사용자 신호만 추출
        
        Y_demod = base_demod(Y,modulation_order);
        num_error(Iter,SNR_index)=biterr(X,Y_demod);     % x,y의 총 error 발생 횟수

    end
end
error_rate=(sum(num_error,1)/(data_size*N_t))/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on