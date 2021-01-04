clc, clear, close all;
%% Parameters
N_t = 4;     % 송신안테나 개수
N_r = 2;     % 수신안테나 개수 / 유저
K = N_t/N_r; % 유저의 수

modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
noise_power = 1;


SNR=0:3:30;
Iteration=10000;

%% channel generate
for i=1:K
    if N_r==1                                                               % 수신 안테나 1개인 경우
        H(i,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    else                                                                    % 수신 안테나 2개 이상인 경우
        idx = N_r*i;                                                        % 인덱싱을 위한 상수
        H(idx-(N_r-1):idx,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    end
end

%% precoding matrix
for i=1:K
    if N_r==1
        H_slash = H;                                                        % H_slash 를 만들기위해 임시 저장
        H_slash(i,:) = [];                                                  % i 번째 채널을 제외한 H_slash_i 생성
        [~,~,V] = svd(H_slash);                                             % singular value decomposition  
        F(:,i) = V(:,N_t);                                                  % precoding matrix (F_j)
    else
        idx = N_r*i; 
        H_slash = H;                                                        % H_slash 를 만들기위해 임시 저장
        H_slash(idx-(N_r-1):idx,:) = [];                                    % i 번째 채널을 제외한 H_slash_i 생성
        [~,~,V] = svd(H_slash);                                             % singular value decomposition  
        F(:,idx-(N_r-1):idx) = V(:,N_t-N_r+1:N_t);                          % precoding matrix (F_j)
    end
end
%% block digonalization (BER performance)
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        X = randi([0 1], [N_t, data_size]);                                 % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, modulation_order);                              % modulation
        X_power = norm(X_mod,'fro');                                        % rho 값 (modulation된 신호의 power, 프로베니우스 놈)
        gamma = trace(F*F')./X_power;                                       % power 정규화 factor
        precoded_X = (F*X_mod)./gamma;                                      % 사용자 신호 전처리
        Y_received = awgn_noise(H*precoded_X,SNR_index);                    % H*전처리된 각 사용자 신호/gamma + noise
        Y = Y_received.*gamma;                                              % 정규화 factor를 곱해주어 신호추출
        H_BD = H*F;                                                         % BD matrix
        Y_ZF = (inv(H_BD'*H_BD)*H_BD')*Y;                                   % ZF를 이용한 single user detection
        
        Y_demod = base_demod(Y_ZF,modulation_order);
        num_error(Iter,SNR_index)=biterr(X,Y_demod);                        % x,y의 총 error 발생 횟수
        
        %% channel capacity
        P = X_power/(noise_power*N_t);
        C(Iter,SNR_index) = log2(det(eye(N_t)+(H_BD*H_BD').*P));
 
    end
end
error_rate=(sum(num_error,1)/(data_size*N_t))/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on