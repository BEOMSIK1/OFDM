clc, clear, close all;
%% Parameters
N_t = 4;     % �۽ž��׳� ����
N_r = 1;     % ���ž��׳� ���� / ����
K = N_t/N_r; % ������ ��

modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;

SNR=0:3:30;
Iteration=5000;

%% channel generate
for i=1:K
    if N_r==1                                                               % ���� ���׳� 1���� ���
        H(i,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    else                                                                    % ���� ���׳� 2�� �̻��� ���
        idx = N_r*i;                                                        % �ε����� ���� ���
        H(idx-(N_r-1):idx,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    end
end

%% channel inversion
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        X = randi([0 1], [N_t, data_size]);              % N_t���� ���� �ٸ� ������ (�۽���)
        X_mod = base_mod(X, modulation_order);           % modulation
        X_power = norm(X_mod,'fro');                     % rho �� (modulation�� ��ȣ�� power, ���κ��Ͽ콺 ��)
        gamma = trace(inv(H*H'))./X_power;               % power ����ȭ factor
        precoded_matrix = (H'*inv(H*H'))./sqrt(gamma);   % ��ó�� ���
        precoded_X = precoded_matrix*X_mod;              % ����� ��ȣ ��ó��
        Y_received = awgn_noise(H*precoded_X,SNR_index); % H*��ó���� �� ����� ��ȣ/sqrt(gamma) + noise
        Y = Y_received.*sqrt(gamma);                     % ����ȭ factor�� �����־� ����� ��ȣ�� ����
        
        Y_demod = base_demod(Y,modulation_order);
        num_error(Iter,SNR_index)=biterr(X,Y_demod);     % x,y�� �� error �߻� Ƚ��

    end
end
error_rate=(sum(num_error,1)/(data_size*N_t))/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on