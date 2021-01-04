clc, clear, close all;
%% Parameters
N_t = 4;     % �۽ž��׳� ����
N_r = 2;     % ���ž��׳� ���� / ����
K = N_t/N_r; % ������ ��

modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
noise_power = 1;


SNR=0:3:30;
Iteration=10000;

%% channel generate
for i=1:K
    if N_r==1                                                               % ���� ���׳� 1���� ���
        H(i,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    else                                                                    % ���� ���׳� 2�� �̻��� ���
        idx = N_r*i;                                                        % �ε����� ���� ���
        H(idx-(N_r-1):idx,:) = (randn(N_r, N_t)+j*randn(N_r, N_t))/sqrt(2);
    end
end

%% precoding matrix
for i=1:K
    if N_r==1
        H_slash = H;                                                        % H_slash �� ��������� �ӽ� ����
        H_slash(i,:) = [];                                                  % i ��° ä���� ������ H_slash_i ����
        [~,~,V] = svd(H_slash);                                             % singular value decomposition  
        F(:,i) = V(:,N_t);                                                  % precoding matrix (F_j)
    else
        idx = N_r*i; 
        H_slash = H;                                                        % H_slash �� ��������� �ӽ� ����
        H_slash(idx-(N_r-1):idx,:) = [];                                    % i ��° ä���� ������ H_slash_i ����
        [~,~,V] = svd(H_slash);                                             % singular value decomposition  
        F(:,idx-(N_r-1):idx) = V(:,N_t-N_r+1:N_t);                          % precoding matrix (F_j)
    end
end
%% block digonalization (BER performance)
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        X = randi([0 1], [N_t, data_size]);                                 % N_t���� ���� �ٸ� ������ (�۽���)
        X_mod = base_mod(X, modulation_order);                              % modulation
        X_power = norm(X_mod,'fro');                                        % rho �� (modulation�� ��ȣ�� power, ���κ��Ͽ콺 ��)
        gamma = trace(F*F')./X_power;                                       % power ����ȭ factor
        precoded_X = (F*X_mod)./gamma;                                      % ����� ��ȣ ��ó��
        Y_received = awgn_noise(H*precoded_X,SNR_index);                    % H*��ó���� �� ����� ��ȣ/gamma + noise
        Y = Y_received.*gamma;                                              % ����ȭ factor�� �����־� ��ȣ����
        H_BD = H*F;                                                         % BD matrix
        Y_ZF = (inv(H_BD'*H_BD)*H_BD')*Y;                                   % ZF�� �̿��� single user detection
        
        Y_demod = base_demod(Y_ZF,modulation_order);
        num_error(Iter,SNR_index)=biterr(X,Y_demod);                        % x,y�� �� error �߻� Ƚ��
        
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