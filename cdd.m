clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % ���� ��� 2:QPSK, 4:QAM
FFT_Size=128;                         % �ݼ��� ����
Data_Size =FFT_Size*Modulation_Order; % ������ ũ��
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % ���߰�� ����
Num_rx_antenna=2;                     % ���� ���׳� ����
Num_tx_antenna=2;                     % �۽� ���׳� ����
cl=3;                                 % constraint length
code_rate=1/2;                        % code rate
SNR=0:3:30;
Iteration=1000;

%% CDD
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (CDD�� ���� �񱳴��)
        
    end
end