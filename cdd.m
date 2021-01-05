clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % 변조 방법 2:QPSK, 4:QAM
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Num_rx_antenna=2;                     % 수신 안테나 갯수
Num_tx_antenna=2;                     % 송신 안테나 갯수
cl=3;                                 % constraint length
code_rate=1/2;                        % code rate
SNR=0:3:30;
Iteration=1000;

%% CDD
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (CDD와 성능 비교대상)
        
    end
end