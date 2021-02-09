clc, clear, close all;
%% Parameters
Modulation_Order=2;                   % 변조 방법 2:QPSK, 4:QAM
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Nt=2;                                  % 송신 안테나 갯수
Nr=1;                                  % 수신 안테나 갯수
SNR=0:3:30;
Iteration=5000;
%% Spatial Phase Coding
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (Uncoded)
        Data = randi([0 1],[Nt Data_Size]);
        mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        h=rayleigh_channel([Multi_path Nt]);                                 % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        Y=MISO_OFDM(mod_data,FFT_Size,GI_Size,h,SNR_index,Nt);
        Y_equalize=Y./H;
        Y_demod=base_demod(Y_equalize,Modulation_Order);                    % equalize된 데이터를 복조
        %% SPC 1-bit
        [pc_x_1bit,H_spc_1bit]=SPC_1bit(FFT_Size,mod_data,H,Nt);
        Y_spc_1bit=MISO_OFDM(pc_x_1bit,FFT_Size,GI_Size,h,SNR_index,Nt);
        X_spc_1bit=sum(Y_spc_1bit)./H_spc_1bit;
        X_demod_spc_1bit=base_demod(X_spc_1bit,Modulation_Order);
        %% SPC 2-bit
        [pc_x_2bit,H_spc_2bit]=SPC_2bit(FFT_Size,mod_data,H,Nt);
        Y_spc_2bit=MISO_OFDM(pc_x_2bit,FFT_Size,GI_Size,h,SNR_index,Nt);
        X_spc_2bit=sum(Y_spc_2bit)./H_spc_2bit;
        X_demod_spc_2bit=base_demod(X_spc_2bit,Modulation_Order);
        %% error
        num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장
        num_error_spc_1bit(Iter,SNR_index)=biterr(Data(1,:),X_demod_spc_1bit);
        num_error_spc_2bit(Iter,SNR_index)=biterr(Data(1,:),X_demod_spc_2bit);
    end
end
error_rate=(sum(num_error,1)/(Data_Size*Nt))/Iteration;
error_rate_spc_1bit=(sum(num_error_spc_1bit,1)/(Data_Size))/Iteration;
error_rate_spc_2bit=(sum(num_error_spc_2bit,1)/(Data_Size))/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')
hold on
semilogy(SNR,error_rate_spc_1bit,'-*')
hold on
semilogy(SNR,error_rate_spc_2bit,'-s')
legend('MISO-OFDM','SPC(1bit)','SPC(2bit)')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on