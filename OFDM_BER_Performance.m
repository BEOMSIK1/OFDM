clc, clear, close all;

%% Parameters
FFT_Size=128;                                      %반송파 개수
GI_Size=FFT_Size/4;                                %CP size
Modulation_Order=4;                                %변조 방법 2:QPSK, 4:16QAM
Data_Size=FFT_Size*Modulation_Order; 
Multi_path=7;                                      %다중경로 개수
SNR=0:3:30;
Iteration=1000;

%% OFDM
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        Data=randi([0 1],[1 Data_Size]);                                    % 0 or 1의 값을 가지는 Data_Size크기만큼의 데이터 생성
        mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);                            % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(FFT_Size-GI_Size+1:end), IFFT_data];         % IFFT연산된 데이터에 CP삽입

        h=rayleigh_channel(Multi_path);                                     % 다중경로 채널
        H=fft(h,FFT_Size);                                                  % 채널 주파수 응답
        hx=conv(Add_CP_data,h);                                             % 전송된 데이터가 다중경로 채널을 통과
        y=awgn_noise(hx,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가
        
        y_remove_CP=y(GI_Size+1:GI_Size+FFT_Size);                          % remove CP
        Y=fft(y_remove_CP,FFT_Size)/sqrt(FFT_Size);                         % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
        Y_equalize=Y./H;                                                    % 등화과정
        Y_demod=base_demod(Y_equalize,Modulation_Order);                    % equalize된 데이터를 복조 
        
        num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장
    end
 
end
error_rate=(sum(num_error,1)/Data_Size)/Iteration;                          %얻은 비트오류 갯수를 합하여 Data_size로 나누어 확률들의 합을 구한뒤 Iteration으로 나누어 평균 값을 구함

%% graph
semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('SISO-OFDM'),axis([0 30 1e-4 1]),grid on