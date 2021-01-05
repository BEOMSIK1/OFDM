clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % 변조 방법 2:QPSK, 4:QAM
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
N=2;                                  % 송신 안테나 갯수
M=1;                                  % 수신 안테나 갯수
cl=3;                                 % constraint length
code_rate=1/2;                        % code rate
                             
SNR=0:3:30;
Iteration=1000;

%% CDD
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (CDD와 성능 비교대상)
        Data = randi([0 1],[N Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);                            % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입

        for Num_channel=1:N                                                 % 송신 안테나 갯수만큼 채널 생성
            
            h(Num_channel,:)=rayleigh_channel(Multi_path);          % time domain channel
            H(Num_channel,:)=fft(h(Num_channel,:),FFT_Size);        % frequency domain channel
            
            hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % 전송된 데이터가 다중경로 채널 통과(h*x)
            y(Num_channel,:)=awgn_noise(hx(Num_channel,:),SNR(SNR_index)); % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
            
            y_remove_CP(Num_channel,:)=y(Num_channel,GI_Size+1:GI_Size+FFT_Size);      %remove CP
            Y(Num_channel,:)=fft(y_remove_CP(Num_channel,:),FFT_Size)/sqrt(FFT_Size);  % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
            Y_equalize(Num_channel,:)=Y(Num_channel,:)./H(Num_channel,:); 
            Y_demod(Num_channel,:)=base_demod(Y_equalize(Num_channel,:),Modulation_Order);                    % equalize된 데이터를 복조 
            

        end
            num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장
    end
end