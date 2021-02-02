clc, clear, close all;
%% Parameters

Modulation_Order=1;                   % 변조 방법 
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Nr=2;
Nt=2;
SNR=0:3:30;
Iteration=5000;
%%
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MIMO-OFDM
        Data = randi([0 1],[Nt Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입
        
        copied_cpdata=repmat(Add_CP_data,Nr,1);
 
        h=rayleigh_channel([Multi_path Nt*Nr]);                                 % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        for k=1:(Nt*Nr)
            hx(k,:)=conv(copied_cpdata(k,:),h(k,:));
        end
        
        for k=1:Nr
            y(k,:)=sum(hx,[k 
         y=awgn_noise(hx,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
         y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
         Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)

         %% reference
         
         
         
         
         %Y_demod=base_demod(Y_equalize,Modulation_Order);                    % equalize된 데이터를 복조 
        
     
        %% error 
         num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장

    end
end
error_rate=(sum(num_error,1)/(Data_Size*N))/Iteration;


%% graph
 semilogy(SNR,error_rate,'-o')

 title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
 grid on