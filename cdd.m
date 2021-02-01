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
Iteration=5000;

%% Cyclic Delay Diversity
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (Uncoded)
        Data = randi([0 1],[N Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입
 
        h=rayleigh_channel([Multi_path N]);                                 % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        
         for Num_channel=1:N
             hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % 전송된 데이터가 다중경로 채널 통과(h*x)
         end
         
         y=awgn_noise(hx,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
         y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
         Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
         Y_equalize=Y./H; 
         Y_demod=base_demod(Y_equalize,Modulation_Order);                    % equalize된 데이터를 복조 
        
        %% Cyclic Delay (Uncoded) 
        
%              if N<=2
%                  CD=(FFT_Size)./(2.^Modulation_Order);
%              else
%                  CD=(FFT_Size/N);
%              end
        
        CD=(FFT_Size)/(2.^Modulation_Order);
        
        for Num_channel=1:N
            cyclic_data(Num_channel,:)=circshift(IFFT_data(1,:),CD*(Num_channel-1));
            for k=1:FFT_Size
                H_eff(Num_channel,k)=H(Num_channel,k).*exp((-j*2*pi*(k-1)*(Num_channel-1)*CD)/FFT_Size);
            end
        end
        cyclic_cp_data=[cyclic_data(:,FFT_Size-GI_Size+1:end), cyclic_data];
        
        for Num_channel=1:N
            cyclic_hx(Num_channel,:)=conv(cyclic_cp_data(Num_channel,:),h(Num_channel,:));  % 전송된 데이터가 다중경로 채널 통과(h*x)
        end

         cyclic_y_combine=awgn_noise(sum(cyclic_hx),SNR(SNR_index));   
         cyclic_y_remove_CP=cyclic_y_combine(:,GI_Size+1:GI_Size+FFT_Size);  %remove CP
         cyclic_Y=fft(cyclic_y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);         % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
         
         cyclic_X=cyclic_Y./sum(H_eff);
         cyclic_X_demod=base_demod(cyclic_X,Modulation_Order);
         
         %% MISO-OFDM (Coded)
         Trellis=poly2trellis(cl,codegenerator(code_rate,cl));
         for k=1:N
             conv_encoded_data(N,:)=convenc(Data(N,:),Trellis);
             interleaver_data(N,:)=interleaver(conv_encoded_data(N,:),1);
         end
        mod_data_cd=base_mod(interleaver_data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
        IFFT_data_cd=ifft(mod_data_cd,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data_cd=[IFFT_data_cd(:,FFT_Size-GI_Size+1:end), IFFT_data_cd];       % IFFT연산된 데이터에 CP삽입

         for Num_channel=1:N
             hx_cd(Num_channel,:)=conv(Add_CP_data_cd(Num_channel,:),h(Num_channel,:));  % 전송된 데이터가 다중경로 채널 통과(h*x)
         end
         
         y_cd=awgn_noise(hx_cd,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
         y_remove_CP_cd=y_cd(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
         Y_cd=fft(y_remove_CP_cd,FFT_Size,2)/sqrt(FFT_Size);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
         Y_equalize_cd=Y_cd./H; 
         Y_demod_cd=base_demod(Y_equalize_cd,Modulation_Order);                    % equalize된 데이터를 복조 
         
        for k=1:N
            deinterleaver_data(N,:)=interleaver(Y_demod_cd(N,:),2);
            conv_decoded_data(N,:)=vitdec(deinterleaver_data(N,:),Trellis,cl*5,'trunc','hard');
        end
        
        %% error 
         num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장
         num_error_cyc(Iter,SNR_index)=biterr(Data(1,:),cyclic_X_demod); 
         num_error_cd(Iter,SNR_index)=biterr(Data,Y_demod_cd);
    end
end
error_rate=(sum(num_error,1)/(Data_Size*N))/Iteration;
error_rate_cyc=(sum(num_error_cyc,1)/(Data_Size))/Iteration;

%% graph
 semilogy(SNR,error_rate,'-o')
 hold on
 semilogy(SNR,error_rate_cyc,'-*')
 legend('OFDM','CDD')
 title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
 grid on