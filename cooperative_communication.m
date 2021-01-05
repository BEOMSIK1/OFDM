clc, clear, close all;

%% Parameters
FFT_Size=256;
GI_Size=FFT_Size/4;
Modulation_Order=4;
Data_Size=FFT_Size*Modulation_Order;
Multi_path=7;
Iteration=10000;
SNR=0:3:30;

%% Cooperative Communication(BER performance)
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% OFDM 신호 송신
        Data=randi([0 1],[1 Data_Size]);          % 0 or 1의 값을 가지는 Data_Size크기만큼의 데이터 생성
        mod_data=base_mod(Data,Modulation_Order); % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);  % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(FFT_Size-GI_Size+1:end), IFFT_data]; % IFFT연산된 데이터에 CP삽입
        
        %% 채널(time domain)
        h=rayleigh_channel(Multi_path);    % 다중경로 채널 (SISO_OFDM)
        h_sr=rayleigh_channel(Multi_path); % source-relay
        h_sd=rayleigh_channel(Multi_path); % source-destination
        h_rd=rayleigh_channel(Multi_path); % relay-destination
        
        hx=conv(Add_CP_data,h);            % 전송된 데이터가 다중경로 채널을 통과 (SISO_OFDM)
        h_sr_x=conv(Add_CP_data,h_sr);     % source-relay 
        h_sd_x=conv(Add_CP_data,h_sd);     % source-destination
        
        y=awgn_noise(hx,SNR(SNR_index));              % 채널 통과된 데이터에 awgn추가 (SISO_OFDM)
        [y_sr,N_0]=awgn_noise(h_sr_x,SNR(SNR_index)); % source-relay 
        y_sd=awgn_noise(h_sd_x,SNR(SNR_index));       % source-destination
        
        
        
        %% 채널(frequency domain)
        H=fft(h,FFT_Size);        % 채널 주파수 응답 (SISO-OFDM)
        H_sr=fft(h_sr,FFT_Size);  % source-relay
        H_sd=fft(h_sd,FFT_Size);  % source-destination
        H_rd=fft(h_rd,FFT_Size);  % relay-destination
        
        %% Relay(time->frequency)
        y_sr_remove_cp=y_sr(GI_Size+1:GI_Size+FFT_Size);     % y_sr 신호 cp제거
        Y_sr=fft(y_sr_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft연산
        Y_sr_equalize=Y_sr./H_sr;
        
        %% Destination(time->frequency)
        y_sd_remove_cp=y_sd(GI_Size+1:GI_Size+FFT_Size);     % y_sd 신호 cp제거
        Y_sd=fft(y_sd_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft연산
        
        
        
        %% AF
        beta_r=sqrt(1./(H_sr.*conj(H_sr)+N_0));                    % amplift factor (송신 전력 p=1)
        Y_sr_amp=beta_r.*Y_sr;                                     % 증폭된 신호(frequency)
        y_sr_amp=ifft(Y_sr_amp,FFT_Size)*sqrt(FFT_Size);           % 증폭된 신호 ifft연산
        y_sr_amp_cp=[y_sr_amp(FFT_Size-GI_Size+1:end), y_sr_amp];  % cp 추가
        
        h_rd_y=conv(y_sr_amp_cp,h_rd);                             % relay-destination 채널 통과
        y_rd_af=awgn_noise(h_rd_y,SNR(SNR_index));                 % awgn 추가
        
        y_rd_af_remove_cp=y_rd_af(GI_Size+1:GI_Size+FFT_Size);     % destination에서 cp제거
        Y_rd_af=fft(y_rd_af_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft 연산 (time->frequency)
       
        a1_af=conj(H_sd)/N_0;                                      % combining factor
        a2_af=(beta_r.*conj(H_sr).*conj(H_rd))./(((beta_r.^2).*(H_rd.*conj(H_rd))+1)*N_0);
        
        Y_af=a1_af.*Y_sd+a2_af.* Y_rd_af;                              % Y_sd, Y_rd 결합
        X_hat_af=(1./(a1_af.*H_sd+a2_af.*(beta_r.*H_rd.*H_sr))).*Y_af; % 최종 추정 X값
        Y_af_demod=base_demod(X_hat_af,Modulation_Order);              % 최종 값 복조
        
        
        %% DF
       
        Y_sr_demod=base_demod(Y_sr_equalize,Modulation_Order);          % source에서 보낸 데이터 추정 값
        
        Y_sr_mod=base_mod(Y_sr_demod,Modulation_Order);                 % decode된 신호 변조
        IFFT_Y_sr=ifft(Y_sr_mod)*sqrt(FFT_Size);                        % IFFT연산
        y_sr_cp=[IFFT_Y_sr(FFT_Size-GI_Size+1:end), IFFT_Y_sr];         % CP추가
        
        h_rd_x_df=conv(y_sr_cp,h_rd);                                   % relay-destination 채널 통과
        y_rd_df=awgn_noise(h_rd_x_df,SNR(SNR_index));                   % awgn 추가
        
        y_rd_df_remove_cp=y_rd_df(GI_Size+1:GI_Size+FFT_Size);          % y_rd cp제거
        Y_rd_df=fft(y_rd_df_remove_cp,FFT_Size)/sqrt(FFT_Size);         % fft연산
       
        
        a1_df=conj(H_sd)/N_0;                                           % combining factor
        a2_df=conj(H_rd)/N_0;
        
        Y_df=a1_df.*Y_sd+a2_df.*Y_rd_df;                                % Y_sd, Y_rd 결합
        X_hat_df=(1./(a1_df.*H_sd+a2_df.*H_rd)).*Y_df;                  % 최종 X추정값
        Y_df_demod=base_demod(X_hat_df,Modulation_Order);               % 최종 값 복조
        
        
        %% SISO OFDM 수신
        y_remove_CP=y(GI_Size+1:GI_Size+FFT_Size);        %remove CP
        Y=fft(y_remove_CP,FFT_Size)/sqrt(FFT_Size);       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
        Y_equalize=Y./H;                                  % 등화과정
        Y_demod=base_demod(Y_equalize,Modulation_Order);  % equalize된 데이터를 복조 
        
        %% 오류 측정
        num_error(Iter,SNR_index)=biterr(Data,Y_demod);          % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장(SISO-OFDM)
        num_error_df(Iter,SNR_index)=biterr(Data,Y_df_demod);    % DF 방식
        num_error_af(Iter,SNR_index)=biterr(Data,Y_af_demod);    % AF 방식
    end
 
end
error_rate=(sum(num_error,1)/Data_Size)/Iteration;       %얻은 비트오류 갯수를 합하여 Data_size로 나누어 확률들의 합을 구한뒤 Iteration으로 나누어 평균 값을 구함
error_rate_df=(sum(num_error_df,1)/Data_Size)/Iteration; % DF 방식 오류
error_rate_af=(sum(num_error_af,1)/Data_Size)/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')                            % SISO-OFDM
hold on
semilogy(SNR,error_rate_df,'-*')                         % DF방식
hold on
semilogy(SNR,error_rate_af,'-+')                         % AF방식
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('SISO-OFDM','Decode-and-Forward','Amplitude-and-Forward'),axis([0 30 1e-4 1]),grid on
