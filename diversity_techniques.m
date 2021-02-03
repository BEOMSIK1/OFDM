clc, clear, close all;

%% Parameters

Modulation_Order=2;                  % 변조 방법 2:QPSK
Data_Size =128;                      % 데이터 크기
FFT_Size=Data_Size/Modulation_Order; % 반송파 갯수
GI_Size=FFT_Size/4;                  % CP size
Multi_path=7;                        % 다중경로 갯수
Num_Antenna=2;                       % 수신 안테나 갯수 (2개 이상)
SNR=0:3:30;
Iteration=1000;



%%  OFDM
for SNR_index=1:length(SNR)
    N_0=10^(-SNR_index/10);                        % noise power
    for Iter=1:Iteration
        %% OFDM 신호 송신
        % SISO-OFDM    -> X_0 사용 antenna(1 X 1)
        % SC, EGC, MRC -> X_0 사용 antenna(1 X N)
        % STBC         -> X_0,X_1 사용 antenna(2 X 1)
        
        X_0=randi([0 1],[1 Data_Size]);            % 0 or 1의 값을 가지는 Data_Size크기만큼의 데이터 생성
        X_1=randi([0 1],[1 Data_Size]); 
        X_0_mod=base_mod(X_0,Modulation_Order);    % 변조 방식에 따른 데이터 변조
        X_1_mod=base_mod(X_1,Modulation_Order);
        x_0=ifft(X_0_mod)*sqrt(FFT_Size);          % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        x_1=ifft(X_1_mod)*sqrt(FFT_Size);
        x_0_cp=[x_0(FFT_Size-GI_Size+1:end), x_0]; % IFFT연산된 데이터에 CP삽입
        x_1_cp=[x_1(FFT_Size-GI_Size+1:end), x_1]; 
        
       % 안테나 갯수만큼 채널 생성(Receiver Diversity)
        h=rayleigh_channel([Multi_path Num_Antenna]);          % time domain channel
        H=fft(h,FFT_Size,2);                                   % frequency domain channel

        %% Tx -> Rx
        for Num_channel=1:Num_Antenna                               % 안테나 갯수만큼 채널 생성(Receiver Diversity)
            % SISO -> h의 1 채널 이용
            % STBC -> h의 1,2 채널 이용

            hx(Num_channel,:)=conv(x_0_cp,h(Num_channel,:));               % 전송된 데이터가 다중경로 채널 통과(h*x)
            y(Num_channel,:)=awgn_noise(hx(Num_channel,:),SNR(SNR_index)); % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
            
            y_remove_CP(Num_channel,:)=y(Num_channel,GI_Size+1:GI_Size+FFT_Size);      %remove CP
            Y(Num_channel,:)=fft(y_remove_CP(Num_channel,:),FFT_Size)/sqrt(FFT_Size);  % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
        end
        %% SC
        
        sc_SNR=(Y.*conj(Y))./N_0;                                               % 전송받은 신호 SNR
        [~,SNR_max_index]=max(sc_SNR);                                          % SNR이 가장큰 안테나 번호
        for Data_num=1:FFT_Size
            Y_sc(Data_num)=Y(SNR_max_index(Data_num),Data_num);                 % Y에서 SNR이 가장 큰 신호 선택
            X_sc(Data_num)=Y_sc(Data_num)./H(SNR_max_index(Data_num),Data_num); % equalize
        end
      
        X_sc_demod=base_demod(X_sc,Modulation_Order);                           % demodualtion
        
        
        
        %% EGC
        
        egc_weight=exp(j*angle(conj(H)));                   % combine weight 설정
        Y_egc_combined=sum((Y.*egc_weight));                % 전송받은 신호 Y와 weight곱하여 결합
        X_egc=Y_egc_combined./sum(abs(H));                  % 채널의 amplitude로 나누어 추정 값 X 구함
        X_egc_demod=base_demod(X_egc,Modulation_Order);     % demodulation
        
       
        
        %% MRC
        
        Y_mrc_combine=sum(conj(H).*Y);                    % 결합된 mrc방식 신호   
        X_mrc=Y_mrc_combine./sum(H.*conj(H));            % equalize      
        X_mrc_demod=base_demod(X_mrc,Modulation_Order);   % demodulation
        
       
        %% STBC
        
        % time=t    -> antenna0: x0         antenna1: x1
        % time=t+T  -> antenna0: -conj(x1)  antenna1: conj(x0)
        % 채널 행렬 h의 1,2채널 이용
       
        %  전송된 데이터가 다중경로 채널 통과
        h0x0=conv(x_0_cp,h(1,:));             % h0*x0
        h1x1=conv(x_1_cp,h(2,:));             % h1*x1
        h0x1=conv(-conj(x_1_cp),h(1,:));      % h0*-conj(x1)
        h1x0=conv(conj(x_0_cp),h(2,:));       % h1*conj(x0)
        
        y_0=awgn_noise(h0x0+h1x1,SNR(SNR_index)); % 수신 신호 time=t
        y_1=awgn_noise(h0x1+h1x0,SNR(SNR_index)); % time=t+T

        y_0_remove_cp=y_0(GI_Size+1:GI_Size+FFT_Size);
        y_1_remove_cp=y_1(GI_Size+1:GI_Size+FFT_Size);
        
        Y_0=fft(y_0_remove_cp,FFT_Size)/sqrt(FFT_Size);   
        Y_1=fft(y_1_remove_cp,FFT_Size)/sqrt(FFT_Size); 
        
        Y_0_combine=conj(H(1,:)).*Y_0+H(2,:).*conj(Y_1);
        Y_1_combine=conj(H(2,:)).*Y_0-H(1,:).*conj(Y_1);
        
        X_stbc_0=Y_0_combine./sum((abs(H).^2));
        X_stbc_1=Y_1_combine./sum((abs(H).^2));
      
        
        X_stbc_demod_0=base_demod(X_stbc_0,Modulation_Order);   % demodulation
        X_stbc_demod_1=base_demod(X_stbc_1,Modulation_Order);   % demodulation
        
        
        %% SISO-OFDM
        % 채널 행렬 h의 1채널 이용
        
        H_siso=H(1,:);                                             % 채널 (frequency domain)
        Y_siso=Y(1,:);                                             % 전송 받은 신호 (HX+n)
        X_siso=Y_siso./H_siso;                                     % equalize
        X_siso_demod=base_demod(X_siso,Modulation_Order);          % demodualtion 
        
        %%  오류측정
        
        num_error_siso(Iter,SNR_index)=biterr(X_0,X_siso_demod);            % 각 SNR당 비트오류갯수를 Iteration마다 배열에 저장
        num_error_sc(Iter,SNR_index)=biterr(X_0,X_sc_demod);
        num_error_egc(Iter,SNR_index)=biterr(X_0,X_egc_demod);
        num_error_mrc(Iter,SNR_index)=biterr(X_0,X_mrc_demod); 
        num_error_stbc0(Iter,SNR_index)=biterr(X_0,X_stbc_demod_0);
        num_error_stbc1(Iter,SNR_index)=biterr(X_1,X_stbc_demod_1);    
    end
 
end

error_rate_siso=(sum(num_error_siso,1)/Data_Size)/Iteration; %얻은 비트오류 갯수를 합하여 Data_size로 나누어 확률들의 합을 구한뒤 Iteration으로 나누어 평균 값을 구함
error_rate_sc=(sum(num_error_sc,1)/Data_Size)/Iteration;
error_rate_egc=(sum(num_error_egc,1)/Data_Size)/Iteration;
error_rate_mrc=(sum(num_error_mrc,1)/Data_Size)/Iteration;
error_rate_stbc0=(sum(num_error_stbc0,1)/Data_Size)/Iteration;
error_rate_stbc1=(sum(num_error_stbc1,1)/Data_Size)/Iteration;

%%  graph

semilogy(SNR,error_rate_siso,'-o')
hold on
semilogy(SNR,error_rate_sc,'-x')
hold on
semilogy(SNR,error_rate_egc,'-+')
hold on
semilogy(SNR,error_rate_mrc,'-*')
hold on
semilogy(SNR,error_rate_stbc0,'s')

title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('SISO-OFDM','sc','egc','mrc','stbc'),axis([0 30 1e-6 1]),grid on
