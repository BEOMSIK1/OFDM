clc, clear, close all;

%% Parameters

Modulation_Order=2;                  % ���� ��� 2:QPSK
Data_Size =128;                      % ������ ũ��
FFT_Size=Data_Size/Modulation_Order; % �ݼ��� ����
GI_Size=FFT_Size/4;                  % CP size
Multi_path=7;                        % ���߰�� ����
Num_Antenna=2;                       % ���� ���׳� ���� (2�� �̻�)
SNR=0:3:30;
Iteration=1000;



%%  OFDM
for SNR_index=1:length(SNR)
    N_0=10^(-SNR_index/10);                        % noise power
    for Iter=1:Iteration
        %% OFDM ��ȣ �۽�
        % SISO-OFDM    -> X_0 ��� antenna(1 X 1)
        % SC, EGC, MRC -> X_0 ��� antenna(1 X N)
        % STBC         -> X_0,X_1 ��� antenna(2 X 1)
        
        X_0=randi([0 1],[1 Data_Size]);            % 0 or 1�� ���� ������ Data_Sizeũ�⸸ŭ�� ������ ����
        X_1=randi([0 1],[1 Data_Size]); 
        X_0_mod=base_mod(X_0,Modulation_Order);    % ���� ��Ŀ� ���� ������ ����
        X_1_mod=base_mod(X_1,Modulation_Order);
        x_0=ifft(X_0_mod)*sqrt(FFT_Size);          % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        x_1=ifft(X_1_mod)*sqrt(FFT_Size);
        x_0_cp=[x_0(FFT_Size-GI_Size+1:end), x_0]; % IFFT����� �����Ϳ� CP����
        x_1_cp=[x_1(FFT_Size-GI_Size+1:end), x_1]; 
        
       % ���׳� ������ŭ ä�� ����(Receiver Diversity)
        h=rayleigh_channel([Multi_path Num_Antenna]);          % time domain channel
        H=fft(h,FFT_Size,2);                                   % frequency domain channel

        %% Tx -> Rx
        for Num_channel=1:Num_Antenna                               % ���׳� ������ŭ ä�� ����(Receiver Diversity)
            % SISO -> h�� 1 ä�� �̿�
            % STBC -> h�� 1,2 ä�� �̿�

            hx(Num_channel,:)=conv(x_0_cp,h(Num_channel,:));               % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
            y(Num_channel,:)=awgn_noise(hx(Num_channel,:),SNR(SNR_index)); % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
            
            y_remove_CP(Num_channel,:)=y(Num_channel,GI_Size+1:GI_Size+FFT_Size);      %remove CP
            Y(Num_channel,:)=fft(y_remove_CP(Num_channel,:),FFT_Size)/sqrt(FFT_Size);  % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
        end
        %% SC
        
        sc_SNR=(Y.*conj(Y))./N_0;                                               % ���۹��� ��ȣ SNR
        [~,SNR_max_index]=max(sc_SNR);                                          % SNR�� ����ū ���׳� ��ȣ
        for Data_num=1:FFT_Size
            Y_sc(Data_num)=Y(SNR_max_index(Data_num),Data_num);                 % Y���� SNR�� ���� ū ��ȣ ����
            X_sc(Data_num)=Y_sc(Data_num)./H(SNR_max_index(Data_num),Data_num); % equalize
        end
      
        X_sc_demod=base_demod(X_sc,Modulation_Order);                           % demodualtion
        
        
        
        %% EGC
        
        egc_weight=exp(j*angle(conj(H)));                   % combine weight ����
        Y_egc_combined=sum((Y.*egc_weight));                % ���۹��� ��ȣ Y�� weight���Ͽ� ����
        X_egc=Y_egc_combined./sum(abs(H));                  % ä���� amplitude�� ������ ���� �� X ����
        X_egc_demod=base_demod(X_egc,Modulation_Order);     % demodulation
        
       
        
        %% MRC
        
        Y_mrc_combine=sum(conj(H).*Y);                    % ���յ� mrc��� ��ȣ   
        X_mrc=Y_mrc_combine./sum(H.*conj(H));            % equalize      
        X_mrc_demod=base_demod(X_mrc,Modulation_Order);   % demodulation
        
       
        %% STBC
        
        % time=t    -> antenna0: x0         antenna1: x1
        % time=t+T  -> antenna0: -conj(x1)  antenna1: conj(x0)
        % ä�� ��� h�� 1,2ä�� �̿�
       
        %  ���۵� �����Ͱ� ���߰�� ä�� ���
        h0x0=conv(x_0_cp,h(1,:));             % h0*x0
        h1x1=conv(x_1_cp,h(2,:));             % h1*x1
        h0x1=conv(-conj(x_1_cp),h(1,:));      % h0*-conj(x1)
        h1x0=conv(conj(x_0_cp),h(2,:));       % h1*conj(x0)
        
        y_0=awgn_noise(h0x0+h1x1,SNR(SNR_index)); % ���� ��ȣ time=t
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
        % ä�� ��� h�� 1ä�� �̿�
        
        H_siso=H(1,:);                                             % ä�� (frequency domain)
        Y_siso=Y(1,:);                                             % ���� ���� ��ȣ (HX+n)
        X_siso=Y_siso./H_siso;                                     % equalize
        X_siso_demod=base_demod(X_siso,Modulation_Order);          % demodualtion 
        
        %%  ��������
        
        num_error_siso(Iter,SNR_index)=biterr(X_0,X_siso_demod);            % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����
        num_error_sc(Iter,SNR_index)=biterr(X_0,X_sc_demod);
        num_error_egc(Iter,SNR_index)=biterr(X_0,X_egc_demod);
        num_error_mrc(Iter,SNR_index)=biterr(X_0,X_mrc_demod); 
        num_error_stbc0(Iter,SNR_index)=biterr(X_0,X_stbc_demod_0);
        num_error_stbc1(Iter,SNR_index)=biterr(X_1,X_stbc_demod_1);    
    end
 
end

error_rate_siso=(sum(num_error_siso,1)/Data_Size)/Iteration; %���� ��Ʈ���� ������ ���Ͽ� Data_size�� ������ Ȯ������ ���� ���ѵ� Iteration���� ������ ��� ���� ����
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
