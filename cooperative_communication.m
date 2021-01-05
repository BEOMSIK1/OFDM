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
        %% OFDM ��ȣ �۽�
        Data=randi([0 1],[1 Data_Size]);          % 0 or 1�� ���� ������ Data_Sizeũ�⸸ŭ�� ������ ����
        mod_data=base_mod(Data,Modulation_Order); % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);  % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(FFT_Size-GI_Size+1:end), IFFT_data]; % IFFT����� �����Ϳ� CP����
        
        %% ä��(time domain)
        h=rayleigh_channel(Multi_path);    % ���߰�� ä�� (SISO_OFDM)
        h_sr=rayleigh_channel(Multi_path); % source-relay
        h_sd=rayleigh_channel(Multi_path); % source-destination
        h_rd=rayleigh_channel(Multi_path); % relay-destination
        
        hx=conv(Add_CP_data,h);            % ���۵� �����Ͱ� ���߰�� ä���� ��� (SISO_OFDM)
        h_sr_x=conv(Add_CP_data,h_sr);     % source-relay 
        h_sd_x=conv(Add_CP_data,h_sd);     % source-destination
        
        y=awgn_noise(hx,SNR(SNR_index));              % ä�� ����� �����Ϳ� awgn�߰� (SISO_OFDM)
        [y_sr,N_0]=awgn_noise(h_sr_x,SNR(SNR_index)); % source-relay 
        y_sd=awgn_noise(h_sd_x,SNR(SNR_index));       % source-destination
        
        
        
        %% ä��(frequency domain)
        H=fft(h,FFT_Size);        % ä�� ���ļ� ���� (SISO-OFDM)
        H_sr=fft(h_sr,FFT_Size);  % source-relay
        H_sd=fft(h_sd,FFT_Size);  % source-destination
        H_rd=fft(h_rd,FFT_Size);  % relay-destination
        
        %% Relay(time->frequency)
        y_sr_remove_cp=y_sr(GI_Size+1:GI_Size+FFT_Size);     % y_sr ��ȣ cp����
        Y_sr=fft(y_sr_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft����
        Y_sr_equalize=Y_sr./H_sr;
        
        %% Destination(time->frequency)
        y_sd_remove_cp=y_sd(GI_Size+1:GI_Size+FFT_Size);     % y_sd ��ȣ cp����
        Y_sd=fft(y_sd_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft����
        
        
        
        %% AF
        beta_r=sqrt(1./(H_sr.*conj(H_sr)+N_0));                    % amplift factor (�۽� ���� p=1)
        Y_sr_amp=beta_r.*Y_sr;                                     % ������ ��ȣ(frequency)
        y_sr_amp=ifft(Y_sr_amp,FFT_Size)*sqrt(FFT_Size);           % ������ ��ȣ ifft����
        y_sr_amp_cp=[y_sr_amp(FFT_Size-GI_Size+1:end), y_sr_amp];  % cp �߰�
        
        h_rd_y=conv(y_sr_amp_cp,h_rd);                             % relay-destination ä�� ���
        y_rd_af=awgn_noise(h_rd_y,SNR(SNR_index));                 % awgn �߰�
        
        y_rd_af_remove_cp=y_rd_af(GI_Size+1:GI_Size+FFT_Size);     % destination���� cp����
        Y_rd_af=fft(y_rd_af_remove_cp,FFT_Size)/sqrt(FFT_Size);    % fft ���� (time->frequency)
       
        a1_af=conj(H_sd)/N_0;                                      % combining factor
        a2_af=(beta_r.*conj(H_sr).*conj(H_rd))./(((beta_r.^2).*(H_rd.*conj(H_rd))+1)*N_0);
        
        Y_af=a1_af.*Y_sd+a2_af.* Y_rd_af;                              % Y_sd, Y_rd ����
        X_hat_af=(1./(a1_af.*H_sd+a2_af.*(beta_r.*H_rd.*H_sr))).*Y_af; % ���� ���� X��
        Y_af_demod=base_demod(X_hat_af,Modulation_Order);              % ���� �� ����
        
        
        %% DF
       
        Y_sr_demod=base_demod(Y_sr_equalize,Modulation_Order);          % source���� ���� ������ ���� ��
        
        Y_sr_mod=base_mod(Y_sr_demod,Modulation_Order);                 % decode�� ��ȣ ����
        IFFT_Y_sr=ifft(Y_sr_mod)*sqrt(FFT_Size);                        % IFFT����
        y_sr_cp=[IFFT_Y_sr(FFT_Size-GI_Size+1:end), IFFT_Y_sr];         % CP�߰�
        
        h_rd_x_df=conv(y_sr_cp,h_rd);                                   % relay-destination ä�� ���
        y_rd_df=awgn_noise(h_rd_x_df,SNR(SNR_index));                   % awgn �߰�
        
        y_rd_df_remove_cp=y_rd_df(GI_Size+1:GI_Size+FFT_Size);          % y_rd cp����
        Y_rd_df=fft(y_rd_df_remove_cp,FFT_Size)/sqrt(FFT_Size);         % fft����
       
        
        a1_df=conj(H_sd)/N_0;                                           % combining factor
        a2_df=conj(H_rd)/N_0;
        
        Y_df=a1_df.*Y_sd+a2_df.*Y_rd_df;                                % Y_sd, Y_rd ����
        X_hat_df=(1./(a1_df.*H_sd+a2_df.*H_rd)).*Y_df;                  % ���� X������
        Y_df_demod=base_demod(X_hat_df,Modulation_Order);               % ���� �� ����
        
        
        %% SISO OFDM ����
        y_remove_CP=y(GI_Size+1:GI_Size+FFT_Size);        %remove CP
        Y=fft(y_remove_CP,FFT_Size)/sqrt(FFT_Size);       % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
        Y_equalize=Y./H;                                  % ��ȭ����
        Y_demod=base_demod(Y_equalize,Modulation_Order);  % equalize�� �����͸� ���� 
        
        %% ���� ����
        num_error(Iter,SNR_index)=biterr(Data,Y_demod);          % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����(SISO-OFDM)
        num_error_df(Iter,SNR_index)=biterr(Data,Y_df_demod);    % DF ���
        num_error_af(Iter,SNR_index)=biterr(Data,Y_af_demod);    % AF ���
    end
 
end
error_rate=(sum(num_error,1)/Data_Size)/Iteration;       %���� ��Ʈ���� ������ ���Ͽ� Data_size�� ������ Ȯ������ ���� ���ѵ� Iteration���� ������ ��� ���� ����
error_rate_df=(sum(num_error_df,1)/Data_Size)/Iteration; % DF ��� ����
error_rate_af=(sum(num_error_af,1)/Data_Size)/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')                            % SISO-OFDM
hold on
semilogy(SNR,error_rate_df,'-*')                         % DF���
hold on
semilogy(SNR,error_rate_af,'-+')                         % AF���
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('SISO-OFDM','Decode-and-Forward','Amplitude-and-Forward'),axis([0 30 1e-4 1]),grid on
