clc, clear, close all;
%% Parameters
Modulation_Order=2;                   % ���� ��� 
FFT_Size=128;                         % �ݼ��� ����
Data_Size =FFT_Size*Modulation_Order; % ������ ũ��
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % ���߰�� ����
Nr=4;
Nt=4;
SNR=0:3:30;
Iteration=1000;
%%
for SNR_index=1:length(SNR)
    N_0=10^(-SNR_index/10); 
    for Iter=1:Iteration
        %% MIMO-OFDM
        Data = randi([0 1],[Nt Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT����� �����Ϳ� CP����
        copied_cpdata=repmat(Add_CP_data,Nr,1);
        h=rayleigh_channel([Multi_path, Nt*Nr]);                             % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        for k=1:(Nt*Nr)
            hx(k,:)=conv(copied_cpdata(k,:),h(k,:));
        end
        for k=1:Nr
            hx_comb(k,:)=sum(hx([(Nt*(k-1)+1):Nt*k],:));
        end
         y=awgn_noise(hx_comb,SNR(SNR_index));                                    % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
         y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
         Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
         %% Zero Forcing-OSIC
         H_rv=reshape(H,Nt,Nr,[]);
         X_hat_zf=ZF_OSIC(FFT_Size,Modulation_Order,Nt,H_rv,Y);
         X_hat_demod_zf=base_demod(X_hat_zf,Modulation_Order);
         num_error_zf(Iter,SNR_index)=biterr(Data,X_hat_demod_zf);
         %% Minimum Mean-Squared Error_OSIC
         X_hat_mmse=MMSE_OSIC(FFT_Size,Modulation_Order,Nt,H_rv,Y,N_0);
         X_hat_demod_mmse=base_demod(X_hat_mmse,Modulation_Order);
         num_error_mmse(Iter,SNR_index)=biterr(Data,X_hat_demod_mmse);
    end
    error_rate_zf=(sum(num_error_zf,1)/(Data_Size*Nt))/Iteration; 
    error_rate_mmse=(sum(num_error_mmse,1)/(Data_Size*Nt))/Iteration;
 end
semilogy(SNR,error_rate_zf,'-o')
hold on
semilogy(SNR,error_rate_mmse,'-*')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('ZF-OSIC','MMSE-OSIC'),axis([0 30 1e-4 1]),grid on