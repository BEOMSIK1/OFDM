clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % ���� ��� 2:QPSK, 4:QAM
FFT_Size=128;                         % �ݼ��� ����
Data_Size =FFT_Size*Modulation_Order; % ������ ũ��
GI_Size=FFT_Size/4;                   % CP size
Multi_path=1;                         % ���߰�� ����
N=2;                                  % �۽� ���׳� ����
M=1;                                  % ���� ���׳� ����
cl=3;                                 % constraint length
code_rate=1/2;                        % code rate
                             
SNR=0:3:30;
Iteration=1000;

%% Cyclic Delay Diversity
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (Uncoded)
        Data = randi([0 1],[N Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT����� �����Ϳ� CP����
 
        h=rayleigh_channel([Multi_path N]);                                 % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        
         for Num_channel=1:N
             hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
         end
         
         y=awgn_noise(hx,SNR(SNR_index));                                    % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
         y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
         Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
         Y_equalize=Y./H; 
         Y_demod=base_demod(Y_equalize,Modulation_Order);                    % equalize�� �����͸� ���� 
        
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
            cyclic_hx(Num_channel,:)=conv(cyclic_cp_data(Num_channel,:),h(Num_channel,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
        end

         cyclic_y_combine=awgn_noise(sum(cyclic_hx),SNR(SNR_index));   
         cyclic_y_remove_CP=cyclic_y_combine(:,GI_Size+1:GI_Size+FFT_Size);  %remove CP
         cyclic_Y=fft(cyclic_y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);         % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
         
         cyclic_X=cyclic_Y./sum(H_eff);
         cyclic_X_demod=base_demod(cyclic_X,Modulation_Order);
         
         %% MISO-OFDM (Coded)
         Trellis=poly2trellis(cl,codegenerator(code_rate,cl));
         for k=1:N
            conv_encoded_data(k,:)=convenc(Data(k,:),Trellis);
            interleaver_data(k,:)=interleaver(conv_encoded_data(k,:),1);
         end
         
         mod_data_cd=base_mod(interleaver_data,Modulation_Order);
         ifft_data_cd=ifft(mod_data_cd,FFT_Size/code_rate,2)*sqrt(FFT_Size/code_rate);
         add_cp_data_cd=[ifft_data_cd(:,(FFT_Size/code_rate)-(GI_Size/code_rate)+1:end), ifft_data_cd];
         H_cd=fft(h,FFT_Size/code_rate,2);  
         for k=1:N
            hx_cd(k,:)=conv(add_cp_data_cd(k,:),h(k,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
         end
         
         
         y_cd=awgn_noise(hx_cd,SNR(SNR_index));                                    % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
         y_remove_CP_cd=y_cd(:,(GI_Size/code_rate)+1:(GI_Size+FFT_Size)/code_rate);                        %remove CP
         Y_cd=fft(y_remove_CP_cd,FFT_Size/code_rate,2)/sqrt(FFT_Size/code_rate);                       % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
         Y_equalize_cd=Y_cd./H_cd; 
         Y_demod_dc=base_demod(Y_equalize_cd,Modulation_Order);                    % equalize�� �����͸� ���� 
        
         for k=1:N
            deinterleaver_data(k,:)=interleaver(Y_demod_dc(k,:),2);
            conv_decoded_data(k,:)=vitdec(deinterleaver_data(k,:),Trellis,cl*5,'trunc','hard');
         end
         
         %% Cyclic Delay (Coded)
         
         CD_cd=(FFT_Size/code_rate)/(2.^Modulation_Order);
         for Num_channel=1:N
            cyclic_data_cd(Num_channel,:)=circshift(ifft_data_cd(1,:),CD_cd*(Num_channel-1));
            for k=1:(FFT_Size/code_rate)
                H_eff_cd(Num_channel,k)=H_cd(Num_channel,k).*exp((-j*2*pi*(k-1)*(Num_channel-1)*CD_cd)/(FFT_Size/code_rate));
            end
        end
        cyclic_cp_data_cd=[cyclic_data_cd(:,(FFT_Size/code_rate)-(GI_Size/code_rate)+1:end), cyclic_data_cd];
        
        for Num_channel=1:N
            cyclic_hx_cd(Num_channel,:)=conv(cyclic_cp_data_cd(Num_channel,:),h(Num_channel,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
        end

         cyclic_y_combine_cd=awgn_noise(sum(cyclic_hx_cd),SNR(SNR_index));   
         cyclic_y_remove_CP_cd=cyclic_y_combine_cd(:,(GI_Size/code_rate)+1:(GI_Size+FFT_Size)/code_rate);  %remove CP
         cyclic_Y_cd=fft(cyclic_y_remove_CP_cd,FFT_Size/code_rate,2)/sqrt(FFT_Size/code_rate);         % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
         
         cyclic_X_cd=cyclic_Y_cd./sum(H_eff_cd);
         cyclic_X_demod_cd=base_demod(cyclic_X_cd,Modulation_Order);
         

         deinterleaver_data_cyc=interleaver(cyclic_X_demod_cd,2);
         conv_decoded_data_cyc=vitdec(deinterleaver_data_cyc,Trellis,cl*5,'trunc','hard');

         
        %% error 
         num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����
         num_error_cyc(Iter,SNR_index)=biterr(Data(1,:),cyclic_X_demod); 
         num_error_cd(Iter,SNR_index)=biterr(Data,conv_decoded_data);
         num_error_cyc_cd(Iter,SNR_index)=biterr(Data(1,:),conv_decoded_data_cyc); 
     end
end
error_rate=(sum(num_error,1)/(Data_Size*N))/Iteration;
error_rate_cyc=(sum(num_error_cyc,1)/(Data_Size))/Iteration;
error_rate_cd=(sum(num_error_cd,1)/(Data_Size*N))/Iteration;
error_rate_cyc_cd=(sum(num_error_cyc_cd,1)/(Data_Size))/Iteration;
%% graph
 semilogy(SNR,error_rate,'-o')
 hold on
 semilogy(SNR,error_rate_cyc,'-*')
 hold on
 semilogy(SNR,error_rate_cd,'-x')
 hold on
 semilogy(SNR,error_rate_cyc_cd,'-s')
 legend('OFDM','CDD','OFDM(coded)','CDD(coded)')
 title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
 grid on