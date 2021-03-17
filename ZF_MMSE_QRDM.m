clc, clear, close all;
%% Parameters
Modulation_Order=2;                   % 변조 방법
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Nr=4;
Nt=4;
M=2;                                  % 후보군

SNR=0:3:30;
Iteration=100;
%%
for SNR_index=1:length(SNR)
    N_0=10^(-SNR(SNR_index)/10);
    for Iter=1:Iteration
        %% MIMO-OFDM
        Data = randi([0 1],[Nt Data_Size]);
        mod_data=base_mod(Data,Modulation_Order)./sqrt(Nt);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입
        copied_cpdata=repmat(Add_CP_data,Nr,1);
        h=rayleigh_channel([Multi_path, Nt*Nr]);                             % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        for k=1:(Nt*Nr)
            hx(k,:)=conv(copied_cpdata(k,:),h(k,:));
        end
        for k=1:Nr
            hx_comb(k,:)=sum(hx([(Nt*(k-1)+1):Nt*k],:));
        end
        y=awgn_noise(hx_comb,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
        y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
        Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size)*sqrt(Nt);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
        
        H_rv=reshape(H,Nt,Nr,[]);                                           % MIMO Channel
        
        %% Reference
        ref=de2bi(0:2^(Modulation_Order)-1,'left-msb');
        mod_ref=base_mod(ref,Modulation_Order);
        
        %% Zero Forcing-
        for K=1:FFT_Size
            H_ch=transpose(H_rv(:,:,K));
            G_zf=inv(H_ch'*H_ch)*H_ch';
            norm_g_zf=sum(G_zf.*conj(G_zf),2);
            [~,norm_idx]=sort(norm_g_zf,'descend');
            P=zeros(Nr,Nr);
            for i=1:Nr
                P(norm_idx(i),i)=1;
            end
            H_s=H_ch*P;
            [Qs,Rs]=qr(H_s);
            Z=Qs'*Y(:,K);
            
            %% Initialize
            idx_temp=[];
            symbol_num=[];
%             e=[];
            selc_ref=[];
            for i=Nt:-1:1
                if i~=Nt
                    for j=1:M
                        rx=Rs(i,i+1:end)*selc_ref(i+1:end,j);
                        e(:,j)=(abs(repmat(Z(i),2^(Modulation_Order),1)-(Rs(i,i)*mod_ref+repmat(rx,2^(Modulation_Order),1))).^2)+e_temp(j);
                    end
                    e=reshape(e,[],1);
                else
                    e=abs(repmat(Z(i),2^(Modulation_Order),1)-Rs(i,i)*mod_ref).^2;
                end
                %% indexing
                [e_temp,idx]=sort(e);
                 e=[];
                idx_temp(Nt+1-i,:)=idx(1:M);
                %% symbol select
                if i~=Nt
                    remainder=rem(idx_temp(Nt+1-i,:),2^(Modulation_Order));
                    symbol_num(i,:)=remainder;
                    n=find(symbol_num(i,:)==0);
                    symbol_num(i,n)=2^(Modulation_Order);
                    remainder=rem(idx_temp(Nt+1-i,:),2^(Modulation_Order));
                    zero_idx=find(remainder==0);
                    idx_temp(Nt+1-i,zero_idx)=idx_temp(Nt+1-i,zero_idx)-1;
                    u=fix(idx_temp(Nt+1-i,:)./2^(Modulation_Order))+1;
                    symbol_num(i+1:end,:)=symbol_num(i+1:end,u);
                    selc_ref(i:end,:)=mod_ref(symbol_num(i:end,:));
                else
                    symbol_num(i,:)=rem(idx_temp,2^(Modulation_Order));
                    n=find(symbol_num(i,:)==0);
                    symbol_num(i,n)=2^(Modulation_Order);
                    selc_ref(i,:)=transpose(mod_ref(symbol_num(i,1:M)));
                end
            end
            X_hat(:,K)=P*selc_ref(:,1);
        end
        X_hat_demod_zf=base_demod(X_hat,Modulation_Order);
        num_error_zf(Iter,SNR_index)=biterr(Data,X_hat_demod_zf);
        %% Minimum Mean-Squared Error_OSIC
        %         for K=1:FFT_Size
        %             H_ch=transpose(H_rv(:,:,K));
        %             H_bar=[H_ch;N_0*eye(Nt)];
        %             Y_bar=[Y;zeros(Nt,FFT_Size)];
        %             G=inv(H_bar'*H_bar)*H_bar';
        %             norm_g_mmse=sum(G.*conj(G),2);
        %             [~,norm_idx]=sort(norm_g_mmse,'descend');
        %             P=zeros(Nr,Nr);
        %             for i=1:Nr
        %                 P(norm_idx(i),i)=1;
        %             end
        %             H_s=H_bar*P;
        %             [Qs,Rs]=qr(H_s);
        %             Z=Qs'*Y_bar(:,K);
        %
        %             %% Initialize
        %             idx_temp=[];
        %             symbol_num=[];
        %             e=[];
        %             selc_ref=[];
        %             for i=Nt:-1:1
        %                 if i~=Nt
        %                     for j=1:M
        %                         rx=Rs(i,i+1:end)*selc_ref(i+1:end,j);
        %                         e(:,j)=(abs(repmat(Z(i),2^(Modulation_Order),1)-(Rs(i,i)*mod_ref+repmat(rx,2^(Modulation_Order),1))).^2)+e_temp(j);             % 누적에러 (Nt-1 ~ 1)layer
        %                     end
        %                     e=reshape(e,[],1);
        %                 else
        %                     e=abs(repmat(Z(i),2^(Modulation_Order),1)-Rs(i,i)*mod_ref).^2;                                                                      % Nt layer 에러
        %                 end
        %                 %% indexing
        %                 [e_temp,idx]=sort(e);
        %                 e=e_temp(1:M);
        %                 idx_temp(Nt+1-i,:)=idx(1:M);
        %                 %% symbol select
        %                 if i~=Nt
        %                     remainder=rem(idx_temp(Nt+1-i,:),2^(Modulation_Order));         % 에러 index를 심볼 개수로 나눔=> 나머지=심볼 번호를 의미
        %                     symbol_num(i,:)=remainder;
        %                     n=find(symbol_num(i,:)==0);                                     % 나머지가 0인경우 마지막 심볼번호를 의미
        %                     symbol_num(i,n)=2^(Modulation_Order);                           % 0인 경우에 마지막 심볼번호 대입
        %                     %             remainder=rem(idx_temp(Nt+1-i,:),2^(Modulation_Order));         %
        %                     zero_idx=find(remainder==0);                                    % 현재 심볼이 이전 어느심볼의 가지에서 나왔는가 추정, 몫이 이전 심볼 번호에 해당
        %                     idx_temp(Nt+1-i,zero_idx)=idx_temp(Nt+1-i,zero_idx)-1;          % 마지막 심볼번호는 몫이 1이 더크므로 -1
        %                     u=fix(idx_temp(Nt+1-i,:)./2^(Modulation_Order))+1;              % 모든 심볼 번호에 +1해주면 이전 심볼의 index를 출력
        %                     symbol_num(i+1:end,:)=symbol_num(i+1:end,u);                    % 이전 심볼과 현재 심볼을 묶어서 열 변환
        %                     selc_ref(i:end,:)=mod_ref(symbol_num(i:end,:));               % 해당 번호의 reference 신호 mapping
        %                 else
        %                     symbol_num(i,:)=rem(idx_temp,2^(Modulation_Order));
        %                     n=find(symbol_num(i,:)==0);
        %                     symbol_num(i,n)=2^(Modulation_Order);
        %                     selc_ref(i,:)=transpose(mod_ref(symbol_num(i,1:M)));
        %                 end
        %             end
        %             X_hat(:,K)=P*selc_ref(:,1);
        %         end
        %         X_hat_demod_zf=base_demod(X_hat,Modulation_Order);
        %         num_error_zf(Iter,SNR_index)=biterr(Data,X_hat_demod_zf);
    end
    error_rate=(sum(num_error_zf,1)/(Data_Size*Nt))/Iteration;
    
end
semilogy(SNR,error_rate,'-o')
% hold on
% semilogy(SNR,error_rate_mmse,'-*')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('MMSE-QRDM'),grid on