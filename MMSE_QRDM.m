function [X_hat] = MMSE_QRDM (FFT_Size,Modulation_Order,M,Nt,Nr,H,Y,N_0,reference)
%%%%%%%%%%% M=후보군 개수 (1~M)%%%%%%%%%%%%%%%%%%%%%%

%% MMSE_QRDM
for K=1:FFT_Size
    H_ch=transpose(H(:,:,K));       % MIMO channel
    H_bar=[H_ch;sqrt(N_0)*eye(Nt)];
    Y_bar=[Y;zeros(Nt,FFT_Size)];
    G=inv(H_bar'*H_bar)*H_bar';
    norm_g_mmse=sum(G.*conj(G),2);
    [~,norm_idx]=sort(norm_g_mmse,'descend');
    P=zeros(Nr,Nr);
    for i=1:Nr
        P(norm_idx(i),i)=1;
    end
    H_s=H_bar*P;
    [Qs,Rs]=qr(H_s);
    Z=Qs'*Y_bar(:,K);
    
    %% Initialize
    idx_temp=[];
    symbol_num=[];
    selc_ref=[];
    %% Error calculate
    for i=Nt:-1:1
        if i~=Nt
            for j=1:M
                rx=Rs(i,i+1:end)*selc_ref(i+1:end,j);
                e(:,j)=(abs(repmat(Z(i),2^(Modulation_Order),1)-(Rs(i,i)*reference+repmat(rx,2^(Modulation_Order),1))).^2)+e_temp(j);             % 누적에러 (Nt-1 ~ 1)layer
            end
            e=reshape(e,[],1);
        else
            e=abs(repmat(Z(i),2^(Modulation_Order),1)-Rs(i,i)*reference).^2;                                                                      % Nt layer 에러
        end
        %% indexing
        [e_temp,idx]=sort(e);
        e=[];
        idx_temp(Nt+1-i,:)=idx(1:M);
        %% symbol select
        if i~=Nt
            remainder=rem(idx_temp(Nt+1-i,:),2^(Modulation_Order));         % 에러 index를 심볼 개수로 나눔=> 나머지=심볼 번호를 의미
            symbol_num(i,:)=remainder;
            n=find(symbol_num(i,:)==0);                                     % 나머지가 0인경우 마지막 심볼번호를 의미
            symbol_num(i,n)=2^(Modulation_Order);                           % 0인 경우에 마지막 심볼번호 대입
            zero_idx=find(remainder==0);                                    % 현재 심볼이 이전 어느심볼의 가지에서 나왔는가 추정, 몫이 이전 심볼 번호에 해당
            idx_temp(Nt+1-i,zero_idx)=idx_temp(Nt+1-i,zero_idx)-1;          % 마지막 심볼번호는 몫이 1이 더크므로 -1
            u=fix(idx_temp(Nt+1-i,:)./2^(Modulation_Order))+1;              % 모든 심볼 번호에 +1해주면 이전 심볼의 index를 출력
            symbol_num(i+1:end,:)=symbol_num(i+1:end,u);                    % 이전 심볼과 현재 심볼을 묶어서 열 변환
            selc_ref(i:end,:)=reference(symbol_num(i:end,:));               % 해당 번호의 reference 신호 mapping
        else
            symbol_num(i,:)=rem(idx_temp,2^(Modulation_Order));
            n=find(symbol_num(i,:)==0);
            symbol_num(i,n)=2^(Modulation_Order);
            selc_ref(i,:)=transpose(reference(symbol_num(i,1:M)));
        end
    end
    X_hat(:,K)=P*selc_ref(:,1);
end
