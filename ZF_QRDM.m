function [X_hat] = ZF_QRDM (FFT_Size,Modulation_Order,M,Nt,Nr,H,Y,reference)
%%%%%%%%%%% M=후보군 개수 (1~M)


for K=1:FFT_Size
    H_ch=transpose(H(:,:,K));
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
    e=[];
    selc_ref=[];
    %% Error calculate
    for i=Nt:-1:1
        if i~=Nt
            for j=1:M
                rx=Rs(i,i+1:end)*selc_ref(i+1:end,j);
                e(:,j)=(abs(repmat(Z(i),2^(Modulation_Order),1)-(Rs(i,i)*reference+repmat(rx,2^(Modulation_Order),1))).^2)+e_temp(j);
            end
            e=reshape(e,[],1);
        else
            e=abs(repmat(Z(i),2^(Modulation_Order),1)-Rs(i,i)*reference).^2;
        end
        %% indexing
        [e_temp,idx]=sort(e);
        e=e_temp(1:M);
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
            selc_ref(i:end,:)=reference(symbol_num(i:end,:));
        else
            symbol_num(i,:)=rem(idx_temp,2^(Modulation_Order));
            n=find(symbol_num(i,:)==0);
            symbol_num(i,n)=2^(Modulation_Order);
            selc_ref(i,:)=transpose(reference(symbol_num(i,1:M)));
        end
    end
    X_hat(:,K)=P*selc_ref(:,1);
end