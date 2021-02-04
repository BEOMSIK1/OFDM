function [X_hat] = ZF_DFE (FFT_Size,Modulation_Order,Nt,Nr,H,Y)

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
    for i=Nt:-1:1
        if i~=Nt
           rx=Rs(i,i+1:end)*Xs_hat(i+1:end,K);
           Xs_slsh=(Z(i)-rx)/Rs(i,i);
           Xs_hat(i,K)=base_mod(base_demod(Xs_slsh,Modulation_Order),Modulation_Order);
         else
           Xs_slsh=Z(i)/Rs(i,i);
           Xs_hat(i,K)=base_mod(base_demod(Xs_slsh,Modulation_Order),Modulation_Order); 
         end
     end
     X_hat(:,K)=P*Xs_hat(:,K);
end