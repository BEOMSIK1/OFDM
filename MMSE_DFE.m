function [X_hat] = MMSE_DFE (FFT_Size,Modulation_Order,Nt,Nr,H,Y,N_0)


for K=1:FFT_Size
    H_ch=transpose(H(:,:,K));
    H_bar=[H_ch;N_0*eye(Nt)];
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