function [X_hat] = ZF_OSIC(FFT_Size,Modulation_Order,Nt,H,Y);
for K=1:FFT_Size
    H_ch=transpose(H(:,:,K));
    temp_idx=[1:Nt];
    for i=1:Nt
        G_zf=inv(H_ch'*H_ch)*H_ch';
        norm_g_zf=sum(G_zf.*conj(G_zf),2);
        [~,min_idx_zf(i)]=min(norm_g_zf);
        min_idx_zf_temp(i)=temp_idx(min_idx_zf(i));
        temp_idx(min_idx_zf(i))=[];
        X_slsh_zf=G_zf(min_idx_zf(i),:)*Y(:,K);
        X_hat(min_idx_zf_temp(i),K)=base_mod(base_demod(X_slsh_zf,Modulation_Order),Modulation_Order);
        Y(:,K)=Y(:,K)-H_ch(:,min_idx_zf(i))*X_hat(min_idx_zf_temp(i),K);
        H_ch(:,min_idx_zf(i))=[];
    end
end