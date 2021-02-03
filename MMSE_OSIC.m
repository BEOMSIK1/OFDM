function [X_hat] = MMSE_OSIC (FFT_Size,Modulation_Order,Nt,H,Y,N_0);

for K=1:FFT_Size
    H_ch=transpose(H(:,:,K));
    H_bar=[H_ch;N_0*eye(Nt)];
    Y_bar=[Y;zeros(Nt,FFT_Size)];
    temp_idx=[1:Nt];
    for i=1:Nt
        G_mmse=inv(H_bar'*H_bar)*H_bar';
        norm_g_mmse=sum(G_mmse.*conj(G_mmse),2);
        [~,min_idx_mmse(i)]=min(norm_g_mmse);
        min_idx_mmse_temp(i)=temp_idx(min_idx_mmse(i));
        temp_idx(min_idx_mmse(i))=[];
        X_slsh_mmse=G_mmse(min_idx_mmse(i),:)*Y_bar(:,K);
        X_hat(min_idx_mmse_temp(i),K)=base_mod(base_demod(X_slsh_mmse,Modulation_Order),Modulation_Order);
        Y_bar(:,K)=Y_bar(:,K)-H_bar(:,min_idx_mmse(i))*X_hat(min_idx_mmse_temp(i),K);
        H_bar(:,min_idx_mmse(i))=[];
    end
 end