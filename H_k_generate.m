function H_k = H_k_generate(K,M,N,h_k,f_k,G,Theta)
H_k=zeros(K,M);
for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);
    H_k(k,:)=h_k_temp+G'*Theta*f_k_temp;
end

end

