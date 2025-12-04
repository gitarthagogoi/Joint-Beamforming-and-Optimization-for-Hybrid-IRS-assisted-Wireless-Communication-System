function [Rsum,gamma_k] = SINR_calculate_SI(K,M,N,H_k,h_k,f_k,G,theta,w_k,sigma2,sigmar2,Pr_max)

gamma_k=zeros(K,1);

U=zeros(N,N);
for k=1:K
    w_k_temp=reshape(w_k(k,:),M,1);
    U=U+diag(G*w_k_temp)*(diag(G*w_k_temp))';
end
U=U+sigmar2*eye(N);

for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);

   theta_temp = theta+diag(theta)*(reshape(H_k(k,:,:),N,N))'*theta;

   theta_temp = min(sqrt(Pr_max/(theta_temp'*U*theta_temp)),1)*theta_temp;

   temp1=h_k_temp+G'*diag(f_k_temp)*theta_temp;
    
    temp2=reshape(w_k(k,:),M,1);
    temp3=temp1'*temp2;
    
    temp4=0;
    for j=1:K
        w_j=reshape(w_k(j,:),M,1);
        temp4=temp4+norm(temp1'*w_j)^2;
    end
    temp4=temp4-norm(temp3)^2+norm(theta_temp'*diag(f_k_temp))^2*sigmar2+sigma2;    
    gamma_k(k)=norm(temp3)^2/temp4;
end
Rsum=0;
for k=1:K
    Rsum=Rsum+log2(1+gamma_k(k));
end
end

