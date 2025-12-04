function [H] = channel_H(N,M,dis,large_fading,f_c)
% N number of receiver
% M number of transmitter
f_c=f_c*10^9;
Lambda=3*10^8/f_c;

H = zeros(N,M);
for aa=1:N
    for bb=1:M
       H(aa,bb) =( randn()+1j*randn())/sqrt(2);
    end
end

H = sqrt(1/2)*exp(1j*2*pi*rand(N,M))+sqrt(1/2)*H;

a = Lambda/4/pi/dis;

a = 10^(-3.73)/(dis^large_fading);
a= sqrt(a);

H = a*H;
end
