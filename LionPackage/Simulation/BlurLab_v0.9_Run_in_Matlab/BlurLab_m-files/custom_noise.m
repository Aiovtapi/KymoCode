function T=custom_noise(P,N,M)

%normalize P
Pnorm=[0 P]/sum(P);

%create cumlative distribution
Pcum=cumsum(Pnorm);

%create random matrix
R=rand(N,M);

%calculate T
T=zeros(N,M);
for i=1:length(P)
    T(and(R>Pcum(i),R<=Pcum(i+1)))=i;
end
