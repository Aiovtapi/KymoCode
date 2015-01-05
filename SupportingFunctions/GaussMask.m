function ar=GaussMask(ar,sigma)
%This function masks an array with a Gaussian window;
%sigma=1 means edge is on one sigma etc.
lar=length(ar);
x0=lar/2;
sig=sigma*x0;
ax=linspace(-x0,x0,lar);
ar=ar.*exp(-(ax/sig).^2);
