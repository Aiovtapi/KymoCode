function [com,comc,moi]=Get_1DCOM(soi);
%Get one-dimensional center of mass (as offset from center)
lp=length(soi);
ax=linspace(-lp/2,lp/2,lp)';
mn=min(min(soi));
if ~isempty(mn),soi2=soi'-mn; else soi2=soi'; end
%background correction
comc=sum(soi2'*ax)/sum(sum(soi2));

com=comc+lp/2;      %(center of mass in array coordinates)

com=max([1 com]); com=min([lp-1 com]); %just to be sure
comc=max([-lp/2 comc]); comc=min([lp/2 comc]); %just to be sure
 moi=sqrt(sum((ax-comc)'.^2*soi2)/sum(sum(soi2)));  %scaled moment of inertia, square root)
 %moi=range(soi2,2)/100;