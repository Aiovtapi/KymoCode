function Iout=poissimnoise(Iin,lambda,basal)

Imin=min(min(Iin));
Imax=max(max(Iin));

if Imin<0
    error('All values in input image must be positive.')
end

if lambda<=0
    error('The intensity multiplier must be greater than zero.')
end

if basal<0
    error('The basal level must be greater than or equal to zero.')
end

Iout=zeros(size(Iin));

for i=Imin:Imax
    list1=find(Iin==i);
    if ~isempty(list1)
        rand1=poissrnd(double(i*lambda+basal),1,length(list1));
        Iout(list1)=rand1;
    end
end