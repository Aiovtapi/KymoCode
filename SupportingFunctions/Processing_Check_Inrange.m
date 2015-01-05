function inrange=Processing_Check_Inrange(repfrs,reppos,frs,lft,rht);
lf=length(repfrs);
overlap=0;
for i=1:lf
    idx=find(frs==repfrs(i));
    if ~isempty(idx);
        cond=reppos(i)>lft(idx)&reppos(i)<rht(idx);
        overlap=overlap+cond;
    end
end
if overlap>1, 
    inrange=1;
else
    inrange=0;
end