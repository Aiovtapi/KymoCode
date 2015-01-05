
function profile=Get_averagedZ_profile_useLocalInfo(Stack,pos,initval)
%This function gets a Z-profile in a stack around a pre-set XY location, by
%averaging over a little square per plane

	hsq=initval.square;  %This is used to find the particel max
    hsqL=20;           %This is used to find and subtract the local background value
	[r,c,p]=size(Stack);
	profile=zeros(p,1);
	
    hix=min([pos(2)+hsq c]); lox=max([pos(2)-hsq 1]);         %define little averaging window; check borders
	hiy=min([pos(1)+hsq r]); loy=max([pos(1)-hsq 1]);
	
    hixL=min([pos(2)+hsqL c]); loxL=max([pos(2)-hsqL 1]);    %define large background window; check borders
	hiyL=min([pos(1)+hsqL r]); loyL=max([pos(1)-hsqL 1]);
    
    for i=1:p
        squ=Stack(lox:hix,loy:hiy,i);                       %max-find window
        squL=Stack(loxL:hixL,loyL:hiyL,i);                  %background-find-window
        profile(i)=max(max(squ))-median(median(squL));
	end
	dum=1;
    lastval=nanmean(profile(p-20:p));
    %profile=(profile-median(profile));
    profile=profile-lastval;
end