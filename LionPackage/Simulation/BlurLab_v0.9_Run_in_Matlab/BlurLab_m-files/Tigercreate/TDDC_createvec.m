function TY = TDDC_createvec(x,y,frames)

rx = round(x*(frames-1));

bends = numel(x);
TY = [];

for i = 1:bends - 1
    thisx = rx(i):(rx(i+1)-1);
    xchange = rx(i+1) - rx(i);
    ychange = y(i+1) - y(i);
    co = ychange/xchange;
    localx = thisx - rx(i);
    thisy = y(i) + localx.*co;
    TY = [TY, thisy];
end

TY = [TY, y(end)];