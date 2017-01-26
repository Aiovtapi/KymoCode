function TitlePos = GetTitPos(Tx,Ty)

XLim = get(gca,'XLim');
YLim = get(gca,'YLim');

Xsize = abs(XLim(2)-XLim(1));
Ysize = abs(YLim(2)-YLim(1));

Xpos = XLim(1)+Tx*Xsize;
Ypos = YLim(2)+Ty*Ysize;

TitlePos=[Xpos,Ypos,0];