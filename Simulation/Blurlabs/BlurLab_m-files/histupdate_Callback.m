function histupdate_Callback(hObject, eventdata, handles, Iout)
if get(handles.checkbox_histogram,'Value')
    axes(handles.axes2)
    
    [Fout,binout]=hist(Iout(Iout>=0),200);
    Fout=Fout/(binout(2)-binout(1));
       
    if get(handles.checkbox_log,'Value')
        Fout=log(Fout);
        Fout(isnan(Fout))=0;
        Fout(isinf(Fout))=0;
    end
    
    minx=round(min(binout));
    maxx=round(max(binout));
    delx=maxx-minx;
    
    if minx==maxx
        maxx=minx+1;
    end
    
    %get noise spectrum
    plot(binout,Fout,'b','Linewidth',1)
    set(handles.axes2,'Xlim',[minx,maxx])
    set(handles.axes2,'XTick',[minx,maxx])
    set(handles.axes2,'XTickLabel',[minx,maxx])
    set(handles.axes2,'YTick',[])
    set(handles.axes2,'Color',[0.9 0.9 0.9])
    set(handles.axes2,'Box','on')
end
