function bx = Removespots(x,thisBacmask)

    frames = size(thisBacmask,2);
    spots = size(x,2);

    for spoti = 1:spots
        
        outbound = 0;
        
        for frami = 1:frames;
            if outbound == 0;
                bacsize = size(thisBacmask{frami});

                xval = round(x{spoti}(frami,2));
                yval = round(x{spoti}(frami,4));

                xlval = xval - 1;
                xhval = xval + 1;
                ylval = yval - 1;
                yhval = yval + 1;

                lbound = xlval < 0;
                rbound = xhval > bacsize(2);
                ubound = ylval < 0;
                dbound = yhval > bacsize(1);

                if ubound && lbound
                    box = buildbox(xval,yval,bacsize(2),1);
                elseif ubound && ~(lbound||rbound)
                    box = buildbox(xval,yval,bacsize(2),2);
                elseif ubound && rbound
                    box = buildbox(xval,yval,bacsize(2),3);
                elseif lbound && ~(ubound||dbound)
                    box = buildbox(xval,yval,bacsize(2),4);
                elseif ~(lbound||rbound||ubound||dbound)
                    box = buildbox(xval,yval,bacsize(2),5);
                elseif rbound && ~(ubound||dbound)
                    box = buildbox(xval,yval,bacsize(2),6);
                elseif dbound && lbound
                    box = buildbox(xval,yval,bacsize(2),7);
                elseif dbound && ~(lbound||rbound)
                    box = buildbox(xval,yval,bacsize(2),8);
                elseif dbound && rbound
                    box = buildbox(xval,yval,bacsize(2),9);
                end

                Bacmaskl = thisBacmask{frami}';
                hit = sum(Bacmaskl(box));

                if hit == 0
                    outbound = 1;
                end
            end
        end
        
        if outbound == 1;
            x{spoti} =[];
        end
    end
    bx = x(~cellfun('isempty',x));
end

function box = buildbox(x,y,l,sit)
    switch sit
        case 1
            box = [1,2,l+1,l+2];
        case 2
            box = [x-1,x,x+1,l+x-1,l+x,l+x+1];
        case 3
            box = [x-1,x,l+x-1,l+x];
        case 4
            box = [(y-2)*l+x,(y-2)*l+x+1,(y-1)*l+x,(y-2)*l+x+1,y*l+x,y*l+x+1];
        case 5
            box = [(y-2)*l+x-1,(y-2)*l+x,(y-2)*l+x+1,(y-1)*l+x-1,(y-1)*l+x,(y-1)*l+x+1,y*l+x-1,y*l+x,y*l+x+1];
        case 6
            box = [(y-2)*l+x-1,(y-2)*l+x,(y-1)*l+x-1,(y-2)*l+x,y*l+x-1,y*l+x];
        case 7
            box = [(y-2)*l+x,(y-2)*l+x+1,(y-1)*l+x,(y-1)*l+x+1];
        case 8
            box = [(y-2)*l+x-1,(y-2)*l+x,(y-2)*l+x+1,(y-1)*l+x-1,(y-1)*l+x,(y-1)*l+x+1];
        case 9
            box = [(y-2)*l+x-1,(y-2)*l+x,y*l-1,y*l];
    end
end
        