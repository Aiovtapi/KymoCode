function Removespotsui(init,celli,Rspot)

    rspots = size(Rspot,1);
    chans = unique(Rspot(:,1));
    chani = 1;
    spoti = 1;
    ubound = 1;
    
    while chani <= numel(chans);
        
        matpath = strcat(init.bacpath,init.flimgname{chans(chani)},...
            init.OSslash,'Results',init.OSslash,'Cell_',num2str(celli,'%03.0f'));
        load(matpath,'bx','ld');

        while Rspot(spoti,1) == chans(chani) && ubound
            bx{Rspot(spoti,2)} = [];
            ld{Rspot(spoti,2)} = [];
            
            if spoti < rspots
                spoti = spoti + 1;
            else
                ubound = 0;
            end
        end
        
        rbx = bx(~cellfun('isempty',bx));
        rld = ld(~cellfun('isempty',bx));
        
        save(matpath,'rbx','rld','-append')
        chani = chani + 1;
    end 
end
