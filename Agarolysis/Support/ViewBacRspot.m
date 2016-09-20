function DataStruct = ViewBacRspot(DataStruct,celli,Rspot)

    rspots = size(Rspot,1);

    for spoti = fliplr(1:rspots);
        thischan = Rspot(spoti,1);
        
        bx = DataStruct(thischan,celli).bx;
        ld = DataStruct(thischan,celli).ld;

        bx{Rspot(spoti,2)} = [];
        ld{Rspot(spoti,2)} = [];
        rbx = bx(~cellfun('isempty',bx));
        rld = ld(~cellfun('isempty',bx));
        
        DataStruct(thischan,celli).bx = rbx;
        DataStruct(thischan,celli).ld = rld;

    end 
end
