function DataStruct = ViewBacNspot(DataStruct,celli,Nspot,Bacmesh)

    nspots = size(Nspot,1);

    for spoti = 1:nspots;
        thischan = Nspot(spoti,9);
        frami = Nspot(spoti,10);
        
        obx = DataStruct(thischan,celli).bx;
        old = DataStruct(thischan,celli).ld;
        
        bx = Nspot(spoti,1:8);
        ld = Nspot(spoti,1:8);
        
        Xval = bx(2);
        Yval = bx(4);
        thisbacmesh = Bacmesh{thischan}{celli,frami};

        [Lval,Dval] = projectToMesh(Xval,Yval,thisbacmesh);
        varval = sqrt(bx(3)^2+bx(5)^2);

        ld(2) = Lval;
        ld(4) = Dval;
        ld(3) = varval;
        ld(5) = varval;
        
        nbx = [obx, bx];
        nld = [old, ld];
        
        DataStruct(thischan,celli).bx = nbx;
        DataStruct(thischan,celli).ld = nld;

    end 
end
