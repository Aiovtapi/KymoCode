function Ecoli=Division_CheckExit(Ecoli,framestate,bacprops,bacno);
%Check if bacterium is too close to exit for proper measurements
    switch bacprops.whenexit
        case 'right'
                    Rightposses=Ecoli(bacno).right;
                    lp=length(Rightposses); 
                    lastpos_R=Rightposses(lp);
                    if lastpos_R>=bacprops.exitpos, 
                        Ecoli(bacno).fate='exit'; 
                    end
        case 'left'
             leftposses=Ecoli(bacno).left;
                    lp=length(leftposses); 
                    lastpos_L=leftposses(lp);
                    if lastpos_L>=bacprops.exitpos, 
                        Ecoli(bacno).fate='exit'; 
                    end
    end
end