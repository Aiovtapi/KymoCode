        function LivingEcoli=Division_GetLivingEcoli(Ecoli);
        NOT USED ANYMORE
        %Find and analyse the living bacteria on this time 
        buf = struct2cell(Ecoli);
        fate=squeeze(buf(2,:,:));   
        pcl=find(strcmp(fate, 'alive'));  %find indices of present bacteria  
        LivingEcoli=Ecoli(pcl) ;                    %...and their names
        