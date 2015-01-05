function Ecoli=Division_CheckEnd(Ecoli,framestate,bacno,kymo);
%Check if measurement ended
[r,~]=size(kymo);
    if r==framestate.frame 
        Ecoli(bacno).fate='endofMeas'; 
       
    end
end