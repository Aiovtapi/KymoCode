function [Bacmask, CBacmask, Bacpics, NMBacpics] = ...
    AgarRotate2(init, Bacmask, CBacmask, Bacpics, NMBacpics, cells, frames, imgflip, chan)

for celli = 1:cells;
    if imgflip(celli) == 1;
        for frami = 1:frames;

            % Rotate Bacpics and masks
            Bacmask{celli,frami} = rot90(Bacmask{celli,frami},2);
            CBacmask{celli,frami} = rot90(CBacmask{celli,frami},2);
            Bacpics{celli,frami} = rot90(Bacpics{celli,frami},2);
            NMBacpics{celli,frami} = rot90(NMBacpics{celli,frami},2);

            % Save rotated bacpic to folder
            bacfolder = strcat(init.bacpath,init.flimgname{chan});
            bacpath=strcat(bacfolder,init.OSslash,'Cell_',num2str(celli,'%03.0f'),init.OSslash);
            thisbacpath = strcat(bacpath,num2str(frami,'%03.0f'),'.tif');
            imwrite(uint16(Bacpics{celli,frami}),thisbacpath)
            
        end
    end
end