function DataStruct = RemoveCells(init, DataStruct, cells,faultycells,fpath)

    goodcells = 1:cells;
    goodcells(faultycells) = [];
    newlabel = 1:numel(goodcells);
    chans = numel(init.flimgname);

    removecells = questdlg('Remove discarded cells?','Menu','Yes','No','No');
    switch removecells
        case 'Yes'
            
            load(fpath,'fremoved')
            if fremoved == 0; % Check if faulty cells are already removed
                renumber = questdlg('Renumber bacpics and results?','Menu','Yes','No','Yes');

                for chan = 1:chans
                    bacfolder = strcat(init.bacpath,init.flimgname{chan},init.OSslash);
%                     resultfolder = strcat(bacfolder,'Results',init.OSslash);
                    
                    disp(['Deleting faulty bacpics and results for ',init.flimgname{chan}]);
                    % Delete faulty bacpics and results
                    for fcelli = faultycells;
                        bacdir = strcat(bacfolder,'Cell_',num2str(fcelli,'%03.0f'));
                        rmdir(bacdir,'s');
                        
                        DataStruct(:,fcelli) = [];

%                         matpath = strcat(resultfolder,'Cell_',num2str(fcelli,'%03.0f'),'.mat');
%                         delete(matpath);
                    end

                    % Relabel folders and results if requested
                    switch renumber
                        case 'Yes'
                            disp('Renaming folders and results...')
                            
                            for celli = newlabel;
                                obacdir = strcat(bacfolder,'Cell_',num2str(goodcells(celli),'%03.0f'));
                                nbacdir = strcat(bacfolder,'Cell_',num2str(celli,'%03.0f'));
%                                 omatpath = strcat(resultfolder,'Cell_',num2str(goodcells(celli),'%03.0f'),'.mat');
%                                 nmatpath = strcat(resultfolder,'Cell_',num2str(celli,'%03.0f'),'.mat');
                                
                                if ~strcmp(obacdir,nbacdir)
                                    movefile(obacdir,nbacdir)
                                end
%                                 if ~strcmp(omatpath,nmatpath)
%                                     movefile(omatpath,nmatpath)
%                                 end
                            end
                        case 'No'
                    end    
                end

                fremoved = 1;
                save(fpath,'fremoved','-append')
            else % Cells are already removed
                msgbox('Faulty cells and results have already been removed','Error','error');
            end
        case 'No'
    end
end