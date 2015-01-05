%ShowData


%ReplicationCluster = 
% 
% struct array with fields:
%     name
%     frames
%     clickpos
%     trackpos
%     fate
%     spot1pos
%     spot2pos
%     firstpos
%     lastpos
%     firstframe
%     lastframe

lR=length(ReplicationCluster);
for i=1:lR
    i
    disp('Name:')
    ReplicationCluster(i).name
    disp('Fate:')
    ReplicationCluster(i).fate
    ReplicationCluster(i).firstframe
    ReplicationCluster(i).lastframe
    reply = input('Do you want more? Y/N [Y]: ', 's');
end

% Bax = struct2cell(Division);
% fate=squeeze(Bax(5,:,:));   
% subset=find(strcmp(fate, 'nonexistent')~=1);  %keep indices of present bacteria
% Division=Division(subset);
% [le,~]=size(Division);
% figure;