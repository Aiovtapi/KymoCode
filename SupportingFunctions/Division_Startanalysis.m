function Division=Division_Startanalysis(Division, bacprops);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate first bacteria (to be replaced by detection, later)

%later to be used in init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single ['1010111'] style numer indicating family tree position
%number length-1= # of prior divisions
            %nth digit: left(1)- or right(0) division in nth division
            %(channel exit is on the far right)
            
%             
% Division = repmat(struct(...
%             'name',0,...            single number indicating family tree position
%             'frames',[1],...         frames            
%             'left',[],...           array containing left edge coordinates from birth until split
%             'right',[],...          array containing right edge coordinates from birth until split
%             'fate','nonexistent',... %'alive' 'exit'  'divided' or 'endofMeas'   
%             'firstfr',[1],...           birth
%             'lastfr',[1],...  %division/end
%             'firstleft',[1],...           birth
%             'firstright',[1],...           birth
%             'lastleft',[1],...  %division/end
%             'lastright',[1]),n_tracks,1);  %division/end
% 
%          %Note: births for this set are put at zero, by definition


for i=1:5
Division(i).fate='alive';
end
for j=6:8, Division(j).fate='exit';  end     %(never seen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

