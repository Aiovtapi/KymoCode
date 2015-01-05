% 2D Gauss fit setup

%Need: division data base; Replication databse; original image stack

%start with timeline selction: pick a division/replication cycle

%2) use  1D division coordinates to define bactrium ROI:
%calcualte back to real image coordinates (start-stop)
%use existing 'image channel map' function using these two coordinates to
%fetch fluorescent area

%option: larger roi, draw ellips for bacterium area

%1)clean 'left spot-right spot': one or two spots
%calc back to ROI image coords, use these for 2x2D fit
