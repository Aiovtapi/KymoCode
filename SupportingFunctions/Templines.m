%A001
zoom: 200


%A050_line 26 
load(Channelpath, 'chanstk_BF','chanstk_FL','endpoints', 'kymo_BF','kymo_FL','presets');
%line 55
save(outname, 'Division', 'DivClicks','-append');