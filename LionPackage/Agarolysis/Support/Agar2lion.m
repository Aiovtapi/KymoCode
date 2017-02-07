function lionval = Agar2lion(init,chan)

lionval.cropindx = init.lioncropindex;
lionval.difchan = init.difchan;
lionval.bacfolder = init.bacpath;
lionval.OSslash = init.OSslash;
lionval.Mainfolder = strcat(init.bacpath,init.OSslash,init.flimgname{chan});
lionval.diffolder = strcat(init.bacpath,init.OSslash,init.difimgname);
lionval.datapath = init.datapath;