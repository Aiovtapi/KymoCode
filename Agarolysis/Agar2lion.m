function lionval = Agar2lion(init)

lionval.viewchan = init.viewchan;
lionval.cropindx = init.lioncropindex;
lionval.difchan = init.difchan;
lionval.bacfolder = init.bacpath;
lionval.OSslash = init.OSslash;
lionval.Mainfolder = strcat(init.bacpath,init.OSslash,init.flimgname);
lionval.diffolder = strcat(init.bacpath,init.OSslash,init.difimgname);