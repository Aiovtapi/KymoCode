%Structure of Databases for bacterial analysis; Jacob Kerssemakers,
%TNW-BN-ND lab 2012; developed for Charl Moolman

%Database 'M' Manual user inputs
M(1).channels
    %       initval: (structure) inititialization values 
    %     endpoints: coordinates of channel extrema in raw data
    %       presets: settings used for making channel images
    %     DivClicks: (structure) manual inputs for division cycles
    %     RepClicks: (structure) manual inputs for replication cycles
    
    %Substructures
M(1).channels.initval            
    %     driftfile: text file containing drift position data
    %                sumHW: ???
    %            spotRoiHW: 2 halfwidth of area using for summing intnesities around spots
    %           basepath: source path data
    %              BFdatapath: subdirctory name brightfield images
    %          BFfiletemplate: brightfield images template name  
    %              FLdatapath: subdirctory name fluorescence images
    %          FLfiletemplate: fluorescence images template name  
            %maxfile: number of images, and thus kymograph lines
    %         channeldistance: (31.1) distance channels for auto-selection
    %               channelno: number of channels (1 for single)
    %              perpadjust: (1) adjust for tilt along channel entrances
    %            kymoangle: ) angle of channels in image
    %           kymolength: length of kymograph line, pixels (140)
    %           kymohwidth: halfwidth of area around kymograph line (12)
    %    
    %       entranceoffset: shift aligned entrance of channel
    %                  nms: channel names, used for lableing channel databases
    %              outname: label used database storage
    %          outname_usr: label used for manual input database (stored separately)
    %            extension: number of frames before initiation and after termination
    %     WorkspaceOutName: name under which data for one channel is stored, e.g. 'Exp2Chan_x103.mat'
  

 M(1).channels.presets %(used only inside code, stored for repeats)
%     twopoints: begin and end coordinates (2x xy) of channel in image
%     type: 'BF' for brightfield, 'FL' for fluorescence
%     adjustxy: 1 if channel start points should be re-adjusted automatically
%     showmap: 1 to show resampling grid of channel
%     storeref: 1 (internal) to store image for later use in code 
%     useref: 1 (internal) to store image for later use in code
%     refpic: example picture of channel

 M(1).channels. DivClicks
%     name (genealogic number; 2^n indicates mother cell) 
%     fate  'divided' ; 'exit' ;'nonexistent' etc.
%     linkedrep: label of associated replication cluster cycle
%     PosClick: clicked coordinates (in time and  x-position along channel, pixels)      
        %       firstfr : estimate of start of bacterial cycle, frames
        %       firstleft: leftmost position at start (pixels)
        %       firstright: leftmost position at start (pixels)
        %       lastfr: estimate of end of bacterial cycle, frames
        %       lastleft: leftmost position at end (pixels)
        %       lastright: leftmost position at right (pixels)
%     accepted  : 1 if this bacterium is considered suitable data

 M(1).channels. RepClicks  %fields:
%   name (genealogic number; 2^n indicates mother cell) 
%   fate  'disassembled' ; 'exit' ;'nonexistent' etc.
%   linkedbac: label of associated division cycle
%   PosClick
%               firstframe: estimate of initiation of replication , frames
%               lastframe:  estimate of termination of replication , frames
%               firstpos: position at start (pixels)
%               lastpos: position at send (pixels)   
S(1).channels
%     initval:  ;   see above
%     DivClicks:  ; see above
%     RepClicks: ;  see above
%     chanstk_BF: image strip of the channel, brightfield
%     chanstk_FL: image strip of the channel, fluorescence
%     kymo_FL: fluorescence kymograph 
%     kymo_BF: brightfield kymograph
%     Autodivision: division data from brightfiled analysis, see below
%     ReplicationCluster: replication data from brightfiled analysis, see below
 S.channels.AutoDivision
%     family
%     Forks
%     divtype
%     divisiontime
%     accepted
%     birthtime
%     birthtype
%     edges



%SubStructures of 'S(1).channels'
S(1).channels.ReplicationCluster
%     PosKyTracCom: (struct) tracking results from first kymograph analysis (see below)
%     PosKyTracGauss: (struct)tracking results from first kymograph analysis (see below)
%     linkedbac: name of corresponding division cycle
%     FluoPropsGen: (struct) various general counts of fluorescence signal(see below)
%     Pos2DPreTrac:(struct) 1D Gauss tracking results (see below)
%     Pos2DFinTrac:(struct) 2D Gauss tracking results (see below)
%     Cycle: (struct) results of automatic refinement of inititiation and termination times   

S(1).channels.ReplicationCluster.PosKyTracCom %tracking results from first kymograph analysis
%       frames: frame numbers
%       clickpos: positions interpolated between manual start and stop
%       trackpos: positions from 1D COM tracking
%       frames_ext: extended frames (few frames before init and after ter)
%       trackpos_ext: positions from 1D COM tracking; (few frames before init and after ter)
S(1).channels.ReplicationCluster.PosKyTracGauss
% PosKyTracGauss
%     spot1pos: Gauss 1 position, pixels
%     spot1amp: Gauss 1 peak height, a.u.
%     spot1sig: common Gauss sigma (user set), pixels
%     spot2pos: Gauss 2 position, pixels
%     spot2amp: Gauss 2 peak height, a.u.
%     spot2sig: common Gauss sigma (user set), pixels 
%     background: background level (subtracted for fit)
S(1).channels.ReplicationCluster.linkedbac : %name of corresponding division cycle
S(1).channels.ReplicationCluster.FluoPropsGen
% general props of fluorescence(FL) signal
%     frs_ext: [1x44 double]  
%     frs: [34x1 double]    frames
%     darklevel: background level just outside bacterium, per frame
%     darklevel_ext: same, extended before init and after ter
%     signalcontent: total FL counts for this bacterium, per
%     frame; minus background
%     signalcontent_ext: same, extended before init and after ter
%     signalpeakval: maximum FL value found in bacterium
%     signalpeakval_ext: same, extended before init and after ter
%     fluospinelevel: median of maxima along bacterium length
%     fluospinelevel_ext: same, extended before init and after ter
%     fluospotscontent: sum of all values above above spine level
%     fluospotscontent_ext: same, extended before init and after ter
%     areasumI1: via Gauss position fits: summed area around spot 1
%     areasumI2: same, spot 2
%     areasumIall: same, all
S(1).channels.ReplicationCluster.Pos2DPreTrac
%Pre-tracking: results of 1D - two-Gauss fit on averages profile
%Note: the 'y' positions are obtained by a 1D single Gauss fit on 
%crosssections taken from the image at the respective X-positions from
%above fit.
%     X0: Spot 1 X-position, pixels
%     X1: Spot 2 X-position, pixels
%     Y0: Spot 1 Y-position, pixels
%     Y1: Spot 2 Y-position, pixels
%     Bck: Background level
%     Pk1: Spot 1 Peak, a.u.
%     Pk2: Spot 2 Peak, a.u.
%     spots: spot number present, according to rejection criteria
S(1).channels.ReplicationCluster.Pos2DFinTrac
%Final tracking: results of 2D, double Gaussian spot fitting using the
%pre-tracking results as first estimates
%     X0: Spot 1 X-position, pixels
%     X1: Spot 2 X-position, pixels
%     Y0: Spot 1 Y-position, pixels
%     Y1: Spot 2 Y-position, pixels
%     Bck: Background level (peak of wide Gaussian) 
%     Pk1: Spot 1 Peak, a.u.
%     Pk2: Spot 2 Peak, a.u.
%     spots: spot number present, according to rejection criteria
S(1).channels.ReplicationCluster.Cycle
    %     Cycle: (struct) results of automatic refinement of inititiation and termination times    
%           normtime_ext: frame numbers normalized to init(=0) and ter(=1) times
%           Init: (stuct)
%                 idx: relative index of initiation (frames)
%                 time:absolute frame number
%           Ter: 
%                 idx: relative index of initiation (frames)
%                 time:absolute frame number


% Fields as of 9-12-2013 -JK; See above for descriptions
 S.channels.RepClicks
%     name
%     fate
%     PosClick
 S.channels.AutoDivision
%     family
%     Forks
%     divtype
%     divisiontime
%     accepted
%     birthtime
%     birthtype
%     edges
 S.channels.ReplicationCluster
%     PosKyTracCom
%     Forks
%     FluoPropsGen
%     Cycle
%     Pos2DPreTrac
%     Pos2DFinTrac
 S(1).channels(1).ReplicationCluster(1).PosKyTracCom
%     frames: [1x34 double]
%     clickpos: [1x34 double]
%     trackpos: [34x1 double]
%     frames_ext: [44x1 double]
%     trackpos_ext: [44x1 double]
 S(1).channels(1).ReplicationCluster(1).Forks
%     left: [19x2 double]
%     right: [17x2 double]
 S(1).channels(1).ReplicationCluster(1).FluoPropsGen
%     frs_ext: [1x44 double]
%     frs: [34x1 double]
%     darklevel: [1x34 double]
%     darklevel_ext: [1x44 double]
%     signalcontent: [1x34 double]
%     signalcontent_ext: [1x44 double]
%     fluosumspinelevel_ext: [1x44 double]
%     signalpeakval: [1x34 double]
%     signalpeakval_ext: [1x44 double]
%     fluospinelevel: [1x34 double]
%     fluospinelevel_ext: [1x44 double]
%     fluospotscontent: [1x34 double]
%     fluospotscontent_ext: [1x44 double]
%     fluospotscontentpadd_ext: [1x44 double]
%     areasumI1: [44x1 double]
%     areasumI2: [44x1 double]
%     areasumIall: [44x1 double]
 S(1).channels(1).ReplicationCluster(1).Cycle
%     normtime_ext: [1x44 double]
%     Init: [1x1 struct]
%     Ter: [1x1 struct]
 S(1).channels(1).ReplicationCluster(1).Pos2DPreTrac
%     X0: [1x44 double]
%     X1: [1x44 double]
%     Y0: [1x44 double]
%     Y1: [1x44 double]
%     contentallspots: [1x44 double]
%     contentspot1: [1x44 double]
%     contentspot2: [1x44 double]
%     spot1OK: [1x44 double]
%     spot2OK: [1x44 double]
 S(1).channels(1).ReplicationCluster(1).Pos2DPreTrac  : Same

