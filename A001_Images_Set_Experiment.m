function initval=A001_Images_Set_Experiment(exp)
%This function stores paths etc. of various bacteria experiments. Kers2012

%common settings---------------------------------------------
user='MarkPC';
switch user
      case 'Jacob', 
          projectpath='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\';
          versionpath='2014_01 DnaNDnaX DualColor';
          toolspath='\SupportingFunctions';          
          initval.skippicturesavingbecauseCharlislogginghislive=0;
      case 'Sriram', 
          projectpath='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\';
          versionpath='2014_01 DnaNDnaX DualColor';
          toolspath='\SupportingFunctions';          
          initval.skippicturesavingbecauseCharlislogginghislive=0;
      case 'Roy' ,
          projectpath='/Users/rleeuw/Work/DataAnalysis/';
          versionpath='201511_TUSdifDnaN_Montage/';
          toolspath='SupportingFunctions';    
      case 'RoyPC' ,
          projectpath='D:\rleeuw\DataAnalysis\';
          versionpath='201511_TUSdifDnaN_Montage\';
          toolspath='SupportingFunctions';
      case 'Peter' ,
          projectpath='D:\peterbrazda\ImageAnalysis\';
          versionpath='KymoCode\';
          toolspath='SupportingFunctions';
      case 'MarkPC'
          projectpath='D:\Users\water_000\Documents\GitHub\';
          versionpath='KymoCode\';
          toolspath='SupportingFunctions';
      case 'Mark'
          projectpath='D:\Users\water\Documents\GitHub\';
          versionpath='KymoCode\';
          toolspath='SupportingFunctions';
end
addpath(strcat(projectpath,versionpath));
addpath(strcat(projectpath,versionpath,toolspath));

initval.plotintermediateresults=0;
% experiment-specific settings-------------------------------------
switch exp    
   case   '001_DnaN_TUS_dif_30122014_DnaNsignal',            initval=Exp001_DnaN_TUS_dif_30122014_DnaNsignal(initval, user);
   case   '001_DnaN_TUS_dif_30122014_TUSsignal',             initval=Exp001_DnaN_TUS_dif_30122014_TUSsignal(initval, user);
   case   '001_DnaN_TUS_dif_30122014_difsignal',             initval=Exp001_DnaN_TUS_dif_30122014_difsignal(initval, user);
   case   '002_oriZ_dif_2304_DnaNsignal',                    initval=Exp002_oriZ_dif_23042015_DnaNsignal(initval, user);
   case   '002_oriZ_dif_2304_R2signal',                      initval=Exp002_oriZ_dif_23042015_R2signal(initval, user);

   case   '005_oriZ_dif_110515_R2signal',                    initval=Exp005_oriZ_dif_11052015_R2signal(initval, user);
   case   '005_oriZ_dif_110515_DnaNsignal',                  initval=Exp005_oriZ_dif_11052015_DnaNsignal(initval, user);
       
   case   '006_oriZ_dif_110515_R2signal',                    initval=Exp006_oriZ_dif_11052015_R2signal(initval, user);
   case   '006_oriZ_dif_110515_DnaNsignal',                  initval=Exp006_oriZ_dif_11052015_DnaNsignal(initval, user);
       
   case   '105_oriZ_dTus_230615_R2signal',                   initval=Exp105_oriZ_dTus_23062015_R2signal(initval, user);
   case   '105_oriZ_dTus_230615_DnaNsignal',                 initval=Exp105_oriZ_dTus_23062015_DnaNsignal(initval, user); 
     
   case   '001_DnaN_TUS_dif_30122014_DnaNsignal_P',          initval=Exp001_DnaN_TUS_dif_30122014_DnaNsignal_P(initval, user);

end

