#voms-proxy-init --voms cms
cd ../../../../../../CMSSW_8_1_0_pre12/src/
cmsenv
cd -
gfal-ls -Hl srm://stormfe1.pi.infn.it/cms/store/user/oiorio/SingleTop/
# latest
#gfal-ls -Hl srm:////stormfe1.pi.infn.it/cms/store/user/oiorio/SingleTop/2018/Apr/8Apr 
#Zurich server:
#gfal-ls -Hl srm://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2018/May28
