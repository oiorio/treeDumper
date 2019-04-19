#for MC 2016 94X
#python submit_all.py -f DAS_Names_2017_bkg.txt -c topplusdmTrees_cfg.py -p version='94X_2016' -d crab_4Oct
#for Data 2016 94X
#python submit_all.py -f DAS_Names_2017_DATA.txt -c topplusdmTrees_cfg.py -p version='94X_2016' -d crab_4Oct --isData

#For Data 2017 94X
#python submit_all.py -f Tprime_LVIGI.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_30Jan
#python submit_all.py -f DAS_Names_MET_2017_94X.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_30Jan --isData
#python submit_all.py -f DAS_Names_MET_part_2017_94X.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_30Jan --isData 

#python submit_all.py -f MC_2016_94X_ST.txt -c topplusdmTrees_cfg.py -p version='94X_2016' -d crab_14March

#python submit_all.py -f Data_Muele_2017_94X.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_30Jan --isData 
python submit_all.py -f DAS_Names_WprimeLep.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_11Apr 
#python submit_all.py -f DAS_Names_MET_part_2017_94X.txt -c topplusdmTrees_cfg.py -p version='94X' -d crab_30Jan --isData
#For a dry run use the -n 1 option
#python submit_all.py -n 1 -f DAS_Names_2017_DATA.txt -c topplusdmTrees_cfg.py -p version='94X_2016' -d crab_4Oct --isData
