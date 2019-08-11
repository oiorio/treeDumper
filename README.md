#Welcome to the tree dumper main page! 

The main code takes entuples and produces custom made trees with weights, event-based systematics, and potentially a baseline selection.

0) To download, from the CMSSW_X_Y_Z/src do:

git clone https://github.com/oiorio/treeDumper.git TreeFWK/treesDumper

Please look below to see what is currently missing!

1) For the impatient user:
The main python code is run with:

cmsRun topplusdmTrees_cfg.py changeJECs=True maxEvts=-1 

Options are:

changeJECs : can be True or False default) to change or not the JECs.
maxEvts : can be an integer or -1 (default)
isData : Data if you are running on data, MC if you are running on MC
version : It considers several features that can be adapted depending on the version : 80X (2016 only),94X (2017 or 2016) , 10X. 

To customize by hand: the content in terms of systematifcs. Not implemented unclustered met.

2) To launch on crab:

python submit_all.py -f Files.txt -c topplusdmTrees_cfg.py -d crab_directory 

options are: 

-f Files.txt : a txt file containing the names of the edm-ntuples to use.
-c topplusdmTrees_cfg.py : chooses the python file to run (as in local).
-d crab_directory : the directory to save info in. Usually it's better to put the date of the launch
--jecVersion: you can choose between Fall17 or Summer16 to get the list of files to include depending on what is needed.  
-n 1 use it for dry run (if you want to make some modifications and see how it goes.

For now you have to modify by hand the folder where you want to save!

3) WHAT IS MISSING

3.1) Unclustered MET 
3.2) High Pt electron triggers: just need to look at the correct names. (need ~half an 
3.3) Resolution scale factors to be added automatically: at the moment there are some dummy ones, need to do it customizable from the outside with the version.
3.4) Custom selection of b-tagging working points.
3.5) Custom selection: a selection that could be customizable.

