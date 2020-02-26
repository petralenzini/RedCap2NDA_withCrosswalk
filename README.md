# RedCap2NDA_withCrosswalk
Use this repository to map HCPA Lifespan data stored in local Redcap Databases (Moises, you'll have to rewrite this stuff to incorporate your yaml and secrets additions...via the config.py library you created in ccf-behavioral and then tell me to destroy my box configuration and redcap config files and create my own secrets file somewhere in my local git clone of this stuff) to validated (but not yet uploaded) NDA structures.

Note the crosswalk is specific to HCP variables encoded in OUR redcap databases, but once you have a crosswalk for YOUR data, you should be able to repurpose the code (start by comparing the overlap between the two apply crosswalk programs...in the end, their only difference should be specialty code that can't be covered by stuff that can be put into a library that is common to all of these behavioral data mapping projects (e.g. box.py, or redcap.py, and maybe a new one called crosswalk.py, perhaps)

THIS IS A WORK IN PROGRESS.

For the time being, applyHCPAcrosswalk.py takes as input HCPA_crosswalk_concatenated_9Jan2020.csv as it's input and produces all the structures in the output path hardcoded in its guts...it assumes you have access to the keys and configurations I've stored elsewhere for Box and Redcap.  This setup is antiquated...if you're not Moises and really want to use this method, create an issue.  If you are MOISES, please replace affected with the stuff necessary to do it your (better) way using config.py from ccf-behavioral.


