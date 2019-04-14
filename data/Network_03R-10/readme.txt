Files used by more than one scenario of this dataset are kept in the root folder. 
Files unique to a scenario are kept in the appropriat scenario_n folder.
The file inputfiles.ini is a text file that indicates where each of the four file types may be found.
The beginning of the file is marked with [INPUTS].
A period (.) after the file type indicates the file is in the root folder.
An x after the file type indicates the file is in the scenario folder(s).
RAW files are, by definition, unique for each scenario
CON files typically vary somewhat for each scenario but are sometimes the same for all scenarios.
ROP files typically are the same for all scenarios, as are the INL files.
The inputfiles.ini file avoids redundant copies of files.

This is a GO Competition Challenge 1 Network_03 dataset with the following characteristics.
			  buses	793
			  loads	568
		   fixed_shunts	 49
		     generators	210
        nontransformer_branches	769
		   transformers	143
    			  areas	  1
	        switched_shunts	 50
 	  generator inl records	210
     generator dispatch records	210
  active power dispatch records	210
piecewise linear cost functions	210
		      scenarios  96

