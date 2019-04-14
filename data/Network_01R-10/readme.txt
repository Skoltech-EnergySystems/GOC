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

This is a GO Competition Challenge 1 Network_01 dataset with the following characteristics.
			  buses	500
			  loads	200
		   fixed_shunts	  0
		     generators	 90
        nontransformer_branches	468
		   transformers	131
    			  areas	  1
	        switched_shunts	 17
 	  generator inl records	 90
     generator dispatch records	 90
  active power dispatch records	 90
piecewise linear cost functions	 90
		      scenarios  10
