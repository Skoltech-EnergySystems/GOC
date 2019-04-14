Files used by more than one scenario of this dataset are kept in the root folder. 
Files unique to a scenario are kept in the appropriat scenario_n folder.
The file inputfiles.ini is a text file that indicates where each of the four file types may be found.
The beginning of the inputfiles.ini file is marked with [INPUTS].
A period (.) after the file type indicates the file is in the root folder.
An x after the file type indicates the file is in the scenario folder(s).
RAW files are, by definition, unique for each scenario
CON files typically vary somewhat for each scenario but are sometimes the same for all scenarios.
ROP files typically are the same for all scenarios, as are the INL files.
The inputfiles.ini file avoids redundant copies of files.

This is a GO Competition Challenge 1 Network_07 dataset with the following characteristics.
			  buses	2312
			  loads	1529
		   fixed_shunts	 121
		     generators	 617
	nontransformer_branches	2156
		   transformers	 857
			  areas	  18
		switched_shunts	 201
	  generator inl records	 617
     generator dispatch records	 617
  active power dispatch records	 617
piecewise linear cost functions	 617
                      scenarios   10
