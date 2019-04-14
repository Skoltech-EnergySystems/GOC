The file inputfiles.ini is a text file that indicates where each of the four file types may be found.
The beginning of the file is marked with [INPUTS].
A period (.) after the file type indicates the file is in the root folder.
An x after the file type indicates the file is in the scenario folder(s).
RAW files are, by definition, unique for each scenario
CON files typically vary somewhat for each scenario but are sometimes the same for all scenarios.
ROP files typically are the same for all scenarios, as are the INL files.
The inputfiles.ini file avoids redundant copies of files.

This is a Network_05 dataset with the following configuration
			  buses	2000
			  loads	1125
		   fixed_shunts	   0
		     generators	 544
	nontransformer_branches	2359
		   transformers	 847
			  areas	   1
		switched_shunts	 153
 	  generator inl records	 544
     generator dispatch records	 544
  active power dispatch records	 544
piecewise linear cost functions	 544
                      scenarios   10
