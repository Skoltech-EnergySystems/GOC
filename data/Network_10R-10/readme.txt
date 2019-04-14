Files used by more than one scenario of this dataset are kept in the root folder. 
Files unique to a scenario are kept in the appropriate scenario_n folder.
The file inputfiles.ini is a text file that indicates where each of the four file types may be found.
The beginning of the file is marked with [INPUTS].
A period (.) after the file type indicates the file is in the root folder.
An x after the file type indicates the file is in the scenario folder(s).
RAW files are, by definition, unique for each scenario.
CON files typically vary somewhat for each scenario but are sometimes the same for all scenarios.
ROP files typically are the same for all scenarios, as are the INL files.
The using inputfiles.ini files avoids redundant copies of files.

This is a Network_10 dataset with the following configuration:
			  buses	7977
			  loads	6671
		   fixed_shunts	   0
		     generators	 610
        nontransformer_branches	8909
		   transformers	2792
    			  areas	   1
	        switched_shunts	 459
 	  generator inl records	 610
     generator dispatch records	 610
  active power dispatch records	 610
piecewise linear cost functions	 610
		      scenarios  10

The following gives the number of contingencies for each scenario along with		      
the Slack Objective and the First Case Objective. 

The Slack Objective is the cost when satisfying all net load by the slack variables directly
(see Eqn. 3 in the Scoring document). No starting point information is used. 

The FirstCase Objective is computed by projecting the starting point onto 
feasible bounds and creating solution files with Delta=0 for each contingency.

All solutions are (hard) feasible but not optimal.

The Real-Time FirstCase Objectives are generally lower (not always) than the corresponding Offline values, 
which are lower than the Slack Objectives. The Network Model Score is the geometric mean of the scenario objectives.

Created On		2019 Mar 20 18:19	2019 Mar 21 13:50	2019 Mar 20 20:19	2019 Mar 20 20:19
Submission ID		30-1553131155		30-1553201448		(30-1553138383)		(30-1553138383)
Network_10*-010-U008k	raw file: any		Offline			Real-Time		Witness
scenario contingencies	Slack Objective		First Case Objective	First Case Objective	First Case Objective
1	 2402	 	270,732,538,720.87 	 249,045,409,083.93 	 46,794,280.91 	 	47,590,664.40 
2	 2289	 	400,820,058,527.55 	 382,569,962,450.09 	 85,008,630.41 	 	82,042,334.55 
3	  629	 	 95,223,785,691.82 	 27,505,726,740.23 	 35,719,475.68 	 	34,459,923.73 
4	  730	 	200,222,350,840.66 	 183,361,120,894.37 	 50,471,477.10 	 	45,892,092.76 
5	 2335	 	280,003,972,215.62 	 260,820,629,091.45 	 92,028,108.90 	 	86,821,678.37 
6	  641	 	442,297,076,242.92 	 286,006,342,969.96 	 33,762,791.32 	 	34,165,061.81 
7	 2319	 	 52,468,593,721.19 	 30,613,588,287.05 	 57,939,878.84 	 	 9,496,295.63 
8	 2242	 	204,108,450,201.21 	 181,403,947,857.46 	 70,551,684.12 	 	66,307,398.88 
9	  655	 	 48,301,178,050.65 	 26,309,102,530.29    1,051,089,962.77 	     1,002,655,366.96 
10	  677	 	 48,182,442,967.70 	 26,313,483,672.35 	 37,349,767.66 	 	36,662,436.88 
Network Model Score	151,538,614,132.99 	 103,259,085,497.96 	 71,670,073.37 	  	69,570,436.29 
