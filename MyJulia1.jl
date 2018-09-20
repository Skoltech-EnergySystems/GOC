function MyJulia1(rawFile, genFile, contFile)

    open("solution1.txt", "w") do f1
        write(f1,"--generation dispatch
bus id,unit id,pg(MW),qg(MVar)
1,'1 ',37.964959681024965,1.8481324992322392
2,'1 ',84.64254846706017,17.528955644536275
3,'1 ',5.801787267635664,26.833484287279358
6,'1 ',110.50327593399243,-13.070044636461493
8,'1 ',0.31982767205661244,7.013410487379017
--end of generation dispatch ")
    end

    open("solution2.txt", "w") do f2
        write(f2,"--contingency generator
conID,genID,busID,unitID,q(MW)
1,l_1,1,'1 ',-21.962758384505367
1,l_2,2,'1 ',36.16681774760686
1,l_3,3,'1 ',14.302337555660847
1,l_4,6,'1 ',-3.298754274417449
1,l_5,8,'1 ',19.237282450762788
--end of contingency generator
--bus
contingency id,bus id,v(pu),theta(deg)
0,1,1.0469439305939694,43.754859634437025
0,2,1.0526258496239331,42.92160215291941
0,3,1.0161219135627484,37.959783015393555
0,4,1.0344596288633165,40.47346725377485
0,5,1.0379725556760226,41.83808734394405
0,6,1.0812385054030371,46.09111843093323
0,7,1.069615995361771,40.09509500729
0,8,1.0999999829050686,40.1225297419414
0,9,1.0569131764900355,39.87806291180732
0,10,1.052216205567539,40.68721372638996
0,11,1.0611586341601618,43.21113775733416
0,12,1.067294026948358,44.963093216263694
0,13,1.0604769704009607,44.51593896270254
0,14,1.0413760446994864,41.03908952564878
1,1,1.0469439155038462,2.435258480538479
1,2,1.0526258640338446,1.600912606334095
1,3,1.0161219135838027,-3.3613157501764173
1,4,1.034266165834819,-0.8543534835409058
1,5,1.037886753488885,0.5221627020144154
1,6,1.0812384734673193,4.847709781022957
1,7,1.069194254156816,-1.2829898996797755
1,8,1.1,-1.2550154995521083
1,9,1.0561156621030308,-1.5269767091876054
1,10,1.051477962666099,-0.6887852563286224
1,11,1.0606944599939436,1.9018024371813997
1,12,1.0385869651306026,2.3702897776312377
1,13,1.0539667493513003,2.759358943520675
1,14,1.038198677974345,-0.5313842931253533
--end of bus
--Delta
contingency id,Delta(MW)
1,0.0020542382015618457
--end of Delta
--line flow
contingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar)
0,'BL',1,2,i_1_2_'BL',21.545839117365425,-19.80714269739382,-21.413180604513258,14.393310993192662
0,'BL',1,5,i_1_5_'BL',16.419120561829605,-2.190055033142641,-16.286105533890176,-2.6076060956840683
0,'BL',2,3,i_2_3_'BL',49.05393129338556,7.3642081490013265,-47.99279717012153,-7.581374150203428
0,'BL',2,4,i_2_4_'BL',27.18659503538017,0.5652888137134382,-26.795824197354722,-3.0824116324611395
0,'BL',2,5,i_2_5_'BL',13.388193204437048,2.6812499703691417,-13.285198689814914,-6.147542843858086
0,'BL',3,4,i_3_4_'BL',-26.867305038186114,-0.4375334519528368,27.33582293540665,0.287663618061643
0,'BL',4,5,i_4_5_'BL',-57.45616395395038,10.30858433534757,57.88126012535054,-8.967700456871198
0,'BL',6,11,i_6_11_'BL',28.129509092733194,-1.7883985417564505,-27.48405497371591,3.140060199711412
0,'BL',6,12,i_6_12_'BL',9.550089891845252,1.3928047518688966,-9.452163519838523,-1.1889926517707108
0,'BL',6,13,i_6_13_'BL',26.32620763138908,4.196419953351376,-25.92408324159438,-3.4045114989513694
0,'BL',7,8,i_7_8_'BL',-0.3198276704775918,-18.449653249736453,0.31982767047764404,18.97389762020335
0,'BL',7,9,i_7_9_'BL',3.892560121739251,12.358193434017556,-3.892560121739188,-12.19676965615112
0,'BL',9,10,i_9_10_'BL',-14.298441790766283,11.388783479101477,14.393595623243638,-11.136017089779262
0,'BL',9,14,i_9_14_'BL',-4.38468415269781,8.218335791795743,4.483415068374046,-8.008321906337867
0,'BL',10,11,i_10_11_'BL',-23.21525217422293,5.582290576002094,23.637752644549128,-4.593263576177942
0,'BL',12,13,i_12_13_'BL',3.8087478890900766,-0.5523381958950645,-3.780022232456147,0.5783280757068677
0,'BL',13,14,i_13_14_'BL',18.03843219663982,-2.4551192688821937,-17.534717193684042,3.480702228741144
0,'BL',4,7,i_4_7_'BL',3.572732451261747,-6.000304547771857,-3.572732451261712,6.091459815718747
0,'BL',4,9,i_4_9_'BL',2.1081196789838095,2.0533561644495837,-2.1081196789837655,-2.011091658933061
0,'BL',5,6,i_5_6_'BL',-35.435084999643585,15.662833246087759,35.43508499964354,-12.613063969479278
1,'BL',1,2,i_1_2_'BL',21.577924962191368,-19.81718891638744,-21.444961695928477,14.404287672651012
1,'BL',1,5,i_1_5_'BL',16.39730590801141,-2.145569468117926,-16.264620575309685,-2.653014431664756
1,'BL',2,3,i_2_3_'BL',49.05764799626011,7.363667602737105,-47.99636374538895,-7.5802011792507535
1,'BL',2,4,i_2_4_'BL',27.286686097692183,0.6509893314504132,-26.892831837795523,-3.158075965823797
1,'BL',2,5,i_2_5_'BL',13.355197056495388,2.7430131128795963,-13.252362199105692,-6.209485369086613
1,'BL',3,4,i_3_4_'BL',-26.76256723149165,-0.366898034769238,27.227463327943457,0.2080404004248452
1,'BL',4,5,i_4_5_'BL',-57.99420032488385,10.228808253987713,58.42700324525414,-8.863614922392218
1,'BL',6,11,i_6_11_'BL',28.769849861607828,-1.808469125981933,-28.09473596321384,3.2222420717657663
1,'BL',6,13,i_6_13_'BL',34.71549391547663,5.5882603490767115,-34.01590355593625,-4.210548540500615
1,'BL',7,8,i_7_8_'BL',-0.32599038508253037,-18.69838093491482,0.32599038508247086,19.237282450762827
1,'BL',7,9,i_7_9_'BL',4.370975153368908,12.7204742353396,-4.370975153368978,-12.546375489800784
1,'BL',9,10,i_9_10_'BL',-14.880898195281675,11.538920379513089,14.98202458658226,-11.270288500704464
1,'BL',9,14,i_9_14_'BL',-3.052661414739417,8.494774950638472,3.1455168798543225,-8.297258935708397
1,'BL',10,11,i_10_11_'BL',-23.803681137561448,5.716561986927835,24.24843363404704,-4.67544544823202
1,'BL',12,13,i_12_13_'BL',-5.643415630748456,-1.7413308476654095,5.714853817220702,1.80596539733051
1,'BL',13,14,i_13_14_'BL',16.635376461304922,-2.876719548956108,-16.19681900516351,3.7696392581136005
1,'BL',4,7,i_4_7_'BL',4.044984768286278,-5.882605732825193,-4.0449847682863185,5.977906699575424
1,'BL',4,9,i_4_9_'BL',2.379270980796775,2.1707209818624227,-2.379270980796826,-2.1200801316168345
1,'BL',5,6,i_5_6_'BL',-36.035149568836914,15.666098572818136,36.03514956883696,-12.52846257913929
--end of line flow
")
    end

end
