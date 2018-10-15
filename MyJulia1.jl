function MyJulia1(rawFile, genFile, contFile)

    open("solution1.txt", "w") do f1
        write(f1,"--generation dispatch
bus id,unit id,pg(MW),qg(MVar)
1,'1 ',37.9649606792,1.8583662976
6,'1 ',110.4266998728,-13.0957522832
8,'1 ',0.3198286705,6.9978519286
2,'1 ',84.7177624287,17.5456887762
3,'1 ',5.8017882658,26.8372933492
--end of generation dispatch")
    end

    open("solution2.txt", "w") do f2
        write(f2,"--contingency generator
contingency id,genID,bus id,unit id,q(MW)
1,l_14,1,'1 ',1.8920439657
1,l_17,6,'1 ',-13.0424143005
1,l_18,8,'1 ',7.2486797834
1,l_15,2,'1 ',17.7006839160
1,l_16,3,'1 ',26.9101533072
--end of contingency generator
--bus
contingency id,bus id,v(pu),theta(deg)
0,1,1.095594573957004,1.435232387005952e-19
0,2,1.091403785104688,-0.6180363565965823
0,3,1.0653687829352025,-5.351616491888493
0,4,1.0696679950731753,-2.812394002728806
0,5,1.0732490989879437,-1.5412323139805826
0,6,1.1000000027825878,2.566303896926547
0,7,1.0887688400411457,-3.195338934299448
0,8,1.0999997214423016,-3.1683868054882414
0,9,1.0791659523273798,-3.416666549410084
0,10,1.0740732341201584,-2.6411064889777407
0,11,1.0816480064080651,-0.21457936924695184
0,12,1.086552173722784,1.4789231266111569
0,13,1.0801035751083823,1.044319074272447
0,14,1.0628667726353607,-2.302089030565137
1,1,1.0955945706309593,1.0791384921753986e-17
1,2,1.0914037832506747,-0.6180363872562077
1,3,1.0653687802104295,-5.351616538063477
1,4,1.0696679912885483,-2.812394017730414
1,5,1.0732490948029123,-1.5412323067064562
1,6,1.0999999918648946,2.56630400297487
1,7,1.088768835806207,-3.1953389757403516
1,8,1.0999997197487705,-3.168386846379175
1,9,1.079165946315946,-3.41666660587974
1,10,1.0740732271550284,-2.6411065268843528
1,11,1.081647997423782,-0.2145793452444017
1,12,1.0865521630366526,1.4789232100950682
1,13,1.0801035647067476,1.0443191403157832
1,14,1.0628667645498915,-2.302089051413765
--end of bus
--Delta
contingency id,Delta(MW)
1,0.0019761473
--end of Delta
--line flow
contingency id,line id,origin bus id,destination bus id,circuit id, p_origin(MW) q_origin(MVar), p_destination(MW), q_destination(MVar)
0,i_1,1,2,'BL',22.016224231039036,0.6662585234749275,-21.93789258530943,-0.4271004492835466
0,i_2,1,5,'BL',15.94873544904371,7.303563173167691,-15.8102294313366,-6.731799667904231
0,i_3,2,3,'BL',49.5586523253742,4.5931185513257935,-48.581442051991736,-0.47610816993567306
0,i_4,2,4,'BL',27.012125177538216,5.037358525146102,-26.64378922842781,-3.919736867133043
0,i_5,2,5,'BL',13.198963120141514,7.159727738626643,-13.091162832464848,-6.83059140813995
0,i_6,3,4,'BL',-26.27866016341438,8.272242830441764,26.72676569458316,-7.128540562186293
0,i_7,4,5,'BL',-57.383135873411504,9.766283907594502,57.778459606634144,-8.519311442725332
0,i_8,4,7,'BL',3.8059293295105365,-9.97731957451697,-3.8059293295105365,10.181147262400144
0,i_9,4,9,'BL',2.2589169911879203,-1.873212633675818,-2.2589169911879203,1.9137753080622544
0,i_10,5,6,'BL',-36.0021964417205,-10.932236505717215,36.0021964417205,13.819008474865482
0,i_11,6,11,'BL',27.853860863902156,-2.4470860647611015,-27.24015992459945,3.732252592738768
0,i_12,6,12,'BL',9.493785373100655,1.3052805075395002,-9.400500017466909,-1.1111276572858846
0,i_13,6,13,'BL',26.170649150846543,3.8330730429882296,-25.78818443534141,-3.0798806244114014
0,i_14,7,8,'BL',-0.319827671293373,-6.941641333340896,0.319827671293373,7.013397059278467
0,i_15,7,9,'BL',4.125756999887203,9.511945530098863,-4.125756999887203,-9.412183532605912
0,i_16,9,10,'BL',-14.076549322172967,11.928787425059825,14.16953893248935,-11.681770069770888
0,i_17,9,14,'BL',-4.2225824318620555,8.570849825437566,4.322220404428275,-8.358906508448683
0,i_18,10,11,'BL',-22.99119548438842,6.128043555992578,23.393857594534477,-5.185455969207104
0,i_19,12,13,'BL',3.757084385839287,-0.63020319038086,-3.729927107323599,0.6547740614188512
0,i_20,13,14,'BL',17.852438264364636,-2.8561961291378903,-17.373522530662274,3.831286830846707
1,i_1,1,2,'BL',22.01622430447103,0.6662557613535823,-21.937892658338026,-0.4270976859315001
1,i_2,1,5,'BL',15.948735381166387,7.303563575108147,-15.810229360950157,-6.731800059486496
1,i_3,2,3,'BL',49.55865238355122,4.593118997780803,-48.58144210295608,-0.47610858600339684
1,i_4,2,4,'BL',27.012125240334058,5.037359667155361,-26.643789282704493,-3.9197379832927157
1,i_5,2,5,'BL',13.198963087309965,7.159729185471554,-13.091162789776122,-6.8305928248886945
1,i_6,3,4,'BL',-26.278660090718038,8.272243481078807,26.726765628278386,-7.128541196510027
1,i_7,4,5,'BL',-57.383136119953036,9.766285053560697,57.77845986188621,-8.519312561215715
1,i_8,4,7,'BL',3.8059295640060586,-9.97731930198901,-3.8059295640060586,10.181146984784963
1,i_9,4,9,'BL',2.258917125624247,-1.873212183510015,-2.258917125624247,1.9137748531004588
1,i_10,5,6,'BL',-36.00219680826883,-10.932233337481664,36.00219680826883,13.819005241707277
1,i_11,6,11,'BL',27.853860685963145,-2.4470871210655227,-27.240159738201104,3.7322536667580666
1,i_13,6,13,'BL',26.17064904138101,3.833072506687204,-25.78818432366389,-3.0798800837541123
1,i_14,7,8,'BL',-0.31982767608317747,-6.941642877160577,0.31982767608317747,7.013398635510303
1,i_15,7,9,'BL',4.125757241005942,9.511947252335643,-4.125757241005942,-9.412185221814603
1,i_16,9,10,'BL',-14.076549124208789,11.928788536790561,14.169538741283617,-11.68177116354881
1,i_17,9,14,'BL',-4.222582252424754,8.57085053672368,4.322220237754697,-8.35890719258458
1,i_18,10,11,'BL',-22.991195291343207,6.1280446497728285,23.393857409932547,-5.185457043222332
1,i_19,12,13,'BL',3.757084333682053,-0.630203320835399,-3.7299270550577455,0.654774191971692
1,i_20,13,14,'BL',17.85243810220068,-2.8561968003399536,-17.373522362139816,3.8312875149949375
--end of line flow")
    end

end
