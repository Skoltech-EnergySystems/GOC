QUESTIONS

since I have complex indeces like g = (i, id) and f = (i, i`, id) a cannot just put my parameters in a list and address them by their relative index in the list
so mb it is better to declare all such parameters as dictionaries to with keys equal to complex indeces (i, id) => "value"


branch zorin_tests

before I had a structure of network where i stored all data in lists, now I create structures for all sources of data (bus, line) and store LISTS with objects of these structures inside network structure

i need to rewrite parsers to call a function for each line of the data to calculate all the parameters and initialize an object of corresponding structure
 but in this case i have a problem storing data and costs inside one type
 since i read all data for generator and then i read costs data for each generator but if the order is different i need to find an appropriate bus to update its costs
 i think to put index of a bus and its index in ListBus to a dictionary so whenever i read costs data for bus i a go to this dictionary and find a position of it in the ListBus which is j for example so it should be faster then just go through the whole list to find it
