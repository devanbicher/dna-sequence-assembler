import sys
import getopt
import string
import os
import re
import glob

outerSet ={1:{11:[111,112,113,114],12:[121,122,123],13:[131]},2:{21:[211,212],22:[221,222,223,224,225]},3:{31:[311,312]},4:{41:[411,412,413],42:[421,422]},5:{51:[511,512,513],52:[521,522],53:[531]}};
print(outerSet);
it = iter(outerSet);
for i in it:
	print(i);
	it1 = iter(outerSet[i]);
	for j in it1:
		print(j);
		print(outerSet[i][j]);

