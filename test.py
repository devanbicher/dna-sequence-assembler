import sys
import getopt
import string
import os
import re
import glob

univar = "this is a universal variable";


def main():
###test function for running individual scripts
	print(univar);
	seq = sys.argv[1];
	hashCode(seq);

###########################################################################
###########################################################################
### Hash Code
###########################################################################
###########################################################################
def hashCode(sequence):
	code = 0;
	num = 0;

	for i in range(15):
		num = pow(4,i)*letterToNumber(sequence[i]);
		code += num;
	
	print code;
	deCode(code);

###########################################################################
###########################################################################
### deCode
###########################################################################
###########################################################################
def letterToNumber(letter):
	num = 0;
		
	if (letter == 'A'):
		num = 0;
	elif (letter == 'C'):
		num = 1;
	elif (letter =='G'):
		num = 2;
	elif (letter == 'T'):
		num = 3;
	##else:
		##Maybe throw some error, could there even be an error?
		
	return num;
	
	
###########################################################################
###########################################################################
### deCode
###########################################################################
###########################################################################
def deCode(code):
	sequence = "";
	
	for i in range(15):
		j = 15 - (i+1);
		if (code == 0):
			while(j >= 0):
				sequence += 'A';
				j-= 1;
			break;	
		elif (code%pow(4,j) > 0):
			current = code - code%pow(4,j);
			if(current == pow(4,j)):
				sequence += 'C';
			elif(current == 2*pow(4,j)):
				sequence += 'G';
			elif(current == 3*pow(4,j)):
				sequence += 'T';
			elif(current == 0):
				sequence += 'A';
			else:
				###if this is ever the case I made a mistake
				print("error");
		else:###when remainder = 0
			if(code/pow(4,j) == 1):
				sequence += 'C';
			elif(code/pow(4,j) == 2):
				sequence += 'G';
			elif(code/pow(4,j) == 3):
				sequence += 'T';
			else:
				print('error');
			
		code = code%pow(4,j);
			
	print sequence[::-1];
	
if __name__ == "__main__":
	main()	