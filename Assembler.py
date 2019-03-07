import sys
import getopt
import string
import os
import re
import glob

sequences = [];
outerSet = {};
allContigs = [];
degree = {};
usedEdges = [];
letterContigs = [];
	
###########################################################################
###########################################################################
### Usage Statement
###########################################################################
###########################################################################
def usage():
	print("This is Devan Bicher's Genome Assembler submission for CSE 308, Project1");
	print("");
	print("The formatting for this is as follows:");
	print("\tpython Assembler.py [inputReads.fasta] [outputfile]");
	print("where: inputReads.fasta is a series of reads(30bp each) in fasta format");
	print("       and outputfile is the file where the complete genome will be written to");
	print("NOTE: [outputfile] is optional, if not supplied the genome will simply print to the screen");
	print("");
	print("");
	print("");

###########################################################################
###########################################################################
### Parse Data
###########################################################################
###########################################################################

################################################
###gobble in a file
################################################
def gobble(fileName):
	#print("Gobbling: "+ fileName);
	
	inputFile = open(fileName, 'rt');
	inputLine = inputFile.readline();
	
	outputList = "".split();  #### list declaration
	
	while inputLine:
		outputList.append(inputLine);
		inputLine = inputFile.readline();

	inputFile.close();
	
	return outputList;

################################################
###parse the fasta file
################################################
def parseFasta(fileName):
	lineList = gobble(fileName);
	seqs = [];
	names = [];

	##first get the mutant sequences.
	tempLine = "";
	for line in lineList:
		if(line[0] == ">"):
			if( tempLine != "" ):
				seqs.append(tempLine);
			tempLine = "";
			names.append(line[1:-1]);
			continue;
		
		tempLine = tempLine + line[0:-1];
	seqs.append(tempLine);
	
	del names;
	
	return seqs;

###########################################################################
###########################################################################
### To and From Hash code
###########################################################################
###########################################################################

################################################
###Convert letters(text) to hashCode
################################################
def hashCode(sequence):
	code = 0;
	num = 0;

	for i in range(15):
		num = (pow(4,i))*letterToNumber(sequence[i]);
		code += num;
	
	return code;

################################################
###convert letters to numbers
################################################
def letterToNumber(letter):
	num = 0;
		
	if(letter == 'A'):
		num = 0;
	elif( letter == 'C'):
		num = 1;
	elif(letter == 'G'):
		num = 2;
	elif(letter == 'T'):
		num = 3;
		
	return num;
	
################################################
###convert the hash code back to text
################################################
def deCode(code):
	sequence = "";
	
	for i in range(15):
		j = 15 - (i+1);
		if(code == 0):
			while(j >= 0):
				sequence += 'A';
				j-=1;
			break;	
		elif(code%pow(4,j) > 0):
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
				print("Converting to code error");
		else:###when remainder = 0
			if(code/pow(4,j) == 1):
				sequence += 'C';
			elif(code/pow(4,j) == 2):
				sequence += 'G';
			elif(code/pow(4,j) == 3):
				sequence += 'T';
			else:
				print('Converting to code error');
			
		code = code%pow(4,j);
			
	return sequence[::-1];

###########################################################################
###########################################################################
### walking the path methods
###########################################################################
###########################################################################
	
################################################
### Lets Do some path walking
################################################
def walking():
	prefIT = iter(outerSet);  ##this is an iterator object, which allows you to run a (for) loop over the iterator
	for i in prefIT:	##iterate for all entries in degree
		while(degree[i]['out'] > degree[i]['in']): ##if outDegree > inDegree
			contig = []; ##reset current contig
			prefix = i;	 ##
			while (prefix != -1):
				contig.append(prefix); ##add this prefix to the current contig
				prefix = genContig(prefix); ##reassign the next input prefix as the suffex of the previous prefix
			allContigs.append(contig);	##add this finished contig to the list of all contigs	
			
################################################
### generate a node to add to the current contig
################################################
def genContig(prefix):
	if(prefix not in outerSet):
		return -1;
	innerIT = iter(outerSet[prefix]);
	for i in innerIT:
		for j in range(len(outerSet[prefix][i])):
			if(not outerSet[prefix][i][j] in usedEdges):
				degree[prefix]['out'] -= 1;  ##de-increment outdegree for this prefix
				degree[i]['in'] -= 1;	##de-increment indegree for this suffix	
				usedEdges.append(outerSet[prefix][i][j]); ##add this edge(read) to the used edges
				return i;
	return -1;
				
###########################################################################
###########################################################################
### Final Checking methods
###########################################################################
###########################################################################
	
################################################
### check to see if there are any loops
################################################
def checkLoops():
	pos = 0;
	j = 0;
	while(j < len(allContigs)-1):
		for i in range(pos,len(allContigs)):
			j = i; ## set placeholder j = i so we know when to end this loop, since we have to continually be breaking.
			if(allContigs[i][0] == allContigs[i][-1]): ## if the first read is the same as the last
				del allContigs[i];  ##delete this contig from the list
				pos = i;  ##since we already checked the previous entries restart the loop at the current position
				break;  ##and break the loop so we don't lose the index that gets shifted down from removing the current one

################################################
### check to see if any sequences are found within another sequence
################################################
def checkOverLaps():
	delete = [];  ##create a delete list with indeces of letterContigs that need to be deleted
	for i in range(len(letterContigs)):  
		if (i in delete):  ##if the index is already in the delete list don't repeat it
			continue;
		for j in range(len(letterContigs)):
			if(i == j):  ## if indeces are the same one will be within the other, don't add this index to delete list
				continue;
			elif(letterContigs[j] in letterContigs[i]): ## otherwise if contig j is within i, or is the same as i add it to the delete list
				if(not j in delete):
					delete.append(j);

	if(len(delete) > 0):
		delete.sort();  ##sorting the list makes more senes in my head
		for i in range(len(delete)): 
			del letterContigs[delete[i]-i]; ##we have to add the -i because the list indeces will change after every deletion

################################################
### check frame shifts
################################################
def checkFrameShifts(nodes):
	if nodes == True:
		contigs = allContigs;
	else:
		contigs = letterContigs
	
	done = False;
	overlap = 0;
	while(done == False):
		for i in range(len(contigs)): ##treat this as the one with suffix
			for j in range(len(contigs)): ##treat this as the one with prefix
				if(i ==j): ## then the two contigs are the same, don't analyze them
					continue;
				elif contigs[i] == contigs[j] and i!=j: ##if the two contigs are the same, but not at the same index
					del contigs[j];	##delete one of the them (we don't want redundancies)
					break;
				elif((contigs[i][-1] in contigs[j]) and (contigs[j][0] in contigs[i])):
				## if the suffix of i is in j, and the prefix of j is in i they will, probably, overlap	
					count = contigs[i].count(contigs[j][0]);
					pos=0;
					while(count > 0): ## find if prefix is only in i once, if it is its because of the overlap
						index = contigs[i][pos:].index(contigs[j][0]) + pos ##find the prefix of j in i
						overlap = len(contigs[i]) - index; ##find length of overlap
						for k in range(overlap):
							if(not contigs[i][index+k] == contigs[j][k]):  ## if the overlap is incorrect at any node: break
								count -= 1;
								pos = index + 1;
								break;
						else: ## if loop ended normally = no breaks = successful overlap
							if nodes == True:
								contigs[i].extend(contigs[j][overlap:]); ## then concatenate j to i (list method)
							else:
								contigs[i] += contigs[j][overlap:];
							del contigs[j];  #delete contig j, which was concatenated to i
							break;  ##break out of the while loop since we have an indexing error from deleting j
					else: ##if while ended normally we didn't have an overlap
						continue; ##so continue the j loop normally
					break;	##if while was broken break through the j loop	
			else:	##if j ended normally
				continue; ##move to the next iteration of i
			break;  ## if j did not end normally break out of i and move back to the while loop
			
		else:  ## if looping for i ends with no breaks we have successfully gone through all of i
			done = True; ## so we can exit the while loop
		
###########################################################################
###########################################################################
### main method
###########################################################################
###########################################################################
def main():
	if(len(sys.argv) == 1):
		usage();
		print("However, you did not enter anything valid so please re-run this script with input as above");
		raise SystemExit, 5
	
	usage();
	finalSequenceFile = "";
	sequences = parseFasta(sys.argv[1]);
	if(len(sys.argv) == 3):
		finalSequenceFile = sys.argv[2];
	
	##making the hash table only needs to happen once so I don't need to make a seperate method for it
	print("generating Hash Table...");
	for i in range(len(sequences)):
		##make suffix and prefix indexes
		Pindex = hashCode(sequences[i][0:15]);
		Sindex = hashCode(sequences[i][15:30]);
		
		##make the hash table from the reads
		if(Pindex in outerSet): ## if the prefix of our read is already in outerSet
			if(Sindex in outerSet[Pindex]):  ##if the suffix of our read is already in the specific prefix table:
				outerSet[Pindex][Sindex].append(i);	##add the seqeunce(read) index to the edge list
				###print(i);
			else: ##read suffix is not in read prefix
				outerSet[Pindex][Sindex] = [i]; ##create a list at suffix containing this read			
		else: ## prefix is not in outerset
			outerSet[Pindex] = {Sindex:[i]}; ##create a dictionary for prefix which is a list at key suffix containing the read
		
		##check if prefix is in degree
		if(Pindex not in degree):
			degree[Pindex] = {'in':0,'out':1}; ## if Pindex is not in degree add it with an outdegree = 1
		else: ##Pindex is in degree already, incremement its outDegree
			degree[Pindex]['out'] += 1;	
		
		##check if the suffix is in degree
		if(Sindex not in degree):
			degree[Sindex] = {'in':1,'out':0}; ## if Sindex is not in degree add it with an indegree = 1
		else: ##Sindex is in degree already, incremement its inDegree
			degree[Sindex]['in'] += 1;	
	
	##I have set up the degrees table and the hashtable, now do the walking
	print("Hash Table and Degree Table Done, walking paths & making contigs...");
	
	counter =0
	while len(sequences) > len(usedEdges) and counter <= 5:
		counter += 1;
		walking();
	
	print("Done walking, performing compression...");
	
	checkLoops();
	print("Loops: Removed!");
	nodes = True;
	checkFrameShifts(nodes);
	print("Node Frameshifts: Concatenated!");
	
	##convert to text to finish the final checks
	print("Converting to text to finish final checks");
	letters = "";
	
	for i in range(len(allContigs)):
		for j in range(len(allContigs[i])):
			letters += deCode(allContigs[i][j]);
		letterContigs.append(letters);
		letters = "";
	
	checkOverLaps();
	print("Overlaps: DONE!");
	nodes = False;
	checkFrameShifts(nodes);
	print("Sequence Frameshifts: Concatenated!");
	
	##print to file if the user specifies that
	if (len(finalSequenceFile) > 0):
		file = open(finalSequenceFile, 'w');
	else:
		print("");
		print("HERE IS THE ASSEMBLED GENOME:");
	
	##final printing statements
	for i in range(len(letterContigs)):
		if (len(finalSequenceFile) > 0):
			file.write(">Contig " + str(i + 1) + "\n" + letterContigs[i] + "\n");
		else:
			print("Contig " + str(i + 1) + ":\n" + letterContigs[i]);
			
	if (len(finalSequenceFile) > 0):	
		file.close();
	
if __name__ == "__main__":
	main()
	