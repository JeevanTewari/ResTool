## compared to V1, this tool takes in input from the console (run with python ResToolV2.py)
## refined to take in motifs like so ('118-120') for MAY motif and returns top 10 frequencies
## compared to V2 now outputs associated residues with each motif -- stored in identification dict
## compared to v3 subsets the sequences by species: use GREP?

import re, sys, csv
from difflib import get_close_matches

def storage_setup(arguments = None, command_line_input = False):
	# taking in arguments
	# 1) species 2) reference sequence 3) indices

	species = arguments[1]
	ref_seq = arguments[2]
	indices = arguments[3]

	print("Species: %s\tReference Sequence: %s\tIndices: %s"
	 % (species, ref_seq, indices))

	with open ('LargeAli.fa', 'r') as f:
		file_contents = f.read()

	fasta_sequences = file_contents.split(">")
	fasta_sequences.pop(0)
	#print(fasta_sequences)

	# split all in to tuples and then assign them as values to a dictionary
	# key of dicitonary should be identifier within | | characters
	seq_storage = {}

	for i in fasta_sequences:
		hold = i.split('\n',1)
		seqname = hold[0]
		sequence = hold[1].replace('\n', '').replace('\r', '') # mysterious carriage returns
		
		code_valid = False
		crds = []
		for i in range(len(seqname)):
			if seqname[i] == "/":
				code_valid = True
				end = i

		if code_valid:
			code = seqname[0:end] # OR identifier code
		else:
			code = seqname

		seq_storage[code] = (seqname, sequence) # adding values to the dictionary

	seq_storage = subset(seq_storage, species)

	# determining length of the sequences
	first_key = list(seq_storage.keys())[0]
	length = len(seq_storage[first_key][1])
	size = len(seq_storage)

	# calculating frequency of residue at particular position
	## INPUTS HERE
	number_of_seqs = len(seq_storage)

	reference_sequence = ref_seq
	index = indices
	reference_sequence = reference_sequence

	if '-' not in index:
		index = int(index)
		motif = False
	else:
		indices = index.split("-")
		index = (int(indices[0]), int(indices[1]))
		motif = True
	#reference_sequence = 'hOR13G1' # use this as a reference sequence


	# identifications are the actual motifs
	# keys of the frequency dict are the sequencce names
	(frequencies, identifications, positions) = calculate_residue_frequency(reference_sequence, index, seq_storage, motif)
	frequencies = {k: v for k, v in sorted(frequencies.items(), key=lambda item: item[1],
		reverse=True)} # sorting the dictionaries

	outputData = []
	print('On the alignment file, the residues aligned to this index are:')
	limit = 10
	for i in frequencies:
		print(' ++ %.2f%% %s --- %s' % (float(frequencies[i]/size * 100), i, identifications[i][0:10]))
		outputData.append([i, round(float(frequencies[i]/size * 100),2), identifications[i]]) ## EDITING CSV
		if limit == 1:
			break
		else:
			limit -= 1


	print('\n')
	for x in positions:
		temp_dict = positions[x]
		temp_dict = {k: v for k, v in sorted(temp_dict.items(), key=lambda item: item[1],
		reverse=True)}

		total = 0
		for i in temp_dict:
			total = total + temp_dict[i]
		tup_list  = []
		for i in temp_dict:
			tup_list.append(i + "=" + str(100*round(float(temp_dict[i]/total),3)))
		print('position' ,index[0]+x, tup_list)
		outputData.append([index[0]+x, tup_list])
 
	print('\n')
	writeCSV(outputData)

def calculate_residue_frequency(reference_sequence, index, storage, motif):
	## return a dictionary with key = residue and value = frequency (ordered)
	## first step  find the character at the index on the reference sequence
	set_of_keys = list(storage.keys())

	if motif:
		length = index[1] - index[0]
		index = index[0]
	
	# imported function that gets the closest string matches and returns as list length n
	match = get_close_matches(reference_sequence, set_of_keys, 1)[0]

	original_index = index
	i = 0
	gaps = 0
	while True:
		if storage[match][1][i] != '-':
			i += 1
		else:
			gaps += 1
			i += 1
			index += 1
		if i == index: # end of loop, now determining if looking at one residue or series
			if motif:
				save_char = storage[match][1][i-1:i+length]
			else:
				save_char = storage[match][1][i-1]
			break
		
	print('\nOn the reference sequence %s, index %d-%d is %s \n' 
		% (match, original_index, original_index+length, save_char))

	## Now calculating the frequency of each residue at that position
	## 

	frequencies = {}
	identifications = {}
	frequenciesIndividuals = {}
	positions = {}
	for x in range(length+1):
			positions[x] = {}

	for i in storage:

		residue = storage[i][1][index-1:index+length] if motif else storage[i][1][index-1]
		frequencies[residue] = frequencies[residue] + 1 if residue in list(frequencies.keys()) else 1 
		# for every unique motif/residue at index i in the alignment create a key for it

		for x in range(len(residue)):
			if residue[x] in list(positions[x].keys()):
				positions[x][residue[x]] = positions[x][residue[x]] + 1
			else:
				positions[x][residue[x]] = 1
			#frequenciesIndividuals[x] = frequenciesIndividuals[x] + 1 if x in list(frequenciesIndividuals[x].keys()) else 1
		
		# this conditional helps secure the sequence codes (in a list) corresponding to every motif
		if residue in list(identifications.keys()): # for every unique motif/residue at index i in the alignment create a key for it
			identifications[residue].append(i) # in the identifications dictionary key is the motif and value is the code for the organism
		else:
			identifications[residue] = [i]

	return (frequencies, identifications, positions)

def subset(seq_storage, species):

	if species.lower() == "no":
		return seq_storage
	species_storage ={}

	## take in species of interest:

	## accout for subsetting for Human ORs
	if species.lower() == 'human':
		for x in seq_storage:
			if x[0:3] == 'hOR':
				species_storage[x] = seq_storage[x]

	elif species.lower() == 'mouse' or species.lower() == 'mice':
		for x in seq_storage:
			if x[0:3] == 'mOR' or x[0:3].lower() == 'ror' or x[0:3].lower() == 'olf':
				species_storage[x] = seq_storage[x]			

	## account for subsetting other ORs
	elif species.upper() != 'NO':
		species = get_close_matches(species, list(seq_storage.keys()), 1)[0]
		for x in seq_storage:
			if species[0:3].lower() in x.lower(): # comparing species to keys in storage dictionary
				species_storage[x] = seq_storage[x]

	seq_storage = species_storage

	return seq_storage # returned modified dictinoary only containing desired species

def writeCSV(data):
	with open('output_data.csv', 'w') as new_file:
		csv_writer = csv.writer(new_file, delimiter = ',')

		csv_writer.writerow(['residues','% conservation', 'associated sequences'])
		for line in data:
			csv_writer.writerow(line)
	return

def getinput():
	species = input('Subset the alignment for a particular species (enter \'no\' or species name)?: ')
	reference_sequence = input('Enter the reference sequence (ie hOR13G1): ')
	index = input('Enter the index of the residue (or indices separated by -): ')

	arguments = ['null', species, reference_sequence, index]
	return arguments

def main():
	print('\n')
	command_line_input = True
	arguments = (sys.argv)
	if len(arguments) != 4: # ensures there is an arg passed from command line length 4
		#arguments = ['null', 'anem_123124', 'anem_118391_153', '215-218'] # for testing
		arguments = getinput()

		command_line_input = False
	storage_setup(arguments, command_line_input)

if __name__ == "__main__":
    main()


## maybe create dictionary with residue type (classifications to provide info on charge
## of motif and what interactions could 
## add function to show what receptors are associated with the output motifs
## return name of OR that corresponds to the position and particular residues

## make folder which is named something unique
## folder will contain all output CSVs
## tool must be contained in while loop

## make it have user promptable OR takes in requests using argument line
## If user wants to use command line, an example input:
## python motiffinderv3.py anem_123124 anem_118391_153 215-218
## python motiffinderv3.py human hOR1A1 215-218
## python python motiffinderv3.py mouse Olfr1500 210-215

## Otherwise user can use prompted interface

## add percentages to single residues instead of numbers - sort by highest percent
## list avaliable species and how to call on them


#


