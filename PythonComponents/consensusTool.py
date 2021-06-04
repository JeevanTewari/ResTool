import re, sys, csv
from difflib import get_close_matches

with open ('anemone.fa', 'r') as f:
	file_contents = f.read()
	
	fasta_sequences = file_contents.split(">")
	fasta_sequences.pop(0)


	sequence_storage = []

	for i in fasta_sequences:
		hold = i.split('\n', 1)[1]
		sequence_storage.append(hold.replace('\n', '').replace('\r', ''))


	length = len(sequence_storage)


	# at this point we have a list of all of the sequences as a string - they are aligned so just need to find
	# what character is most common at each point

	residue_dictionary = {} # key = index on the sequence, value = dictionary containing frequencies of residues

	for x in sequence_storage: # keeps track of the frequence of amino acids at a particular index
		for index, amino_acid in enumerate(x):

			if index in residue_dictionary:
				if amino_acid in residue_dictionary[index]:
					residue_dictionary[index][amino_acid] = residue_dictionary[index][amino_acid] + 1
				else:
					residue_dictionary[index][amino_acid] = 1
			else:
				residue_dictionary[index] = dict()
				residue_dictionary[index][amino_acid] = 1

	#for x in residue_dictionary:
		#print(x, residue_dictionary[x])

	consensus = ""

	for index in range(len(sequence_storage[0])):
		maximum = 0
		character = ''
		for internal_dict in residue_dictionary[index]:
			if residue_dictionary[index][internal_dict] > maximum:
				if internal_dict == "-" and residue_dictionary[index][internal_dict] < length/2:
					continue
				else:
					maximum = residue_dictionary[index][internal_dict] 
					character = internal_dict
		consensus += character

	print(consensus, )


	# find number of sequences: if frequency of "-" is > 50%, use that character
	# adjust if in determining max to exclude "-"s