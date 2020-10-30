import re
from difflib import get_close_matches

def storage_setup():
	with open ('ORS-ALI-COMP-13-091818.fa', 'r') as f:
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
		#print(sequence)
		
		code_valid = False
		crds = []
		for i in range(len(seqname)):
			if seqname[i] == "|":
				crds.append(i)
				code_valid = True

		if code_valid:
			code = seqname[crds[0]+1:crds[1]] # OR identifier code
		else:
			code = seqname[0:8]

		seq_storage[code] = (seqname, sequence) # adding values to the dictionary

	# determining length of the sequences
	first_key = list(seq_storage.keys())[0]
	length = len(seq_storage[first_key][1])
	size = len(seq_storage)

	# calculating frequency of residue at particular position
	## INPUTS HERE
	number_of_seqs = len(seq_storage)
	index = 118
	reference_sequence = 'hOR1G1' # use this as a reference sequence
	reference_sequence = reference_sequence.upper()

	frequencies = calculate_residue_frequency(reference_sequence, index, seq_storage)
	frequencies = (sorted(frequencies.items(), key = 
             lambda kv:(kv[1], kv[0]), reverse = True)) # sorting the dictionaries

	print('On the alignment file, the residues aligned to this index are:')
	for i in frequencies:
		print(' ++ %.4f%% %s' % ((float(i[1])/size) * 100, i[0]))

def calculate_residue_frequency(reference_sequence, index, storage):
	## return a dictionary with key = residue and value = frequency (ordered)
	## first step find the character at the index on the reference sequence
	set_of_keys = list(storage.keys())

	# imported function that gets the closest string matches and returns as list length n
	match = get_close_matches(reference_sequence, set_of_keys, 1)[0]
	#print(match)

	original_index = index
	i = 0
	while True:
		if storage[match][1][i] != '-':
			i += 1
		elif storage[match][1][i] == '-':
			i += 1
			index += 1
		if i == index:
			save_char = storage[match][1][index-1]
			break

		
	print('On the reference sequence %s, index %d is %s \n' 
		% (match, original_index, save_char))

	## Now calculating the frequency of each residue at that position
	## 

	frequencies= {}
	for i in storage:

		residue = storage[i][1][index-1]
		if residue in list(frequencies.keys()):
			frequencies[residue] = frequencies[residue] + 1
		else:
			frequencies[residue] = 1

	return frequencies

def main():
    storage_setup()

if __name__ == "__main__":
    main()

## NOTES - since jalview can do exactly what this code does with
## annotations --> autocalculated annotation --> show consensus logo,
## maybe try to show how well a particular motif is conserved in the grouping
## EG frequency of MAYDRY vs MAYERY vs MAYDRA etc.
## not sure if this will help us realize anything substantial but potentially
## could help group OR's based on particular mechanisms

## table with all possible amino acids all positions
## Rows: amino acid types
## Cols: positions