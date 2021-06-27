import sys
import json

# An Annot object represents an annotation on a genome. There are three
# attributes: the start and end of the annotation, and its name
class Annot:
    
    def __init__(self, start, end, name=''):        
        self.start = start
        self.end = end
        self.name = name

# A Seq object represents a DNA sequence. It has three attributes:
# the sequence (as a String), a name, and an array of annotations.
class Seq:
    
    def __init__(self, sequence, name='', annotations=[]):        
        self.sequence = sequence
        self.name = name
        self.annotations = annotations

def read_fasta(file_path):
    sequences = []
    names = []
    
    with open(file_path, 'r') as file_in:        
        line = file_in.readline()
        while line != '':
            # The name of every sequence starts with '>'
            if line.startswith('>'):
                # Split the name to extract the first part only
                entry = line.split()
                name = entry[0][1:] # Remove ">"
                names.append(name)
                
                # Read the next lines until the next '>' is reached
                seq = ""
                line = file_in.readline()
                while (not line.startswith('>')) and (line != ''):
                    # Append the line to sequence
                    seq = seq + line.strip()
                    line = file_in.readline()
                
                # Append the sequence to the sequences array
                sequences.append(seq)     
                
    return sequences, names

def read_gff3(file_path):
    annotations = []
    
    with open(file_path, 'r') as file_in:
        for line in file_in:
            if ('CDS' in line) and ('+' in line):
                entry = line.split()
                
                # The name is the first string in the entry
                # The start is the 4th element and the end is the 5th
                annot = Annot(int(entry[3]),int(entry[4]), name=entry[0])
                
                # Append the annotation to the list
                annotations.append(annot)
    
    return annotations

# Read both the fasta and gff3 and combine it into a single array of Seqs. This
# array is the genome.
def parse_inputs(path_fasta, path_gff3):
    sequences, names = read_fasta(path_fasta)
    annotations = read_gff3(path_gff3)

    genome = []

    # Go through the list of names and put the appropriate
    # sequences, names and annatations into a Seq object
    for i in range(0, len(names)):
        annotation_seq = []
        
        # Create an annotation array for sequence[i]
        for j in range(0,len(annotations)):
            if (annotations[j].name == names[i]):
                annotation_seq.append(annotations[j])
        
        seq = Seq(sequences[i], names[i], annotation_seq)        
        genome.append(seq)
    
    return genome

# Get the intergenic sequences. The output is stored in an array.
def get_intergenic(genome):
    intergenic_seqs = []
    
    for seq in genome:
        # If there are no annotations, add the entire sequence in
        if seq.annotations == []:
            intergenic_seqs.append(seq)
        else:
            # Extract the section of the sequence before the first annotation
            first_intergenic = seq.sequence[0:seq.annotations[0].start-1]
            # Only add it if the sequence is not an empty string
            if first_intergenic != '': 
                intergenic_seqs.append(Seq(first_intergenic))

            # Loop through the annotations of the sequence
            for i in range(0, len(seq.annotations)-1):
                curr_annot = seq.annotations[i]
                next_annot = seq.annotations[i+1]
                if curr_annot.end < next_annot.start:
                    intergenic = seq.sequence[curr_annot.end: next_annot.start-1]
                    intergenic_seqs.append(Seq(intergenic))
            
            # Extract the section of the sequence after the last annotation
            last_annot = seq.annotations[-1]
            # Only add it if the sequence is not an empty string
            if last_annot.end < len(seq.sequence):
                last_intergenic = seq.sequence[last_annot.end:]    
                intergenic_seqs.append(Seq(last_intergenic))

    return intergenic_seqs

# Get the genic sequences. The output is stored in an array.
def get_genic(genome):
    genic_seqs = []

    # Loop through every sequence in the genome    
    for seq in genome:               
        # Go through all the annotations in a sequence
        for annot in seq.annotations:            
            # Add the sequence that falls between the start and end
            # of an annotation
            genic = seq.sequence[annot.start-1: annot.end]
            genic_seqs.append(Seq(genic))

    return genic_seqs

def avg_seq_len(seqs):
    lengths = []
    
    # Store the lengths of all the intergenic sequences
    for seq in seqs:
        lengths.append(len(seq.sequence))
    
    # Find and return the average
    avg = sum(lengths)/len(lengths)
    return avg

# Get the intergenic frequencies.
def intergenic_freqs(intergenic_seqs):
    intergenic_freq = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }

    for seq in intergenic_seqs:
        for s in seq.sequence:
            if s == 'A':
                intergenic_freq['A'] +=1
            elif s == 'C':
                intergenic_freq['C'] +=1
            elif s == 'G':
                intergenic_freq['G'] +=1
            elif s == 'T':
                intergenic_freq['T'] +=1
            else:
                sys.exit("Error: invalid nucleotide")
    
    total = sum(intergenic_freq.values())
    for nucleotide in intergenic_freq:
        intergenic_freq[nucleotide] /= total

    return intergenic_freq

# Helper function that converts integers to codons.
def ints2codon(ints):

    def int2base(i):
        if i == 0:
            return 'A'
        elif i == 1:
            return 'C'
        elif i == 2:
            return 'G'
        else:
            return 'T'

    codon = int2base(ints[0])+int2base(ints[1])+int2base(ints[2])
    return codon

# Get the frequencies of each codon.
def codon_freqs(genic_seqs):
    
    # Intialize start, middle, stop frequencies
    start_freq = {}
    for i in range(0,4):
        for j in range(0,4):
            for k in range(0,4):
                codon = ints2codon([i,j,k])
                start_freq[codon] = 0
    
    mid_freq = {}
    for i in range(0,4):
        for j in range(0,4):
            for k in range(0,4):
                codon = ints2codon([i,j,k])
                mid_freq[codon] = 0

    stop_freq = {}
    for i in range(0,4):
        for j in range(0,4):
            for k in range(0,4):
                codon = ints2codon([i,j,k])
                stop_freq[codon] = 0

    # Go through every genic sequence
    for seq in genic_seqs:
        # Extract the first codon (start codon)
        start_codon =seq.sequence[0] + seq.sequence[1] + seq.sequence[2]
        start_freq[start_codon] += 1
        
        # Loop through the sequence in chunks of 3 and extract the codon
        for s in range(3, len(seq.sequence)-3, 3):
            codon =seq.sequence[s]+ seq.sequence[s+1] + seq.sequence[s+2]            
            mid_freq[codon] += 1

        # Extract the last codon (stop codon)
        stop_codon = seq.sequence[-3] + seq.sequence[-2] + seq.sequence[-1]
        stop_freq[stop_codon] += 1
        
    # Divide every occurences by the total number of codons
    total = sum(start_freq.values())
    for codon in start_freq:
        start_freq[codon] /= total

    total = sum(mid_freq.values())
    for codon in mid_freq:
        mid_freq[codon] /= total
    
    total = sum(stop_freq.values())
    for codon in stop_freq:
        stop_freq[codon] /= total
    
    return start_freq, mid_freq, stop_freq

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    gff3_path = sys.argv[2]

    genome = parse_inputs(fasta_path, gff3_path)
    intergenic_seqs = get_intergenic(genome)
    genic_seqs = get_genic(genome)

    avg_intergenic_len = avg_seq_len(intergenic_seqs)
    avg_genic_len = avg_seq_len(genic_seqs)
    intergenic_freqs = intergenic_freqs(intergenic_seqs)    
    start_freq, middle_freq, stop_freq = codon_freqs(genic_seqs)

    data = {
        "avg_intergenic_length": avg_intergenic_len,
        "avg_genic_length": avg_genic_len,
        "intergenic_frequencies": intergenic_freqs,
        "start_frequencies": start_freq,
        "middle_frequencies": middle_freq,
        "stop_frequencies": stop_freq 
    }

    with open('config.json', 'a') as file_out:
        json.dump(data, file_out, indent=4)
        
