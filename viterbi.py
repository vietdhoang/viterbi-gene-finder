import sys
import json
import numpy as np
import pandas as pd

# Usage:
# python viterbi.py PATH/TO/config.json PATH/TO/fasta_file
# config.json is produced by config.py

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

#  Read the fasta file. The function returns an array of sequences, called 
# the genome
def read_fasta(file_path):
    genome = []
    
    with open(file_path, 'r') as file_in:
        line = file_in.readline()
        while line != '':
            # The name of every sequence starts with '>'
            if line.startswith('>'):
                # Split the name to extract the first part only
                entry = line.split()
                name = entry[0][1:] # Remove ">"
                
                # Read the next lines until the next '>' is reached
                seq = ""
                line = file_in.readline()
                while (not line.startswith('>')) and (line != ''):
                    # Append the line to sequence
                    seq = seq + line.strip()
                    line = file_in.readline()
                
                # Append the sequence to the sequences array
                sequence = Seq(seq, name=name)
                genome.append(sequence)     
                
    return genome

# Get the log transition probabilities of HMM. A log transformation had to be
# applied because the probabilities can get very small.
def get_log_transition_probs(avg_intergenic_len, avg_genic_len):
    trans_prob = [[-float('inf') for j in range(0, 4)] for i in range(0,4)]
    
    # trans_prob[i][j]= probability of going from state i to state j
    # 0:= Intergenic, 1:=Start, 2:=Middle, 3:=Stop
    trans_prob[0][0] = np.log((avg_intergenic_len-1)/avg_intergenic_len)
    trans_prob[0][1] = np.log(1/avg_intergenic_len)
    trans_prob[1][2] = np.log(1)
    trans_prob[2][2] = np.log((avg_genic_len-1)/avg_genic_len)
    trans_prob[2][3] = np.log(1/avg_genic_len)
    trans_prob[3][0] = np.log(1)

    return trans_prob

# Get the log emmision probabilities of HMM. A log transformation had to be
# applied because the probabilities can get very small.
def get_log_emmission_probs(freq):
    
    for key in freq:
        if freq[key] == 0:
            freq[key] = -float('inf')
        else:
            freq[key] = np.log(freq[key])   
    
    return freq

# Implementation of the viterbi algorithm
def viterbi(sequence, log_trans_probs, log_intergenic_probs,
            log_start_probs, log_mid_probs, log_stop_probs):    
    
    V = np.zeros((4,len(sequence)))
    V.fill(-float('inf'))    
    back_point = [[[-1,-1] for j in range(0, len(sequence))] for i in range(0,4)]    

    # Base cases

    # Beta states: 0:=intergenic, 1:= start, 2:=middle, 3:=stop

    # Always start with an intergenic nucleotide
    # Beta = 0, i=0
    V[0][0] = 0

    # Beta = 0, i=1 
    V[0][1] = V[0][0] + log_intergenic_probs[sequence[1]] + log_trans_probs[0][0]
    back_point[0][1] = [0,0]

    # Beta = 0, i=2 
    V[0][2] = V[0][1] + log_intergenic_probs[sequence[2]] + log_trans_probs[0][0]
    back_point[0][2] = [0,1]

    # Beta = 0, i=3 
    V[0][3] = V[0][2] + log_intergenic_probs[sequence[2]] + log_trans_probs[0][0]
    back_point[0][3] = [0,2]

    # Beta = 1, i=3
    V[1][3] = V[0][0] + log_start_probs[sequence[1:4]] + log_trans_probs[0][1]
    back_point[1][3] = [0,0]

    # Inductive step        
    for i in range(4, len(sequence)):
        
        for beta in range(0, 4):
            
            if beta == 0:
                prev_col = np.array([
                    # Gamma = 0
                    V[0][i-1] + log_intergenic_probs[sequence[i]] + log_trans_probs[0][beta],
                    # Gamma = 1
                    V[1][i-1] + log_intergenic_probs[sequence[i]] + log_trans_probs[1][beta],
                    # Gamma = 2
                    V[2][i-1] + log_intergenic_probs[sequence[i]] + log_trans_probs[2][beta],
                    # Gamma = 3
                    V[3][i-1] + log_intergenic_probs[sequence[i]] + log_trans_probs[3][beta]
                ])

                V[beta][i] = np.max(prev_col)
                back_point[beta][i] = [np.argmax(prev_col), i-1]

            elif beta == 1:
                prev_col = np.array([
                    # Gamma = 0
                    V[0][i-3] + log_start_probs[sequence[i-2:i+1]] + log_trans_probs[0][beta],
                    # Gamma = 1
                    V[1][i-3] + log_start_probs[sequence[i-2:i+1]] + log_trans_probs[1][beta],
                    # Gamma = 2
                    V[2][i-3] + log_start_probs[sequence[i-2:i+1]] + log_trans_probs[2][beta],
                    # Gamma = 3
                    V[3][i-3] + log_start_probs[sequence[i-2:i+1]] + log_trans_probs[3][beta]
                ])
                
                V[beta][i] = np.max(prev_col)
                back_point[beta][i] = [np.argmax(prev_col), i-3]

            elif beta == 2:
                prev_col = np.array([
                    # Gamma = 0
                    V[0][i-3] + log_mid_probs[sequence[i-2:i+1]] + log_trans_probs[0][beta],
                    # Gamma = 1
                    V[1][i-3] + log_mid_probs[sequence[i-2:i+1]] + log_trans_probs[1][beta],
                    # Gamma = 2
                    V[2][i-3] + log_mid_probs[sequence[i-2:i+1]] + log_trans_probs[2][beta],
                    # Gamma = 3
                    V[3][i-3] + log_mid_probs[sequence[i-2:i+1]] + log_trans_probs[3][beta]
                ])
                
                V[beta][i] = np.max(prev_col)
                back_point[beta][i] = [np.argmax(prev_col), i-3]
            
            else: # beta == 3
                prev_col = np.array([
                    # Gamma = 0
                    V[0][i-3] + log_stop_probs[sequence[i-2:i+1]] + log_trans_probs[0][beta],
                    # Gamma = 1
                    V[1][i-3] + log_stop_probs[sequence[i-2:i+1]] + log_trans_probs[1][beta],
                    # Gamma = 2
                    V[2][i-3] + log_stop_probs[sequence[i-2:i+1]] + log_trans_probs[2][beta],
                    # Gamma = 3
                    V[3][i-3] + log_stop_probs[sequence[i-2:i+1]] + log_trans_probs[3][beta]
                ])
                
                V[beta][i] = np.max(prev_col)
                back_point[beta][i] = [np.argmax(prev_col), i-3]
    
    return V, back_point

# Backtrack portion of the algorithm
def backtrack(sequence, V, back_point):
    path = ''
    
    index = len(sequence)-1
    state = np.argmax(V[:,index])
    if state == 0:
        path = 'I' + path
        index -= 1
    else:
        path = 'GGG' + path
        index -= 3

    pointer = back_point[state][len(sequence)-1]
    while index >= 0:
        if pointer[0] == 0:
            path = 'I' + path
            index -= 1
        else:
            path = 'GGG' + path
            index -= 3
        pointer = back_point[pointer[0]][pointer[1]]

    return path

def get_annotations(path):
    annotations = []

    index = 0
    while index < len(path):
        if path[index] == 'G':
            start = index + 1
            index += 1
            while (index<len(path)) and (path[index] == 'G'):
                index += 1
            end = index 
            annot = Annot(start, end)
            annotations.append(annot)
        else:
            index +=1
    
    return annotations

# Make the gene predictions using the Viterbi algorithm
def predict_annotations(genome, log_trans_probs, log_intergenic_probs,
                        log_start_probs, log_mid_probs, log_stop_probs):
    
    for seq in genome:
        V, back_point = viterbi(seq.sequence, log_trans_probs, log_intergenic_probs,
                                log_start_probs, log_mid_probs, log_stop_probs)
    
        path = backtrack(seq.sequence, V, back_point)
          
        annotations = get_annotations(path)        
        seq.annotations = annotations
    
    return genome

# Output the gene predictions as a gff3 file
def write_gff3(genome):
    
    with open('viterbi_out.gff3', 'a') as file_out:
        file_out.write('##gff-version 3\n')        
        for seq in genome:
            for annot in seq.annotations:
                entry = (seq.name + "\tena\tCDS\t" + 
                        str(annot.start) + "\t" + str(annot.end) +
                        "\t.\t+\t0\t.\n"
                        )
                
                file_out.write(entry)

if __name__ == "__main__":
    config_path = sys.argv[1]
    fasta_path = sys.argv[2]

    genome = read_fasta(fasta_path)

    # Open the config.json file to extract the necessary data for the 
    # Viterbi algorithm
    data = {}
    with open(config_path, 'r') as file_in:
        data = json.load(file_in)

    avg_intergenic_len = data['avg_intergenic_length']
    avg_genic_len = data['avg_genic_length']
    intergenic_freq = data['intergenic_frequencies']
    start_freq = data['start_frequencies']
    mid_freq = data['middle_frequencies']
    stop_freq = data['stop_frequencies']

    log_trans_probs = get_log_transition_probs(avg_intergenic_len, avg_genic_len)

    log_intergenic_probs = get_log_emmission_probs(intergenic_freq)
    log_start_probs = get_log_emmission_probs(start_freq)
    log_mid_probs = get_log_emmission_probs(mid_freq)
    log_stop_probs = get_log_emmission_probs(stop_freq)

    # Make predictions
    genome = predict_annotations(genome, log_trans_probs, log_intergenic_probs,
                                log_start_probs, log_mid_probs, log_stop_probs)
    
    write_gff3(genome)
