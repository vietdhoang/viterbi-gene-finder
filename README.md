# Viterbi Gene Finder
Implemention of the Viterbi algorithm and its application to the gene finding hidden Markov model to predict genes in a bacterial genome

## Pipeline
1. Run `config.py`, which produces `config.json`. This JSON file contains the necessary genic and intergenic frequencies for the Viterbi algorithm.
```bash
python config.py PATH/TO/input_fasta PATH/TO/input_gff3
```
2. Run `viterbi.py`, which will use the Viterbi algorithm to make gene predictions.
```bash
python viterbi.py PATH/TO/config.json PATH/TO/input_fasta
```
