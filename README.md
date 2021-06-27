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
In this repo, example output predictions are provided for *Vibrio cholerae* and *Vibrio vulnificus*. *V. cholerae* and *V. vulnificus* are bacteria that cause cholera [1,2].

## References
1. Heidelberg, J., Eisen, J., Nelson, W. et al. DNA sequence of both chromosomes of the cholera pathogen Vibrio cholerae . Nature 406, 477–483 (2000). https://doi.org/10.1038/35020000
2. Centers for Disease Control and Prevention (CDC). “Infectious disease and dermatologic conditions in evacuees and rescue workers after Hurricane Katrina--multiple states, August-September, 2005.” MMWR. Morbidity and mortality weekly report vol. 54,38 (2005): 961-4.
