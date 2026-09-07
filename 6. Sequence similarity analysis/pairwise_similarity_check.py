import sys
import os
import numpy as np
from pyfaidx import Fasta

base_dir = sys.argv[1]

def get_similarity(data, sub1, sub2):
    seq1 = [data[g][:].seq for g in data.keys() if sub1 in g][0]
    seq2 = [data[g][:].seq for g in data.keys() if sub2 in g][0]

    same = 0
    valid = 0

    for x, y in zip(seq1, seq2):
        if x.upper() == 'N' or y.upper() == 'N':
            continue
        valid += 1
        if x.upper() == y.upper():
            same += 1
    
    return float(same) / valid if valid > 0 else 0

# Triplet numbers
all_files = os.listdir(base_dir)
triplet_numbers = set()
for f in all_files:
    if f.endswith('.aln') and ('AB_pair_' in f or 'BR_pair_' in f or 'AR_pair_' in f):
        num = f.split('_')[-1].replace('.aln', '')
        triplet_numbers.add(int(num))

# Comparisons in dict
comparisons = [('AB_pair', 'Av', 'Bv'), ('AR_pair', 'Av', 'SECCE'), ('BR_pair', 'Bv', 'SECCE')]


# Calculation for each triplet
for num in triplet_numbers:
    scores = []

    for prefix, sub1, sub2 in comparisons:
        filename = f'{prefix}_{num}.aln'
        filepath = os.path.join(base_dir, filename)

        if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
            data = Fasta(filepath)
            score = get_similarity(data, sub1, sub2)
            scores.append(score)
    
    if len(scores) == 3:
        mean_score = np.mean(scores)
        print(f'triplet_{num} {mean_score}')