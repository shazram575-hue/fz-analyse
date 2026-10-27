"""Equivalent evaluation of DeepPurpose's native 8420 AAC features, not a new encoding."""
from collections import Counter,defaultdict
from itertools import product
import numpy as np
from DeepPurpose.pybiomed_helper import AALetter

AA=list(AALetter)
PAIRS=["".join(words) for words in product(AA,repeat=2)]
TRIPLES=["".join(words) for words in product(AA,repeat=3)]
PAIR_INDEX={word:index for index,word in enumerate(PAIRS)}
TRIPLE_INDEX={word:index for index,word in enumerate(TRIPLES)}

def nonoverlapping_counts(sequence,k):
    counts=defaultdict(int); next_start={}
    for position in range(len(sequence)-k+1):
        word=sequence[position:position+k]
        if position>=next_start.get(word,0):
            counts[word]+=1;next_start[word]=position+k
    return counts

def fast_aac(sequence):
    n=len(sequence)
    result=np.zeros(8420,dtype=np.float32)
    if n<2:return result
    counts=Counter(sequence)
    for index,letter in enumerate(AA):result[index]=round(float(counts[letter])/n*100,3)
    for word,count in nonoverlapping_counts(sequence,2).items():
        index=PAIR_INDEX.get(word)
        if index is not None:result[20+index]=round(float(count)/(n-1)*100,2)
    for word,count in nonoverlapping_counts(sequence,3).items():
        index=TRIPLE_INDEX.get(word)
        if index is not None:result[420+index]=count
    return result
