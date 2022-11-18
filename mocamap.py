import numpy as np
import screed 
from collections import defaultdict
import sys

'''
MoCaMap
By Alejandro Castellanos and Santiago Morales
Required arguments: 
    Reference Genome
    Long genomic reads sequencing data
    Epsilon
'''

# Invertible Hash function
def InvertibleHash(x, p):
    m = (2 ** p) - 1
    x = ((~x) + (x << 21)) & m
    x = x ^ (x >> 24)
    x = (x + (x << 3) + (x << 8)) & m
    x = x ^ (x >> 14)
    x = (x + (x << 2) + (x << 4)) & m
    x = x ^ (x >> 28)
    x = (x + (x << 31)) & m
    return x

# Natural hash
def NaturalHash(kmer, k):
    values =  {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    x = 0
    for i in range(k):
        x += values[kmer[k - 1 - i]]  * (4 ** i)
    return x

# Composition of string and integer hash to avoid errors with Poly-A's
def Phi(kmer, k):
    return InvertibleHash(NaturalHash(kmer, k), 2 * k)

# Compute minimizers
def LocalMinimizers(queue):
    i = 1
    ranked = sorted(queue)
    t = ranked[0][0]
    if t < np.Inf:
        for M in ranked[1:]:
            if M[0] == t: 
                i += 1
            else:
                break 
        return ranked[:i]
    return []


def MinimizerSketch(s, w, k): 
    queue = [] 
    M = []
    for i in range(w):
        kmer = s[i: i + k]
        rckmer = screed.rc(kmer)
        u = Phi(kmer, k)
        v = Phi(rckmer, k)
        if u < v:
            queue.append((u, i, 0))
        if u == v: 
            queue.append((np.Inf, -1, -1))
        if u > v:
            queue.append((v, i, 1))
    M.extend(LocalMinimizers(queue))


    for i in range(w, len(s) - k + 1):
        kmer = s[i: i + k]
        rckmer = screed.rc(kmer)
        u = Phi(kmer, k)
        v = Phi(rckmer, k)
        if u < v:
            queue.append((u, i, 0))
        if u == v: 
            queue.append((np.Inf, -1, -1))
        if u > v:
            queue.append((v, i, 1))
        
        queue.pop(0)

        lastm = M[-1][0]
        lasti = M[-1][1]
        m = queue[-1][0]
        i = queue[-1][1]
        if m <= lastm:
            M.append(queue[-1])
        elif i - lasti >= w:
            M.extend(LocalMinimizers(queue))
    
    return M

# Index target sequences
def Index(T, w, k):
    A = []
    for t in range(len(T)):
        M = MinimizerSketch(T[t], w, k)
        for minimizer in M:
            h, i, r = minimizer
            seqminimizer = (h, t, i, r)
            A.append(seqminimizer)
    A.sort()
    H = defaultdict(list)
    for a in A:
        H[a[0]] = []
    for a in A:
        H[a[0]].append((a[1], a[2], a[3])) 
    return H

# Map a query sequence

def Map(H, q, w, k, epsilon):
    A = []
    M = MinimizerSketch(q, w, k)
    for minimizer in M:
        h, i, r = minimizer
        if h in H.keys():
            for hminimizer in H[h]:
                t, i_h, r_h = hminimizer
                if r == r_h:
                    A.append((t, 0, i - i_h, i_h)) 
                else:
                    A.append((t, 1, i + i_h, i_h))
    
    A.sort()
    b = 0
    maxchain = []
    for e in range(len(A)):
        if e == len(A) - 1 or A[e + 1][0] != A[e][0] or (
            A[e + 1][1] != A[e][1] or A[e + 1][2] - A[e][2] >= epsilon):
            chain = A[b: (e + 1)]
            if len(chain) >= 4 and len(chain) > len(maxchain):
                    maxchain = chain
            b = e + 1
    return maxchain

# Get sequenes from file
def GetSequencesFromFile(fileName):
    sequences = []
    for record in screed.open(fileName):
        sequences.append(record.sequence)
    return sequences

target_sequence_file_name = sys.argv[1]
sequence_file_name = sys.argv[2]
epsilon = sys.argv[3] 

target_sequence = [GetSequencesFromFile(target_sequence_file_name)[0]]
w = 5
k = 15


H = Index(target_sequence, w, k)

# Output sam_file
with open(sequence_file_name.rstrip('.fasta') + '_mocamap.sam', 'w') as f:
    for record in screed.open(sequence_file_name):
        Seq = record.sequence
        maxchain = Map(H, Seq, w, k, epsilon)
        Seq_ID = record.name
        if not maxchain:
            flag = 4
        elif maxchain[0][1] == 0:
            flag = 0
        else:
            flag = 16
        if maxchain:
            pos = maxchain[0][3] + 1
        else: 
            pos = 0
        f.write(f'{Seq_ID}\t{flag}\t*\t{pos}\t0\t*\t*\t0\t0\t{Seq}\t*\n')


print('Mapping completed. Thank you for using MoCaMap.')