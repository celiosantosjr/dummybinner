import sys
import gzip
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm
from Bio import SeqIO
from itertools import product
from sklearn.cluster import OPTICS
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import is_valid_y
from multiprocessing import Pool, cpu_count

def dummy_binner(mode, sequences, coverages=None, minlen=1000, k=4, threshold=0.1):
    """
    A simple sequence binning function based on tetranucleotide frequencies
    
    Args:
    mode (str): Mode of binner targets (contigs/reads)
    sequences (list): A list of sequences (strings)
    minlen (int): Minimum sequence length
    k (int): The length of the k-mers to use for frequency calculations (default: 4)
    threshold (float): The threshold above which sequences are considered to belong to the same bin (default: 0.1)
    coverages (dataframe): Dataframe with sequences in the indexes and coverage in the columns (optional)
    
    Returns:
    list: A list of lists, where each sub-list contains the indices of the sequences in that bin
    """
    
    print('# Calculate tetranucleotide frequencies, GC% and coverage for each sequence (if available for contigs)')
    freqs = []
    for i, sequence in tqdm(enumerate(sequences), total=len(sequences)):
        if len(sequence) > minlen:
            kmerlist = {''.join(x): 0 for x in product('ACTGN', repeat=k)}
            for idx in range(len(sequence) - k + 1):
                kmerlist[sequence[idx:idx+k]] += 1
            gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
            freq = [kmerlist.get(x, 0)/(len(sequence) - k + 1) for x in kmerlist]
            freq.append(gc)
            if isinstance(coverages, pd.DataFrame):
                freq += coverages.iloc[i].tolist()
            freqs.append(freq)
        
    if mode == 'contigs':
        print('# Calculate pairwise distances between sequences based on tetranucleotide frequencies')
        dists = squareform(pdist(freqs))  
        dists = squareform(dists)     
        valid = is_valid_y(dists)
        
        if not valid:
            print("# dists matrix is not properly condensed!")
        
        print('# Perform hierarchical clustering on the distances')
        linkage = hierarchy.linkage(dists, method='average')
        clusters = hierarchy.fcluster(linkage, threshold, criterion='distance')
    
    if mode == 'reads':
        print('# Calculate clusters using OPTICS algorithm')
        clustering = OPTICS(min_samples=10).fit(freqs)
        clusters = clustering.labels_            

    print('# Organize sequences into bins based on clustering')
    bins = [[] for _ in range(max(clusters))]
    for i, cluster in enumerate(clusters):
         bins[cluster-1].append(i)
           
    return bins
    
def merge_reads(read_pair):
    r1, r2 = read_pair
    seq = r1.seq + r2.seq.reverse_complement()
    return (r1.id, str(seq))

def load_seqs(infile, infile2=None, mode='contigs'):
    print('# Loading sequences into memory')
    hs = []
    if mode == 'contigs':
        if infile.endswith('gz'):
            infile = gzip.open(infile, 'rt')
        for record in SeqIO.parse(infile, 'fasta'):
            hs.append((record.id, str(record.seq)))

    elif mode == 'reads' and infile2:
        if infile.endswith('gz'):
            infile = gzip.open(infile, 'rt')
        if infile2.endswith('gz'):
            infile2 = gzip.open(infile2, 'rt')

        infile = SeqIO.parse(infile, 'fastq')
        infile2 = SeqIO.parse(infile2, 'fastq')    

        print('# Merge the paired reads using multiple processes')
        with Pool(processes=cpu_count()) as pool:
            hs = pool.map(merge_reads, zip(infile, infile2))
            print(f'It was obtained {len(hs)} merged pairs')

    elif mode == 'reads':
        if infile.endswith('gz'):
            infile = gzip.open(infile, 'rt')
        for record in SeqIO.parse(infile, 'fastq'):
            hs.append((record.id, str(record.seq)))
                
    return hs

def main(args):
    if (args.mode == 'contigs') and args.f:
        exit('ERROR: Mode set to contigs and reads file provided')

    if (args.mode == 'reads') and (args.infile):
        exit('ERROR: Mode set to reads and contigs file provided')
            
    if (args.mode == 'reads') and (args.minlen == 1000):
        args.minlen = 50
        
    if args.mode == 'contigs':
        hs = load_seqs(args.infile, mode=args.mode)
    elif args.mode == 'reads':
        if args.f and args.r:
            hs = load_seqs(infile=args.f, infile2=args.r, mode=args.mode)
        elif args.f:
            hs = load_seqs(infile1=args.f, mode=args.mode)
        
    if args.cov:
        cov = pd.read_table(args.cov, header=None)
        cov = cov.drop(0, axis=1)
        bins = dummy_binner(args.mode, [s for _, s in hs], cov, args.minlen, args.k, args.threshold)
    else:    
        bins = dummy_binner(args.mode, [s for _, s in hs], None, args.minlen, args.k, args.threshold)

    print('# Exporting clusters')
    if args.mode == 'contigs':
        for idx, sublist in enumerate(bins):
            if len(sublist) >= args.mincontigs:
                with open(f'{args.otag}__{idx}.fa', 'w') as handle:
                    for cidx in sublist:
                        h, s = hs[cidx]
                        handle.write(f'>{h}\n{s}\n')
    
    if args.mode == 'reads':
        for idx, sublist in enumerate(bins):
            with open(f'{args.otag}__reads_clusters.fasta', 'w') as handle:
                for cidx in sublist:
                    h, s = hs[cidx]
                    handle.write(f'>{h} | cluster_{cidx}\n{s}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple contig binning function based on tetranucleotide frequencies')
    parser.add_argument('--r', type=str, required=False, help='R2 reads file (fastq or fastq.gz)')
    parser.add_argument('--mode', type=str, required=True, help='Mode of operation, contigs or reads')
    parser.add_argument('--otag', type=str, default='binout', help='Output file tag (defaul: binout)')
    parser.add_argument('--minlen', type=int, default=1000, help='Minimum contig length (default: 1000bp)')
    parser.add_argument('--k', type=int, default=4, help='K-mer length to measure the contigs similarity (default: 4)')
    parser.add_argument('--threshold', type=int, default=0.1, help='Maximum distance to consider two contigs in the same bin (default: 0.1)')
    parser.add_argument('--mincontigs', type=int, default=2, help='Minimum number of contigs to print a bin (default: 2)')
    parser.add_argument('--cov', type=str, help='Coverage table without header and 2 columns (contig name and coverage) - optional') 

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--infile', type=str, help='Input FASTA file containing contigs')
    group.add_argument('--f', type=str, help='R1 reads file (fastq or fastq.gz)')

    args = parser.parse_args()
   
    main(args)
    
