# dummybinner

Dummy binner of contigs and reads that uses kmer frequencies, GC content and coverage.
Few dependencies, basic results.

## Citation

Dummy binner is GNU publicly licensed. Enjoy it.

In case you cite it:

```
Santos-Junior, C.D. (2023) Dummy binner - Using tetranucleotide frequencies, GC and
coverage to cluster contigs and reads. Software in GitHub Repository available
in <https://github.com/celiosantosjr/dummybinner>. Accessed in Month/Year.
```

## Installing

Dummy binner uses python>3.9 and a set of packages easilly installed using 

```
$ pip install biopython numpy pandas scipy scikit-learn tqdm
```

The versions required follow:

| Requirement | Version |
| :---: | :---:|
| biopython | 1.79 |
| numpy | 1.22.4 |
| pandas | 1.4.4 |
| scipy | 1.9.1 |
| sklearn | 0.0 |
| tqdm | 4.64.0 |

## Usage

```
usage: dummy_binner2.py [-h] [--r R] --mode MODE [--otag OTAG]
                        [--minlen MINLEN] [--k K] [--threshold THRESHOLD]
                        [--mincontigs MINCONTIGS] [--cov COV]
                        (--infile INFILE | --f F)
```

Dummy binner can be used in the reads mode, where it accepts the minimum sequence
length (--minlen), the R1 and R2 files (--f and --r, respectively), and you can
also specify the output tag for the output filename (--otag). The basic usage cases
are:

```
python3 dummy_binner2.py --f <R1.fq.gz> --r <R2.fq.gz> --mode reads
```


For contigs, the binning is a bit more parameterized, you should input a contigs
file (--infile), can also specify the minimum contig length to be accepted (--minlen)
and the minimum number of contigs in a bin to be considered (--mincontigs). You also
can alternatively give the kmer size (--k) and the maximum distance between two
contigs (--threshold). The basic usage cases are:

```
python3 dummy_binner2.py --infile <contigs.fa> --mode contigs
```

## Contact

In case any issues, please contact [me](mailto:celio.diasjunior@gmail.com).
