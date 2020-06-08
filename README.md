[![Anaconda-Server Badge](https://anaconda.org/bioconda/shark/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

# Shark-sbt
## Shark implementation with another data structure
Fast tool for mapping-free gene separation of reads, using Bloom filter.

## Dependencies
Shark-sbt requires Intel TBB (tested with version 2017 update 7)

## Download and Installation
To install the tool, run the following steps.

First, clone the repository and move into it.
```shell
git clone https://github.com/DavidePizzoliCoso/shark.git
make
```

## Usage
```
Usage: shark -r <references> -1 <sample1> [OPTIONAL ARGUMENTS]

Arguments:
      -r, --reference                   reference sequences in FASTA format (can be gzipped)
      -1, --sample1                     sample in FASTQ (can be gzipped)

Optional arguments:
      -h, --help                        display this help and exit
      -2, --sample2                     second sample in FASTQ (optional, can be gzipped)
      -o, --out1                        first output sample in FASTQ (default: sharked_sample.1)
      -p, --out2                        second output sample in FASTQ (default: sharked_sample.2)
      -k, --kmer-size                   size of the kmers to index (default:17, max:31)
      -c, --confidence                  confidence for associating a read to a gene (default:0.6)
      -b, --bf-size                     bloom filter size in Kb (default:1024)
      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)
      -s, --single                      report an association only if a single gene is found
      -t, --threads                     number of threads (default:1)
      -m, --method                      subject of the condition [base / kmer] (default: base)
      -x, --xxhash                      number of hash functions
      -v, --verbose                     verbose mode
```

## Output format

`shark-sbt` outputs to `stdout` a ssv file reporting associations between reads and genes.
The first element of each line is the name of the read whereas the second element is the gene identifier.

Reads in the samples that pass the filter step are stored in the files passed as argument to `-o` and `-p`.

## Example

A small example is provided in the example directory.
* `ENSG00000277117.fa` is the gene sequence of gene ENSG00000277117 in FASTA format
* `sample_1.fq` and `sample_2.fq` are a paired-end RNA-Seq sample simulated from genes ENSG00000277117 and ENSG00000275464

To filter out the reads sequenced from gene ENSG00000277117, run `shark-sbt` as follows:

```
./shark -r example/ENSG00000277117.fa -1 example/sample_1.fq \
                                      -2 example/sample_2.fq \
                                      -o example/sharked.sample_1.fq \
                                      -p example/sharked.sample_2.fq > example/ENSG00000277117.ssv
```

The results should be equal to: `example/ENSG00000277117.truth.ssv`, `example/sharked.sample_1.truth.fq`, and  `example/sharked.sample_2.truth.fq`.
