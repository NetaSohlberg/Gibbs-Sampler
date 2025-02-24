# Motif Finding - Gibbs Sampler Algorithm

## Overview

This project implements the Gibbs Sampler Algorithm for motif finding in DNA sequences. The Gibbs Sampler is a Markov Chain Monte Carlo (MCMC) algorithm used to obtain a sequence of approximated observations from a probability distribution when direct sampling is difficult.

## Authors

Yeynit Asraf

Neta Sohlberg

## Description

The script uses the Gibbs Sampler algorithm to identify recurring motifs in multiple DNA sequences. Given a set of DNA sequences and a motif length (k), the algorithm finds the most probable motifs present in all sequences using iterative sampling and probability calculations.


## Usage

### Input

The script reads DNA sequences from a text file where each line represents a DNA sequence.

#### Running the Script

python gibbs_sampler.py

Upon execution, the script will prompt the user for the motif length (k). It then runs the Gibbs Sampler algorithm multiple times and outputs the best set of motifs found.

## Example

Input File (dna_seq.txt):
```python
TTACCTTAAC
GATGTCTGTC
CCGGCGTTAG
CACTAACGAG
CGTCAGAGGT
```
Sample Output:
```python
k: 4
['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']
```

## Functions
```python
FindMotifs(dna, k)
```
Selects a random k-mer from each DNA sequence.
```python
CountOccurances(motifs)
```
Generates a position-specific count matrix for motifs.
```python
CalcScore(motifs)
```
Calculates the score of mismatches in motifs.
```python
FindProfile(motifs)
```
Generates a probability matrix (PSSM) for motifs.
```python
ProfileProb(sequence, profile)
```
Computes the probabilities of k-mers in a sequence based on a given profile.
```python
GibbsSampler(dna, k, N)
```
Runs the Gibbs Sampling algorithm for N iterations.
```python
GetDNA(filename)
```
Reads DNA sequences from a file.
```python
repeatGibbsSampler(dna, k, N, repeats)
```
Runs the Gibbs Sampler multiple times to find the best motifs.

## Parameters

```dna```: List of DNA sequences.

```k```: Length of the motif to find.

```N```: Number of iterations in Gibbs Sampling.

```repeats```: Number of times the Gibbs Sampler is repeated for convergence.
