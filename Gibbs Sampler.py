'''
Motive Finding- Gibbs Sampler Algorithm
/ Yeynit Asraf and Neta Sohlberg

"In statistics, Gibbs sampling or a Gibbs sampler is a Markov chain Monte Carlo (MCMC) algorithm for obtaining a sequence of observations which are approximated from a specified multivariate probability distribution,
 when direct sampling is difficult." (Wikipedia)

I used Gibbs Sampler Algorithm to find similar motifs in several DNA sequences. 
'''

import random

"""params: list, int
function:chose a randomal kmer from every sequence in the list in length 'k' 
return: list of kmers
"""
def FindMotifs(dna, k):
    length=len(dna[0])
    kmers=[]
    for s in dna:
        index=random.randint(0,length-k)
        kmers.append(s[index:index+k])
    return kmers

"""params: list
function: generates a position specific counting matrix without adding psuidocount.
return: a list of frenquency lists for every index in the motifs
"""
def CountOccurances(motifs):
    counter=[]
    length=len(motifs[0])
    for i in range(length):
        freq=[0] * 4
        dic={'A':0,'T':1,'C':2,'G':3}
        for j in motifs:
            index=dic[j[i]]
            freq[index]+=1
        counter.append(freq)
    return counter

"""param: list
function: counts the number of mismatches between nucleotides  in the list "motifs" and returns the sum of the mismatch counts as the score.
return: the score of the mismatches""" 
def CalcScore(motifs):
    score=0
    counter=CountOccurances(motifs)
    length=len(motifs)
    for i in counter:
        score+=length-max(i)
    return score

"""param: list
function: generates a 4xk matrix containing the probabilities of PSSM
return: a list of profile lists
"""
def FindProfile(motifs):
    length=len(motifs)*2
    counter=CountOccurances(motifs)
    profile=[]
    for i in counter:
         profile+=[[(x+1)*1.0/length for x in i]]
    return profile

"""param: string,list of lists
function: computes the probabilities of all k-mers in the removed sequence
return: a list of probabilities
"""
def ProfileProb(sequence, profile):
    kmers=[]
    length=len(profile)
    for i in range(len(sequence)-length+1):
        kmers+=[sequence[i:i+length]]
    dic={'A':0,'T':1,'C':2,'G':3}
    prob=[]
    for kmer in kmers:
        compute=1.0
        for i in range(len(kmer)):
            compute*=profile[i][dic[kmer[i]]]
        prob.append(compute)
    return prob

"""param:list,int,int
function: The function implements the Gibbs sampling algorithm
return:a list of Best Motifs and the best score
"""
def GibbsSampler(dna, k, N):
    Motifs=FindMotifs(dna,k)
    BestMotifs = Motifs
    for i in range(N):
        index=random.randint(0,len(Motifs)-1)
        outMotifs=Motifs[index]
        Motifs.remove(outMotifs)
        profile=FindProfile(Motifs)
        score=ProfileProb(dna[index],profile)
        Motifs.insert(index,dna[index][score.index(max(score)):score.index(max(score))+k])
        if CalcScore(Motifs)<CalcScore(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs ,CalcScore(BestMotifs)        

"""param: string
function: The function open the file and make a list of the dns sequence
return: a list of DNA sequences
"""
def GetDNA(filename):
    f=open(filename)
    dna=[]
    for line in f:
        dna.append(line.strip())
    f.close()
    return dna

"""params: list, int, int, int
 function: This function repeats the whole process in order to reach convergence.
return: a list of kmers """

def repeatGibbsSampler(dna, k, N, repeats):
    bestMotifs = FindMotifs(dna,k)    
    bestScore = CalcScore(bestMotifs)    
    for i in range(repeats):
        (motifs, score) = GibbsSampler(dna, k, N)
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs


#main:
dna= GetDNA("dna_seq.txt")
N=20
repeats=40
k=int(input ("What is the length of the motif you are interested in? "))
print ("k: ",k)
print (repeatGibbsSampler(dna,k,N,repeats))

#examples:
"""

file:
AACC
ATCT
GGAC
TCAC

output:
['AC', 'AT', 'AC', 'AC']

########################################################

file:
TTACCTTAAC
GATGTCTGTC
CCGGCGTTAG
CACTAACGAG
CGTCAGAGGT

output:
['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']

########################################################

file:
TTTTCTGAAC
TTTTTCTGTC
TTTTCGTTAG
TTTTAACGAG
TTTTAGAGGT

output:
['TTTT', 'TTTT', 'TTTT', 'TTTT', 'TTTT']
"""
