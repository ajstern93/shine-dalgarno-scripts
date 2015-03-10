import sys
import os
import re
import glob

from collections import OrderedDict

import Bio.Alphabet as AB
from Bio.Seq import Seq
import Bio.SeqIO as SeqIO
import Bio.SeqUtils as SeqUtils
import numpy as np
from random import *
import math
from scipy.optimize import minimize, show_options
from scipy import stats
from sklearn import mixture
import bimodal_test


def create_temp_file(sequences, antiSequences):
    
    '''
    stores ViennaRNA data in a temporary file for parsing
    '''
    
    idx = 0
    with open('temp.seq', 'w') as temp:
        for site in sequences:
            temp.write('>temporary RNA cofold file\n')
            temp.write('%s&%s\n' %(sequences[idx], antiSequences[idx]))
            idx+=1
    return 

def get_binding_energy(sequence, antiSequence, directory = ''):
    
    '''
    runs ViennaRNA for a list of sequences and antisequences
    
    directory can be modified if desired
    '''
    
    create_temp_file(sequence,antiSequence)
    f = os.popen('RNAcofold -p <'+directory+'temp.seq')
    output = f.readlines()
    idx = 3
    energies = []
    while idx<len(output):
        energies.append(float(re.search('[0-9\-]+\.[0-9]*', output[idx]).group()))
        idx+=5
    f.close()
    return energies

def get_min_energy(genome,feature,bounds,width,random=False):
    
    '''
    returns the lowest binding energy of a subsequence
    of a defined width in a given UTR
    '''
    
    if feature.location.strand==1:
        utrStart = int(feature.location.start)+bounds[1]
        utrEnd = int(feature.location.start)+bounds[0]
        utr = genome.seq.__str__()[utrStart:utrEnd]

    elif feature.location.strand==-1:
        utrStart = int(feature.location.end)-bounds[0]
        utrEnd = int(feature.location.end)-bounds[1]
        utrRev = Seq(genome.seq.__str__()[utrStart:utrEnd],AB.generic_dna)
        utr = utrRev.reverse_complement().__str__()

    # cycle through kmers
    oldKmers = []
    newKmers = []
    newKmerASDs = []
    
    for i in range(0,len(utr)-width+1):
        kmer = utr[i:(i+width)]
        if kmer not in stored:
            newKmers.append(kmer)
            newKmerASDs.append(rRNA[width])
        else:
            oldKmers.append(kmer)
            
    if len(newKmers)>=1:
        energies = get_binding_energy(newKmers,newKmerASDs)
        for (idx,energy) in list(enumerate(energies)):
            stored[newKmers[idx]]=energy
    else:
        energies = []
        
    for kmer in oldKmers:
        energies.append(stored[kmer])
    
    if len(energies)>0:    
        return min(energies)
    else:
        return None


def get_rand_energy(genome,feature,bounds,width,trials):
    
    '''
    randomizes a UTR and returns the lowest 
    binding energy for a defined width in that
    sequence. Trials and size of UTR may be specified
    '''
    
    UTR_SIZE = 30 # This is an arbitrary number of nt before the start codon.
    
    if feature.location.strand==1:
        utrStart = int(feature.location.start)-UTR_SIZE
        utrEnd = int(feature.location.start)
        utr = genome.seq.__str__()[utrStart:utrEnd]
    elif feature.location.strand==-1:
        utrStart = int(feature.location.end)
        utrEnd = int(feature.location.end)+UTR_SIZE
        utrRev = Seq(genome.seq.__str__()[utrStart:utrEnd],AB.generic_dna)
        utr = utrRev.reverse_complement().__str__()

    minEnergies = []
    for i in range(0,trials):
        startFrame = 0
        normal = list(utr)
        shuffle(normal)
        shuffled = ''.join(normal)

        energies = []
        # cycle through kmers 
        oldKmers = []
        newKmers = []
        newKmerASDs = []
    
        for i in range(0,len(utr)-width+1):
            kmer = shuffled[i:(i+width)]
            if kmer not in stored:
                newKmers.append(kmer)
                newKmerASDs.append(rRNA[width])
            else:
                oldKmers.append(kmer)
            
        if len(newKmers)>=1:
            energies = get_binding_energy(newKmers,newKmerASDs)
            for (idx,energy) in list(enumerate(energies)):
                stored[newKmers[idx]]=energy
        
        for kmer in oldKmers:
            energies.append(stored[kmer])
        if len(energies)>0:    
            minEnergies.append(min(energies))
        else:
            return None 
    return minEnergies
    
def bind(records,bounds,width,trials):
    
    '''
    simulates binding of a given UTR region (width and bounds)
    for a whole genome and trials# of randomized genomes
    '''
    storedLen.append(len(stored))
    data = {}
    randoms = {}
    n = 0
    
    for record in records:
        if record.description.lower().find('plasmid') == -1:
            for feature in record.features:
                if feature.type=='CDS':
                    n+=1	
                    gene = 'gene_#'+str(n) # feature.qualifiers['gene'][0]
                    data[gene]=get_min_energy(record,feature,bounds,width)
                    randoms[gene]=get_rand_energy(record,feature,bounds,width,trials)
                    if randoms[gene]==None or data[gene]==None:
                        del randoms[gene]
                        del data[gene]
    randNumList = []
    
    randAvgList = []
    randStdList = []   
    for idx in range(0,trials):
        numRand = 0
        trialList = []

        for k in randoms:
            trialList.append(randoms[k][idx])
            if randoms[k][idx]<-4: #arbitrary cutoff
                numRand+=1
        randAvgList.append(np.mean(trialList))
        randStdList.append(np.std(trialList))
        randNumList.append(numRand)

    numMu = np.mean(randNumList)
    numStd = np.std(randNumList)
        
    avgMu = np.mean(randAvgList)
    avgStd = np.std(randAvgList)

    stdMu = np.mean(randStdList)
    stdStd = np.std(randStdList)
        
    numReal = 0

    realEnergies = np.abs(np.array(data.values()))# mean-center
    trimmedValues = np.array(filter(lambda x:x>0.5,np.abs(np.array(realEnergies))));

    for energy in realEnergies:
        if energy<-4:
            numReal+=1

    # Likelihood testing of uni and bimodal distributions

    results_unimodal, results_bimodal   = bimodal_test.get_gaussian_fits(trimmedValues)
    aic1, aic2, p_1_is_better           = bimodal_test.compare_models(results_unimodal,results_bimodal)


    zscores = [(numReal-numMu)/numStd,(np.mean(realEnergies)-avgMu)/avgStd,(np.std(realEnergies)-stdMu)/stdStd]

    return OrderedDict([	('z-score(num<thresh)',		zscores[0]),
		                    ('z-score(avgE)',   		zscores[1]),
		                    ('z-score(stdE)',   		zscores[2]),
		                    ('likelihood unimodal',		results_unimodal.fun),
		                    ('likelihood bimodal fit',	results_bimodal.fun),
		                    ('AIC unimodal',			aic1),
		                    ('AIC bimodal',			    aic2),
		                    ('p(model1 fits better)',	p_1_is_better)       ])

		




rRNA = {4: 'CCUC',
        5: 'CCUCC',
        6: 'UCCUCC',
        7: 'UUCCUCC',
        8: 'AUUCCUCC',
        9: 'AUUCCUCCA',
        10: 'AUUCCUCCAC',
        11: 'AUUCCUCCACU',
        12: 'AUUCCUCCACUA'}

stored = {}
storedLen = []
rComp = {'A':'U','U':'A','G':'C','C':'G'}

pathToGenbanks = sys.argv[1]
bounds = [int(sys.argv[2]),int(sys.argv[3])] 
width = 6
trials = int(sys.argv[4])
fileRates = sys.argv[5]

growthRates = open(fileRates,'rU').readlines()
results = open('results.csv','w')

first = True
for line in growthRates:

    columns = line.split(',')
    gbk = columns[2]
    if 'GCA' in gbk:
        
        records = list(SeqIO.parse(pathToGenbanks+gbk+'.gbk','genbank'))
        output = bind(records,bounds,width,trials)

    # if first pass through, print out the column labels
        if first:
            firstLine = ''
            for label in ['name1','name2','filename','growth temp','growth rate']:
                firstLine = firstLine+label+'\t'
            for label in output:
                firstLine = firstLine+label+'\t'
            results.write(firstLine.strip('\t')+'\n')
            first = False

	
        flat1 = '\t'.join(line.strip('\n').split(','))
        vals = output.values()
        flat2 = '\t'.join(str(x) for x in vals)
        outLine = '\t'.join([flat1,flat2])
        
        print outLine
        results.write(outLine+'\n')




