#!/usr/bin/python
# -*- coding: ascii -*-
#import warnings

#warnings.simplefilter("error")
import sys
import os
import re
import glob
import Bio.Alphabet as AB
from Bio.Seq import Seq
import Bio.SeqIO as SeqIO
import Bio.SeqUtils as SeqUtils
import numpy as np
import math
from random import *
import matplotlib.pyplot as plt


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
    energy = []
    while idx<len(output):
        energy.append(float(re.search('[0-9\-]+\.[0-9]*', output[idx]).group()))
        idx+=5
    f.close()
    
    idx = 0
    for seq in sequence:
        if idx>len(energy)-1:
            print idx
            print sequence[0]
            print antiSequence[0]
        stored[seq]=energy[idx]
        idx+=1
    return energy

def bind(records,bounds,width,trials):
    
    '''
    simulates binding of a given UTR region (width and bounds)
    for a whole genome and trials# of randomized genomes
    '''
    print len(stored)
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
    return [data,randoms]

def get_min_energy(genome,feature,bounds,width,random=False):
    
    '''
    returns the lowest binding energy of a subsequence
    of a defined width in a given UTR
    '''
    
    if feature.location.strand==1:
        utrStart = int(feature.location.start)+bounds[1]
        utrEnd = int(feature.location.start)+bounds[0]
        utr = genome.seq.__str__()[utrStart:utrEnd]
       # print feature.qualifiers['gene'], feature.location.strand, utr
    elif feature.location.strand==-1:
        utrStart = int(feature.location.end)-bounds[0]
        utrEnd = int(feature.location.end)-bounds[1]
        utrRev = Seq(genome.seq.__str__()[utrStart:utrEnd],AB.generic_dna)
        utr = utrRev.reverse_complement().__str__()
       # print feature.qualifiers['gene'], utr

#    rRNAList = []
#    seqList = []
#    startFrame = 0
#    energies = []
#    while startFrame+width<abs(bounds[1]-bounds[0])+1:
#        currentFrame = utr[startFrame:startFrame+width]
#        if currentFrame not in stored:
#            seqList.append(currentFrame)
#            rRNAList.append(rRNA[width])
#        else: 
#            energies.append(stored[currentFrame])
#        startFrame+=1
#    newEnergies = get_binding_energy(seqList,rRNAList)
#    for e in newEnergies:
#        energies.append(e)
#    return min(energies)
    
    ###### We check to see if binding of the kmer to the conserved
    ###### antimer is already booked:
    # print len(utr)
    if utr not in stored:
        stored[utr] = get_binding_energy([utr],[rRNA[width]])[0]
    return stored[utr]


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

   # rRNAList = []
   # seqList = []
    energies = []
    for i in range(0,trials):
        startFrame = 0
        normal = list(utr)
        shuffle(normal)
        shuffled = ''.join(normal)
#        trialLen = 0
        randSD = shuffled[-bounds[0]:-bounds[0]+width]
#        if i == 0: print randSD	
#	print len(randSD)	# sanity check
        if randSD not in stored:
            newEnergy = get_binding_energy([randSD],[rRNA[width]])
            stored[randSD] = newEnergy[0]
            energies.append(newEnergy[0])
        else:
            energies.append(stored[randSD])
    return energies    

#        while startFrame+width<abs(bounds[1]-bounds[0])+1:
#            currentFrame = shuffled[-bounds[0]+startFrame:-bounds[0]+startFrame+width]
#            if currentFrame not in stored:
#                seqList.append(currentFrame)
#                rRNAList.append(rRNA[width])
#                energies.append(1)
#            else: energies.append(stored[currentFrame])
#            startFrame+=1
#            trialLen = max(trialLen,startFrame)
#    
#    newEnergies = get_binding_energy(seqList,rRNAList)
#    for e in newEnergies:
#        idx = 0
#        while energies[idx]!=1:
#            idx+=1
#        energies[idx] = e
#    minEnergies = []
#    for i in range(0,trials):
#        minEnergies.append(min(energies[trialLen*i:trialLen*(i+1)]))
#    return minEnergies

def make_models(records,bounds,trials):    # bounds = window (e.g. [0,-20])

    '''
    finds z-score of SD genes in whole genome 
    at a given definition (width, energy) of SD
    '''
    
    widthSpace = np.linspace(4,11,8)
    rsSpace = np.linspace(bounds[0],bounds[1],abs(bounds[1]-bounds[0])+1)
    mydata = np.empty((np.size(widthSpace),np.size(rsSpace)))
    i = 0

    for width in widthSpace:
        print width
        j=0
        
        for rs in rsSpace:
            [reals,randoms] = bind(records,[bounds[0]-j,bounds[0]-j-int(width)],int(width),trials)
            
          
            
            
            #print randoms
            #randSD = []
            #for idx in range(0,trials):
            #    randSD.append(0)
            #    for gene in randoms:
            #        if randoms[gene][idx]<minBind[width][0]*thresh: 
            #            randSD[idx]+=1 #PDF vs CDF??
            randAvgList = []
            trialList = []
            #randTrials1 = []
	        #randTrials2 = []
            #randTrials3 = []
        
            for idx in range(0,trials):
                for k in randoms:
                    energy = randoms[k][idx]
                    trialList.append(energy)
                    #if idx==1:
			            #randTrials1.append(energy)
		            #if idx==2:
			            #randTrials2.append(energy)
		            #if idx==3:
			            #randTrials3.append(energy)
                randAvgList.append(np.mean(trialList))
                trialList = []
            
            mu = np.mean(randAvgList)
            std = np.std(randAvgList)
            
            #num = 0
            #for gene in reals:
            #    if reals[gene]<minBind[width][0]*thresh: 
            #        num+=1
            realList = []
            for k in reals:
                realList.append(reals[k])
            avg = np.mean(realList)

            if math.isnan((float(avg)-mu)/std) or std==0:
                mydata[i,j]=0 ## code to special colors?
            else:
                mydata[i,j] = (float(avg)-mu)/std
            
            #if i==0 and j==6:
	    	    #plt.hist(realList,bins=40,alpha=0.4,label='actual')
		        #plt.hist(randTrials1,bins=40,alpha=0.4,label='random 1')
		        #plt.hist(randTrials2,bins=40,alpha=0.4,label='random 2')
		        #plt.hist(randTrials3,bins=40,alpha=0.4,label='random 3')
		        #plt.ylabel('Number of genes')
    		    #plt.xlabel('Energy')
    		    #ax.set_xticklabels(range(bounds[0],bounds[1]-1,-1))# whatever you want
    		    #ax.set_yticklabels(range(4,13,1))
    		    #plt.savefig('MyHist.png')
    	
            j+=1
        i+=1
    return mydata

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
# minBind = {}
#for width in rRNA:
#    maxComp = ''
#    for nt in rRNA[width]:
#       maxComp+=rComp[nt]
#    minBind[width] = get_binding_energy([maxComp[::-1]],[rRNA[width]])

# takes path input of form 'RELATIVE_PATH/TO/DIR_CONTAINING_GB/'

pathToGenbanks = sys.argv[1]
bounds = [int(sys.argv[2]),int(sys.argv[3])]
trials = int(sys.argv[4])
fileRates = sys.argv[5]

results = open('results.csv','w')
storedData = open('stored_status.txt','w')

#ascii = sorted(glob.glob(pathToGenbanks+'*.gbk')) 
#for count,i in enumerate(ascii):
#    print i
#genome = SeqIO.read(i,'genbank')

growthRates = open(fileRates,'rU').readlines()
for line in growthRates:

    columns = line.split(',')
    gbk = columns[2]
    if 'GCA' in gbk:
        records = list(SeqIO.parse(pathToGenbanks+gbk+'.gbk','genbank'))
        mydata = make_models(records,bounds,trials)
        ### if matching genbank, get growth rate data and store with:
        ########################### 1. min zscore
        ########################### 2. width at min zscore
        ########################### 3. relative spacing at min zscore
        minZ = np.min(mydata)
        minIndices = np.argwhere(minZ==mydata)[0]
        widthIndex = minIndices[0]
        spacingIndex = minIndices[1]
        outLine = ','.join([line.strip('\n'),str(minZ),str(widthIndex+4),str(bounds[0]-spacingIndex),'\n'])
        print outLine
        results.write(outLine)
      
storedData.write(str(storedLen))     
    
    
    
    
    
    #fig,ax=plt.subplots()
    #im=ax.pcolormesh(mydata)
    #fig.colorbar(im)
    #ax.axis('tight')
    #plt.ylabel('SD width')
    #plt.xlabel('Relative Spacing')
    #ax.set_xticklabels(range(bounds[0],bounds[1]-1,-1))# whatever you want
    #ax.set_yticklabels(range(4,13,1))
    #plt.savefig(i.split('/')[2].strip('.gb')+'.png')

    

















