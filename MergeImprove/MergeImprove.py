'''
Created on Oct 26, 2012

@author: James Holt
'''

import pysam #@UnresolvedImport
#import Synchronize
import Pileup
import random
import logging
import argparse as ap
import util
import os
import sys
import re
import multiprocessing
from multiprocessing.managers import SyncManager

import pylab

DESC = "A merger of genome alignments."
VERSION = '0.2.1'
PKG_VERSION = '0.2.1'

#constant flags from the sam specs
MULTIPLE_SEGMENT_FLAG = 1 << 0 #0x01
#BOTH_ALIGNED_FLAG = 1 << 1#0x02
SEGMENT_UNMAPPED_FLAG = 1 << 2#0x04
OTHER_UNMAPPED_FLAG = 1 << 3#0x08
FIRST_SEGMENT_FLAG = 1 << 6#0x40
#SECOND_SEGMENT_FLAG = 1 << 7#0x80

#mergerTypes
RANDOM_MERGE = 0
UNION_MERGE = 1
PILEUP_MERGE = 2

#Save types
SAVE_UNION = 0
SAVE_KEPT_ALL = 1
SAVE_QUALITY = 2
SAVE_PILEUP = 3
SAVE_RANDOM = 4

#multithread stages
QUALITY_STAGE = 0
MERGE_PILEUP_STAGE = 1
CONDENSE_PILEUP_STAGE = 2
PILEUP_STAGE = 3
CLEANUP_STAGE = 4
FINISHED_STAGE = 5
END_OF_QUEUE = 'EOQ'

#tags from either alignment or lapels
OLD_CIGAR_TAG = 'OC'
OLD_MISMATCH_TAG = 'OM'
SNP_TAG = 's0'
INSERTION_TAG = 'i0'
DELETION_TAG = 'd0'
YA_TAG = 'YA'

#Tags added by me
PARENT_OF_ORIGIN_TAG = 'po'
CHOICE_TYPE_TAG = 'ct'
M_CIGAR_TAG = 'mc'
M_MISMATCH_TAG = 'mm'
P_CIGAR_TAG = 'pc'
P_MISMATCH_TAG = 'pm'
PILEUP_HI_TAG = 'pH'

#maternal side
M_SNP_TAG = 'ms'
M_IN_TAG = 'mi'
M_DEL_TAG = 'md'

#paternal side
P_SNP_TAG = 'ps'
P_IN_TAG = 'pi'
P_DEL_TAG = 'pd'

#constants for the clusters
MAX_CLUSTER_SIZE = 100

verbosity = False

class MergeWorker(multiprocessing.Process):
    
    def __init__(self, readsQueue, resultsQueue, isMaster, firstFilename, secondFilename, outputFilename, mergerType, isRandomFilter, outHeader,
                 managerAddress, numProcs, workerID, max_queue_size, pileupDict, keepAll):
        '''
        @param readsQueue - the shared queue where the master will place indicators for the workers to decide what they should process
        @param resultsQueue - the shared queue for storing statistics on the alignments from each worker
        @param isMaster - boolean, indicates if this process is the master or not
        @param firstFilename - the first fn to be merged
        @param secondFilename - the second fn to be merged
        @param outputFilename - the overall output filename where the final merge will be, worker won't use this unless a single process scenario
        @param mergerType - union, quality, pileup are the types offered right now
        @param isRandomFilter - boolean determining whether we random or union at the end
        @param outHeader - the header to use in any output files
        @param managerAddress - address so that we can connect and access synchronization primitives
        @param numProcs - number of processes
        @param workerID - integer ID, 0 means master
        @param max_queue_size - mostly for the master to tell it how large the queue can get up to
        '''
        multiprocessing.Process.__init__(self)
        
        #always get these stats
        self.statistics = {'K':{'1':0,'2':0,'3':0},
                           'U':{'1':0,'2':0,'3':0},
                           'Q':{'1':0,'2':0,'3':0},
                           'P':{'1':0,'2':0,'3':0},
                           'R':{'1':0,'2':0,'3':0},
                           'tot':{'1':0,'2':0,'3':0}}
        self.percentageChoice = [0 for x in range(0, 101)]
        
        #save the inputs from the init
        self.readsQueue = readsQueue
        self.resultsQueue = resultsQueue
        self.isMaster = isMaster
        self.firstFilename = firstFilename
        self.secondFilename = secondFilename
        self.baseOutputFN = outputFilename
        self.outHeader = outHeader
        self.isRandomFilter = isRandomFilter
        self.keepAll = keepAll
        
        if numProcs > 1:
            self.outputFilename = outputFilename+'.tmp'+str(workerID)+'.bam'
            self.isSingleProcess = False
        elif isMaster:
            self.outputFilename = outputFilename+'.tmp0.bam'
            self.isSingleProcess = True
        else:
            logger.error('Created non-master process with only 1 process allowed')
        
        self.pileupTempPrefix = outputFilename+'.tmp'+str(workerID)+'.pileup'
        self.numProcs = numProcs
        self.workerID = workerID
        self.maxQueueSize = max_queue_size
        
        if mergerType == UNION_MERGE:
            self.mergerType = UNION_MERGE
        elif mergerType == PILEUP_MERGE:
            self.mergerType = PILEUP_MERGE
        else:
            self.mergerType = RANDOM_MERGE
        
        self.stage = QUALITY_STAGE
        
        #finally start our own manager
        class LocalSynchronyManager(SyncManager):
            pass
        
        LocalSynchronyManager.register('numPastStage')
        LocalSynchronyManager.register('finishedStage')
        
        self.ltm = LocalSynchronyManager(managerAddress)
        self.ltm.connect()
        self.pileupDict = pileupDict

    '''
    def finishAndBarrierSync(self, stage):
        self.ltm.finishedStage(stage)
        while int(str(self.ltm.numPastStage(stage))) < self.numProcs:
            time.sleep(1)
        
        #once we get here we know everything is at or past this stage
    '''
        
    def run(self):
        '''
        required since it's a process
        '''
        if self.isMaster:
            self.fillMergeQueue()
        else:
            self.emptyMergeQueue()
        
        #return the same thing regardless of what happened earlier
        self.statistics['percentageChoice'] = self.percentageChoice
        self.statistics['outputFilename'] = self.outputFilename
        self.resultsQueue.put(self.statistics)
        
    def fillMergeQueue(self):
        '''
        This is the main function that should be called by an outside entity to perform a merge
        @param firstFilename - the name of the file that holds the first parent's alignments
        @param secondFilename - the name of the file that holds the second parent's alignments
        @param outputFilename - the name of the destination file for saving the merge results
        '''
        random.seed()
        #import pydevd;pydevd.settrace()
        #attempt to open both bam files for merging
        firstFile = pysam.Samfile(self.firstFilename, 'rb')
        secondFile = pysam.Samfile(self.secondFilename, 'rb')
        
        #open the output file as well
        self.outputFile = pysam.Samfile(self.outputFilename, 'wb', header=self.outHeader, referencenames=firstFile.references)
        
        if self.mergerType == PILEUP_MERGE:
            #this one stores pieces counted in the pileup
            self.pileupTempFile = pysam.Samfile(self.pileupTempPrefix+'.bam', 'wb', template=firstFile)
            
            #this one stores those to be resolved by pileup
            self.tempFN = self.outputFilename+'.pileup_tmp.bam'
            self.tempFile = pysam.Samfile(self.tempFN, 'wb', template=firstFile)
        
        #start with no qname
        currentQname = None
        currClusterSize = 0
        clusterStart = None
        clusterEnd = None
        
        #get the first read ID from each file
        firstTell = firstFile.tell()
        firstRead = firstFile.next()
        firstQname = firstRead.qname
        
        secondTell = secondFile.tell()
        secondRead = secondFile.next()
        secondQname = secondRead.qname
        
        clusterStarts = [firstTell, secondTell]
        
        #pick the 'lowest' or first qname
        if isQnameBefore(firstQname, secondQname):
            currentQname = firstQname
        else:
            currentQname = secondQname
        
        firstReads = []
        secondReads = []
        
        #used for output and debugging purposes
        count = 0
        maxReads = None
        isMatch = True
        
        #while we have a read left to process OR if there are matched values from before, this needs to loop at least once more
        while firstRead != None or secondRead != None or isMatch:
            #used for breaking if we need to, strictly for debug purposes
            if maxReads != None and count > maxReads:
                break
            
            #default there is no match
            isMatch = False
            
            #check if the first one has a match, if so increment it and get the next read
            while firstQname == currentQname and firstQname != None:
                #always append the read
                firstReads.append(firstRead)
                
                #get the next one for analysis
                try:
                    firstTell = firstFile.tell()
                    firstRead = firstFile.next()
                    firstQname = firstRead.qname
                except:
                    firstRead = None
                    firstQname = None
                    
                isMatch = True
                
            #check the second also and get the next read if needed
            while secondQname == currentQname and secondQname != None:
                #always append the read
                secondReads.append(secondRead)
                
                #get the next one for analysis
                try:
                    secondTell = secondFile.tell()
                    secondRead = secondFile.next()
                    secondQname = secondRead.qname
                except:
                    secondRead = None
                    secondQname = None
                    
                isMatch = True
                
            #check if there are no reads matching the current qname
            if not isMatch:
                
                #debug statement every 100k counts where a count is a distinct read name
                if count % 100000 == 0:
                    logger.info('[Master] Processed '+str(count)+' read names...')
                    
                #increment the count always when we update the qname
                count += 1
                
                if self.isSingleProcess:
                    #only a single process, us, handle it here
                    self.handleSingleIdMerge(firstReads, secondReads)
                else:
                    #handling multi-process 
                    #if we have no start, make the start this name
                    if clusterStart == None:
                        clusterStart = currentQname
                        
                    #if we are within the boundaries, extend the cluster
                    if currClusterSize < MAX_CLUSTER_SIZE and currentQname != None:
                        clusterEnd = currentQname
                        currClusterSize += 1
                    
                    if currClusterSize == MAX_CLUSTER_SIZE:
                        #the cluster is ready for submission
                        #check if the queue has space
                        if self.readsQueue.qsize() >= self.maxQueueSize.value:
                            #no space in the queue
                            if clusterEnd == currentQname:
                                #we just added this to the very end AND we know the cluster is full, so we shouldn't process it
                                pass
                            else:
                                #cluster is full, BUT this wasn't added to the end
                                self.handleSingleIdMerge(firstReads, secondReads)
                        else:
                            #full cluster and space to submit, so do so
                            #self.readsQueue.put([clusterStart, clusterEnd])
                            self.readsQueue.put(clusterStarts)
                            clusterStarts = [firstTell, secondTell]
                            
                            #there is space in the queue
                            if clusterEnd == currentQname:
                                #we just added this to the end of the cluster AND we know the cluster is full, send it to be processed
                                #reset
                                clusterStart = None
                                clusterEnd = None
                                currClusterSize = 0
                            else:
                                self.handleSingleIdMerge(firstReads, secondReads)
                                
                                clusterStart = None
                                clusterEnd = None
                                currClusterSize = 0
                                
                                '''
                                #we should restart the cluster with this extra read
                                clusterStart = currentQname
                                clusterEnd = currentQname
                                currClusterSize = 1
                                '''
                #reset the pairing
                firstReads = []
                secondReads = []
                
                #get the next qname
                if firstQname != None and (secondQname == None or isQnameBefore(firstQname, secondQname)):
                    currentQname = firstQname
                else:
                    currentQname = secondQname
        
        if clusterStart != None:
            #remainder cluster, put it on 
            #self.readsQueue.put([clusterStart, clusterEnd])
            self.readsQueue.put(clusterStarts)
            
            #reset, just in case
            clusterStart = None
            clusterEnd = None
            currClusterSize = 0
            
        #tell the other processes it's the end of the queue
        for x in range(0, self.numProcs-1):
            self.readsQueue.put(END_OF_QUEUE)
        
        logger.info('[Master] All reads scanned.  Waiting on workers to finish processing...')
        
        #this handles the sorting of this processes pileup file
        if self.mergerType == PILEUP_MERGE:
            self.prepareForPileup()
        
        #close the files
        firstFile.close()
        secondFile.close()
        self.outputFile.close()
        
        logger.info('[0] Completed first pass')
        
        
        
        #self.resultsQueue.put(self.percentageChoice)
        
    def emptyMergeQueue(self):
        '''
        This is executed by non-master worker processes
        '''
        random.seed()
        
        #attempt to open both bam files for merging
        firstFile = pysam.Samfile(self.firstFilename, 'rb')
        secondFile = pysam.Samfile(self.secondFilename, 'rb')
        self.outputFile = pysam.Samfile(self.outputFilename, 'wb', header=self.outHeader, referencenames=firstFile.references)
        
        firstRead = firstFile.next()
        firstQname = firstRead.qname
        
        secondRead = secondFile.next()
        secondQname = secondRead.qname
        
        #open the temp pileup file if necessary
        if self.mergerType == PILEUP_MERGE:
            #this one stores pieces counted in the pileup
            self.pileupTempFile = pysam.Samfile(self.pileupTempPrefix+'.bam', 'wb', template=firstFile)
            
            #this one stores those to be resolved by pileup
            self.tempFN = self.outputFilename+'.pileup_tmp.bam'
            self.tempFile = pysam.Samfile(self.tempFN, 'wb', template=firstFile)
            
        while self.stage == QUALITY_STAGE:
            firstReads = None
            secondReads = None
            
            try:
                qnames = self.readsQueue.get(True, 1)
                
            except:
                self.maxQueueSize.value += 10
                logger.info('['+str(self.workerID)+'] Increasing queue size to '+str(self.maxQueueSize.value))
                qnames = None
            
            if qnames == None:
                pass
            elif qnames == END_OF_QUEUE:
                self.stage = PILEUP_STAGE
            else:
                
                #move to the relevant points in our files
                firstFile.seek(qnames[0])
                secondFile.seek(qnames[1])
                
                try:
                    firstRead = firstFile.next()
                    firstQname = firstRead.qname
                except:
                    firstRead = None
                    firstQname = None
                    
                try:
                    secondRead = secondFile.next()
                    secondQname = secondRead.qname
                except:
                    secondRead = None
                    secondQname = None
                
                if firstRead == None:
                    currentQname = secondQname
                elif secondRead == None:
                    currentQname = firstQname
                elif isQnameBefore(firstQname, secondQname):
                    currentQname = firstQname
                else:
                    currentQname = secondQname
                
                numProcessed = 0
                while currentQname != None and numProcessed < MAX_CLUSTER_SIZE:
                    #each round through this loop is a new set of first and second reads
                    firstReads = []
                    secondReads = []
                
                    #get and save all reads equal to this
                    while firstRead != None and firstQname == currentQname:
                        firstReads.append(firstRead)
                        try:
                            firstRead = firstFile.next()
                            firstQname = firstRead.qname
                        except:
                            firstRead = None
                            firstQname = None
                    
                    
                    #get and save all reads equal to this
                    while secondRead != None and secondQname == currentQname:
                        secondReads.append(secondRead)
                        try:
                            secondRead = secondFile.next()
                            secondQname = secondRead.qname
                        except:
                            secondRead = None
                            secondQname = None
                
                    #merge the two sets
                    self.handleSingleIdMerge(firstReads, secondReads)
                    numProcessed += 1
                    
                    #get the next qname
                    if firstRead == None:
                        currentQname = secondQname
                    elif secondRead == None:
                        currentQname = firstQname
                    elif isQnameBefore(firstQname, secondQname):
                        currentQname = firstQname
                    else:
                        currentQname = secondQname
        
        
        firstFile.close()
        secondFile.close()
        
        #this handles the sorting and indexing of pileup files
        if self.mergerType == PILEUP_MERGE:
            self.prepareForPileup()
        
        firstFile.close()
        secondFile.close()
        self.outputFile.close()
        
        logger.info('['+str(self.workerID)+'] Completed first pass')
        
        #self.resultsQueue.put(self.statistics)
        #self.resultsQueue.put(self.percentageChoice)
        
    def handleSingleIdMerge(self, firstReads, secondReads):
        '''
        This function takes two sets of reads and attempts to pick one of the best alignments to save using region based likelihood
        @param firstReads - the set of reads coming from the first parent
        @param secondReads - the set of reads coming from the second parent
        '''
        if verbosity:
            dumpReads(firstReads, secondReads)
        
        #import pydevd;pydevd.settrace()
        
        #pair the reads using their HI tags
        [firstPairs, firstSingles] = pairReads(firstReads, 'HI')
        [secondPairs, secondSingles] = pairReads(secondReads, 'HI')
        
        '''
        #remove this later
        if (len(firstPairs) > 0 and len(firstSingles) > 0) or (len(secondPairs) > 0 and len(secondSingles) > 0):
            dumpReads(firstReads, secondReads)
        '''
        
        #combine reads but keep all of them regardless of score
        possiblePairs = combinePairs(firstPairs, secondPairs)
        possibleSingles = combineSingles(firstSingles, secondSingles)
        
        #stats re-done
        pairLen = len(possiblePairs)
        if pairLen == 0:
            #this is when the two ends may be unique but there isn't a unique pair
            for seq in possibleSingles:
                #only true if A) there are no pairs and B) there is only one alignment for this sequence
                if len(possibleSingles[seq]) == 1:
                    #single end, mark union and save it
                    setTag(possibleSingles[seq][0], CHOICE_TYPE_TAG, 'U')
                    self.saveRead(possibleSingles[seq][0])
                    
                    if self.mergerType == PILEUP_MERGE:
                        self.parseReadForTree(possibleSingles[seq][0])
                    
                    #clear the possible singles here so we cannot save it twice by accident
                    possibleSingles[seq] = []
                    
        elif pairLen == 1 and len(possibleSingles) == 0:
            #single pair, no singles
            setTag(possiblePairs[0][0], CHOICE_TYPE_TAG, 'U')
            setTag(possiblePairs[0][1], CHOICE_TYPE_TAG, 'U')
            
            #this can be completely solved using the union
            self.saveRead(possiblePairs[0][0])
            self.saveRead(possiblePairs[0][1])
            
            if self.mergerType == PILEUP_MERGE:
                self.parseReadForTree(possiblePairs[0][0])
                self.parseReadForTree(possiblePairs[0][1])
            
            #clear this so we cannot save it twice by accident
            possiblePairs = []
            
        else:
            #nothing unique from straight union
            pass
            
        if self.mergerType == UNION_MERGE:
            #save all pairs
            for pair in possiblePairs:
                setTag(pair[0], CHOICE_TYPE_TAG, 'K')
                setTag(pair[1], CHOICE_TYPE_TAG, 'K')
                
                self.saveRead(pair[0])
                self.saveRead(pair[1])
                
            #save all singles
            for seq in possibleSingles:
                for single in possibleSingles[seq]:
                    setTag(single, CHOICE_TYPE_TAG, 'K')
                    self.saveRead(single)
                    
        else:
            #we now need to reduce the remaining reads based on quality
            #the following functions weed out lower quality reads
            possiblePairs = reducePairsByScore(possiblePairs)
            possibleSingles = reduceSinglesByScore(possibleSingles)
            
            if len(possiblePairs) > 0:
                #if there's only one possible pair remaining, then it's a quality based match
                if len(possiblePairs) == 1:
                    setTag(possiblePairs[0][0], CHOICE_TYPE_TAG, 'Q')
                    setTag(possiblePairs[0][1], CHOICE_TYPE_TAG, 'Q')
                    
                    self.saveRead(possiblePairs[0][0])
                    self.saveRead(possiblePairs[0][1])
                    
                    if self.mergerType == PILEUP_MERGE:
                        self.parseReadForTree(possiblePairs[0][0])
                        self.parseReadForTree(possiblePairs[0][1])
                        
                else:
                    if self.mergerType == RANDOM_MERGE:
                        #we randomly choose a pair to keep if this is the case
                        self.saveRandomPair(possiblePairs)
                    
                    elif self.mergerType == PILEUP_MERGE:
                        #store all of the pairs remaining for sorting through later
                        pHI = 0
                        for pair in possiblePairs:
                            setTag(pair[0], PILEUP_HI_TAG, pHI)
                            setTag(pair[1], PILEUP_HI_TAG, pHI)
                            self.tempFile.write(pair[0])
                            self.tempFile.write(pair[1])
                            pHI += 1
            else:
                #we randomly choose for each sequence a single alignment to keep
                pHI = 0
                for seq in possibleSingles:
                    if len(possibleSingles[seq]) == 1:
                        #only one option for this end of the read, save it
                        setTag(possibleSingles[seq][0], CHOICE_TYPE_TAG, 'Q')
                        self.saveRead(possibleSingles[seq][0])
                        
                        #parse the read if it's a pileup merge
                        if self.mergerType == PILEUP_MERGE:
                            self.parseReadForTree(possibleSingles[seq][0])
                        
                    elif len(possibleSingles[seq]) > 1:
                        #multiple options for this end of the read
                        if self.mergerType == RANDOM_MERGE:
                            #random merge, random save
                            self.saveRandomSingle(possibleSingles[seq])
                        elif self.mergerType == PILEUP_MERGE:
                            #store these alignments for pileup checking later
                            for single in possibleSingles[seq]:
                                setTag(single, PILEUP_HI_TAG, pHI)
                                self.tempFile.write(single)
                                pHI += 1
    
    def prepareForPileup(self):
        '''
        Executed when done processing the first pass of reads from quality scores, basically sort our data so we can merge and 
        build the pileup heights
        '''
        #close the file that has the reads for the pileup calculations
        self.tempFile.close()
        self.pileupTempFile.close()
        
        logger.info('['+str(self.workerID)+'] Sorting the pileup data...')
        
        #sort the current output file into normal positional sort
        pysam.sort(self.pileupTempPrefix+'.bam', self.pileupTempPrefix+'.sorted')
        
    def cleanUpFiles(self):
        '''
        this function removes the temporary files from doing the pileup calculations
        '''
        if self.mergerType == PILEUP_MERGE and not self.keepAll:
            #remove all the extra pileup files
            os.remove(self.baseOutputFN+'.tmp'+str(self.workerID)+'.pileup.bam')
            os.remove(self.baseOutputFN+'.tmp'+str(self.workerID)+'.pileup.sorted.bam')
            #os.remove(self.baseOutputFN+'.tmp'+str(self.workerID)+'.pileup.sorted.bam.bai')
    
            if self.isMaster:
                os.remove(self.baseOutputFN+'.tmp.pileup_all.bam')
                os.remove(self.baseOutputFN+'.tmp.pileup_all.bam.bai')
                
    def saveRead(self, readToSave):
        '''
        @param readToSave - the read that should be added to the output
        @param parent - 1, 2, 3, or SET; 1, 2, and 3 indicate the parent type, SET means it's already added, don't add again
        @param out - whether you want to see output from this save or not
        '''
        if readToSave == None:
            return
        
        #save the read
        self.outputFile.write(readToSave)
        
        choice = getTag(readToSave, CHOICE_TYPE_TAG)
        parent = getTag(readToSave, PARENT_OF_ORIGIN_TAG)
        
        if parent == None or choice == None:
            logger.error('Missing parent and/or choice tags:'+str(readToSave))
        else:
            if self.statistics.has_key(choice) and self.statistics[choice].has_key(parent):
                self.statistics[choice][parent] += 1
                self.statistics['tot'][parent] += 1
            else:
                print 'Poorly structured tag:'
                print readToSave
            
        if verbosity:
            logger.info('Saving from '+getTag(readToSave, PARENT_OF_ORIGIN_TAG)+': '+str(readToSave))
        
    def saveRandomPair(self, possiblePairs):
        '''
        @param possiblePairs - a list of possible pairs that could be saved
        @param parent - the parent of origin for these pairs, can be SET indicating a mix of pairs that already have it set
        '''
        if len(possiblePairs) > 1:
            if self.isRandomFilter:
                #pick a random pair
                rv = random.randint(0, len(possiblePairs)-1)
            
                #random choice
                setTag(possiblePairs[rv][0], CHOICE_TYPE_TAG, 'R')
                setTag(possiblePairs[rv][1], CHOICE_TYPE_TAG, 'R')
                
                #save the pair
                self.saveRead(possiblePairs[rv][0])
                self.saveRead(possiblePairs[rv][1])
            else:
                for pair in possiblePairs:
                    setTag(pair[0], CHOICE_TYPE_TAG, 'K')
                    setTag(pair[1], CHOICE_TYPE_TAG, 'K')
                    
                    self.saveRead(pair[0])
                    self.saveRead(pair[1])
        else:
            #save the pair
            self.saveRead(possiblePairs[0][0])
            self.saveRead(possiblePairs[0][1])
        
        
    def saveRandomSingle(self, possibleSingles):
        '''
        @param possibleSingles - a list of possible singles that could be saved
        @param parent - the parent of origin for these singles, can be SET indicating a mix of pairs that already have it set
        '''
        if len(possibleSingles) > 1:
            if self.isRandomFilter:
                #pick a random single and save it
                rv = random.randint(0, len(possibleSingles)-1)
                setTag(possibleSingles[rv], CHOICE_TYPE_TAG, 'R')
                self.saveRead(possibleSingles[rv])
            else:
                for single in possibleSingles:
                    setTag(single, CHOICE_TYPE_TAG, 'K')
                    self.saveRead(single)
        else:
            self.saveRead(possibleSingles[0])
    
    def parseReadForTree(self, readToParse):
        '''
        @param readToParse - the read we want to save as part of the pileup calculations later
        '''
        if readToParse == None:
            return
        
        #save the read
        self.pileupTempFile.write(readToParse)
        
def saveChoiceChart(data, fn, titleFn):
    #import pydevd;pydevd.settrace()
    logger.info('Generating pileup dominance percentage chart')
    pylab.figure(1)
    xaxis = [x for x in range(0, len(data))]
    pylab.bar(xaxis, data)
    pylab.title(titleFn)
    pylab.xlabel('Pileup dominance (%)')
    pylab.xlim(xmin=0, xmax=101)
    pylab.ylabel('Number of read ends')
    pylab.grid(True)
    pylab.savefig(fn)
    logger.info('Chart data: '+str(data))
    logger.info('Chart saved to "'+fn+'".')

def calculatePairScore(readPair):
    '''
    @param readPair - the paired end read that needs a score calculated
    @return - the score associated with the paired end read
    '''
    #get the tags from the pair
    oc1 = getTag(readPair[0], OLD_CIGAR_TAG)
    oc2 = getTag(readPair[1], OLD_CIGAR_TAG)
    
    ed1 = getTag(readPair[0], 'OM')
    ed2 = getTag(readPair[1], 'OM')

    #calculate and return the score
    score = calculateScore(oc1, ed1) + calculateScore(oc2, ed2)
    return score

def calculateReadScore(read):
    '''
    @param read - the unpaired read that needs a score calculated
    @return - the score for this unpaired read
    '''
    #get tags
    oc = getTag(read, OLD_CIGAR_TAG)
    ed = getTag(read, 'OM')

    #calculate the score
    score = calculateScore(oc, ed)
    return score
    
def calculateScore(cigar, editDistance):
    '''
    @param cigar - the cigar string from the read
    @param editDistance - the 'OM' tag from the read, basically the edit distance (SNPs and indel count)
    @return - the score for this read as calculated by Bowtie
    '''
    
    #make sure we have both values
    if cigar == None or editDistance == None:
        return 0
    
    #scoring constants as gathered from the bowtie website
    #TODO: make these user-configurable
    #MATCH = 2
    MATCH = 0
    MISMATCH = -6
    GAP_OPEN = -5
    EXTENSION = -3
    
    #init
    number = ''
    cigType = ''
    
    score = 0
    foundEdits = 0
    
    #parse through the cigar string
    for symbol in cigar:
        
        #if it's a digit, it's part of the number
        if symbol.isdigit():
            number += symbol
            
        #we're looking at an identifier for scoring
        else:
            #get the symbol and parse the number
            cigType = symbol
            intNum = int(number)
            
            #match type
            if cigType == 'M':
                score += MATCH*intNum
            
            #indel type
            elif cigType == 'I' or cigType == 'D':
                score += GAP_OPEN + intNum*EXTENSION
                foundEdits += intNum
            
            #gap type
            elif cigType == 'N':
                #ignore these, they don't effect scoring
                pass
            
            #unhandled type
            else:
                logger.warning('Unhandled cigar string type:'+symbol)
            
            #reset these values
            number = ''
            cigType = ''
    
    #remove indels and all that's left is mismatches, so factor that in
    #NOTE: we already counted this as a MATCH, so we need to remove that AND add the mismatch penalty
    score += (editDistance-foundEdits)*(MISMATCH-MATCH)
    
    #return the final score
    return score


def getTag(read, tag):
    '''
    @param read - the read we want to get the tag from
    @param tag - the tag type, 'OM', 'OC', etc.
    @return - the tag from read or None if that tag isn't found
    '''
    #make sure we have a read
    if read == None:
        return None
    
    #get all the tags
    allTags = read.tags
    
    #search for our tag and return the associated value
    for t in allTags:
        if t[0] == tag:
            return t[1]
        
    #wasn't found, return nothing
    return None

def setTag(read, tag, value):
    '''
    @param read - the read to modify
    @param tag - the tag to add
    @param value - the value of the tag
    '''
    read.tags = read.tags + [(tag, value)]
    
def dumpReads(firstReads, secondReads):
    '''
    This function is mainly for debugging
    @param firstReads - the first set of reads
    @param secondReads - the second set of reads
    '''
    #old debug junk
    print 'First:'+str(len(firstReads))
    for read in firstReads:
        print read
    
    print 'Second:'+str(len(secondReads))
    for read in secondReads:
        print read
        
    print
    
def pairReads(reads, pairTag):
    '''
    @param reads - the reads we want to pair up if able
    @return - [the list of paired reads after processing, the list of unpaired reads after processing]
    '''
    #init
    pairs = []
    singles = []
    readsLookup = {}
    
    #TODO: if THIS segment is unmapped, we should remove it from any return value
    
    #iterate through the reads
    for read in reads:
        if isFlagSet(read.flag, MULTIPLE_SEGMENT_FLAG):
            #if isFlagSet(read.flag, BOTH_ALIGNED_FLAG):
            if not isFlagSet(read.flag, OTHER_UNMAPPED_FLAG):
                #this is a read that is paired and both ends aligned
                #get the hiTag
                hiTag = getTag(read, pairTag)
                
                #check if we have first side of this pair
                if readsLookup.has_key(hiTag):
                    #other side of pair exists, pair them up
                    otherRead = readsLookup[hiTag]
                    if isFlagSet(read.flag, FIRST_SEGMENT_FLAG):
                        singlePair = [read, otherRead]
                    else:
                        singlePair = [otherRead, read]
                    
                    #add the pair
                    pairs.append(singlePair)
                    
                    #remove this
                    readsLookup.pop(hiTag)
                    
                else:
                    #other side of pair isn't found yet, so save this end
                    readsLookup[hiTag] = read
                    
                    if hiTag == None and int(getTag(read, 'NH')) > 1:
                        logger.error('NH > 1 but no \''+pairTag+'\' tag: '+str(read))
            else:
                #only one aligned, so it's not paired
                singles.append(read)
        else:
            #there are not multiple segments, so it's not paired
            singles.append(read)
    
    if len(readsLookup) != 0:
        print pairTag
        print readsLookup
        
        dumpReads(reads, [])
        
        for tagValue in readsLookup:
            singles.append(readsLookup[tagValue])
    
    #return both lists
    return [pairs, singles]

def combinePairs(pairs1, pairs2):
    '''
    @param pairs1 - the first set of pairs
    @param pairs2 - the second set of pairs
    @return returns a merge set of pairs with the appropriate PO and FC tags already set
    '''
    #unique pairs for potential saving
    uniquePairs = []
    
    #iterate through trying to match up pairs
    for p1 in pairs1:
        foundSame = False
        for p2 in pairs2:
            #get the original p2 and store it
            origP2 = p2
            
            #check if this pair is is the same position
            if isPositionSame(p1[0], p2[0]) and isPositionSame(p1[1], p2[1]):
                foundSame = True
            elif isPositionSame(p1[0], p2[1]) and isPositionSame(p1[1], p2[0]):
                foundSame = True
                
                #for this case, we need to swap them so they are in the same order when cigar strings are compared
                p2 = [p2[1], p2[0]]
            else:
                pass
            
            #break out if we get a match
            if foundSame:
                break
        
        if foundSame:
            #check scores and keep the best ones
            score1 = calculatePairScore(p1)
            score2 = calculatePairScore(p2)
            
            #same position, therefor YA is 3
            setTag(p1[0], YA_TAG, '3')
            setTag(p1[1], YA_TAG, '3')
            setTag(p2[0], YA_TAG, '3')
            setTag(p2[1], YA_TAG, '3')
            
            #get the first pair's first end s0 tag
            ms1 = getTag(p1[0], SNP_TAG)
            mi1 = getTag(p1[0], INSERTION_TAG)
            md1 = getTag(p1[0], DELETION_TAG)
            setTag(p1[0], M_SNP_TAG, ms1)
            setTag(p2[0], M_SNP_TAG, ms1)
            setTag(p1[0], M_IN_TAG, mi1)
            setTag(p2[0], M_IN_TAG, mi1)
            setTag(p1[0], M_DEL_TAG, md1)
            setTag(p2[0], M_DEL_TAG, md1)
            
            #get the first pair's second end s0 tag
            ms2 = getTag(p1[1], SNP_TAG)
            mi2 = getTag(p1[1], INSERTION_TAG)
            md2 = getTag(p1[1], DELETION_TAG)
            setTag(p1[1], M_SNP_TAG, ms2)
            setTag(p2[1], M_SNP_TAG, ms2)
            setTag(p1[1], M_IN_TAG, mi2)
            setTag(p2[1], M_IN_TAG, mi2)
            setTag(p1[1], M_DEL_TAG, md2)
            setTag(p2[1], M_DEL_TAG, md2)
            
            #get the second pair's first end s0 tag
            ps1 = getTag(p2[0], SNP_TAG)
            pi1 = getTag(p2[0], INSERTION_TAG)
            pd1 = getTag(p2[0], DELETION_TAG)
            setTag(p1[0], P_SNP_TAG, ps1)
            setTag(p2[0], P_SNP_TAG, ps1)
            setTag(p1[0], P_IN_TAG, pi1)
            setTag(p2[0], P_IN_TAG, pi1)
            setTag(p1[0], P_DEL_TAG, pd1)
            setTag(p2[0], P_DEL_TAG, pd1)
            
            #get the second pair's second end s0 tag
            ps2 = getTag(p2[1], SNP_TAG)
            pi2 = getTag(p2[1], INSERTION_TAG)
            pd2 = getTag(p2[1], DELETION_TAG)
            setTag(p1[1], P_SNP_TAG, ps2)
            setTag(p2[1], P_SNP_TAG, ps2)
            setTag(p1[1], P_IN_TAG, pi2)
            setTag(p2[1], P_IN_TAG, pi2)
            setTag(p1[1], P_DEL_TAG, pd2)
            setTag(p2[1], P_DEL_TAG, pd2)
            
            #if cigars are different, save the second
            oc1 = getTag(p1[0], OLD_CIGAR_TAG)
            oc2 = getTag(p2[0], OLD_CIGAR_TAG)
            om1 = getTag(p1[0], OLD_MISMATCH_TAG)
            om2 = getTag(p2[0], OLD_MISMATCH_TAG)
            if oc1 != oc2 or om1 != om2:
                setTag(p1[0], P_CIGAR_TAG, oc2)
                setTag(p1[0], P_MISMATCH_TAG, om2)
                setTag(p2[0], M_CIGAR_TAG, oc1)
                setTag(p2[0], M_MISMATCH_TAG, om1)
                
            oc1 = getTag(p1[1], OLD_CIGAR_TAG)
            oc2 = getTag(p2[1], OLD_CIGAR_TAG)
            om1 = getTag(p1[1], OLD_MISMATCH_TAG)
            om2 = getTag(p2[1], OLD_MISMATCH_TAG)
            if oc1 != oc2 or om1 != om2:
                setTag(p1[1], P_CIGAR_TAG, oc2)
                setTag(p1[1], P_MISMATCH_TAG, om2)
                setTag(p2[1], M_CIGAR_TAG, oc1)
                setTag(p2[1], M_MISMATCH_TAG, om1)
            
            if score1 == score2:
                #set the PO tag
                setTag(p1[0], PARENT_OF_ORIGIN_TAG, '3')
                setTag(p1[1], PARENT_OF_ORIGIN_TAG, '3')
                
                #add the modified pair to our list
                uniquePairs.append(p1)
                
            else:
                #differing score, keep both as different
                setTag(p1[0], PARENT_OF_ORIGIN_TAG, '1')
                setTag(p1[1], PARENT_OF_ORIGIN_TAG, '1')
                uniquePairs.append(p1)
                
                setTag(p2[0], PARENT_OF_ORIGIN_TAG, '2')
                setTag(p2[1], PARENT_OF_ORIGIN_TAG, '2')
                uniquePairs.append(p2)
            
            #remove this from pairs2 so we don't mark it as unmatched down there
            pairs2.remove(origP2)
            
        else:
            #no match, mark it as such
            #set the PO tag
            setTag(p1[0], PARENT_OF_ORIGIN_TAG, '1')
            setTag(p1[1], PARENT_OF_ORIGIN_TAG, '1')
            setTag(p1[0], YA_TAG, '1')
            setTag(p1[1], YA_TAG, '1')
            uniquePairs.append(p1)
            
    for p2 in pairs2:
        #add anything leftover here
        #set the PO tag
        setTag(p2[0], PARENT_OF_ORIGIN_TAG, '2')
        setTag(p2[1], PARENT_OF_ORIGIN_TAG, '2')
        setTag(p2[0], YA_TAG, '2')
        setTag(p2[1], YA_TAG, '2')
        uniquePairs.append(p2)
    
    #at this point uniquePairs contains each unique pair and the FC has been set if necessary
    return uniquePairs
    
def combineSingles(singles1, singles2):
    '''
    @param singles1 - the list of unpaired alignments from the first parent
    @param singles2 - the list of unpaired alignments from the second parent
    @return returns a map with two lists of possible alignments for each endpoint
    '''
    #init scoring
    uniqueSingles = {}
    
    #now we loop through trying to pair up the singles
    for s1 in singles1:
        foundSame = False
        
        #iterate through each second single, comparing the start position
        for s2 in singles2:
            #check if this pair is is the same position
            if isPositionSame(s1, s2):
                foundSame = True
                break
        
        if foundSame:
            #check scores and keep the best ones
            score1 = calculateReadScore(s1)
            score2 = calculateReadScore(s2)
            
            #YA tag is 3
            setTag(s1, YA_TAG, '3')
            setTag(s2, YA_TAG, '3')
            
            ms = getTag(s1, SNP_TAG)
            mi = getTag(s1, INSERTION_TAG)
            md = getTag(s1, DELETION_TAG)
            setTag(s1, M_SNP_TAG, ms)
            setTag(s2, M_SNP_TAG, ms)
            setTag(s1, M_IN_TAG, mi)
            setTag(s2, M_IN_TAG, mi)
            setTag(s1, M_DEL_TAG, md)
            setTag(s2, M_DEL_TAG, md)
            
            ps = getTag(s2, SNP_TAG)
            pi = getTag(s2, INSERTION_TAG)
            pd = getTag(s2, DELETION_TAG)
            setTag(s1, P_SNP_TAG, ps)
            setTag(s2, P_SNP_TAG, ps)
            setTag(s1, P_IN_TAG, pi)
            setTag(s2, P_IN_TAG, pi)
            setTag(s1, P_DEL_TAG, pd)
            setTag(s2, P_DEL_TAG, pd)
            
            #if cigars are different, save the second
            oc1 = getTag(s1, OLD_CIGAR_TAG)
            oc2 = getTag(s2, OLD_CIGAR_TAG)
            om1 = getTag(s1, OLD_MISMATCH_TAG)
            om2 = getTag(s2, OLD_MISMATCH_TAG)
            if oc1 != oc2 or om1 != om2:
                setTag(s1, P_CIGAR_TAG, oc2)
                setTag(s1, P_MISMATCH_TAG, om2)
                setTag(s2, M_CIGAR_TAG, oc1)
                setTag(s2, M_MISMATCH_TAG, om1)
                
            #unique singles has two empty arrays from the sequence keys
            if score1 == score2:
                #set the PO tag
                setTag(s1, PARENT_OF_ORIGIN_TAG, '3')
                isFirstSegment = isFlagSet(s1.flag, FIRST_SEGMENT_FLAG)
                
                if not uniqueSingles.has_key(isFirstSegment):
                    uniqueSingles[isFirstSegment] = []
                uniqueSingles[isFirstSegment].append(s1)
                
            else:
                isFirstSegment = isFlagSet(s1.flag, FIRST_SEGMENT_FLAG)
                if not uniqueSingles.has_key(isFirstSegment):
                    uniqueSingles[isFirstSegment] = []
                
                #different so keep both around
                setTag(s1, PARENT_OF_ORIGIN_TAG, '1')
                uniqueSingles[isFirstSegment].append(s1)
                
                setTag(s2, PARENT_OF_ORIGIN_TAG, '2')
                uniqueSingles[isFirstSegment].append(s2)
            
            singles2.remove(s2)
            
        else:
            isFirstSegment = isFlagSet(s1.flag, FIRST_SEGMENT_FLAG)
            if not uniqueSingles.has_key(isFirstSegment):
                uniqueSingles[isFirstSegment] = []
            
            #no match, mark it as such
            setTag(s1, PARENT_OF_ORIGIN_TAG, '1')
            setTag(s1, YA_TAG, '1')
            uniqueSingles[isFirstSegment].append(s1)
            
            
    for s2 in singles2:
        isFirstSegment = isFlagSet(s2.flag, FIRST_SEGMENT_FLAG)
        if not uniqueSingles.has_key(isFirstSegment):
            uniqueSingles[isFirstSegment] = []

        #add anything leftover here
        setTag(s2, PARENT_OF_ORIGIN_TAG, '2')
        setTag(s2, YA_TAG, '2')
        uniqueSingles[isFirstSegment].append(s2)
    
    return uniqueSingles
    
def reducePairsByScore(pairs):
    '''
    Iterate through a list of paired-end reads and only return those with the highest scores
    @param pairs - the pairs to scan
    '''
    #get the best overall paired score for anything that may be leftover
    bestScore = None
    bestPairs = []
    for pair in pairs:
        score = calculatePairScore(pair)
        if bestScore == None or score > bestScore:
            bestScore = score
            bestPairs = []
            
        if score == bestScore:
            bestPairs.append(pair)
            
    #after the above loops, we have the best score between all detected pairs in either alignment
    return bestPairs

def reduceSinglesByScore(singles):
    '''
    Iterate through a list of single-end reads and only return those with the highest scores
    @param singles - the list of single-end reads to scan
    '''
    bestSingles = {}
    for seq in singles:
        bestScore = None
        bests = []
        
        for single in singles[seq]:
            score = calculateReadScore(single)
            if bestScore == None or score > bestScore:
                bestScore = score
                bests = []
            
            if score == bestScore:
                bests.append(single)
                
        bestSingles[seq] = bests
    
    return bestSingles
    
def isFlagSet(value, FLAG):
    '''
    @param value - the flag value to be checked
    @param FLAG - a constant representing which flag we want to check (see top of file for some flag constants)
    @return - True if 'value' has the bit for FLAG set
    '''
    
    if (value & FLAG) == 0:
        return False
    else:
        return True

def isQnameBefore(qname1, qname2):
    '''
    @param qname1 - the first qname
    @param qname2 - the second qname
    @return - True if qname1 comes before qname2
    '''
    ret = None
    
    if qname2 == None:
        ret = False
    
    if qname1 == None:
        ret = True
    
    pos1 = 0
    pos2 = 0
    
    parts1 = re.split('([0-9]+)', qname1)
    parts2 = re.split('([0-9]+)', qname2)
    
    isNumerical = False

    while ret == None and pos1 < len(parts1) and pos2 < len(parts2):
        if isNumerical:
            num1 = int(parts1[pos1])
            num2 = int(parts2[pos2])

            #compare those numbers
            if num1 < num2:
                ret = True
            elif num1 > num2:
                ret = False
            else:
                #identical
                pass
        else:
            l1 = len(parts1[pos1])
            l2 = len(parts2[pos2])

            comp1 = parts1[pos1]
            comp2 = parts2[pos2]
        
            if l1 == l2:
                #looking at normal letters/symbols
                pass
            elif l1 < l2:
                if pos1+1 < l1:
                    comp1 = comp1 + parts1[pos1+1][0]
            else:
                if pos2+1 < l2:
                    comp2 = comp2 + parts2[pos2+1][0]

            if comp1 < comp2:
                ret = True
            elif comp1 > comp2:
                ret = False
            else:
                #identical
                pass

        pos1 += 1
        pos2 += 1
        isNumerical = not isNumerical
    
    if ret == None:
        if len(parts1) < len(parts2):
            ret = True
        else:
            #identical or greater
            ret = False
    
    #print ret
    return ret

def isPositionSame(read1, read2):
    '''
    @param read1 - the first read to check
    @param read2 - the second read to check
    @return - True if they have the same starting position and cigar string OR if they are both None; False otherwise
    '''
    
    #if both none, True
    if read1 == None and read2 == None:
        return True
    
    #if only one is None, False
    if read1 == None or read2 == None:
        return False
    
    #finally, compare
    #if read1.seq == read2.seq and read1.rname == read2.rname and read1.pos == read2.pos:
    #same segment, same chrom, same pos
    if isFlagSet(read1.flag, FIRST_SEGMENT_FLAG) == isFlagSet(read2.flag, FIRST_SEGMENT_FLAG) and read1.rname == read2.rname and read1.pos == read2.pos:
        #now check the cigars
        cig1 = read1.cigar
        cig2 = read2.cigar
        
        #go through each position in the cigar
        for cigLoc in range(0, len(cig1)):
            if cig1[cigLoc][0] == cig2[cigLoc][0] and cig1[cigLoc][1] == cig2[cigLoc][1]:
                #same type and same position, we're still fine
                pass
            else:
                #different type or position, return false
                return False
        
        #all values matched, good to go
        return True
    else:
        return False

def initLogger():
    '''
    This code taken from Shunping's Lapels for initializing a logger
    '''
    global logger
    logger = logging.getLogger('root')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
def getLogger():
    return logging.getLogger('root')

def performInputChecks(fn1, fn2):
    bam1 = pysam.Samfile(fn1, 'rb')
    bam2 = pysam.Samfile(fn2, 'rb')
    
    retHeader = dict(bam1.header.items())
    h2 = dict(bam2.header.items())
    
    if retHeader['HD']['SO'] != 'query_name':
        logger.error('File "'+fn1+'" is not sorted by name.')
        return None
    
    if h2['HD']['SO'] != 'query_name':
        logger.error('File "'+fn2+'" is not sorted by name.')
        return None
    
    #first combine the header if h2 has a header to combine with
    if h2.has_key('PG'):
        #set h2 to different keys before combining
        for pg in h2['PG']:
            pg['ID'] = pg['ID'] + ' #2'
            if pg.has_key('PP'):
                pg['PP'] = pg['PP'] + ' #2'
                
        if retHeader.has_key('PG'):
            #finally combine them
            retHeader['PG'] = h2['PG'] + retHeader['PG']
            
        else:
            #no data was in the original PG, so copy the one from h2
            retHeader['PG'] = h2['PG']
    
        
    #now try adding my header so final result should be suspenders, h2, h1
    try:
        retHeader['PG'] = [{'ID': 'Suspenders', 'VN': PKG_VERSION,
                            'PP': retHeader['PG'][0]['ID'],
                            'CL': ' '.join(sys.argv)}] + retHeader['PG']
    except KeyError:
        retHeader['PG'] = [{'ID': 'Suspenders', 'VN': PKG_VERSION,
                            'CL': ' '.join(sys.argv)}]
    bam1.close()
    bam2.close()
    
    return retHeader

def mainRun():
    #start up the logger
    initLogger()
    
    #attempt to parse the arguments
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #version data
    p.add_argument('-V', '--version', action='version', version='%(prog)s' + \
                   ' %s in Suspenders %s' % (VERSION, PKG_VERSION))
    
    group = p.add_mutually_exclusive_group()
    group.add_argument('-u', '--union', dest='mergeType', action='store_const', const=UNION_MERGE, help='merge the files using the union filter', default=RANDOM_MERGE)
    group.add_argument('-s', '--quality', dest='mergeType', action='store_const', const=RANDOM_MERGE, help='merge the files using the quality filter (default)', default=RANDOM_MERGE)
    group.add_argument('-t', '--pileup', dest='mergeType', action='store_const', const=PILEUP_MERGE, help='merge the files using the quality filter, then the pileup height filter', default=RANDOM_MERGE)
    
    #optional arguments
    group2 = p.add_mutually_exclusive_group()
    group2.add_argument('-k', '--union-filter', dest='isRandomFilter', action='store_false', help='use union filter if other filters fail (always active for --union)', default=True)
    group2.add_argument('-r', '--random-filter', dest='isRandomFilter', action='store_true', help='use random filter if other filters fail (default for --quality and --pileup)', default=True)
    
    #number of processes is next in importance
    p.add_argument('-p', metavar='numProcesses', dest='numProcesses', type=int, default=1, help='number of processes to run (default: 1)')
    
    #superfluous/debugging arguments
    p.add_argument('-c', metavar='chartFilename', dest='chartFilename', type=str, help='save pileup chart to an image (default: none)', default=None)
    #p.add_argument('-v', dest='verbose', action='store_true', help='verbose (default: no)')
    
    p.add_argument('-e', '--keep-temp', dest='keepTemp', action='store_true', help='keep the temporary merge files', default=False)
    
    #required main arguments
    p.add_argument('motherSortedBam', type=util.readableFile, help='the bam file of the sorted by name alignment to the mother pseudogenome')
    p.add_argument('fatherSortedBam', type=util.readableFile, help='the bam file of the sorted by name alignment to the father pseudogenome')
    p.add_argument('outMergedBam', type=util.writableFile, help='the output bam file merged and sorted by name')
    
    #parse the arguments
    args = p.parse_args()
    
    #TODO: check if we should be wordy
    #if args.verbose:
    #    verbosity = True
    
    if args.mergeType == UNION_MERGE:
        mt = 'Union'
        ft = 'Unique->Kept-All'
    elif args.mergeType == RANDOM_MERGE:
        mt = 'Quality'
        ft = 'Unique->Quality->'
    elif args.mergeType == PILEUP_MERGE:
        mt = 'Pileup'
        ft = 'Unique->Quality->Pileup->'
    else:
        mt = 'ERROR: Unknown'
        ft = 'ERROR: Unknown'
        
    if args.mergeType != UNION_MERGE:
        if args.isRandomFilter:
            ft += 'Random'
        else:
            ft += 'Kept-All'

    logger.info('Merge Type:'+mt)
    logger.info('Filter Type: '+ft)
    logger.info('Mother: '+args.motherSortedBam)
    logger.info('Father: '+args.fatherSortedBam)
    logger.info('Output: '+args.outMergedBam)
    logger.info('Number of processes: '+str(args.numProcesses))
    if args.keepTemp:
        logger.info('Keep temporary files: True')
    
    outHeader = performInputChecks(args.motherSortedBam, args.fatherSortedBam)
    if outHeader == None:
        #there was an error reported, end program
        return
    
    readsQueue = multiprocessing.Queue()
    resultsQueue = multiprocessing.Queue()
    
    max_queue_size = multiprocessing.Value('i', 100*args.numProcesses)
    #pileupTrees = {}
    
    #create the manager for the pileups
    class SynchronyManager(SyncManager):
        pass

    #SynchronyManager.register('numPastStage', Synchronize.numPastStage)
    #SynchronyManager.register('finishedStage', Synchronize.finishedStage)
    
    manager = SynchronyManager()
    manager.start()
    myPileupDict = manager.dict()
    
    try:
        #create any extra workers
        workerThreads = []
        for i in range(1, args.numProcesses):
            p = MergeWorker(readsQueue, resultsQueue, False, args.motherSortedBam, args.fatherSortedBam, args.outMergedBam, args.mergeType, args.isRandomFilter, outHeader, 
                            manager.address, args.numProcesses, i, max_queue_size, myPileupDict, args.keepTemp)
            p.start()
            workerThreads.append(p)
        
        #create the master
        masterProcess = MergeWorker(readsQueue, resultsQueue, True, args.motherSortedBam, args.fatherSortedBam, args.outMergedBam, args.mergeType, args.isRandomFilter, outHeader, 
                                    manager.address, args.numProcesses, 0, max_queue_size, myPileupDict, args.keepTemp)
        masterProcess.start()
        
        masterProcess.join()
        for p in workerThreads:
            p.join()
        
        normalEnd = True
    except KeyboardInterrupt:
        logger.info('Terminating program')
        masterProcess.terminate()
        for p in workerThreads:
            p.terminate()
        
        normalEnd = False
        
    if normalEnd and args.mergeType == PILEUP_MERGE:
        #master needs to merge the pileup files, sort, and index them
        logger.info('[Master] Merging '+str(args.numProcesses)+' pileup files...')
        mergeArgs = ['-f', args.outMergedBam+'.tmp.pileup_all.bam']
            
        for i in range(0, args.numProcesses):
            mergeArgs.append(args.outMergedBam+'.tmp'+str(i)+'.pileup.sorted.bam')
        
        if args.numProcesses > 1:
            pysam.merge(*mergeArgs)
        else:
            pysam.sort(args.outMergedBam+'.tmp0.pileup.sorted.bam', args.outMergedBam+'.tmp.pileup_all')
            
        logger.info('[Master] Creating index for pileup file...')
        pysam.index(args.outMergedBam+'.tmp.pileup_all.bam')
    
        #create the pileup stuff using pools
        #in queue a bunch of chromosome names, sorted descending
        firstFile = pysam.Samfile(args.motherSortedBam)
        sqHead = firstFile.header['SQ']
        sqHead = sorted(sqHead, key=lambda k: k['LN'], reverse=True)
        firstFile.close()
        
        constructedArgs = []
        for x in sqHead:
            constructedArgs.append((args.outMergedBam+'.tmp.pileup_all.bam', x))
        
        #begin processing
        myPool = multiprocessing.Pool(args.numProcesses)
        results = myPool.map(Pileup.createPH, constructedArgs)
        
        #use the same semaphore lock for everything
        #NOTE: this isn't really locking because we allow all processes access at the same time, it's read only so should be fine
        sem = multiprocessing.Semaphore(args.numProcesses)
        
        sharedDict = {}
        for index in range(0, len(results)):
            res = results[index]
            chrom = res.keys()[0]
            values = res[chrom]
            sharedDict[chrom] = multiprocessing.Array('i', len(values), lock=sem)
            for x in range(0, len(values)):
                sharedDict[chrom][x] = values[x]
            results[index] = None
        
        try:
            pws = []
            
            for x in range(0, args.numProcesses):
                #create more pileup files that we can merge in
                pw = Pileup.PileupWorker(sharedDict, resultsQueue, args.outMergedBam, args.isRandomFilter, outHeader, x, args.keepTemp)
                pw.start()
                pws.append(pw)
            
            for pw in pws:
                pw.join()
                
            normalEnd = True
        except KeyboardInterrupt:
            logger.info('Terminating program')
            for pw in pws:
                pw.terminate()
            
            normalEnd = False
        
        if normalEnd:
            os.remove(args.outMergedBam+'.tmp.pileup_all.bam')
            os.remove(args.outMergedBam+'.tmp.pileup_all.bam.bai')
    
    if normalEnd:
        #first parse the results
        res = []
        filenameToResults = {}
        while True:
            try:
                data = resultsQueue.get_nowait()
            except:
                break
            
            pcData = data['percentageChoice']
            filenameToResults[data['outputFilename']] = data
            res.append(pcData)
        
        if args.numProcesses > 1 or args.mergeType == PILEUP_MERGE:
            mergeArgs = ['-fn', args.outMergedBam]
            
            #this loop makes sure that there is something in the file we're trying to merge, it crashes samtools merge if 
            #we try to merge something with an empty file
            for i in range(0, args.numProcesses):
                totVals = filenameToResults[args.outMergedBam+'.tmp'+str(i)+'.bam']['tot']
                if totVals['1'] > 0 or totVals['2'] > 0 or totVals['3'] > 0:
                    mergeArgs.append(args.outMergedBam+'.tmp'+str(i)+'.bam')
                
                if args.mergeType == PILEUP_MERGE:
                    totVals = filenameToResults[args.outMergedBam+'.tmp'+str(i)+'.bam.pileup_complete.bam']['tot']
                    if totVals['1'] > 0 or totVals['2'] > 0 or totVals['3'] > 0:
                        mergeArgs.append(args.outMergedBam+'.tmp'+str(i)+'.bam.pileup_complete.bam')
            
            #print mergeArgs
            logger.info('Merging '+str(len(mergeArgs)-2)+' results files...')
            if len(mergeArgs) == 2:
                logger.warning('No files to merge, output not generated.  Check inputs for valid data.')
            elif len(mergeArgs) == 3:
                #only one file to "merge", just rename it
                os.rename(mergeArgs[2], mergeArgs[1])
            else:
                pysam.merge(*mergeArgs)
            
            if not args.keepTemp:
                logger.info('Cleaning up...')
                for i in range(0, args.numProcesses):
                    try:
                        os.remove(args.outMergedBam+'.tmp'+str(i)+'.bam')
                    except:
                        logger.info('Failed to remove '+args.outMergedBam+'.tmp'+str(i)+'.bam from the file system.')
                        
                    if args.mergeType == PILEUP_MERGE:
                        try:
                            os.remove(args.outMergedBam+'.tmp'+str(i)+'.bam.pileup_complete.bam')
                        except:
                            logger.info('Failed to remove '+args.outMergedBam+'.tmp'+str(i)+'.bam.pileup_complete.bam from the file system.')
        
        logger.info('Merge complete!')
        
        if args.chartFilename != None and args.mergeType == PILEUP_MERGE:
            sumData = [sum(a) for a in zip(*res)]
            saveChoiceChart(sumData, args.chartFilename, args.outMergedBam)
        
if __name__ == '__main__':
    mainRun()