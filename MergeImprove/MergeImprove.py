'''
Created on Oct 26, 2012

@author: holtjma@cs.unc.edu
'''

import pysam #@UnresolvedImport
import re
import Synchronize
import random
import logging
import argparse as ap
import util
import os
import time
import multiprocessing
from multiprocessing.managers import SyncManager

DESC = "A merger of genome alignments."
VERSION = '0.1.1'
PKG_VERSION = '0.1.1'

#constant flags from the sam specs
MULTIPLE_SEGMENT_FLAG = 1 << 0 #0x01
BOTH_ALIGNED_FLAG = 1 << 1#0x02
#SEGMENT_UNMAPPED_FLAG = 1 << 2#through testing, this flag was never used
OTHER_UNMAPPED_FLAG = 1 << 3#0x08, for all our reads, this SHOULD be (not BOTH_ALIGNED_FLAG)
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
PRE_PILEUP_STAGE = 1
PILEUP_STAGE = 2
CLEANUP_STAGE = 3
FINISHED_STAGE = 4
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

verbosity = False

#pileupTrees = {}
#treeDictLock = threading.Lock()

#manager = AvlTree.TreeManager()
#manager.start()
#treeDictLock = manager.Lock()
#treeLocks = manager.dict()
#pileupTrees = manager.dict() 

class MergeWorker(multiprocessing.Process):
    
    def __init__(self, readsQueue, resultsQueue, isMaster, firstFilename, secondFilename, outputFilename, mergerType, managerAddress, numProcs, workerID, max_queue_size):
        multiprocessing.Process.__init__(self)
        
        #always get these stats
        self.statistics = {'K':{'1':0,'2':0,'3':0},'U':{'1':0,'2':0,'3':0},'Q':{'1':0,'2':0,'3':0},'P':{'1':0,'2':0,'3':0},'R':{'1':0,'2':0,'3':0}, 'tot':{'1':0,'2':0,'3':0}}
        
        #save the inputs from the init
        self.readsQueue = readsQueue
        self.resultsQueue = resultsQueue
        self.isMaster = isMaster
        self.firstFilename = firstFilename
        self.secondFilename = secondFilename
        self.baseOutputFN = outputFilename
        self.outputFilename = outputFilename+'.tmp'+str(workerID)+'.bam'
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
        class LocalTreeManager(SyncManager):
            pass
        
        LocalTreeManager.register('updateTree')
        LocalTreeManager.register('inorder')
        LocalTreeManager.register('getPileupHeight')
        LocalTreeManager.register('getTotalPileupHeight')
        LocalTreeManager.register('getTreeKeys')
        LocalTreeManager.register('numPastStage')
        LocalTreeManager.register('finishedStage')
        
        self.ltm = LocalTreeManager(managerAddress)
        self.ltm.connect()

    def finishAndBarrierSync(self, stage):
        '''
        @param stage - the stage we are waiting to be complete
        '''
        self.ltm.finishedStage(stage)
        while int(str(self.ltm.numPastStage(stage))) < self.numProcs:
            time.sleep(1)
        
        #once we get here we know everything is at or past this stage
        
    def run(self):
        if self.isMaster:
            self.fillMergeQueue()
        else:
            self.emptyMergeQueue()
    
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
        
        #TODO: samfile open checks
        
        #open the output file as well
        self.outputFile = pysam.Samfile(self.outputFilename, 'wb', template=firstFile)
        
        if self.mergerType == PILEUP_MERGE:
            #this one stores pieces counted in the pileup
            self.pileupTempFile = pysam.Samfile(self.pileupTempPrefix+'.bam', 'wb', template=firstFile)
            
            #this one stores those to be resolved by pileup
            self.tempFN = self.outputFilename+'.pileup_tmp.bam'
            self.tempFile = pysam.Samfile(self.tempFN, 'wb', template=firstFile)
        
        #TODO: update headers?
        #print firstFile.header
        
        #start with no qname
        currentQname = None
        MAX_CLUSTER_SIZE = 100
        currClusterSize = 0
        clusterStart = None
        clusterEnd = None
        
        #get the first read ID from each file
        firstRead = firstFile.next()
        firstQname = firstRead.qname
        firstFlags = firstRead.flag
        
        secondRead = secondFile.next()
        secondQname = secondRead.qname
        secondFlags = secondRead.flag
        
        #pick the 'lowest' or first qname
        if isQnameBefore(re.split(':|#', firstQname), re.split(':|#', secondQname)):
            currentQname = firstQname
        else:
            currentQname = secondQname
        
        #this is stat gathering items
        pairingsDictionary = {}
        pairingRep = [0, 0]
        pairingRep2 = [0, 0]
        firstPairingSeq = None
        
        firstReads = []
        secondReads = []
        
        #used for output and debugging purposes
        count = 0
        maxReads = None#4000000 #TODO: change before upload
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
                #TODO: handle cases of unmapped reads, there is a flag for this
                #check if we should 'count' it
                if isFlagSet(firstFlags, MULTIPLE_SEGMENT_FLAG) and \
                   (isFlagSet(firstFlags, BOTH_ALIGNED_FLAG) or not isFlagSet(firstFlags, OTHER_UNMAPPED_FLAG)) and \
                   not isFlagSet(firstFlags, FIRST_SEGMENT_FLAG):
                    #if it has multiple segments, both are aligned, and this isn't the first segment, don't count it
                    pass
                else:
                    #count it
                    if firstPairingSeq == None:
                        firstPairingSeq = firstRead.seq
                        
                    if firstRead.seq == firstPairingSeq:
                        pairingRep[0] += 1
                    else:
                        pairingRep2[0] += 1
                        
                #get the next one for analysis
                try:
                    firstRead = firstFile.next()
                    firstQname = firstRead.qname
                    firstFlags = firstRead.flag
                    
                except:
                    firstRead = None
                    firstQname = None
                    firstFlags = None
                    
                isMatch = True
                
            #check the second also and get the next read if needed
            while secondQname == currentQname and secondQname != None:
                #always append the read
                secondReads.append(secondRead)
                
                #check if we should 'count' it
                if isFlagSet(secondFlags, MULTIPLE_SEGMENT_FLAG) and \
                   (isFlagSet(secondFlags, BOTH_ALIGNED_FLAG) or not isFlagSet(secondFlags, OTHER_UNMAPPED_FLAG)) and \
                   not isFlagSet(secondFlags, FIRST_SEGMENT_FLAG):
                    #if it has multiple segments, both are aligned, and this isn't the first segment, don't count it
                    pass
                else:
                    #count it
                    if firstPairingSeq == None:
                        firstPairingSeq = secondRead.seq
                        
                    if secondRead.seq == firstPairingSeq:
                        pairingRep[1] += 1
                    else:
                        pairingRep2[1] += 1
                    
                #get the next one for analysis
                try:
                    secondRead = secondFile.next()
                    secondQname = secondRead.qname
                    secondFlags = secondRead.flag
                except:
                    secondRead = None
                    secondQname = None
                    secondFlags = None
                    
                isMatch = True
                
            #check if there are no reads matching the current qname
            if not isMatch:
                
                #debug statement every 100k counts where a count is a distinct read name
                if count % 100000 == 0:
                    logger.info('[Master] Processed '+str(count)+' read names...')
                    
                #increment the count always when we update the qname
                count += 1
                
                #no match, time for the next qname after storing the data
                if not pairingsDictionary.has_key(pairingRep[0]):
                    pairingsDictionary[pairingRep[0]] = {}
                
                if not pairingsDictionary.has_key(pairingRep2[0]):
                    pairingsDictionary[pairingRep2[0]] = {}
                
                if not pairingsDictionary[pairingRep[0]].has_key(pairingRep[1]):
                    pairingsDictionary[pairingRep[0]][pairingRep[1]] = 0
                
                if not pairingsDictionary[pairingRep2[0]].has_key(pairingRep2[1]):
                    pairingsDictionary[pairingRep2[0]][pairingRep2[1]] = 0
                
                #increment the number of values with that pair type by 1, stats gathering info 
                pairingsDictionary[pairingRep[0]][pairingRep[1]] += 1
                
                if pairingRep2[0] != 0 or pairingRep2[1] != 0:
                    pairingsDictionary[pairingRep2[0]][pairingRep2[1]] += 1
                    
                #pick from the reads which ones make the most sense
                #self.handleSingleIdMerge(firstReads, secondReads)
                #self.readsQueue.put([firstReads, secondReads])
                #self.readsQueue.put([firstReads, ["c", "d"]])
                
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
                        self.readsQueue.put([clusterStart, clusterEnd])
                            
                        #there is space in the queue
                        if clusterEnd == currentQname:
                            #we just added this to the end of the cluster AND we know the cluster is full, send it to be processed
                            #reset
                            clusterStart = None
                            clusterEnd = None
                            currClusterSize = 0
                        else:
                            #we should restart the cluster with this extra read
                            clusterStart = currentQname
                            clusterEnd = currentQname
                            currClusterSize = 1
                            
                    
                #else:
                #    self.readsQueue.put(currentQname)
                
                #reset the pairing
                pairingRep = [0, 0]
                pairingRep2 = [0, 0]
                firstPairingSeq = None
                firstReads = []
                secondReads = []
                
                #get the next qname
                if firstQname != None and (secondQname == None or isQnameBefore(re.split(':|#', firstQname), re.split(':|#', secondQname))):
                    currentQname = firstQname
                else:
                    currentQname = secondQname
        
        if clusterStart != None:
            #remainder cluster, put it on 
            self.readsQueue.put([clusterStart, clusterEnd])
            
            #reset, just in case
            clusterStart = None
            clusterEnd = None
            currClusterSize = 0
            
        #tell the other processes it's the end of the queue
        for x in range(0, self.numProcs-1):
            self.readsQueue.put(END_OF_QUEUE)
        
        logger.info('[Master] All reads scanned.  Waiting on workers to finish processing...')
        
        #this handles the sorting and indexing of pileup files
        if self.mergerType == PILEUP_MERGE:
            self.prepareForPileup()
        
        #mark that we're done and barrier sync
        self.finishAndBarrierSync(QUALITY_STAGE)
        
        if self.mergerType == PILEUP_MERGE:
            #master needs to merge the pileup files, sort, and index them
            logger.info('[Master] Merging '+str(self.numProcs)+' pileup files...')
            
            mergeArgs = ['-f', self.baseOutputFN+'.tmp.pileup_all.bam']
                
            for i in range(0, self.numProcs):
                mergeArgs.append(self.baseOutputFN+'.tmp'+str(i)+'.pileup.sorted.bam')
            
            pysam.merge(*mergeArgs)
            
            #logger.info('[Master] Sorting pileup file for searching...')
            
            #pysam.sort(self.baseOutputFN+'.tmp.pileup_all.bam', self.baseOutputFN+'.tmp.pileup_all.sorted')
            
            logger.info('[Master] Creating index for pileup file...')
            
            pysam.index(self.baseOutputFN+'.tmp.pileup_all.bam')
        
        #mark that we're done and barrier sync
        self.finishAndBarrierSync(PRE_PILEUP_STAGE)
        
        #do the pileup analysis if necessary
        if self.mergerType == PILEUP_MERGE:
            self.startPileupAnalysis()
        
        #mark that pileup is done and barrier sync
        self.finishAndBarrierSync(PILEUP_STAGE)
        
        #clean up excess files
        self.cleanUpFiles()
        
        #sync and finish
        self.finishAndBarrierSync(CLEANUP_STAGE)
        
        #close the files
        firstFile.close()
        secondFile.close()
        self.outputFile.close()
        
        '''
        treekeys = self.ltm.getTreeKeys()
        for treeId in treekeys:
            print str(treeId)+':'+str(self.ltm.inorder(treeId))
        '''
        
        logger.info('[0] Finished!')
        self.resultsQueue.put(self.statistics)
        
        return pairingsDictionary
    
    def emptyMergeQueue(self):
        
        random.seed()
        
        #attempt to open both bam files for merging
        firstFile = pysam.Samfile(self.firstFilename, 'rb')
        secondFile = pysam.Samfile(self.secondFilename, 'rb')
        self.outputFile = pysam.Samfile(self.outputFilename, 'wb', template=firstFile)
        
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
            
        #while (not self.isMaster and self.stage == QUALITY_STAGE) or (self.isMaster and ((self.readingComplete and self.stage == QUALITY_STAGE) or self.readsQueue.qsize() > MAX_QUEUE_SIZE / 2)):
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
                
                startQname = qnames[0]
                endQname = qnames[1]
                
                #skip all reads before the start qname in the first file
                while firstRead != None and isQnameBefore(re.split(':|#', firstQname), re.split(':|#', startQname)):
                    try:
                        firstRead = firstFile.next()
                        firstQname = firstRead.qname
                    except:
                        firstRead = None
                        firstQname = None
                
                #skip all reads before the start qname in the second file
                while secondRead != None and isQnameBefore(re.split(':|#', secondQname), re.split(':|#', startQname)):
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
                elif isQnameBefore(re.split(':|#', firstQname), re.split(':|#', secondQname)):
                    currentQname = firstQname
                else:
                    currentQname = secondQname
                
                #while we have a current qname AND (that qname is before the end of the block or is at the end of the block)
                while currentQname != None and (isQnameBefore(re.split(':|#', currentQname), re.split(':|#', endQname)) or currentQname == endQname):
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
                    
                    #get the next qname
                    if firstRead == None:
                        currentQname = secondQname
                    elif secondRead == None:
                        currentQname = firstQname
                    elif isQnameBefore(re.split(':|#', firstQname), re.split(':|#', secondQname)):
                        currentQname = firstQname
                    else:
                        currentQname = secondQname
        
        
        firstFile.close()
        secondFile.close()
        
        #this handles the sorting and indexing of pileup files
        if self.mergerType == PILEUP_MERGE:
            self.prepareForPileup()
        
        #mark that we're done and barrier sync
        self.finishAndBarrierSync(QUALITY_STAGE)
        
        #we don't have to do anything in the pre-pileup, it's all on the master
        self.finishAndBarrierSync(PRE_PILEUP_STAGE)
        
        if self.mergerType == PILEUP_MERGE:
            self.startPileupAnalysis()
        
        #mark that pileup is done and barrier sync
        self.finishAndBarrierSync(PILEUP_STAGE)
        
        #clean up excess files
        self.cleanUpFiles()
        
        #sync and finish
        self.finishAndBarrierSync(CLEANUP_STAGE)
        
        self.outputFile.close()
        
        '''
        treekeys = self.ltm.getTreeKeys()
        for treeId in treekeys:
            print str(treeId)+':'+str(self.ltm.inorder(treeId))
        '''
        
        logger.info('['+str(self.workerID)+'] Finished!')
        
        self.resultsQueue.put(self.statistics)
    
    def handleSingleIdMerge(self, firstReads, secondReads):
        '''
        This function takes two sets of reads and attempts to pick one of the best alignments to save using region based likelihood
        @param firstReads - the set of reads coming from the first parent
        @param secondReads - the set of reads coming from the second parent
        '''
        
        #if len(firstReads) > 0 and firstReads[0].qname == 'UNC11-SN627_0160:4:1101:2486:133373#TGACCA':
        #dumpReads(firstReads, secondReads)
        
        #pair the reads using their HI tags
        [firstPairs, firstSingles] = pairReads(firstReads, 'HI')
        [secondPairs, secondSingles] = pairReads(secondReads, 'HI')
        
        #TODO: remove this later
        if (len(firstPairs) > 0 and len(firstSingles) > 0) or (len(secondPairs) > 0 and len(secondSingles) > 0):
            dumpReads(firstReads, secondReads)
        
        interestList = []
        '''
        interestList = ['UNC11-SN627_0160:4:1101:1181:46840#TAGCTT', 'UNC11-SN627_0160:4:1101:1181:97059#TAGCTT', 'UNC11-SN627_0160:4:1101:1181:140573#TAGCTT', 
                        'UNC11-SN627_0160:4:1101:1181:172414#TAGCTT', 'UNC11-SN627_0160:4:1101:1182:3400#TAGCTT', 'UNC11-SN627_0160:4:1101:1184:123367#TAGCTT']
        '''
        
        if len(firstReads) > 0:
            for value in interestList:
                if firstReads[0].qname == value:
                    dumpReads(firstReads, secondReads)
        
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
                        #import pydevd;pydevd.settrace()
        
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
        #close the file that has the reads for the pileup calculations
        self.pileupTempFile.close()
        
        logger.info('['+str(self.workerID)+'] Sorting the pileup data...')
        
        #sort the current output file into normal positional sort
        pysam.sort(self.pileupTempPrefix+'.bam', self.pileupTempPrefix+'.sorted')
        
        #logger.info('['+str(self.workerID)+'] Building the pileup index...')
        
        #build the index for the file
        #pysam.index(self.pileupTempPrefix+'.sorted.bam')
        
        #logger.info('['+str(self.workerID)+'] Pileup data complete.  Waiting for other processes...')
        
    def startPileupAnalysis(self):
        logger.info('['+str(self.workerID)+'] Beginning post pileup merge process...')
        
        #close the file, any discrepancies are added already
        self.tempFile.close()
        
        #open the file for reading now
        #TODO: remove this?
        self.percentageChoice = [0 for x in range(0, 101)]
        self.tempFile = pysam.Samfile(self.tempFN, 'rb')
        
        #self.pileupFiles = []
        #for wid in range(0, self.numProcs):
        #    self.pileupFiles.append(pysam.Samfile(self.baseOutputFN+'.tmp'+str(wid)+'.pileup.sorted.bam', 'rb'))
        
        self.pileupFile = pysam.Samfile(self.baseOutputFN+'.tmp.pileup_all.bam', 'rb')
        
        #run through the extras in the temp file
        try:
            read = self.tempFile.next()
            readQname = read.qname
        except:
            read = None
            readQname = None
        
        #init
        currentQname = readQname
        isMatch = False
        readsInSet = []
        
        #main loop to get all the reads
        while read != None or isMatch:
            isMatch = False
            
            while readQname == currentQname and readQname != None:
                #append the read
                readsInSet.append(read)
                isMatch = True
                
                #see if you can get another
                try:
                    read = self.tempFile.next()
                    readQname = read.qname
                except:
                    read = None
                    readQname = None
            
            #now we have all the reads with the same name
            if not isMatch:
                #magic to handle the set
                self.handlePostPileupMerge(readsInSet)
                
                #reset these
                readsInSet = []
                currentQname = readQname
        
        #for fileToClose in self.pileupFiles:
        #    fileToClose.close()
        self.pileupFile.close()
        
        #close and remove file from system
        self.tempFile.close()
        os.remove(self.tempFN)
        
        logger.info('['+str(self.workerID)+'] Finished post pileup merge process.  Waiting for other processes...')
    
    def handlePostPileupMerge(self, reads):
        #dumpReads(reads, [])
        
        avgSum = 0
        [pairs, singles] = pairReads(reads, PILEUP_HI_TAG)
        
        if len(pairs) != 0:
            bestAvgPileup = -1
            bestPairs = []
            
            for pair in pairs:
                [tot1, bases1] = self.calcPileupStats(pair[0])
                [tot2, bases2] = self.calcPileupStats(pair[1])
                
                avgPileup = float(tot1+tot2)/(bases1+bases2)
                avgSum += avgPileup
                
                if avgPileup > bestAvgPileup:
                    bestAvgPileup = avgPileup
                    bestPairs = []
                
                if avgPileup == bestAvgPileup:
                    bestPairs.append(pair)
            
            #stats
            if len(bestPairs) == 1:
                setTag(bestPairs[0][0], CHOICE_TYPE_TAG, 'P')
                setTag(bestPairs[0][1], CHOICE_TYPE_TAG, 'P')
                
            #save one of the best pileup pairs
            self.saveRandomPair(bestPairs)
            
            if(bestAvgPileup == 0):
                self.percentageChoice[0] += 1
            else:
                self.percentageChoice[int(100*bestAvgPileup/avgSum)] += 1
            
        else:
            #do this over singles
            bestAvgPileup = {}
            bestReads = {}
            avgSum = 0
            
            for read in singles:
                #if there's nothing yet for this sequence, set it's best as -1 so it gets overwritten below
                if not bestAvgPileup.has_key(read.seq):
                    bestAvgPileup[read.seq] = -1
                    bestReads[read.seq] = []
                
                #get the pileup calculation
                [tot, bases] = self.calcPileupStats(read)
                avgPileup = float(tot)/bases
                avgSum += avgPileup
                
                #if it's better, keep it
                if avgPileup > bestAvgPileup[read.seq]:
                    bestAvgPileup[read.seq] = avgPileup
                    bestReads[read.seq] = []
                
                if avgPileup == bestAvgPileup[read.seq]:    
                    bestReads[read.seq].append(read)
            
            #save the best from each end
            finAvg = 0
            for seq in bestReads:
                brs = bestReads[seq]
                
                if len(brs) == 1:
                    setTag(brs[0], CHOICE_TYPE_TAG, 'P')
                    
                self.saveRandomSingle(brs)
                finAvg += bestAvgPileup[seq]
                
            if finAvg == 0:
                self.percentageChoice[0] += 1
            else:
                self.percentageChoice[int(100*finAvg/avgSum)] += 1

    def calcPileupStats(self, read):
        '''
        @param read - the read we are calculating
        @return a pair [total coverage, bases] representing the sum of all coverages over 'bases' number of base pairs
        '''
        total = 0
        bases = 0
        
        cig = read.cigar
        pos = read.pos
        chrom = read.rname
        
        for cigTypePair in cig:
            cigType = cigTypePair[0]
            cigLen = cigTypePair[1]
            
            if cigType == 0:
                '''
                total += int(str(self.ltm.getTotalPileupHeight(chrom, pos, pos+cigLen-1)))
                '''
                #iterate through each pileup file
                #for pf in self.pileupFiles:
                
                #get the pileups from this file for that range
                iterList = self.pileupFile.pileup(self.pileupFile.getrname(chrom), pos, pos+cigLen-1)
                for pc in iterList:
                    #add the number of pileups at this base
                    total += pc.n
            
                #modify the bases and position
                bases += cigLen
                pos += cigLen
                
            elif cigType == 1:
                #insertion to the reference
                #no change to position, these will be ignored bases though
                pass
            elif cigType == 2:
                #deletion to the reference
                #modify position by the second value
                pos += cigLen
            elif cigType == 3:
                #skipped region
                #modify the position by the second value
                pos += cigLen
            else:
                print 'UNHANDLED CIG TYPE:'+cigType
                
        return [total, bases]
        
    def cleanUpFiles(self):
        if self.mergerType == PILEUP_MERGE:
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
        
        #TODO: remove this print error
        if parent == None or choice == None:
            print 'ERROR:'+str(readToSave)
        else:
            self.statistics[choice][parent] += 1
            self.statistics['tot'][parent] += 1
            
            #TODO: remove this
            #if choice == 'Q' and parent == '3':
            #    print readToSave
            
        if verbosity:
            logger.info('Saving from '+getTag(readToSave, PARENT_OF_ORIGIN_TAG)+': '+str(readToSave))
        
    def saveRandomPair(self, possiblePairs):
        '''
        @param possiblePairs - a list of possible pairs that could be saved
        @param parent - the parent of origin for these pairs, can be SET indicating a mix of pairs that already have it set
        '''
        #pick a random pair
        rv = random.randint(0, len(possiblePairs)-1)
        
        if len(possiblePairs) > 1:
            #random choice
            setTag(possiblePairs[rv][0], CHOICE_TYPE_TAG, 'R')
            setTag(possiblePairs[rv][1], CHOICE_TYPE_TAG, 'R')
            
        #save the pair
        self.saveRead(possiblePairs[rv][0])
        self.saveRead(possiblePairs[rv][1])
        
        
    def saveRandomSingle(self, possibleSingles):
        '''
        @param possibleSingles - a list of possible singles that could be saved
        @param parent - the parent of origin for these singles, can be SET indicating a mix of pairs that already have it set
        '''
        #pick a random single and save it
        rv = random.randint(0, len(possibleSingles)-1)
        
        if len(possibleSingles) > 1:
            setTag(possibleSingles[rv], CHOICE_TYPE_TAG, 'R')
        
        self.saveRead(possibleSingles[rv])
    
    def parseReadForTree(self, readToParse):
        if readToParse == None:
            return
        
        #save the read
        self.pileupTempFile.write(readToParse)
        
        '''
        #extract the chromosome ID so we know which tree to modify
        chromId = readToParse.rname
        
        #now we have to parse the reads into regions
        cig = readToParse.cigar
        pos = readToParse.pos
        
        for pair in cig:
            if pair[0] == 0 or pair[0] == 2:
                #this is a match/mismatch, add the starting position, the shift and add the ending position
                success = self.ltm.updateTree(chromId, pos, True)
                pos += pair[1]
                success = self.ltm.updateTree(chromId, pos, False)
            elif pair[0] == 1:
                #do nothing, this is an insertion TO the ref
                pass
            elif pair[0] == 3:
                #skip this region
                pos += pair[1]
            else:
                logger.warning('Unhandled cigar type:'+pair[0])
        '''
            
def splitSingles(singles):
    '''
    @param singles - the set of unpaired reads from a parent
    @return - two separate lists of values such that the reads are grouped based on the sequence they aligned two
    '''
    #init
    firstSeq = None
    firstSet = []
    
    secondSeq = None
    secondSet = []
    
    #time to split the singles up
    for read in singles:
        #if there isn't a first sequence, mark this one as the first
        if firstSeq == None:
            firstSeq = read.seq
            firstSet.append(read)
        
        #if they are equal append to list
        elif firstSeq == read.seq:
            firstSet.append(read)
            
        #not part of the first, so check if we have a second, then add it to second
        elif secondSeq == None:
            secondSeq = read.seq
            secondSet.append(read)
        else:
            secondSet.append(read)
    
    #return the first list and the second list
    return [firstSet, secondSet]


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
    MATCH = 2
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
                logger.warning('Unhandled Cigar string type:'+symbol)
            
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
    #import pydevd;pydevd.settrace()
    #old debug junk
    print 'First:'+str(len(firstReads))
    for read in firstReads:
        print read
        #print hex(read.flag)
    
    print 'Second:'+str(len(secondReads))
    for read in secondReads:
        print read
        #print hex(read.flag)
        
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
    
    #iterate through the reads
    for read in reads:
        if isFlagSet(read.flag, MULTIPLE_SEGMENT_FLAG):
            if isFlagSet(read.flag, BOTH_ALIGNED_FLAG):
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
                setTag(p1[1], P_MISMATCH_TAG, oc2)
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
                
                if not uniqueSingles.has_key(s1.seq):
                    uniqueSingles[s1.seq] = []
                uniqueSingles[s1.seq].append(s1)
                
            else:
                if not uniqueSingles.has_key(s1.seq):
                    uniqueSingles[s1.seq] = []
                
                #different so keep both around
                setTag(s1, PARENT_OF_ORIGIN_TAG, '1')
                uniqueSingles[s1.seq].append(s1)
                
                setTag(s2, PARENT_OF_ORIGIN_TAG, '2')
                uniqueSingles[s2.seq].append(s2)
            
            singles2.remove(s2)
            
        else:
            if not uniqueSingles.has_key(s1.seq):
                uniqueSingles[s1.seq] = []
            
            #no match, mark it as such
            setTag(s1, PARENT_OF_ORIGIN_TAG, '1')
            setTag(s1, YA_TAG, '1')
            uniqueSingles[s1.seq].append(s1)
        
    for s2 in singles2:
        if not uniqueSingles.has_key(s2.seq):
            uniqueSingles[s2.seq] = []

        #add anything leftover here
        setTag(s2, PARENT_OF_ORIGIN_TAG, '2')
        setTag(s2, YA_TAG, '2')
        uniqueSingles[s2.seq].append(s2)
    
    return uniqueSingles
    
def reducePairsByScore(pairs):
    #get the best overall paired score for anything that may be leftover
    bestScore = 0
    bestPairs = []
    for pair in pairs:
        score = calculatePairScore(pair)
        if score > bestScore:
            bestScore = score
            bestPairs = []
            
        if score == bestScore:
            bestPairs.append(pair)
            
    #after the above loops, we have the best score between all detected pairs in either alignment
    return bestPairs

def reduceSinglesByScore(singles):
    bestSingles = {}
    for seq in singles:
        bestScore = 0
        bests = []
        
        for single in singles[seq]:
            score = calculateReadScore(single)
            if score > bestScore:
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

def isQnameBefore(qnameArray1, qnameArray2):
    '''
    @param qnameArray1 - the array version of the first qname split along the ':'
    @param qnameArray2 - the array version of the second qname split along the ':'
    @return - True if qnameArray1 comes before qnameArray2
    '''
    
    if len(qnameArray2) == 0:
        return False
    
    if len(qnameArray1) == 0:
        return True
    
    #if they are identical, recurse
    if qnameArray1[0] == qnameArray2[0]:
        return isQnameBefore(qnameArray1[1:], qnameArray2[1:])
    
    #they are not identical
    if qnameArray1[0].isdigit() and qnameArray2[0].isdigit():
        return int(qnameArray1[0]) < int(qnameArray2[0])
    else:
        return qnameArray1[0] < qnameArray2[0]

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
    if read1.seq == read2.seq and read1.rname == read2.rname and read1.pos == read2.pos:
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
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    #global multilogger
    #multilogger = multiprocessing.log_to_stderr()
    #multilogger.setLevel(multiprocessing.SUBDEBUG)

#if __name__ == '__main__':
def mainRun():
    #start up the logger
    initLogger()
    
    #attempt to parse the arguments
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #version data
    p.add_argument('-V', '--version', action='version', version='%(prog)s' + \
                   ' %s in Suspenders %s' % (VERSION, PKG_VERSION))
    
    group = p.add_mutually_exclusive_group()
    group.add_argument('-u', '--union', dest='mergeType', action='store_const', const=UNION_MERGE, help='merge the files using the union technique', default=RANDOM_MERGE)
    group.add_argument('-s', '--quality', dest='mergeType', action='store_const', const=RANDOM_MERGE, help='merge the files using the quality score, then randomly choosing (default)', default=RANDOM_MERGE)
    group.add_argument('-t', '--pileup', dest='mergeType', action='store_const', const=PILEUP_MERGE, help='merge the files using the quality score, then the pileup heights, then randomly choosing', default=RANDOM_MERGE)
    
    #optional argument
    #p.add_argument('-c', dest='showChart', action='store_true', help='show distribution chart (default: no)')
    #p.add_argument('-v', dest='verbose', action='store_true', help='verbose (default: no)')
    
    p.add_argument('-p', metavar='numProcesses', dest='numProcesses', type=int, default=1, help='number of processes to run (default: 1)')
    
    #required main arguments
    p.add_argument('motherSortedBam', type=util.readableFile, help='the bam file of the sorted by name alignment to the mother pseudogenome')
    p.add_argument('fatherSortedBam', type=util.readableFile, help='the bam file of the sorted by name alignment to the father pseudogenome')
    p.add_argument('outMergedBam', type=util.writableFile, help='the output bam file merged and sorted by name')
    
    #parse the arguments
    args = p.parse_args()
    
    #check if we should be wordy
    #if args.verbose:
    #    verbosity = True
    
    if args.mergeType == UNION_MERGE:
        mt = 'Union'
    elif args.mergeType == RANDOM_MERGE:
        mt = 'Quality->Random'
    elif args.mergeType == PILEUP_MERGE:
        mt = 'Quality->Pileup->Random'
    else:
        mt = 'ERROR: Unknown'
    
    logger.info('Merge Type: '+mt)
    logger.info('Mother: '+args.motherSortedBam)
    logger.info('Father: '+args.fatherSortedBam)
    logger.info('Output: '+args.outMergedBam)
    logger.info('Number of processes: '+str(args.numProcesses))
    
    readsQueue = multiprocessing.Queue()
    resultsQueue = multiprocessing.Queue()
    
    max_queue_size = multiprocessing.Value('i', 100)
    #pileupTrees = {}
    
    #create a merger
    #merger = MergeImprove(args.mergeType, args.numProcesses)
    
    #create the manager for the pileups
    class TreeManager(SyncManager):
        pass

    #TreeManager.register('Tree', AvlTree)
    #TreeManager.register('AvlNode', AvlNode)
    #TreeManager.register('updateTree', AvlTree.updateTree)
    #TreeManager.register('inorder', AvlTree.inorder)
    #TreeManager.register('getPileupHeight', AvlTree.getPH)
    #TreeManager.register('getTotalPileupHeight', AvlTree.getTotalPH)
    #TreeManager.register('getTreeKeys', AvlTree.getTreeKeys, ListProxy)
    TreeManager.register('numPastStage', Synchronize.numPastStage)
    TreeManager.register('finishedStage', Synchronize.finishedStage)
    
    manager = TreeManager()
    manager.start()
    
    try:
        #create any extra workers
        workerThreads = []
        for i in range(1, args.numProcesses):
            p = MergeWorker(readsQueue, resultsQueue, False, args.motherSortedBam, args.fatherSortedBam, args.outMergedBam, args.mergeType, manager.address, args.numProcesses, i, max_queue_size)
            p.start()
            workerThreads.append(p)
        
        #create the master
        masterProcess = MergeWorker(readsQueue, resultsQueue, True, args.motherSortedBam, args.fatherSortedBam, args.outMergedBam, args.mergeType, manager.address, args.numProcesses, 0, max_queue_size)
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
    
    if normalEnd:
        if args.numProcesses > 1:
            logger.info('Merging '+str(args.numProcesses)+' result files...')
            
            mergeArgs = ['-f', args.outMergedBam]
            
            for i in range(0, args.numProcesses):
                mergeArgs.append(args.outMergedBam+'.tmp'+str(i)+'.bam')
            
            pysam.merge(*mergeArgs)
            
            logger.info('Cleaning up...')
            for i in range(0, args.numProcesses):
                os.remove(args.outMergedBam+'.tmp'+str(i)+'.bam')
            
            logger.info('Merge complete!')
            
        #TODO: after these all join, I have several files that need merging
        #TODO: get the merge statistics
        
        #merge the files
        #pairingsDict = merger.mergeAnnotatedBamFiles(args.motherSortedBam, args.fatherSortedBam, args.outMergedBam)
        
        #get some stats dumped for us
        #merger.dumpStats()
        
        #check for validity
        '''
        if pairingsDict != None:
            #create a chart of the output
            logger.info('Finished!')
            
            if args.showChart:
                gc = GradientChart.GradientChart(['Mother Potential Alignments', 'Father Potential Alignments'], pairingsDict)
                gc.drawChart()
                
                merger.showChoiceChart()
        else:
            logger.info('Failure during merge process!')
        '''
if __name__ == '__main__':
    mainRun()