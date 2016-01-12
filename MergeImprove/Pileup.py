'''
Created on Mar 12, 2013

@author: holtjma
'''

import MergeImprove
import pysam#@UnresolvedImport
import multiprocessing
import logging
import os
import random

class PileupWorker(multiprocessing.Process):
    def __init__(self, sharedDict, resultsQueue, outputFilename, isRandomFilter, outHeader, workerID, keepAll, numInputs):
        '''
        @param sharedDict - a dictionary full of multiprocessing.Arrays
        @param resultsQueue - the shared queue for storing statistics on the alignments from each worker
        @param outputFilename - the overall output filename where the final merge will be, worker won't use this unless a single process scenario
        @param isRandomFilter - boolean determining whether we random or union at the end
        @param outHeader - the header to use in any output files
        @param workerID - integer ID, 0 means master
        '''
        global logger
        logger = MergeImprove.getLogger()
        
        #this first
        multiprocessing.Process.__init__(self)
        
        #my vars
        self.sharedDict = sharedDict
        
        #always get these stats
        '''
        self.statistics = {'K':{'1':0,'2':0,'3':0},
                           'U':{'1':0,'2':0,'3':0},
                           'Q':{'1':0,'2':0,'3':0},
                           'P':{'1':0,'2':0,'3':0},
                           'R':{'1':0,'2':0,'3':0},
                           'tot':{'1':0,'2':0,'3':0}}
        '''
        cts = ['K', 'U', 'Q', 'P', 'R', 'tot']
        self.statistics = {}
        maxPO = 2**numInputs
        for ct in cts:
            self.statistics[ct] = [0 for x in range(0, maxPO)]
        
        self.percentageChoice = [0 for x in range(0, 101)]
        
        #save the inputs from the init
        self.resultsQueue = resultsQueue
        self.baseOutputFN = outputFilename
        self.outHeader = outHeader
        self.isRandomFilter = isRandomFilter
        
        self.pileupTempFN = outputFilename+'.tmp'+str(workerID)+'.bam'+'.pileup_tmp.bam'
        self.outputFN = outputFilename+'.tmp'+str(workerID)+'.bam'+'.pileup_complete.bam'
        self.workerID = workerID
        self.keepAll = keepAll
        
    def run(self):
        self.inputFile = pysam.Samfile(self.pileupTempFN, 'rb')
        self.outputFile = pysam.Samfile(self.outputFN, 'wb', header=self.outHeader, referencenames=self.inputFile.references)
        
        self.startPileupAnalysis()
            
        #clean up excess files
        self.cleanUpFiles()
        
        #close the files
        self.inputFile.close()
        self.outputFile.close()
        
        #self.resultsQueue.put(self.statistics)
        #self.resultsQueue.put(self.percentageChoice)
        self.statistics['percentageChoice'] = self.percentageChoice
        self.statistics['outputFilename'] = self.outputFN
        self.resultsQueue.put(self.statistics)
        
    def startPileupAnalysis(self):
        '''
        Begin handling pileup data
        TODO: this is mostly incomplete or too slow right now, needs a rework
        '''
        logger.info('['+str(self.workerID)+'] Beginning post pileup merge process...')
        
        #run through the extras in the temp file
        try:
            read = self.inputFile.next()
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
                    read = self.inputFile.next()
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
        
        #close and remove file from system
        logger.info('['+str(self.workerID)+'] Finished post pileup merge process.  Waiting for other processes...')
    
    def handlePostPileupMerge(self, reads):
        '''
        This function compares a group of alignments and decides which one to keep
        @param reads - a set of reads with the same name to be compared using pileup
        '''
        avgSum = 0
        [pairs, singles] = MergeImprove.pairReads(reads, MergeImprove.PILEUP_HI_TAG)
        
        if len(pairs) != 0:
            bestAvgPileup = -1
            bestPairs = []
            
            for pair in pairs:
                [tot1, bases1] = self.calcPileupStats(pair[0])
                [tot2, bases2] = self.calcPileupStats(pair[1])
                if bases1+bases2 == 0:
                    bases1 = 1
                
                avgPileup = float(tot1+tot2)/(bases1+bases2)
                avgSum += avgPileup
                
                if avgPileup > bestAvgPileup:
                    bestAvgPileup = avgPileup
                    bestPairs = []
                
                if avgPileup == bestAvgPileup:
                    bestPairs.append(pair)
            
            #stats
            if len(bestPairs) == 1:
                MergeImprove.setTag(bestPairs[0][0], MergeImprove.CHOICE_TYPE_TAG, 'P')
                MergeImprove.setTag(bestPairs[0][1], MergeImprove.CHOICE_TYPE_TAG, 'P')
            
                #save pileup stats
                if(bestAvgPileup == 0):
                    self.percentageChoice[0] += 2
                else:
                    self.percentageChoice[int(100*bestAvgPileup/avgSum)] += 2
                
            #save one of the best pileup pairs
            self.saveRandomPair(bestPairs)
            
        else:
            #do this over singles
            bestAvgPileup = {}
            bestReads = {}
            avgSum = {False: 0, True: 0}
            
            for read in singles:
                #if there's nothing yet for this sequence, set it's best as -1 so it gets overwritten below
                isFirst = MergeImprove.isFlagSet(read.flag, MergeImprove.FIRST_SEGMENT_FLAG)
                if not bestAvgPileup.has_key(isFirst):
                    bestAvgPileup[isFirst] = -1
                    bestReads[isFirst] = []
                
                #get the pileup calculation
                [tot, bases] = self.calcPileupStats(read)
                if bases == 0:
                    bases = 1
                avgPileup = float(tot)/bases
                avgSum[isFirst] += avgPileup
                
                #if it's better, keep it
                if avgPileup > bestAvgPileup[isFirst]:
                    bestAvgPileup[isFirst] = avgPileup
                    bestReads[isFirst] = []
                
                if avgPileup == bestAvgPileup[isFirst]:    
                    bestReads[isFirst].append(read)
            
            #save the best from each end
            for end in bestReads:
                brs = bestReads[end]
                
                if len(brs) == 1:
                    MergeImprove.setTag(brs[0], MergeImprove.CHOICE_TYPE_TAG, 'P')
                    if bestAvgPileup[end] == 0:
                        self.percentageChoice[0] += 1
                    else:
                        self.percentageChoice[int(100*bestAvgPileup[end]/avgSum[end])] += 1

                self.saveRandomSingle(brs)
    
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
                
                #get the pileups from this file for that range
                #iterList = self.pileupFile.pileup(self.pileupFile.getrname(chrom), pos, pos+cigLen-1)
                #for pc in iterList:
                    #add the number of pileups at this base
                #    total += pc.n
                total += getTotalPileup(self.sharedDict[chrom], pos, pos+cigLen-1)
                
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
        '''
        this function removes the temporary files from doing the pileup calculations
        '''
        if not self.keepAll:
            #remove all the extra pileup files
            os.remove(self.baseOutputFN+'.tmp'+str(self.workerID)+'.pileup.bam')
            os.remove(self.baseOutputFN+'.tmp'+str(self.workerID)+'.pileup.sorted.bam')
            os.remove(self.pileupTempFN)
        
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
        
        choice = MergeImprove.getTag(readToSave, MergeImprove.CHOICE_TYPE_TAG)
        parent = MergeImprove.getTag(readToSave, MergeImprove.PARENT_OF_ORIGIN_TAG)
        
        #TODO: remove this print error
        if parent == None or choice == None:
            print 'ERROR:'+str(readToSave)
        else:
            if self.statistics.has_key(choice) and parent < len(self.statistics[choice]):
                self.statistics[choice][parent] += 1
                self.statistics['tot'][parent] += 1
            else:
                print 'Poorly structured tag:'
                print self.statistics
                print readToSave
            
        if MergeImprove.verbosity:
            logger.info('Saving from '+MergeImprove.getTag(readToSave, MergeImprove.PARENT_OF_ORIGIN_TAG)+': '+str(readToSave))
        
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
                MergeImprove.setTag(possiblePairs[rv][0], MergeImprove.CHOICE_TYPE_TAG, 'R')
                MergeImprove.setTag(possiblePairs[rv][1], MergeImprove.CHOICE_TYPE_TAG, 'R')
                
                #save the pair
                self.saveRead(possiblePairs[rv][0])
                self.saveRead(possiblePairs[rv][1])
            else:
                for pair in possiblePairs:
                    MergeImprove.setTag(pair[0], MergeImprove.CHOICE_TYPE_TAG, 'K')
                    MergeImprove.setTag(pair[1], MergeImprove.CHOICE_TYPE_TAG, 'K')
                    
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
                MergeImprove.setTag(possibleSingles[rv], MergeImprove.CHOICE_TYPE_TAG, 'R')
                self.saveRead(possibleSingles[rv])
            else:
                for single in possibleSingles:
                    MergeImprove.setTag(single, MergeImprove.CHOICE_TYPE_TAG, 'K')
                    self.saveRead(single)
        else:
            self.saveRead(possibleSingles[0])
        
def createPH(tup):
    fn = tup[0]
    chrom = tup[1]['SN']
    logger = logging.getLogger('root')
    #MergeImprove.initLogger()
    #logger = MergeImprove.getLogger()
    
    myFile = pysam.Samfile(fn, 'rb')
    logger.info('['+chrom+'] Creating Pileup for '+chrom+'...')
    chromNum = myFile.gettid(chrom)
    
    curPos = -1
    curCount = 0
    
    ph = []
    
    for pc in myFile.pileup(chrom):
        #print str(pc.pos)+': '+str(pc.n)
        if curPos + 1 < pc.pos:
            #we skipped a region, add a stop block and change the curCount to zero
            ph.append(curPos+1)
            ph.append(0)
            curCount = 0
            
        #update the current position
        curPos = pc.pos
        
        #check if the counts are different
        if pc.n != curCount:
            #they are, so add a block telling us so
            curCount = pc.n
            ph.append(curPos)
            ph.append(curCount)
        
    ph.append(curPos+1)
    ph.append(0)
    #resultDict[chromNum] = ph
    
    myFile.close()
    logger.info('['+chrom+'] Finished Pileup for '+chrom+'.')
    return {chromNum: ph}

def getPHIndex(pileups, location):
    minIndex = 0
    maxIndex = len(pileups)-2
    
    while True:
        if minIndex >= maxIndex:
            if pileups[minIndex] <= location:
                return minIndex
            else:
                return minIndex-2
        
        midIndex = (minIndex+maxIndex)/2
        if midIndex % 2 == 1:
            midIndex -= 1
        mid = pileups[midIndex]
        
        if mid == location:
            #return mid[1]
            return midIndex
        elif mid > location:
            maxIndex = midIndex-2
        elif mid < location:
            minIndex = midIndex+2

def getTotalPileup(pileup, start, end):
    '''
    Note, this function is inclusive, so [0, 10] covers 11 bases
    @param pileup - the pileup array
    @param start - the first base to be counted
    @param end - the last base to be counted
    '''
    botIndex = getPHIndex(pileup, start)
    phMax = len(pileup)-1
    
    currIndex = botIndex
    currBase = start
    ret = 0
    
    while currIndex < phMax and currBase <= end:
        if currIndex+2 >= phMax:
            ret += (end-currBase+1) * pileup[currIndex+1]
        elif pileup[currIndex+2] <= end:
            ret += (pileup[currIndex+2]-currBase) * pileup[currIndex+1]
        else:
            # >=
            ret += (end-currBase+1) * pileup[currIndex+1]
            
        currIndex += 2
        if currIndex < phMax:
            currBase = pileup[currIndex]
        
    return ret
