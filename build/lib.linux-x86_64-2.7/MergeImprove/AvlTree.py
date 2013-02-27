'''
Created on Jan 14, 2013

@author: holtjma
'''

from multiprocessing.managers import SyncManager, ListProxy
from multiprocessing import Process

class AvlNode:
    def __init__(self, value):
        self.value = value
        self.height = 1
        self.delta = 0
        self.left = None
        self.right = None
        
    def getHeight(self):
        return self.height
    
    def incrementDelta(self, change=1):
        self.delta += change
        
    def decrementDelta(self, change=1):
        self.delta -= change
        
    def getLeftHeight(self):
        if self.left == None:
            return 0
        else:
            return self.left.getHeight()
        
    def getRightHeight(self):
        if self.right == None:
            return 0
        else:
            return self.right.getHeight()
    
    def evaluateHeight(self):
        lh = self.getLeftHeight()
        rh = self.getRightHeight()
        self.height = max(lh, rh)+1
    
    def inorder(self):
        ret = '['
        
        if self.left != None:
            ret += self.left.inorder()+','
            
        ret += ' '+str(self.value)+':'+str(self.delta)+' h'+str(self.height)
        
        if self.right != None:
            ret += ','+self.right.inorder()
            
        ret += ']'
        return ret
    
class AvlTree:
    def __init__(self):
        self.root = None
        #self.treeLock = threading.Lock()

    def addStart(self, location):
        #self.treeLock.acquire()
        if self.root == None:
            self.root = AvlNode(location)
            self.root.incrementDelta()
            #self.treeLock.release()
            return True
            
        else:
            return self.updateNode(location, True)
        
    def addEnd(self, location):
        #self.treeLock.acquire()
        if self.root == None:
            self.root = AvlNode(location)
            self.root.decrementDelta()
            #self.treeLock.release()
            return True
        
        else:
            return self.updateNode(location, False)
        
    def getPileupHeight(self, value):
        '''
        This is the main search structure
        @param value - the location on the genome to get the pileup height for
        @return positive int value denoting how many reads put into the tree overlap that position
        '''
        #init
        #self.treeLock.acquire()
        pileupSum = 0
        currentNode = self.root
        isValueFound = False
        
        while currentNode != None and not isValueFound:
            if currentNode.value == value:
                pileupSum += currentNode.delta
                isValueFound = True
                
            elif currentNode.value < value:
                pileupSum += currentNode.delta
                currentNode = currentNode.right
                
            else:
                currentNode = currentNode.left
        
        #self.treeLock.release()
        return pileupSum
    
    def getTotalPileupHeight(self, start, end):
        '''
        This is hopefully a faster way then doing each piece of a pileup sequentially
        @param start - the first base to be considered as part of the pileup
        @param end - the last base to be considered as part of the pileup
        @return positive integer value denoting the summation of the pileup counts for each base pair (to get average, divide by length later)
        '''
        #self.treeLock.acquire()
        pileupSum = 0
        currentNode = self.root
        
        totals = self.getSumPileups(currentNode, pileupSum, start, end)
        #self.treeLock.release()
        return totals
        
    def getSumPileups(self, currentNode, delta, start, end):
        '''
        This is a recursive function for when we want to do an average pileup height
        NOTE: this function alone is not thread-safe, call getTotalPileupHeight(...) instead
        @param currentNode - the node we're looking at right now
        @param delta - the delta so far
        @param start - the first base to be considered of the pileup
        @param end - the last base to be considered as part of the pileup
        @return a sum of the pileup (not average)
        '''
        #general break condition
        if start > end:
            return 0
        
        pileupSum = delta
        isInRange = False
        while currentNode != None and not isInRange:
            if currentNode.value < start:
                pileupSum += currentNode.delta
                currentNode = currentNode.right
            elif currentNode.value > end:
                currentNode = currentNode.left
            else:
                isInRange = True
                
        #at this point we're within the range from before
        if currentNode == None:
            #all the values have this pileup size
            total = pileupSum * (end-start+1)
        else:
            total = pileupSum+currentNode.delta
            total += self.getSumPileups(currentNode, pileupSum, start, currentNode.value-1)
            total += self.getSumPileups(currentNode, pileupSum, currentNode.value+1, end)
        
        return total
        
    def updateNode(self, value, isStartOfRead):
        '''
        This is the main body of the balancing structure
        @param value - the location on the genome
        @param isStartOfRead - indicates if this the start or end position of a read
        @return - True if it's successfully updated
        '''
        
        path = []
        position = -1
        currentNode = self.root
        foundValue = False
        newNode = False
        
        while currentNode != None and not foundValue:
            #first add this to the path
            path.append(currentNode)
            position += 1
            
            if currentNode.value == value:
                #increment or decrement appropriately
                if isStartOfRead:
                    currentNode.incrementDelta()
                else:
                    currentNode.decrementDelta()
                    
                foundValue = True
            
            else:
                if currentNode.value > value:
                    #if there's no node, make it
                    if currentNode.left == None:
                        currentNode.left = AvlNode(value)
                        newNode = True
                    
                    #modify our value appropriately since it's to the left
                    if isStartOfRead:
                        currentNode.incrementDelta()
                    else:
                        currentNode.decrementDelta()
                        
                    #move to the left and loop back up
                    currentNode = currentNode.left
                
                else:
                    #if there's no node, make it
                    if currentNode.right == None:
                        currentNode.right = AvlNode(value)
                        newNode = True
                        
                    #shift right and loop back up
                    currentNode = currentNode.right
        
        #now check if we need to do restructuring from a new node
        if newNode:
            #decrement the position once
            position -= 1
            
            while position >= 0:
                #get the node and its child heights
                currentNode = path[position]
                lh = currentNode.getLeftHeight()
                rh = currentNode.getRightHeight()
                
                #at this point we have both heights, so we should check for a rebalance
                diff = lh - rh
                if diff > 1:
                    if currentNode.left.getRightHeight() > currentNode.left.getLeftHeight():
                        #perform a minor rotation first
                        lTree = currentNode.left
                        currentNode.left = lTree.right
                        
                        #reassign edges
                        lTree.right = currentNode.left.left
                        currentNode.left.left = lTree
                        
                        #reweight the edges
                        currentNode.left.left.evaluateHeight()
                        currentNode.left.evaluateHeight()
                        
                        #modify the delta values since we shifted
                        currentNode.left.incrementDelta(currentNode.left.left.delta)
                        
                    #now do the larger rotation
                    if position == 0:
                        #replace the root
                        self.root = currentNode.left
                        parent = None
                    else:
                        #replace a point in our path
                        parent = path[position-1]
                        if parent.right == currentNode:
                            parent.right = currentNode.left
                        else:
                            parent.left = currentNode.left
                        
                    #get the new top and modify the old delta by removing the counts
                    newTop = currentNode.left
                    currentNode.decrementDelta(newTop.delta)
                    
                    #get the left branch and move it over to the old top
                    branch = currentNode.left.right
                    currentNode.left.right = currentNode
                    currentNode.left = branch
                    
                    #update heights
                    currentNode.evaluateHeight()
                    newTop.evaluateHeight()
                
                elif diff < -1:
                    #first, check if we need to minor rotate first
                    if currentNode.right.getLeftHeight() > currentNode.right.getRightHeight():
                        #perform a minor rotation
                        rTree = currentNode.right
                        currentNode.right = rTree.left
                        
                        #reassign edges
                        rTree.left = currentNode.right.right
                        currentNode.right.right = rTree
                        
                        #reweight the heights
                        currentNode.right.right.evaluateHeight()
                        currentNode.right.evaluateHeight()
                        
                        #modify the delta value since it's no longer represented
                        rTree.decrementDelta(currentNode.right.delta)
                        
                    #now do the larger rotation
                    #parent needs to point to new child
                    if position == 0:
                        #replace the root
                        self.root = currentNode.right
                        parent = None
                    else:
                        #replace a point in our path
                        parent = path[position-1]
                        if parent.right == currentNode:
                            parent.right = currentNode.right
                        else:
                            parent.left = currentNode.right
                    
                    #get the new top and update its delta by adding its parent's delta in
                    newTop = currentNode.right
                    newTop.incrementDelta(currentNode.delta)
                    
                    #get the left branch and move it over to the old top
                    branch = currentNode.right.left
                    currentNode.right.left = currentNode
                    currentNode.right = branch
                    
                    #update heights
                    currentNode.evaluateHeight()
                    newTop.evaluateHeight()
                
                else:
                    #relatively balanced, double check heights
                    currentNode.evaluateHeight()
                    pass
                
                
                position -= 1
            
        #self.treeLock.release()
        return True
                
    def inorder(self):
        if self.root == None:
            return None
        else:
            return self.root.inorder()



def testProcess(isMaster, serverAddress):
    print serverAddress
    
    class LocalTreeManager(SyncManager):
        pass
    
    LocalTreeManager.register('updateTree')
    LocalTreeManager.register('inorder')
    LocalTreeManager.register('getTreeKeys')
    
    ltm = LocalTreeManager(serverAddress)
    ltm.connect()
    
    ltm.updateTree(0, 2, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 5, False)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 3, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 7, False)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 1, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 3, False)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 6, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 9, False)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 6, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 11, False)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 4, True)
    print ltm.inorder(0)
    
    ltm.updateTree(0, 10, False)
    print ltm.inorder(0)
    
    keys = ltm.getTreeKeys()
    print keys
    print keys[0]
    for k in keys:
        print k
    
    '''
    myDictLock.acquire()
    if myTreeDict.has_key(0):
        pass
    else:
        myTreeDict[0] = manager.Tree()
    myDictLock.release()
    
    myDictLock.acquire()
    myTree = myTreeDict[0]
    myTree.addStart(2)
    myTreeDict[0] = myTree
    print myTreeDict[0].inorder()
    myDictLock.release()
    
    myDictLock.acquire()
    myTree = myTreeDict[0]
    myTree.addEnd(5)
    myTreeDict[0] = myTree
    print myTreeDict[0].inorder()
    myDictLock.release()
    '''
    '''
    myTreeDict[0].addStart(3)
    print myTreeDict[0].inorder()
    myTreeDict[0].addEnd(7)
    print myTreeDict[0].inorder()
    
    myTreeDict[0].addStart(1)
    print myTreeDict[0].inorder()
    myTreeDict[0].addEnd(3)
    print myTreeDict[0].inorder()
    
    myTreeDict[0].addStart(6)
    print myTreeDict[0].inorder()
    myTreeDict[0].addEnd(9)
    print myTreeDict[0].inorder()
    
    myTreeDict[0].addStart(6)
    print myTreeDict[0].inorder()
    myTreeDict[0].addEnd(11)
    print myTreeDict[0].inorder()
    
    myTreeDict[0].addStart(4)
    print myTreeDict[0].inorder()
    myTreeDict[0].addEnd(10)
    print myTreeDict[0].inorder()
    '''
    
myTreeDict = {}
def updateTree(chrom, position, isStart):
    if not myTreeDict.has_key(chrom):
        myTreeDict[chrom] = AvlTree()
        
    if isStart:
        myTreeDict[chrom].addStart(position)
    else:
        myTreeDict[chrom].addEnd(position)
        
    return True
    
def inorder(chrom):
    return myTreeDict[chrom].inorder()
    
def getPH(chrom, pos):
    return myTreeDict[chrom].getPileupHeight(pos)

def getTotalPH(chrom, start, end):
    if myTreeDict.has_key(chrom):
        return myTreeDict[chrom].getTotalPileupHeight(start, end)
    else:
        return 0
def getTreeKeys():
    return myTreeDict.keys()

stageCounter = {}
def numPastStage(stage):
    if stageCounter.has_key(stage):
        return stageCounter[stage]
    else:
        return 0
    
def finishedStage(stage):
    if stageCounter.has_key(stage):
        stageCounter[stage] += 1
    else:
        stageCounter[stage] = 1
    
if __name__ == '__main__':
    class TreeManager(SyncManager):
        pass

    #TreeManager.register('Tree', AvlTree)
    #TreeManager.register('AvlNode', AvlNode)
    TreeManager.register('updateTree', updateTree)
    TreeManager.register('inorder', inorder)
    TreeManager.register('getPileupHeight', getPH)
    TreeManager.register('getTreeKeys', getTreeKeys, ListProxy)
    
    manager = TreeManager()
    manager.start()
    print manager.address
    
    #myTreeDict = manager.dict()
    #myDictLock = manager.Lock()
    
    p = Process(target=testProcess, args=(True, manager.address))
    p2 = Process(target=testProcess, args=(False, manager.address))
    
    p.start()
    p2.start()
    p.join()
    p2.join()
    
    '''
    myTestTree = AvlTree()
    
    myTestTree.addStart(2)
    print myTestTree.inorder()
    myTestTree.addEnd(5)
    print myTestTree.inorder()
    
    myTestTree.addStart(3)
    print myTestTree.inorder()
    myTestTree.addEnd(7)
    print myTestTree.inorder()
    
    myTestTree.addStart(1)
    print myTestTree.inorder()
    myTestTree.addEnd(3)
    print myTestTree.inorder()
    
    myTestTree.addStart(6)
    print myTestTree.inorder()
    myTestTree.addEnd(9)
    print myTestTree.inorder()
    
    myTestTree.addStart(6)
    print myTestTree.inorder()
    myTestTree.addEnd(11)
    print myTestTree.inorder()
    
    myTestTree.addStart(4)
    print myTestTree.inorder()
    myTestTree.addEnd(10)
    print myTestTree.inorder()
    '''
    #myTestTree = myTreeDict[0]
    
    print 'Pileup counts:'
    for x in range(0, 12):
        print str(x)+': '+str(manager.getPileupHeight(0, x))#str(myTestTree.getPileupHeight(x))
    
    '''
    print
    print 'Total pileups in range:'
    for x in range(3, 8):
        print '4-'+str(x)+': '+str(myTestTree.getTotalPileupHeight(4, x))
    '''