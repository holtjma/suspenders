'''
Created on Feb 27, 2013

@author: James Holt
'''

stageCounter = {}
def numPastStage(stage):
    '''
    Shared function, returns the number of processes that have reported completion of the given stage
    @param stage - stage value
    '''
    if stageCounter.has_key(stage):
        return stageCounter[stage]
    else:
        return 0
    
def finishedStage(stage):
    '''
    Shared function, increment the number of processes that have reported completion of a stage
    @param stage - stage value
    '''
    if stageCounter.has_key(stage):
        stageCounter[stage] += 1
    else:
        stageCounter[stage] = 1