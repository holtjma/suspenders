'''
Created on Feb 27, 2013

@author: holtjma
'''

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