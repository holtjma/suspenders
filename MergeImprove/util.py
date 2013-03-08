'''
Created on Jan 11, 2013

@author: Shunping Huang

These functions were copied from Shunping Huang's Lapels
Just used for argparse to check that a file exists and/or is read/writable
'''

import os
import argparse as ap

def readableFile(fileName):
    if os.path.isfile(fileName) and os.access(fileName, os.R_OK):
        return fileName
    else:
        raise ap.ArgumentTypeError("Cannot read file '%s'." % fileName)


def writableFile(fileName):
    if os.access(os.path.dirname(fileName), os.W_OK):
        return fileName
    else:        
        raise ap.ArgumentTypeError("Cannot write file '%s'." % fileName)
