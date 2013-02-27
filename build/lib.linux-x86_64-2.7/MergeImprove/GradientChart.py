'''
Created on Oct 30, 2012

@author: holtjma
'''

import Image, ImageDraw
import math

class GradientChart(object):
    '''
    classdocs
    '''


    def __init__(self, pairing, data):
        '''
        Constructor
        '''
        
        self.data = data
        
        self.xMax = 1000
        self.yMax = 1000
        
        self.labels = pairing
        
        self.image = Image.new('RGBA', (self.xMax, self.yMax), "#ffffff")
        
    def drawChart(self):
        self.xAxisYLoc = int(self.image.size[1]*.95)
        self.yAxisXLoc = int(self.image.size[0]*.05)
        
        #this is for the data gathering output
        self.maxRow = max(self.data)
        self.maxColumn = 0
        
        self.logBase = 6
        
        #find the max
        for row in self.data:
            maxForRow = max(self.data[row])
            if maxForRow > self.maxColumn:
                self.maxColumn = maxForRow
        
        #set the piece width and height
        self.pieceWidth = int((self.xMax-(self.yAxisXLoc+1))/(self.maxRow+1))
        self.pieceHeight = int(self.xAxisYLoc-1)/(self.maxColumn+1)
        
        
        incX = int(self.maxRow / 4)
        if(incX == 0):
            incX = 1
        numIncX = int(self.maxRow/incX)
        
        incY = int(self.maxColumn / 4)
        if(incY == 0):
            incY = 1
        numIncY = int(self.maxColumn/incY)
        
        #data time
        self.intensity = {}
        draw = ImageDraw.Draw(self.image)
        
        #output the actual pairing counts
        for row in range(0, self.maxRow+2):
            self.intensity[row] = {}
            for column in range(0, self.maxColumn+2):
                if self.data.has_key(row) and self.data[row].has_key(column):
                    value = self.data[row][column]
                else:
                    value = 0
                
                #now we do something with the value
                if value == 0:
                    self.intensity[row][column] = value
                else:
                    if value == 1:
                        value = 2 #make sure it shows
                    self.intensity[row][column] = math.log(value, self.logBase)
                    
                if(row > 0 and column > 0):
                    self.drawArea(row-1, column-1, draw)
                    
        
        #draw some axes
        draw.line([self.yAxisXLoc, self.xAxisYLoc, self.xMax, self.xAxisYLoc], "#000000")
        draw.line([self.yAxisXLoc, 0, self.yAxisXLoc, self.xAxisYLoc], "#000000")
        
        for x in range(0, numIncX+1):
            marker = x*incX
            tcoordinates = [self.yAxisXLoc+(marker+.5)*self.pieceWidth-2, self.xAxisYLoc]
            draw.text(tcoordinates, str(marker), "#000000")
            
            lineCoor = [self.yAxisXLoc+marker*self.pieceWidth, 0, self.yAxisXLoc+marker*self.pieceWidth, self.xAxisYLoc]
            draw.line(lineCoor, "#000000")
        
        for y in range(0, numIncY+1):
            marker = y*incY
            tcoordinates = [self.yAxisXLoc-10, self.xAxisYLoc-(marker+.5)*self.pieceHeight-5]
            draw.text(tcoordinates, str(marker), "#000000")
            
            lineCoor = [self.yAxisXLoc, self.xAxisYLoc-marker*self.pieceHeight, self.xMax, self.xAxisYLoc-marker*self.pieceHeight]
            draw.line(lineCoor, "#000000")
        
        (labelWidth, labelHeight) = draw.textsize(self.labels[0])
        labelCoordinates = [self.yAxisXLoc + (self.xMax-self.yAxisXLoc)/2 - labelWidth/2, self.xAxisYLoc+(self.yMax-self.xAxisYLoc)/2]
        draw.text(labelCoordinates, self.labels[0], "#000000")
    
        #have to create a new image, rotate 90 degrees and paste
        temp = Image.new('RGB', (self.yMax, self.yAxisXLoc/2), "#ffffff")
        tempDraw = ImageDraw.Draw(temp)
        (tempX, tempY) = tempDraw.textsize(self.labels[1])
        tempCoor = (self.yMax - self.xAxisYLoc/2 - tempX/2, 0)
        tempDraw.text(tempCoor, self.labels[1], "#000000")
        temp = temp.rotate(90)
        self.image.paste(temp, (0, 0))
        del tempDraw
    
        del draw
        
        self.image.show()
        
    def drawArea(self, x, y, draw):
        color = 16
        i = int(255 - self.intensity[x][y]*color)
        
        if self.data.has_key(x) and self.data[x].has_key(y):
            value = self.data[x][y]
        else:
            value = 0
        
        
        
        startXPixel = int(self.pieceWidth*x + (self.yAxisXLoc + 1))
        startYPixel = int((self.xAxisYLoc - 1) - self.pieceHeight*y)
        endXPixel = int(self.pieceWidth*(x+1) + (self.yAxisXLoc + 1))
        endYPixel = int((self.xAxisYLoc - 1) - self.pieceHeight*(y+1))
        
        draw.rectangle([startXPixel, startYPixel, endXPixel, endYPixel], (i, i, 255))
        
        #only show text for non-zero regions
        if value != 0:
            suffixes = ['K', 'M', 'B']    
            suffixIndex = 0
            
            displayValue = str(value)
            while len(displayValue) > 4 and suffixIndex < 3:
                value = int(value/1000)
                displayValue = str(value)+suffixes[suffixIndex]
                suffixIndex += 1
            (labelWidth, labelHeight) = draw.textsize(displayValue)
            draw.text([(startXPixel+endXPixel)/2-labelWidth/2, (startYPixel+endYPixel)/2-labelHeight/2], displayValue, (0, 0, 0))
        
        
    def drawGradientArea(self, x, y):
        
        color = 32
        
        i00 = 255 - self.intensity[x][y]*color
        i10 = 255 - self.intensity[x+1][y]*color
        i01 = 255 - self.intensity[x][y+1]*color
        i11 = 255 - self.intensity[x+1][y+1]*color
        
        startXPixel = self.pieceWidth*x + (self.yAxisXLoc + 1)
        startYPixel = (self.xAxisYLoc - 1) - self.pieceHeight*y
        endXPixel = self.pieceWidth*(x+1) + (self.yAxisXLoc + 1)
        endYPixel = (self.xAxisYLoc - 1) - self.pieceHeight*(y+1)
        
        pixels = self.image.load()
        for xLoc in range(startXPixel, endXPixel):
            for yLoc in range(startYPixel, endYPixel, -1):
                #time to do some calculations
                ax = endXPixel - startXPixel
                ay = startYPixel - endYPixel
                p1 = int((i00*float(yLoc-endYPixel)/ay + i01*(1-float(yLoc-endYPixel)/ay))*(1-float(xLoc-startXPixel)/ax))
                p2 = int((i10*float(yLoc-endYPixel)/ay + i11*(1-float(yLoc-endYPixel)/ay))*float(xLoc-startXPixel)/ax)
                pc = p1+p2
                pixels[xLoc, yLoc] = (pc, pc, 255)
                
            
            