# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:36:13 2022

@author: antoine


TODO:
    - change the reference to correct image drift
    - keep stars in order of brightness in plotDif rescaling.
    - add an option to let users to fine tune tresholling:                Done
    - improve asteroid detection
    - change method to compute fwhm
    - add delete star passages algorythm to photometry
    - add methode to creat a log file
    - center time if needed
    
PROBLEMS:
    - for hist upper to 32000... tresholing problem with yGauss and x
    - pltFid problem when binning != 1

"""

import numpy as np
import os



import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arrow
from matplotlib.backend_bases import MouseButton
import matplotlib.animation as animation


from utils import *
from exceptions import *
from detector import Detector
from occultation import *
from data_structurs import Apertures
from singleton import Singleton

from astropy.time import Time

import glob

class Photometry:
    def __init__(self, detector = None):
        self.detector = detector
        self.Appertur = []
        self.stars = []
        self.results = []
        self.fwhm = []
        self.nbOfAst = 0
        self.nbOfStar = 0
        self.singleton = Singleton()
        self.starPassages = None
  
        
        
    def FWHM(self, eps = 10):
        
        for i in range(len(self.detector.seqManager)):
            img = self.detector.getData(i)
          
            s = self.detector.correctStars(self.stars, i, 1)
            s = Utils.centred(s , img, eps)
            
            
            fwhmOfOneImage = np.zeros((len(s), 2))
            
            for j in range(len(s)):
                imstar = img[int(s[j,1] - eps) : int(s[j,1] + eps), int(s[j,0] - eps) : int(s[j,0] + eps)]
                imstar = imstar - np.nanmedian(imstar)
                # Utils.imshow(imstar)
                mmax = np.max(imstar)
                idxOfMax = np.asarray(np.where(imstar == mmax))
                idxOfTopOfFWHM = np.asarray(np.where(imstar >= mmax / 2.0))
                
                idxOfSameX = idxOfTopOfFWHM.T[idxOfTopOfFWHM.T[:,1] == idxOfMax[1,0]]
                idxOfSameY = idxOfTopOfFWHM.T[idxOfTopOfFWHM.T[:,0] == idxOfMax[0,0]]
                
                # print('i',i , 'j',j, 'x', idxOfSameY , 'y', idxOfSameX)
                
                xFwhm = idxOfSameY.shape[0]
                yFwhm = idxOfSameX.shape[0]
                
                fwhmOfOneImage[j,0] = xFwhm
                fwhmOfOneImage[j,1] = yFwhm
                
            self.fwhm.append(np.mean(fwhmOfOneImage, axis = 0))
                
    def getFwhm(self, idx = 0):
        return np.max(self.fwhm[idx])
            
    def isAstToFast(self, ast, speed):
        astCorrected = np.asarray(ast) - (self.detector.getImg(-1).getTime(self.detector.key, self.detector.format) - self.detector.getImg().getTime(self.detector.key, self.detector.format)).to_value('sec')*np.asarray(speed)
        return (astCorrected < 0).any() or (astCorrected > self.detector.getImgShape(0)).any()
    
    def isOutOfBoundaries(self, pos, re):
        extrema = self.detector.findDriftExtrema()
        return (pos + extrema + re > self.detector.shape).any() or (pos - extrema - re < np.zeros((2))).any()
    
    def isOutOfCircle(self, pos, center, delt = 0):

        pos = pos - center
        
        minn = np.min(center) - delt
        return (Utils.dist(pos)+4 > minn).any()
    
    def isStarsToClose(self, s, re):
        for i in self.detector.getStarsListOfImg(0):
            for j in i:
                if Utils.dist(s, j) < re and Utils.dist(s, j) != 0:
                    return True

    def findGoodStars(self, nbOfStars, maxVal, re):
        img = self.detector.getData()
        bx = 5

        for i, s in enumerate(self.detector.getStarsListOfImg(0)):
           
            if self.isOutOfBoundaries(s, re) or self.isOutOfCircle(s, self.detector.getImgCenter(), re):
                continue
            
            if np.max(img[s[1] - bx:s[1] + bx, s[0] - bx:s[0] + bx]) > maxVal:
                continue

            if self.isStarsToClose(s, re):
                continue
            
            self.stars.append(s)
            
            if len(self.stars) == nbOfStars:
                break
            
        if len(self.stars) < nbOfStars:
            print("Warning: not enougth stars. requiere:", nbOfStars, "find: ", len(self.stars))
    
    def showAp(self, idx):
        
        img = self.detector.getReducedData(idx)
        r, ri, re = self.getAppSize(idx)
        
        fig,ax = plt.subplots(1)
        
        aps = self.Appertur[idx].positions
        
        for i,pos in enumerate(aps):
            col = 'r'
            if i < len(self.detector.asteroidsPositionFirstImage):
                col = 'black'
            c = Circle(pos, radius = r, fill=False, color = col)
            c1 = Circle(pos, radius = ri, fill=False, color = col)
            c2 = Circle(pos, radius = re, fill=False, color = col)
            ax.add_patch(c)
            ax.add_patch(c1)
            ax.add_patch(c2)
            
        ax.imshow(img, cmap="Greys", vmin = np.mean(img[img>0])*0.5, vmax = np.mean(img[img>0])/0.5)
        # plt.savefig(str(idx) + '.png')
        
    def selectOnManuel(self, arr, im_idx = 0, r = 5):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        isClose = [False]
        img = self.detector.getData(im_idx)
        
        ax.imshow(img, cmap='Greys', vmin = np.median(img)*0.5, vmax = np.median(img)/0.5)
        
        if len(self.detector.asteroidsPositionFirstImage) != 0 :
            for i in self.detector.asteroidsPositionFirstImage:
                
                c = Circle(i, radius = 50, fill=False)
                ax.add_patch(c)
        
        plt.connect('button_press_event', lambda event: Utils.getPosFromImage(event, arr, isClose))
        plt.show(block=False)
            
        
        while not isClose[0]:
            plt.pause(0.01)
        
        if len(arr) != 0 :
            img1 = img[int(arr[0][1]) - r: int(arr[0][1]) + r, int(arr[0][0]) - r: int(arr[0][0]) + r]
            arr = Utils.centred(arr, img, r)
            print("centré", arr, 'image', im_idx)
            
        
    def resetLists(self):
        self.stars = []
        self.Appertur = []
        self.results = []
        self.fwhm = []
        
    def getAppSize(self, idx):
        r = self.getFwhm(idx)*2.1
        ri = 1.5*r
        re = 2*r
        
        return r, ri, re
        
    
    def startByTracking(self, r, ri, re, eps = 10):
        print('select stars:')
        self.selectOnManuel(self.stars)
        print('select asteroid:')
        firstPos = []


        self.selectOnManuel(firstPos, 0)
        
        for i in range(len(self.detector)):
            img = self.detector.getData(i)
            s = Utils.centred(firstPos+ self.stars , img, eps)
            
            
            self.Appertur.append(Appertur(s, len(firstPos), r, ri, re))
            self.results.append(self.Appertur[-1].Photom(self.detector.getImg(i), self.detector.key, self.detector.format, self.detector.isTimeCentred, self.detector.seqManager.getExpo(i, self.detector.exposurKey)))
            
            print("image: ", i, '\n', self.Appertur[-1].Photom(self.detector.getImg(i), self.detector.key, self.detector.format), '\n', type(self.Appertur[-1].Photom(self.detector.getImg(i), self.detector.key, self.detector.format)))
        
    def start(self, nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000):
        
        self.resetLists()
        self.singleton.clearDataStarsPassages()
        
        if self.detector == None:
            raise NoneDectector()
       
        ismanual = str(input("MANUAL STARS SELECTION? (YES/NO)")).replace(' ','')
        isManualSelection = ismanual.lower() == 'yes' or ismanual.lower() == 'y'
        
        if not isManualSelection:
            self.findGoodStars(nbOfStars, maxVal, 20)    
        else:
            self.selectOnManuel(self.stars)
            
            
        self.FWHM()
        ismanual = str(input("MANUAL ASTEROIDS SELECTION? (YES/NO)")).replace(' ','')
        isManualSelection = ismanual.lower() == 'yes' or ismanual.lower() == 'y'
        
        if isManualSelection:
            
            print("select asteroid on the first image")
            firstPos = []
            lastPos = []
            idxF, idxE = self.detector.findBestIdx()
            print("first image:", self.detector.seqManager.seq[idxF])
            print("last image:", self.detector.seqManager.seq[idxE])
            self.selectOnManuel(firstPos, idxF)
            firstPos = self.detector.correctStarsFromRot(firstPos, idxF, -1) - self.detector.avgDrif(idxF)
            print("first", firstPos)
            print("select asteroid on the second image IN THE SAME ORDER than on the first image")
            self.selectOnManuel(lastPos, idxE)
            lastPos = self.detector.correctStarsFromRot(lastPos, idxE, -1) - self.detector.avgDrif(idxE)
            
            
            time1 = self.detector.getImg(idxF).getTime(self.detector.key, self.detector.format)
            time2 = self.detector.getImg(idxE).getTime(self.detector.key, self.detector.format)
            
            self.detector.asteroidsSpeed = -((lastPos-firstPos)/(time2 - time1).to_value('sec'))
            print('speed', self.detector.asteroidsSpeed, self.detector.asteroidsSpeed * (self.detector.getImg(0).getTime(self.detector.key,  self.detector.format) - time1).to_value('sec')  )
            self.detector.asteroidsPositionFirstImage = firstPos - self.detector.asteroidsSpeed * (self.detector.getImg(0).getTime(self.detector.key,  self.detector.format) - time1).to_value('sec')  
            

            # img2 = img[int(self.stars[0][1]) - r: int(self.stars[0][1]) + r, int(self.stars[0][0]) - r: int(self.stars[0][0]) + r]

        self.nbOfStar = len(self.stars)
        self.nbOfAst = len(self.detector.asteroidsPositionFirstImage)
        
        self.starPassages = np.zeros((self.nbOfAst, len(self.detector)))
        
        if not isManualSelection:
            self.eliminateUnconsistanteAst(np.max(self.fwhm)*2)
        
        correction = np.zeros((2))
        
        listOfStarsPassages = np.concatenate(self.StarPassages(starPassageOfs))
        
        for i in range(len(self.detector)):
            
            r, ri, re = self.getAppSize(i)
         
            apertures = np.concatenate((self.detector.getAstPositionAtImg(i), np.asarray(self.stars))) + self.detector.avgDrif(i) + correction
            apertures = self.detector.correctStarsFromRot(apertures, i, 1)#correct drift and rot
           
            
            if center :
                ap = Utils.centred(apertures , self.detector.getData(i), r)
                
            else:
                ap = apertures
            
            
            correction = apertures - ap
            
            starPassagesCorrection = self.detector.correctStars(listOfStarsPassages, i, 1)
           
            for ast in range(self.nbOfAst):
                for star in starPassagesCorrection:
                    if Utils.dist(star, ap[ast, :]) < 3*self.getFwhm(i):
                        self.singleton.appDataStarsPassages(self.detector.seqManager.getFileName(i) + '\n')
                        self.starPassages[ast, i] = True
                        break
            
            self.Appertur.append(Apertures(ap, len(self.detector.asteroidsPositionFirstImage), r, ri, re))
            self.results.append(self.Appertur[-1].photom(self.detector.getImg(i), self.detector.key, self.detector.format, self.detector.isTimeCentred, self.detector.seqManager.getExpo(i, self.detector.exposurKey)))
            
            
            print("image: ", i, '\n', self.Appertur[-1].photom(self.detector.getImg(i), self.detector.key, self.detector.format), '\n', type(self.Appertur[-1].photom(self.detector.getImg(i), self.detector.key, self.detector.format)))
            

            
            if i == 0 or i == (len(self.detector) - 1) / 2 or i == len(self.detector) - 1: 
                # pass   
            
                self.showAp(i)
               # for j in ap:
               #     Utils.imshowap(self.detector.getData(i), r, j)
                        
    def eliminateUnconsistanteAst(self, re):
       
        cpt = 0
    
        while cpt != len(self.detector.asteroidsPositionFirstImage):
            ast = self.detector.asteroidsPositionFirstImage[cpt]
           
            if self.isOutOfBoundaries(ast, re) or self.isOutOfCircle(ast - (self.detector.getImg(-1).getTime(self.detector.key,  self.detector.format) - self.detector.getImg().getTime(self.detector.key,  self.detector.format)).to_value('sec')*np.asarray(self.detector.asteroidsSpeed[cpt]), self.detector.getImgCenter(), re) or self.isAstToFast(ast, self.detector.asteroidsSpeed[cpt]):
                self.detector.asteroidsPositionFirstImage.pop(cpt)
                self.detector.asteroidsSpeed.pop(cpt)
            else:
                cpt += 1
    
    def getAstPhot(self, binning, inMag=True, withoutStarPassages = True):
        
        data = np.asarray(self.astPhot())
   
        if withoutStarPassages:
            data[self.starPassages.T == True] = np.nan
        
        data = np.asarray(Utils.binn(binning, data))
            
        if inMag:
            print("result ast in Mag")
            return -2.5*np.log10(data)
        print("result ast in ADU")
        return data
    
    def getStarPhot(self, binning, inMag=True):
        
        data = np.asarray(Utils.binn(binning, self.starsPhot()))
        
        if inMag:
            print("result star in Mag")
            return -2.5*np.log10(data)
        print("result star in ADU")
        return data
    
    def getTime(self, binning, forma= "jd"):
         t = Utils.binn(binning, self.extractTime())
        
         precision = len(str(t[0]).split('.')[-1])
         if precision > 9:
             precision = 9
         
         t = Time(t, format = "jd", precision= precision)
         return self.timeFormat(t, forma)
    
    
      
    
    def plotDif(self, refS = 0, ast = -1, yRange = None, binning = 1, resc = True, forma = 'jd', xtick = None, inMag = True, rmExtremPoint = False, cStd = 2, deg = 4, displayRmFit = False, starPassage = False, markerSize = 100, lineWidths = 5):
        
        
        if binning == -1: # Auto binning
        
            t = self.getTime(1)
            
            dth = (t[-1] - t[0])*24
            NdthIn12H = 12 / dth
            
            nbOfPoint = NdthIn12H * t.shape[0]
            
            if (nbOfPoint <= 240 and nbOfPoint >= 100) or nbOfPoint < 100 :
                binning = 1
            
            else:
                binning = int(nbOfPoint / 240 + 0.5)
                
            print("binning: ", binning)
        
        idxOfExtremPoints = self.idxOfDeletedPoint(refS, binning, d = deg, cStd = cStd, plot = displayRmFit) 
        
        asts = self.getAstPhot(binning, inMag)
        stars =self.getStarPhot(binning, inMag)

        
        x = self.getTime(binning, forma)
        
        if ast == -1:
            plt.figure()
            for i in range(asts.shape[1]):
                plt.scatter(x, asts[:,i] - stars[:, refS], label='ast - s' + str(refS+1), linewidths=lineWidths, marker='x', s=markerSize)
                
                if rmExtremPoint:
                    plt.scatter(x[idxOfExtremPoints[i]], (asts[:,i] - stars[:, refS])[idxOfExtremPoints[i]], linewidths=lineWidths, color = 'r')
        else:
            plt.figure()
            plt.scatter(x, asts[:,ast] - stars[:, refS], label='ast - s' + str(refS+1), linewidths=lineWidths, marker='x', s = markerSize)
            
            if rmExtremPoint:
                plt.scatter(x[idxOfExtremPoints[ast]], (asts[:,ast] - stars[:, refS])[idxOfExtremPoints[ast]], linewidths=lineWidths, color = 'r')
        
        of = np.zeros(stars.shape[1])
        
        if resc:
            of = Utils.rescalMulti((asts.T  - stars[:, refS]).T, (stars.T - stars[:, refS]).T)
        
        for j in range(stars.shape[1]):
            if j != refS:
                if ast != -1:
                    plt.scatter(x, (stars[:,j] - stars[:, refS]) + of[j, ast], label='s' + str(j+1) +' - s' + str(refS+1), linewidths=lineWidths, marker='x', s = markerSize)
                else:
                    plt.scatter(x, (stars[:,j] - stars[:, refS]) + np.min(of[j]), label='s' + str(j+1) +' - s' + str(refS+1), linewidths=lineWidths, marker='x', s = markerSize)
                
        plt.legend()        
        if yRange != None:
            plt.ylim(yRange)

        if xtick != None:
            plt.xticks(np.linspace(plt.xlim()[0], plt.xlim()[-1], xtick))
        
        if inMag:
            plt.gca().invert_yaxis()
            
        plt.show()
        
    def plot(self, yRange = None, binning = 1, selection = None, inMag = True, forma = 'jd', xtick = None, markerSize = 100):
        
        if inMag:
            print("result in Mag")
            asts = -2.5*np.log10(Utils.binn(binning, self.astPhot()))
            stars = -2.5*np.log10(Utils.binn(binning, self.starsPhot()))
        else:
            print("result in ADU")
            asts = Utils.binn(binning, self.astPhot())
            stars = Utils.binn(binning, self.starsPhot())
        
        print(stars.shape)
        
        x = self.toTime(binning)
        x = self.timeFormat(x, forma)
        plt.figure()
        if selection == None:
            for i in range(asts.shape[1]):
                plt.scatter(x, asts[:, i], label='ast ' + str(i+1), linewidths=1, marker='x', s = markerSize)
            for j in range(stars.shape[1]):
                plt.scatter(x, stars[:, j], label='s' + str(j+1), linewidths=1, marker='x', s = markerSize)
        else:
            for i in range(len(selection)):
                if selection[i] < asts.shape[1]:
                    plt.scatter(x, asts[:, i], label='ast ' + str(i+1), linewidths=1, marker='x', s = markerSize)
                else:
                    j = i - asts.shape[1]
                    plt.scatter(x, stars[:, j], label='s' + str(j+1), linewidths=1, marker='x', s = markerSize)
                
                
        plt.legend()   
        if yRange != None:
           plt.ylim(yRange)
           
        if xtick != None:
           plt.xticks(np.linspace(plt.xlim()[0], plt.xlim()[-1], xtick))
        
        if inMag:
            plt.gca().invert_yaxis()
            
            
    def plotMeanDif(self, yRange = None, binning = 1, forma = 'jd', xtick = None, inMag = True):
        
        if inMag:
            print("result in Mag")
            asts = -2.5*np.log10(Utils.binn(binning, self.astPhot()))
            stars = -2.5*np.log10(Utils.binn(binning, self.starsPhot()))
        else:
            print("result in ADU")
            asts = Utils.binn(binning, self.astPhot())
            stars = Utils.binn(binning, self.starsPhot())
        
        x = self.toTime(binning)
        x = self.timeFormat(x, forma)
        
        for i in range(asts.shape[1]):
            plt.figure()
            plt.scatter(x, asts[:,i] - stars[:, 0], label='ast - s1', linewidths=1, marker='x')
           
        plt.legend()

        if yRange != None:
           plt.ylim(yRange)
           
        if xtick != None:
            plt.xticks(np.linspace(plt.xlim()[0], plt.xlim()[-1], xtick))
            
        if inMag:
            plt.gca().invert_yaxis()
        
    def extractTime(self):
        x = []
        for i in range(len(self.results)):
            x.append(float(self.results[i][0][4]))
        return np.asarray(x)
    
    # def extractTimeb(self):
    #     x = []
    #     for i in range(len(self.detector)):
    #         x.append(self.detector.getImg(i).getTime(self.detector.key,  self.detector.format).to_value('jd'))
    #     return np.asarray(x)
    
    def starsPhot(self):
        
        stars = []
        for res in self.results:
            star = []
            for i in range(self.nbOfStar):
                star.append(res[self.nbOfAst + i][7])
            stars.append(star)
        return np.asarray(stars)
                
        
    def astPhot(self):
        
        asts = []
        
        for res in self.results:
            ast = []
            for i in range(self.nbOfAst):
                ast.append(res[i][7])
               
            asts.append(ast)
        return np.asarray(asts)
    
    def toTime(self, binning):
        t = Utils.binn(binning, self.extractTime())
        
        precision = len(str(t[0]).split('.')[-1])
        
        if precision > 9:
            precision = 9
        
        return Time(t, format = 'jd', precision = precision)
        
    def timeFormat(self, tArr, forma):
        newArr = []
        #for i in tArr:
        newArr = tArr.to_value(forma)
        return newArr
    
    def toCsv(self, path):
    
        header = self.buildCsvHeader()
        
        csv = []
        for im, table in enumerate(self.results):
            csvRow = np.zeros(( len(table) * len(table[0]) + self.starPassages.shape[0]))
            for i, row in enumerate(table):
                if i < len(self.detector.asteroidsPositionFirstImage):
                    csvRow[i] = 1
                else:
                    csvRow[i] = 0
                
                for j in range(len(row)-1):
                    if j == 0 or j == 1:
                        csvRow[(j+1)*len(table) + i] = row[1+j].value
                    else:
                        csvRow[(j+1)*len(table) + i] = row[1+j]
            
            for i in range(self.starPassages.shape[0]):
                csvRow[len(table) * len(table[0]) + i] = self.starPassages[i, im]
            
            csv.append(csvRow)
            
        df = pd.DataFrame(csv, columns= header)
        df.to_csv(path, index=False)
            
        return df
    
    def buildCsvHeader(self):
        
        NOfAp = len(self.results[0])
        NbOfRowElements = len(self.results[0][0])
        header = [0]*(NOfAp*NbOfRowElements + self.starPassages.shape[0])
        
        for i in range(NOfAp):
            header[i] = 'isAst (' + str(i+1) + ')'
    
        for i in range(NbOfRowElements-1):
            for j in range(NOfAp):
                if i == 0:
                    header[NOfAp + i*NOfAp + j] = "xcenter (" + str(j+1) + ')'
                elif i == 1:
                    header[NOfAp + i*NOfAp + j] = "ycenter (" + str(j+1) + ')'
                elif i == 2:
                    header[NOfAp + i*NOfAp + j] = "aperture_sum (" + str(j+1) + ')'
                elif i == 3:
                    header[NOfAp + i*NOfAp + j] = "Time (" + str(j+1) + ')'
                elif i == 4:
                    header[NOfAp + i*NOfAp + j] = "annulus_median (" + str(j+1) + ')'
                elif i == 5:
                    header[NOfAp + i*NOfAp + j] = "aper_bkg (" + str(j+1) + ')'
                elif i == 6:
                    header[NOfAp + i*NOfAp + j] = "aper_sum_bkgsub (" + str(j+1) + ')'
            
        
        for i in range(self.starPassages.shape[0]):
            header[NOfAp*NbOfRowElements + i] = "star_passages_ast (" + str(i+1) + ')'
        
        return header 
    
    def timeLimit(self):
        """
        return the time needed for the slowest asteroid to travel 1 fwhm on the frame

        Returns
        -------
        float
            DESCRIPTION. Time needed for the slowest asteroid to travel 1 fwhm

        """
        return 2*np.min(self.fwhm) / self.detector.astSpeed(self.detector.slowestAst())
    
    
    def idxOfImages(self):
        lt = self.timeLimit()
        idx = [0]
        d = self.detector
        
        for i in range(len(d)):
           print((d.seqManager.getTime(d.key, d.format, i) - d.seqManager.getTime(d.key, d.format, idx[-1])).to_value('sec'), lt)
           if (d.seqManager.getTime(d.key, d.format, i) - d.seqManager.getTime(d.key, d.format, idx[-1])).to_value('sec') >= lt:
               idx.append(i)
        
        return idx
    
    
    def astPathLinearCoef(self):
        firstPos = self.detector.getAstPositionAtImg(0)
        lastPos = self.detector.getAstPositionAtImg(-1)
        
        return Utils.linearBtwTwoPoints(firstPos, lastPos)
    
    
    def findCorner(self, a, idx_Img):
        """
        according to the first and last position of asteroids (information given by the slope a), find two points
        on the normal at a distance + 2fwhm and - 2fwhm on the image idx_Img

        Parameters
        ----------
        a : array
            Slopes of the line done by asteroids between their first and last position.
        idx_Img : int
            idx of image to produce the two point.

        Returns
        -------
        p1 : array
            coordinate (x,y) of the first point.
        p2 : array
            coordinate (x, y) of the second point.

        """
        
        c = 3.5
        cfwhm = c*self.getFwhm(idx_Img)
        firstPos = self.detector.getAstPositionAtImg(idx_Img)
        
        p1 = np.zeros(firstPos.shape)
        p2 = np.zeros(firstPos.shape)
        
        
        
        p1[:, 0] = 2*np.sqrt( a**2/(a**2 + 1) * cfwhm ) + firstPos[:, 0]
        p1[:, 1] = (-2/a)*np.sqrt( a**2/(a**2 + 1) * cfwhm ) + firstPos[:, 1]
        
        p2[:, 0] = -2*np.sqrt( a**2/(a**2 + 1) * cfwhm ) + firstPos[:, 0]
        p2[:, 1] = (2/a)*np.sqrt( a**2/(a**2 + 1) * cfwhm ) + firstPos[:, 1]
        
        return p1, p2
        
        
    
    def buildBoxAreaArrountAstPath(self, a):
        """
        according to the slope done by two point(asteroids first position and last position), 
        return corners in clockwise order of the box of thickness of 2fwhm
    
        o------------------------------------o
        |                                    |
        *------------------------------------*
        |                                    |
        o------------------------------------o
        
        o : points
        * : two point used for slope (asteroid first and last position)
        | : 2fwhm

        Parameters
        ----------
        a : array
            array of slope.

        Returns
        -------
        list
            list of points of corner of the box in clockwise order.

        """
        
        
        p1 = self.findCorner(a, 0)
        p2 = self.findCorner(a, -1)

        return [p1[0], p2[0], p2[1], p1[1]]
    
    def isStarsInBox(self, star):
        """

        Parameters
            star: array representing the point to check

        Returns
            list of boleans. True if the point is inside the rectangle. False otherwise.
        """
        
        a, b = self.astPathLinearCoef()
        vertices = np.asarray(self.buildBoxAreaArrountAstPath(a))
        
        sol = []
        
        for i in range(vertices[0].shape[0]):
            sol.append(all(Utils.dotProdWithSharedStart(vertices[j - 1, i, :], v, star) > 0 for j, v in enumerate(vertices[:,i,:])))
         
        return sol
    
    def sisa(self, s, ast, eps = 2):
        return s[0] < ast[0] + eps and s[0] > ast[0] - eps and s[1] < ast[1] + eps and s[1] > ast[1] - eps
               
    def StarPassages(self, ofs = 15000):
        
        idxF, idxE = self.detector.findBestIdx()
        
        s1 = self.detector.getImg(idxF).findStars(self.detector.getImg(idxF).getTresh() + ofs)

        s2 = self.detector.getImg(idxE).findStars(self.detector.getImg(idxE).getTresh() + ofs)
        
        s1 = self.detector.correctStars(s1, idxF)
        s2 = self.detector.correctStars(s2, idxE)
                
        astsFirstPos = self.detector.getAstPositionAtImg(idxF)
        astsLastPos = self.detector.getAstPositionAtImg(idxE)
        
        listOfStars = [[] for _ in range(self.detector.nofa())] # create an empty list of empty lists [[],[]] of len = N of asteroids. [[]]*N dont works because each list inside are reference to the same object
        
        
        for s in s1:
            isInBoxs = self.isStarsInBox(s)
            for i, box in enumerate(isInBoxs):
                if box and not self.sisa(s, astsFirstPos[i]):
                    listOfStars[i].append(s)
             
                    
        for s in s2:
            isInBoxs = self.isStarsInBox(s)
            for i, box in enumerate(isInBoxs):
                if box and not self.sisa(s, astsLastPos[i]):
                    listOfStars[i].append(s)

        return listOfStars
        
    def checkBox(self, ofs):
        
        a, b = self.astPathLinearCoef()
        vertices = np.asarray(self.buildBoxAreaArrountAstPath(a))
        idxF, idxE = self.detector.findBestIdx()
        print("first image idx: ", idxF, "last image idx: ", idxE)
        
        
        pot = self.detector.asteroidsPositionFirstImage[-1]
        astsFirstPos = self.detector.getAstPositionAtImg(idxF)
        astsLastPos = self.detector.getAstPositionAtImg(idxE)
        
        im = self.detector.getData(idxF)
        
        plt.imshow(im, vmin = np.median(im) - np.std(im), vmax = np.median(im) + np.std(im), cmap = "Greys")
        
        starsInBox = np.asarray(self.StarPassages(ofs))
        
        
        s1 = self.detector.getImg(idxF).findStars(self.detector.getImg(idxF).getTresh() + ofs)

        s2 = self.detector.getImg(idxE).findStars(self.detector.getImg(idxE).getTresh() + ofs)
        
        s1 = self.detector.correctStars(s1, idxF)
        s2 = self.detector.correctStars(s2, idxE)

        
        
        for ver in vertices:
            plt.scatter(ver[:, 0], ver[:, 1], color = 'blue')
        
        for ss in starsInBox:
            plt.scatter(ss[:,0], ss[:, 1], marker = "x", linewidths=1, color = 'red')
        
        plt.scatter([astsFirstPos[0,0], astsLastPos[0,0]], [astsFirstPos[0,1], astsLastPos[0,1]], color= "orange", s = 10 )
        plt.scatter(pot[0], pot[1], color = "green", s=10)
        
  
        
        
                
    def toGif(self, path):
        
        snapshots = [ self.detector.getReducedData(i) for i in range(len(self.detector)) ]

        fig,ax = plt.subplots(1)

        a = snapshots[0]
        im = ax.imshow(a, interpolation='none', cmap='Greys', vmin = np.median(a[a>0])*0.5, vmax = np.median(a[a>0]) / 0.5)
        
        idxF, idxL = self.detector.findBestIdx()
       
        
        def animate_func(i):
            ax.clear()
            im = ax.imshow(snapshots[i], interpolation='none', cmap='Greys', vmin = np.median(snapshots[i][snapshots[i]>0])*0.5, vmax = np.median(snapshots[i][snapshots[i]>0]) / 0.5)
            
            if len(self.Appertur) == 0:
                return [im]
            
            
            for j, pos in enumerate(self.Appertur[i].positions):
               
                col = 'black'
                
                if j >= self.Appertur[i].idxOfStars:
                    col = 'black'
                    c = Circle(pos, radius = self.Appertur[i].aperture.r, fill=False, color = col, linewidth = 0.5)
                    c1 = Circle(pos, radius = self.Appertur[i].annulus_aperture.r_in, fill=False, color = col, linewidth = 0.5)
                    c2 = Circle(pos, radius = self.Appertur[i].annulus_aperture.r_out, fill=False, color = col, linewidth = 0.5)
                    ax.add_patch(c)
                    ax.add_patch(c1)
                    ax.add_patch(c2)
                    continue
                
                
                if self.starPassages[j, i]:
                    col = "red"
                            
                c = Circle(pos, radius = self.Appertur[i].aperture.r, fill=False, color = col, linewidth = 0.5)
                c1 = Circle(pos, radius = self.Appertur[i].annulus_aperture.r_in, fill=False, color = col, linewidth = 0.5)
                c2 = Circle(pos, radius = self.Appertur[i].annulus_aperture.r_out, fill=False, color = col, linewidth = 0.5)
                ax.add_patch(c)
                ax.add_patch(c1)
                ax.add_patch(c2)
                    
            return [im]

        anim = animation.FuncAnimation(fig, animate_func, frames = len(self.detector)) #interval = 1000 / 30)
        writer = animation.PillowWriter(fps=25) 
        anim.save(path, writer=writer)

        print('Done!')
        
        
    def readCsv(self, path):
        data = pd.read_csv(path)
        
        header = data.columns.values
        qTableHeader = ['id','xcenter','ycenter','aperture_sum','Time','annulus_median','aper_bkg','aper_sum_bkgsub']

        NofAp = len([i for i in header if 'isAst' in i])
        Nast = len([i for i in data.iloc[0,:NofAp] if i == 1])
        
        self.nbOfStar = NofAp - Nast
        self.nbOfAst = Nast
        
        self.results = []
        
        self.starPassages = np.zeros((self.nbOfAst, data.shape[0]))
        
        for i in range(data.shape[0]):
            qTArr = np.zeros((NofAp, len(qTableHeader)))
            
            for row in range(qTArr.shape[0]):
                for col in range(qTArr.shape[1]):
                    
                  
                    if col == 0:
                        qTArr[row, col] = row + 1
                    else:
                        
                        qTArr[row, col] = data.iloc[i, col*NofAp+row] 
                        
            self.results.append(QTable(qTArr, names = qTableHeader))
            
            if data.shape[1] > 32: # previously csv didnt get any star passages flag, so we have to be able to read all format
                for astIdx in range(self.nbOfAst):
                    self.starPassages[astIdx, i] = data.iloc[i, -1-astIdx]
        
    def toDat(self, path, filename, binning = 1, forma = 'mjd', refS = -1, deg = 4, cStd = 2, displayRmFit = False):
        """
        produce 4 files which corresponds to aperture mesurement in ADU, in Mag, differential photometry and mean
        differential photometry
        
        INPUT:
            path: string (path where files should be save)
            filename: string (res_filename.[1-4])
            binning = int (default = 1, -1 for automatic selection) 
            forma = string (default = 'mjd' for other formats, refer to Time.FORMATS from astropy.time)
        
        """
        
        
        if binning == -1: # Auto binning
        
            t = self.toTime(1)
            t = self.timeFormat(t, 'jd')
            
            dth = (t[-1] - t[0])*24
            NdthIn12H = 12 / dth
            
            nbOfPoint = NdthIn12H * t.shape[0]
            
            if (nbOfPoint <= 240 and nbOfPoint >= 100) or nbOfPoint < 100 :
                binning = 1
            
            else:
                binning = int(nbOfPoint / 240 + 0.5)
                
            print("binning: ", binning)
            
        idxOfExtremPoints = self.idxOfDeletedPoint(refS, binning, d = deg, cStd = cStd, plot = displayRmFit) 
        print(idxOfExtremPoints)
        
        asts = self.getAstPhot(binning, False)
        stars =self.getStarPhot(binning, False)
        linesByAsts = []
        
        for i in range(asts.shape[1]):
            if os.path.exists(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.1'):
                os.remove(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.1')
            if os.path.exists(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.2'):
                os.remove(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) +'.2')
            if os.path.exists(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.3'):
                os.remove(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.3')
            if os.path.exists(path + "res_" + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) +'.4'):
                os.remove(path + "res_" + filename + "_ast_" + str(i+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.4')
            linesByAsts.append([ [], [], [], [] ])
            
            
        t = self.toTime(binning)
        t = self.timeFormat(t, forma)
        

        
        if refS == -1:
            refS = np.where(stars[0] == np.max(stars[0]))[0][0]
            print("star ref:", refS + 1)
        

    
        for j in range(asts.shape[1]):
            for i in range(len(t)):
                if i not in idxOfExtremPoints[j]:
                    print(i)
                    strg1 = str(round(t[i] + 0.5, 6)) + ' '
                    strg2 = strg1
                    strg3 = strg1
                    strg4 = strg1
                    
                    #if not self.starPassages[j,i]: #delete stars passages 
                    if not np.isnan(asts[i][j]): 
                            
                         strg1 += str(round(asts[i][j], 4) ) + ' '
                         strg2 += str(round(-2.5*np.log10(asts[i][j]), 4)) + ' '
                         strg3 += str(round(-2.5*np.log10(asts[i][j]) - -2.5*np.log10(stars[i][refS]), 4)) + ' '
                         strg4 += str(round(-2.5*np.log10(asts[i][j]) - np.median(-2.5*np.log10(stars[i])), 4)) + ' '
                   
                    # else:
                    #      print(i, j)
                    #      strg1 += ' ' * (len(str(round(asts[i][j]+1, 4))))
                    #      strg2 += ' ' * (len(str(round(-2.5*np.log10(asts[i][j]+1), 4))))
                    #      strg3 += ' ' * (len(str(round(-2.5*np.log10(asts[i][j]) - -2.5*np.log10(stars[i][refS]) + 1, 4))))
                    #      strg4 += ' ' * (len(str(round(-2.5*np.log10(asts[i][j]) - np.median(-2.5*np.log10(stars[i])) + 1, 4))))
                        
                        
                
                    for k in range(stars.shape[1]):
                        strg1 += str(round(stars[i][k], 4)) + ' '
                        strg2 += str(round(-2.5*np.log10(stars[i][k]), 4)) + ' '
                        strg3 += str(round(-2.5*np.log10(stars[i][k]) - -2.5*np.log10(stars[i][refS]), 4)) + ' '
                    
                    strg1 = strg1[:-1] + '\n'
                    strg2 = strg2[:-1] + '\n'
                    strg3 = strg3[:-1] + '\n'
                    strg4 = strg4[:-1] + '\n'  
                
                    linesByAsts[j][0].append(strg1)
                    linesByAsts[j][1].append(strg2)
                    linesByAsts[j][2].append(strg3)
                    linesByAsts[j][3].append(strg4)
        
            with open(path + 'res_' + filename + "_ast_" + str(j+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.1', 'w') as f:
                for i in linesByAsts[j][0]:
                    f.write(i)
                
            with open(path + 'res_' + filename + "_ast_" + str(j+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) +'.2', 'w') as f:
                for i in linesByAsts[j][1]:
                    f.write(i)
                
            with open(path + 'res_' + filename + "_ast_" + str(j+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.3', 'w') as f:
                for i in linesByAsts[j][2]:
                    f.write(i)
                
            with open(path + 'res_' + filename + "_ast_" + str(j+1) + '_bin_' + str(binning) + "_refs_" + str(refS+1) + '.4', 'w') as f:
                for i in linesByAsts[j][3]:
                    f.write(i)
                
                
                
    def idxOfDeletedPointBis(self, refS, binning, win=5, cStd = 2.5, plot=False):
        

        asts = -2.5*np.log10(Utils.binn(binning, self.astPhot()))
        stars = -2.5*np.log10(Utils.binn(binning, self.starsPhot()))
        
        asts = asts - stars[:, refS].reshape(stars[:, refS].shape[0], 1)
        
        slidingAvg = Utils.slidingAvg(asts, win)
        
        normLc = (asts[:-win + 1] - slidingAvg)**2

        med = np.nanmedian(normLc, axis = 0)
        std = np.nanstd(normLc, axis = 0)
        
        idxByAst = []
        
        for i in range(asts.shape[1]):
            idxByAst.append(np.asarray(np.where(np.invert((normLc < med + cStd*std) & (normLc > med - cStd*std))))[0])
            
        if plot:
            plt.figure()
            plt.scatter(np.arange(asts.shape[0]), asts, linewidths=1, marker='x', label = 'ast')
            plt.scatter(np.arange(slidingAvg.shape[0]), slidingAvg, linewidths=1, marker='x', label = "sliding average")
            plt.scatter(np.arange(normLc.shape[0]), normLc, linewidths=1, marker='x', label = 'normalisation')
            plt.plot(np.ones(normLc.shape)*med)
            plt.plot(np.ones(normLc.shape)*(med + cStd * std))
            plt.plot(np.ones(normLc.shape)*med - cStd * std)
            plt.legend()

            plt.show()
            
            print("med", med, "std", std)
            
        return idxByAst
        
    def idxOfDeletedPoint(self, refS, binning, d = 5, cStd = 1.5, plot=False):
        
        print("poly degree:", d, "c*std:", cStd)

        asts = self.getAstPhot(binning, True, True)
        stars = self.getStarPhot(binning, True)
        
        asts = asts - stars[:, refS].reshape(stars[:, refS].shape[0], 1)
        
        x = self.getTime(binning)
        
        fits = []
       
        normLc = np.zeros(asts.shape)
        
        for i in range(asts.shape[1]):
            fits.append(np.polynomial.polynomial.Polynomial.fit(x[~(np.isnan(asts).any(axis=1))], asts[~(np.isnan(asts).any(axis=1))][:,i], d, window= [0,1]))
            
            normLc[:,i] = asts[:,i] / fits[i](x)

        med = np.nanmedian(normLc, axis = 0)
        std = np.nanstd(normLc, axis = 0)

        
        idxByAst = []
        
        for i in range(asts.shape[1]):
            idxByAst.append(np.asarray(np.where(np.invert((normLc[:,i] < med[i] + cStd*std[i]) & (normLc[:,i] > med[i] - cStd*std[i]))))[0])
            
        
        
        
        if plot:
            for i in range(asts.shape[1]):
                plt.figure()
                plt.scatter(x, asts[:,i], linewidths=1, marker='x', label = 'ast' + str(i))
                plt.scatter(x[~(np.isnan(asts).any(axis=1))], fits[i](x[~(np.isnan(asts).any(axis=1))]), linewidths=1, marker='x', label = "fit" + str(i))
                plt.scatter(x, normLc[:,i], linewidths=1, marker='x', label = 'normalisation' + str(i))
                plt.plot(x, np.ones(normLc.shape[0])*med[i])
                plt.plot(x, np.ones(normLc.shape[0])*(med[i] + cStd * std[i]))
                plt.plot(x, np.ones(normLc.shape[0])*med[i] - cStd * std[i])
                plt.legend()
                plt.gca().invert_yaxis()
                plt.show()
            
                print("med", med, "std", std)
            
        return idxByAst


    def log(self, path, name = "log.txt"):
        """
        

        Parameters
        ----------
        path : STRING
            path where to save the file.
        name : STRING, optional
            name of the file. should end by .txt. The default is "log.txt".

        Returns
        -------
        None.

        """
        
        f = open(os.path.join(path, name), "w")
        

        
        f.write("data rejected during the correction phase:\n\n")
        f.writelines(self.singleton.getDataDeleted())
        f.write("\n")
        f.write("-----------------------------------------\n")
        f.write("star passages:\n\n")
        f.writelines(self.singleton.getDataStarsPassages())
        f.write("\n")
        f.write("-----------------------------------------\n")
        f.write("FWHM of images\n")
        for i in range(len(self.fwhm)):
            f.write(self.detector.seqManager.getFileName(i) + ' fwhm: ' + str(self.getFwhm(i)) + '\n')
        f.close()
        
                
        
if __name__ == '__main__':
  
    
    
    path = r"C:\Users\antoi\OneDrive\Documents\PHD\lightcurve\entrainement\Suh\20-10-22Suh\2020-10-22_gunlod/"
    res = r"C:\Users\antoi\OneDrive\Documents\PHD\lightcurve\entrainement\Suh\res\20-10-22Suh\pres/"
    

    seq = list(np.asarray(glob.glob(path + "*gun*.f*t*"))[:])
    # seq = [seq[0], seq[73],seq[-1]]
    dark = glob.glob(path + "*dark*.f*t*")
    flat = glob.glob(path + "*flat*.f*t*")
    bias = glob.glob(path + "*bias*.f*t*")
    
    #------------ex photometry---------------

    if len(dark) == 0:
        dark = None
        print('DARK EMPTY')
    if len(bias) == 0:
        bias = None
        print('BIAS EMPTY')
    if len(flat) == 0:
        flat = None
        print('FLAT EMPTY')

    d = Detector(seq, flatSeq = flat, biasSeq = bias , darkSeq = dark)
    d.computeImagesCorrection(1500, True)
    print('done')    

    d.findAsteroid(1000, True)
    
    phot = Photometry(d)
    phot.start(3, False, starPassageOfs = 500)
    
    # ----------------ex occult-----------------
    # path = 'C:\\Users\\antoine\\OneDrive\\Documents\\PHD\\lightcurve\\entrainement\\Adorea\\reduced\\target-c1\\'
    # path = 'C:\\Users\\antoine\\OneDrive\\Documents\\PHD\\lightcurve\\entrainement\\21-12-25PST3c\\Imprinetta_occult\\raw\\'
    # seq = glob.glob(path + "*.fit*")
    # occult = Occult(seq)
    # r = 4
    # occult.start(r, r*1.5, r*2)
    
    
    
        