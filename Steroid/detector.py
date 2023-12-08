# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:32:53 2022

@author: antoine
"""
from corrector import Corrector

import numpy as np
from tensorflow import keras

from data_structurs import SeqManager, Fit
from utils import *
from singleton import Singleton
import os


class Detector(Corrector):
    def __init__(self, imageSeq, flatSeq = None, biasSeq = None, darkSeq = None):
        
        self.seqManager = SeqManager(imageSeq)
        
        # self.corrector = Corrector(self.seqManager)
        print(self.getHeader())
        self.key = str(input("ENTER THE KEY OF A TIME REF WITH | SEPARATOR IF NEEDED (exemple JD | jd or DATE-OBS | isot. Refere to Time.FORMATS from astropy.time): "))#(exemple JD | jd or DATE-OBS | %Y-%m-%dT%H:%M:%S.%f): "))
        self.isTimeCentred = str(input("DO YOU WANT TO CENTRED TIME ? (Yes/No) "))
        
        if self.isTimeCentred.lower() == 'Yes' or self.isTimeCentred.lower() == 'y':
            self.isTimeCentred = True
        else:
            self.isTimeCentred = False
        
        exposurKey = str(input("ENTER THE KEY OF EXPOSURE (no sperator): "))
        
        self.format = self.key.split("|")[-1].replace(' ', '')
        self.key =  self.key.split("|")[0].replace(' ', '')
        
        self.model = keras.models.load_model(os.getcwd()+"/best_model_alexnet_good3_99.h5")
        super().__init__(self.seqManager, flatSeq, biasSeq, darkSeq, exposurKey)
        
        self.asteroidsPositionFirstImage = []
        self.shape = self.getData().shape
        self.asteroidsSpeed = []

    def __len__(self):
        return len(self.seqManager)

    
    def computeImagesCorrection(self, offsetTreshStarsDetection = 0, treshOnReduced = False):
        """
        
        Parameters
        ----------
        offsetTreshStarsDetection : int, optional
            offset to add to fine tune the sensitivity of stars detection of frames.
            The default is 0.
        
        treshOnReduced : bool, optional 
            choose if tresholing should be done on reduced frame.
            

        """
        
        
        print("start computing drift")
        self.imagesCorrection(offsetTreshStarsDetection, treshOnReduced)
        print("done \n star rejected bad data")
        self.rejectBadData()
        print('done')

    def findDifStar(self, stars1, stars2, eps):
        """
        After image correction, try to find which object are present at a position on a frame and not present at the same position in an other frame.

        Parameters
        ----------
        stars1 : list
            DESCRIPTION. list of object detected on the frame 1
        
        stars2 : list
            DESCRIPTION. list of object detected on the frame 2
        
        eps: int, optional
            DESCRIPTION. represent the minimum distance (in pixels) that two object should be from each other on the first and the last frame to be detected as star or moving object 
        
        Returns
        -------
        list of all object present in stars1 and not in stars2 at the same position 

        """
        
        
        dif = []
        isStar = False
        for i in stars1:
           for j in stars2:
               if (np.abs(i - j) <= eps).all():
                   isStar = True
                   continue
           if not isStar:
               dif.append(i)
           isStar = False
           
        return dif
    
    def isStarOutOfShape(self, star, idxOfImage):
        return star[0] < 0 or star[0] > self.getImgShape(idxOfImage)[1] or star[1] < 0 or star[1] > self.getImgShape(idxOfImage)[0]
        
    def isStar(self, starPosition, idxOfFirst, idxOfSecond):
        
        dr = self.avgDrif(-1)*-1
        
        img1 = Utils.rotate_image(self.getData(idxOfFirst), -Utils.rtod(self.avgAng(idxOfFirst))).astype(float)
        img2 = np.roll(Utils.rotate_image(self.getData(idxOfSecond), -Utils.rtod(self.avgAng(idxOfSecond))).astype(float), int(dr[0]+0.5), axis=1)
        img2 = np.roll(img2, int(dr[1]+0.5), axis=0)
        
        delt = 5
       
        if self.isStarOutOfShape(starPosition - delt, idxOfFirst) or self.isStarOutOfShape(starPosition - delt, idxOfSecond) or self.isStarOutOfShape(starPosition + delt, idxOfFirst) or self.isStarOutOfShape(starPosition + delt, idxOfSecond):
            return True
        
        img1 = img1[int(starPosition[1] - delt) : int(starPosition[1] + delt), int(starPosition[0] - delt) : int(starPosition[0] + delt)]
        img2 = img2[int(starPosition[1] - delt) : int(starPosition[1] + delt), int(starPosition[0] - delt) : int(starPosition[0] + delt)]
        
        img1 = np.pad(img1, 10 - delt, 'constant', constant_values=(0))
        img2 = np.pad(img2, 10 - delt, 'constant', constant_values=(0))
        
        max1 = np.max(img1)
        max2 = np.max(img2)
        
        avg1 = np.nanmean(img1)
        avg2 = np.nanmean(img2)

        
        if max1 == 0 or max2 == 0:
            return True
        
        elif max1 < avg1 or max2 < avg2:
            return False
        
        modelPredictionOnImage1 = self.model.predict(img1.reshape(1,20,20,1))[0,0]
        modelPredictionOnImage2 = self.model.predict(img2.reshape(1,20,20,1))[0,0]
        
        print("IA:  im1 shape: ", img1.shape, " im1 prediction: ", modelPredictionOnImage1, " im2 shape: ", img2.shape , " im2 prediction: ", modelPredictionOnImage2  )
        
        return modelPredictionOnImage2 > 0.5 and  modelPredictionOnImage1 > 0.5
        
    def isOutOfBoundaries(self, pos, eps = 20):
        
        extrema = self.findDriftExtrema()
        
        return (pos + extrema + eps > self.shape).any() or (pos - extrema - eps < np.zeros((2))).any()
    
    def extractAst(self, dif, idxofFirstImg, idxOfSecondIng, distMax = 400):
        ast = []
        idx = []
        time0 = self.getImg(0).getTime(self.key, self.format)
        time1 = self.getImg(idxofFirstImg).getTime(self.key, self.format)
        time2 = self.getImg(idxOfSecondIng).getTime(self.key, self.format)
        print("img1: ", idxofFirstImg, "img2: ", idxOfSecondIng)
        print("dif: ", dif)
        for i, e in enumerate(dif):
            for j in range(i, len(dif)):
                f = dif[j]
                
                if Utils.dist(e,f) < distMax and not self.isOutOfBoundaries(e) and not self.isOutOfBoundaries(f) and not Utils.isPresent(idx, i) and i != j and not self.isStar(e, idxofFirstImg, idxOfSecondIng) and not self.isStar(f, idxofFirstImg, idxOfSecondIng):
                   
                    print("e: ", e, "f: ", f,"dist: ", Utils.dist(e,f) < distMax, "e out boundarie: ", not self.isOutOfBoundaries(e), "f out boundarie: ", not self.isOutOfBoundaries(f), "i is present: ", not Utils.isPresent(idx, i), "i!=j: ", i != j, "e not isstar: ", not self.isStar(e, idxofFirstImg, idxOfSecondIng), "f not isstar: ", not self.isStar(f, idxofFirstImg, idxOfSecondIng))
                    self.asteroidsSpeed.append((e-f)/(time2 - time1).to_value('sec'))
                    ast.append(e - self.asteroidsSpeed[-1] * (time0 - time1).to_value('sec'))
                    
                    idx.append(i)
                    
                    
                    print('-------------------------------------------------------------------------\natsteroid Position first img   asteroid position last img   deltImg   speed by img,             estimation pos of ast on last\n', e, '                   ', f, '              ' , time2 - time1, '       ', (e-f)/(time2 - time1), e-(time2 - time1)*(e-f)/(time2 - time1), '\n-------------------------------------------------------------------------')
                
        return ast

    def findBestIdx(self, plot = False):
        
        nbOfStarsDetected = np.zeros(( len(self.starsPosition) ))
        
        for i in range(len(self.starsPosition)):
            nbOfStarsDetected[i] = len(self.starsPosition[i])
            
        med = np.nanmedian(nbOfStarsDetected)
        std = np.nanstd(nbOfStarsDetected)
        
        idxFirst = None
        idxLast = None
        
        for i in range(len(nbOfStarsDetected)):
            if nbOfStarsDetected[i] > med*0.9 and nbOfStarsDetected[i] < med*1.1 and idxFirst == None:
                idxFirst = i
            if nbOfStarsDetected[-i-1] > med*0.9 and nbOfStarsDetected[-i-1] < med*1.1 and idxLast == None:
                idxLast = -i-1
        
        if plot:
            plt.figure()
            plt.plot(nbOfStarsDetected)
            plt.plot(np.ones(nbOfStarsDetected.shape)*med)
            plt.plot(np.ones(nbOfStarsDetected.shape)*( med*0.9 ), color = 'green')
            plt.plot(np.ones(nbOfStarsDetected.shape)*( med*1.1 ), color = 'green')
            
            plt.xlabel("images of the sequence")
            plt.ylabel("number of stars detected")

        return idxFirst, idxLast
    
    def findAsteroid(self, offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2):
        """
        find asteoird based on what object are present on first frame and not on the last one and reciprocaly. 

        Parameters
        ----------
        offsetTreshAstsDetection : float, optional
            DESCRIPTION. Offset add to the treshold value .The default is 0.
        
        treshOnReduced: boolean, optional
            DESCRIPTION. if true, will estimate the threshold value on reduced frames. The default is False
        
        eps: int, optional
            DESCRIPTION. represent the minimum distance (in pixels) that two object should be from each other on the first and the last frame to be detected as star or moving object 

        """
        
        
        idxOfFirst, idxOfSecond = self.findBestIdx()
    
        print('idx: ', idxOfFirst, idxOfSecond)
        dr = self.avgDrif(idxOfFirst)
        dr2 = self.avgDrif(idxOfSecond)

        stars1 = np.asarray(self.getImg(idxOfFirst).findStars(self.getImg(idxOfFirst).getTresh(treshOnReduced) + offsetTreshAstsDetection, treshOnReduced))  #np.median(self.getData()) + np.std(self.getData())*0.4
        stars1 = self.correctStarsFromRot(stars1, idxOfFirst)
        stars1 = stars1 - dr.astype(float)
        
        stars2 = np.asarray(self.getImg(idxOfSecond).findStars(self.getImg(idxOfFirst).getTresh(treshOnReduced) + offsetTreshAstsDetection, treshOnReduced)) #- dr  #+ np.std(self.getData(idxOfSecond))
        stars2 = self.correctStarsFromRot(stars2, idxOfSecond)
        stars2 = stars2 - dr2.astype(float)
        
        
        print("dif1:", self.findDifStar(stars1, stars2, eps), "dif2:",self.findDifStar(stars2, stars1, eps))
        dif = self.findDifStar(stars1, stars2, eps) + self.findDifStar(stars2, stars1, eps)
        self.asteroidsPositionFirstImage = self.extractAst(dif, idxOfFirst, idxOfSecond)
        
        
    def findDriftExtrema(self):
        maxx = -100000
       
        for i in self.drifts:
            val = np.max(np.abs(i))
            
            if val > maxx:
                maxx = val
                
        return np.ones((2), dtype=int)*maxx
    
        
    def astSpeed(self, idx = 0):
        """
        return the total speed of the asteroid number (idx)

        Parameters
        ----------
        idx : int, optional
            DESCRIPTION. index of the asteroid to get the speed .The default is 0.

        Returns
        -------
        int. the total speed of the asteroid

        """
        
        return np.sqrt( np.sum( self.asteroidsSpeed[idx]**2 ) )
    

    def fasterAst(self):
        """
        return the index of the fasted asteroid

        Returns
        -------
        idx : int
            DESCRIPTION. index of the fasted asteroid

        """
        
        
        idx = 0
        speed = 0
        
        for i in range(len(self.asteroidsSpeed)):
            if self.astSpeed(i) > speed:
                speed = self.astSpeed(i)
                idx = i
        
        return idx
            
                       
    def slowestAst(self):
        """
        return the index of the fasted asteroid
    
        Returns
        -------
        idx : int
            DESCRIPTION. index of the fasted asteroid

        """
            
        idx = 0
        speed = 1000000
            
        for i in range(len(self.asteroidsSpeed)):
            if self.astSpeed(i) < speed:
                speed = self.astSpeed(i)
                idx = i
            
        return idx
        
    def getAstPositionAtImg(self, idx): #pos0 - (timei - time0) * speed
        return np.asarray(self.asteroidsPositionFirstImage) - ( self.getImg(idx).getTime(self.key, self.format) - self.getImg().getTime(self.key,  self.format)).to_value("sec")*np.asarray(self.asteroidsSpeed)
    
    def nofa(self):
        return np.asarray(self.asteroidsPositionFirstImage).shape[0]