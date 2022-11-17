# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:29:08 2022

@author: antoine
"""

import numpy as np
from cmath import rect, phase
import gc

from tqdm import tqdm

from utils import *
from data_structurs import *


class Corrector:
    def __init__(self, seqManager, flatSeq = None, biasSeq = None, darkSeq = None, exposurKey = None):
        self.seqManager = seqManager
        self.starsPosition = []
        self.paterns = []
        self.drifts = [[np.zeros((2))]]
        self.angles = [np.zeros((2))]
        
        self.exposurKey = exposurKey
      
        self.darkExp = None
        if darkSeq != None:
            if self.exposurKey != None:
                self.darkExp = float(Fit(darkSeq[0]).getHeader()[self.exposurKey])
            else:
                self.darkExp = None
        
        self.biasMaster = SeqManager(biasSeq).getAvg()
        self.darkMaster = SeqManager(darkSeq).getAvg(biasM=self.biasMaster)
        self.flatMaster = SeqManager(flatSeq).getAvg(1, biasM=self.biasMaster)
        print(self.flatMaster)
        
        self.seqManager.setBiasDarkFlat(self.biasMaster, self.darkMaster, self.flatMaster, self.darkExp, self.exposurKey)
        
    def __len__(self):
        return len(self.seqManager)

        
    def getImgShape(self, idx = 0, idx_HDU = 0):
        return self.seqManager.getData(idx, idx_HDU).shape
    
    def getImg(self, idx = 0):
        return self.seqManager.getImg(idx)
    
    def getData(self, idx = 0, HDU = 0):
        return self.seqManager.getData(idx, HDU)
    
    def getReducedData(self, idx = 0, HDU = 0):
        return self.getImg(idx).getReducedData(HDU)
    
    def getInfo(self, idx = 0):
        return self.seqManager.getInfo(idx)
    
    def getHeader(self, idx = 0, HDU = 0):
        return self.seqManager.getHeader(idx)
    
    def histogram(self, idx = 0 , idx_HDU = 0):
        return self.seqManager.histogram(idx, idx_HDU)
    
    def getStarsListOfImg(self, idx):
        return self.starsPosition[idx]
    
    def detectStars(self, offsetTreshStarsDetection = 0, treshOnReduced = False):
        for i in tqdm(range(len(self.seqManager))):
            img = self.seqManager.getImg(i)
            self.starsPosition.append(img.findStars(img.getTresh(treshOnReduced) + offsetTreshStarsDetection, treshOnReduced))#np.median(img.getData()) + 1*np.std(img.getData())
            del img
            gc.collect()
        
    def getImgCenter(self, idx_img = 0, idx_HDU = 0):
        return self.seqManager.getCenter(idx_img, idx_HDU)
        
    def findCloseststars(self, starList, idxOfStar):
        stars = np.zeros((5), dtype = int)
        distance = np.ones((5))*1000000
        
        stars[0] = idxOfStar
        distance[0] = 0
        
        for i,s in enumerate(starList):
            if i == idxOfStar:
                continue
            
            idx, distance = self.insertInArraySorted(distance,  Utils.dist(s, starList[idxOfStar]))

            if idx != -1:
                stars = Utils.insert(stars, idx, i)

        return stars
    
    def buildTriangles(self, stars):
        
        t1 = Triangle(stars[0], stars[1], stars[2])
        t2 = Triangle(stars[0], stars[1], stars[3])
        t3 = Triangle(stars[0], stars[2], stars[4])
        t4 = Triangle(stars[0], stars[2], stars[3])
        t5 = Triangle(stars[0], stars[2], stars[4])
        
        return t1, t2, t3, t4, t5
        
    def buildPatern(self, starsList):
        patern = []
        for i in range(len(starsList)): #iterate over all stars of an image
            stars = self.findCloseststars(starsList, i)
            t1, t2, t3, t4, t5 = self.buildTriangles(starsList[stars.astype(int)])
            patern.append(Patern(t1, t2, t3, t4, t5))
        return patern
    
    def buildPaterns(self):
        for i in tqdm(self.starsPosition):
            self.paterns.append(self.buildPatern(i))
            
    def pop(self, idx = -1):
        self.seqManager.pop(idx)
        self.drifts.pop(idx)
        self.angles.pop(idx)
        self.paterns.pop(idx)
        self.starsPosition.pop(idx)
        
    
    def computeAngles(self, paternsOfImage1, paternsOfImage2):
        
        angle = []
        
        for i in paternsOfImage1:
            for j in paternsOfImage2:
                if i ==j:
                    angle.append(i.computeAngle(j))   
                    break
        return angle
    
    def computeDrift(self, paternsOfImage1, paternsOfImage2):
        
        drift = []
        for i in paternsOfImage1:
            for j in paternsOfImage2:
                
                if i == j:   
                    # print("pat 1: \n",i, "\npat 2: \n", j, "\ndist: ", i.computeDistance(j))
                    drift.append(i.computeDistance(j))
                    break
                    
        return drift
    
    def correctStarsFromRot(self, arrayToCorrect, idx, coefMultAngle = -1):
        
        angle = coefMultAngle*self.avgAng(idx)
        if np.isnan(angle):
            return
        
        arrayToCorrect = np.asarray(arrayToCorrect)

        if arrayToCorrect.shape[0] == 0:
            return arrayToCorrect
        
        center = np.flip(self.getImgCenter(idx))
 
        stars = arrayToCorrect - center
        rot = np.asarray([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])
     
        return np.dot(stars, rot) + center
     
    def correctStarsFromRotTest(self, arrayToCorrect, idx, coefMultAngle = -1):
        
        arrayToCorrect = np.asarray(arrayToCorrect)
        print("before", arrayToCorrect[0])
        if arrayToCorrect.shape[0] == 0:
            return arrayToCorrect
        
        angle = coefMultAngle*self.avgAng(idx)
        
        if np.isnan(angle):
            return
        
        newArry = np.column_stack((arrayToCorrect, np.ones((arrayToCorrect.shape[0]))))
        center = self.getImgCenter(idx)

        alpha = np.cos(angle)
        beta = np.sin(angle)
        
        rot = np.asarray([[alpha, beta, (1-alpha)*center[0] - beta*center[1]],[-beta, alpha, beta * center[0] + (1 - alpha)*center[1]]])
        print("after", np.dot(rot, newArry.T).T[0])
        
        return np.dot(rot, newArry.T).T
    
    def correctStars(self, arrayToCorrect, idx, coefMultAngle = -1):
        
        arrayToCorrect = np.asarray(arrayToCorrect)
        
        if arrayToCorrect.shape[0] == 0:
            return arrayToCorrect
        
        return self.correctStarsFromRot(arrayToCorrect+ coefMultAngle*self.avgDrif(idx), idx, coefMultAngle) 
    
      
        
    def correctPaternFromRot(self, paterns, idx):
        for p in paterns:
            p.correctRot(-self.avgAng(idx), np.flip(self.getImgCenter(idx)))
        
    def imagesDrift(self, offsetTreshStarsDetection = 0, treshOnReduced = False):
        self.detectStars(offsetTreshStarsDetection, treshOnReduced)
        print("stars detection done")
        print("begin to build patern")
        self.buildPaterns()
        print("patern building done")
        print("start partn comparaison")
    
        for i, p1 in enumerate(tqdm(self.paterns)): # iterate along images
            if i == 0:
                continue
            
            
            
            a = self.computeAngles(self.paterns[0], p1)
            self.angles.append(np.asarray(a))
            
            if len(self.angles[-1]) != 0:
                self.starsPosition[i] = self.correctStarsFromRot(self.starsPosition[i], i)
                self.correctPaternFromRot(p1, i)
        
            d = self.computeDrift(p1, self.paterns[0])
            self.drifts.append(d)
      
        print("patern comparaison done")
        
    def insertInArraySorted(self, arr, val):
        for i, e in enumerate(arr):
            if e > val:
                arr = Utils.insert(arr, i, val)
                return i, arr
        return -1, arr
    
    def medDrif(self, idx):
        return np.nanmedian(self.drifts[idx], axis = 0)
    
    def avgDrif(self, idx):
        return np.nanmean(self.drifts[idx], axis = 0)
    
    def medAng(self, idx):
        return np.nanmedian(self.angles[idx][np.logical_not(np.isnan(self.angles[idx]))])
    
    def avgAng(self, idx):
        angs = self.angles[idx][np.logical_not(np.isnan(self.angles[idx]))]
        return phase(sum(rect(1, d) for d in angs)/len(angs))
    
    def correctedImg(self, idx = 0, HDU_idx = 0):
        """
        

        Parameters
        ----------
        idx : int, optional
            DESCRIPTION. Index of the images to correct. The default is 0.
        HDU_idx : int, optional
            DESCRIPTION. HDU idex in fits file. The default is 0.

        Returns
        -------
        img : numpy arrat 
            Return the image corrected from drift and rotation.

        """
        
        drift = self.avgDrif(idx)
        angle = self.avgAng(idx)
        Utils.imshow(self.getData(idx, HDU_idx))
        img = Utils.rotate_image(self.getData(idx, HDU_idx), -Utils.rtod(angle))
        
        Utils.imshow(img)
        img = np.roll(img, -int(drift[0]+0.5), axis=1)
        Utils.imshow(img)
        img = np.roll(img, -int(drift[1]+0.5), axis=0)
        
        return img
    
    def getSuperImg(self, idx_ims = None, HDU_idx = 0):
        """
        Parameters
        ----------
        idx_ims: numpy array or list of int
            DESCRIPTION. array of index of images to combine
        HDU_idx : int, optional
            DESCRIPTION. Index of the hdu in fits file.The default is 0.
        

        Returns
        -------
        numpy array of shape (number of image, images raw, images col)

        """
        if idx_ims == None:
            idx_ims = np.arange(len(self))
            
        imgs = np.zeros(self.getData(0).shape, dtype = float)
       
        for i in idx_ims:
            imgs = imgs + self.correctedImg(i)
        
        return imgs / len(self)