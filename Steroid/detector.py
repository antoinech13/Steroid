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

    
    def rejectBadData(self):
        """
        reject bad data. (if they are no drift value or angle which was found
        """
        singleton = Singleton()
        
        cpt = 0
        N = len(self)
        
        
        while len(self.drifts) != cpt:
            element = self.drifts[cpt]
   
            if len(element) == 0 or len(self.angles[cpt]) == 0 or np.isnan(self.avgAng(cpt)):
                
                if type(self.starsPosition[cpt]) != type(None):
                   starDetected =  len(self.starsPosition[cpt])
                else:
                    starDetected = 0
                
                
                print('pop', 'drift:', len(element), 'angle:', len(self.angles[cpt]),"nb of stars", starDetected)
                singleton.appDataDeleted( self.seqManager.getFileName(cpt) + ', drift: ' + str(len(element)) + ', angle: ' + str( len(self.angles[cpt])) +  ", nb of stars " + str(starDetected)  + '\n')
                
                self.pop(cpt)
            else:
                cpt+=1
                
        print(N - cpt, "data rejected")
    
    def computeImagesDrift(self, offsetTreshStarsDetection = 0, treshOnReduced = False):
        """
        
        Parameters
        ----------
        offsetTreshStarsDetection : int, optional
            offset to add to fine tune the sensitivity of stars detection of frames.
            The default is 0.
        
        treshOnReduced : bool, optional 
            choose if tresholing should be done on reduced frame.
            

        Returns
        -------
        None.

        """
        
        
        print("start computing drift")
        self.imagesDrift(offsetTreshStarsDetection, treshOnReduced)
        print("done \n star rejected bad data")
        self.rejectBadData()
        print('done')

    def findDifStar(self, stars1, stars2, eps):
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

        
        if (starPosition == np.asarray([1422., 381])).all():
            print('max1', max1, "max2", max2, "avg1", avg1, "avg2", avg2, "med1", np.nanmedian(img1), "med2", np.nanmedian(img2))
        
        if max1 == 0 or max2 == 0:
            return True
        
        elif max1 < avg1 or max2 < avg2:
            return False
        print("IA", img1.shape, img1.shape ,self.model.predict(img2.reshape(1,20,20,1))[0,0] ,  self.model.predict(img1.reshape(1,20,20,1))[0,0])
        return self.model.predict(img2.reshape(1,20,20,1))[0,0] > 0.5 and  self.model.predict(img1.reshape(1,20,20,1))[0,0] > 0.5
        
    def isOutOfBoundaries(self, pos, eps = 20):
        extrema = self.findDriftExtrema()
        return (pos + extrema + eps > self.shape).any() or (pos - extrema - eps < np.zeros((2))).any()
    
    def extractAst(self, dif, idxofFirstImg, idxOfSecondIng, distMax = 200):
        ast = []
        idx = []
        time1 = self.getImg(idxofFirstImg).getTime(self.key, self.format)
        time2 = self.getImg(idxOfSecondIng).getTime(self.key, self.format)
        print("img1: ", idxofFirstImg, "img2: ", idxOfSecondIng)
        print("dif: ", dif)
        for i, e in enumerate(dif):
            for j in range(i, len(dif)):
                f = dif[j]
                # print("e: ", e, "f: ", f,"dist: ", Utils.dist(e,f) < distMax, "e out boundarie: ", not self.isOutOfBoundaries(e), "f out boundarie: ", not self.isOutOfBoundaries(f), "i is present: ", not Utils.isPresent(idx, i), "i!=j: ", i != j, "e not isstar: ", not self.isStar(e, idxofFirstImg, idxOfSecondIng), "f not isstar: ", not self.isStar(f, idxofFirstImg, idxOfSecondIng))
                if Utils.dist(e,f) < distMax and not self.isOutOfBoundaries(e) and not self.isOutOfBoundaries(f) and not Utils.isPresent(idx, i) and i != j and not self.isStar(e, idxofFirstImg, idxOfSecondIng) and not self.isStar(f, idxofFirstImg, idxOfSecondIng):
                    # print("idx: ", Utils.isPresent(idx, i),idx, i, j, 'eeee: ', e, ' fff: ', f)
                    print("e: ", e, "f: ", f,"dist: ", Utils.dist(e,f) < distMax, "e out boundarie: ", not self.isOutOfBoundaries(e), "f out boundarie: ", not self.isOutOfBoundaries(f), "i is present: ", not Utils.isPresent(idx, i), "i!=j: ", i != j, "e not isstar: ", not self.isStar(e, idxofFirstImg, idxOfSecondIng), "f not isstar: ", not self.isStar(f, idxofFirstImg, idxOfSecondIng))
                    ast.append(e)
                    idx.append(i)
                    
                    self.asteroidsSpeed.append((e-f)/(time2 - time1).to_value('sec'))
                    print('-------------------------------------------------------------------------\natsteroid Position first img   asteroid position last img   deltImg   speed by img,             estimation pos of ast on last\n', e, '                   ', f, '              ' , time2 - time1, '       ', (e-f)/(time2 - time1), e-(time2 - time1)*(e-f)/(time2 - time1), '\n-------------------------------------------------------------------------')
                
        return ast

    def findBestIdx(self, plot = False):
        nbOfStarsDetected = np.zeros((len(self.starsPosition)))
        
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
        
        idxOfFirst, idxOfSecond = self.findBestIdx()
    
        print('idx: ', idxOfFirst, idxOfSecond)
        dr = self.avgDrif(idxOfFirst)
        dr2 = self.avgDrif(idxOfSecond)

        stars1 = np.asarray(self.getImg(idxOfFirst).findStars(self.getImg(idxOfFirst).getTresh(treshOnReduced) + offsetTreshAstsDetection, treshOnReduced))  #np.median(self.getData()) + np.std(self.getData())*0.4
        stars1 = self.correctStarsFromRot(stars1, idxOfFirst)
        stars1 = stars1 - dr.astype(float)
        stars2 = np.asarray(self.getImg(idxOfSecond).findStars(self.getImg(idxOfFirst).getTresh(treshOnReduced) + offsetTreshAstsDetection, treshOnReduced)) #- dr  #+ np.std(self.getData(idxOfSecond))
        print(stars2)
        stars2 = self.correctStarsFromRot(stars2, idxOfSecond)
        stars2 = stars2 - dr2.astype(float)
        print(stars2)
        
        Utils.imshowstar(self.getData(idxOfFirst), stars1, -self.avgAng(idxOfFirst), dr)
        Utils.imshowstar(self.getData(idxOfSecond), stars2, -self.avgAng(idxOfSecond), dr2)
        Utils.imshowstar(self.getData(idxOfFirst), stars2)
        
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
    
    def checkPaterns(self, idxOfImage = 0, patidx = None):
        def patToArw(pat):
            s1 = pat.t1.s1
            s2 = pat.t1.s2
            s3 = pat.t1.s3
            s4 = pat.t2.s3
            s5 = pat.t3.s3
            
            d2 = s2 - s1
            d3 = s3 - s1
            d4 = s4 - s1
            d5 = s5 - s1
            
            c = [random.random(), random.random(), random.random()]
            a2 = Arrow(s1[0], s1[1], d2[0], d2[1], color=c)
            a3 = Arrow(s1[0], s1[1], d3[0], d3[1], color=c)
            a4 = Arrow(s1[0], s1[1], d4[0], d4[1], color=c)
            a5 = Arrow(s1[0], s1[1], d5[0], d5[1], color=c)
            
            return a2, a3, a4, a5
        
        fig,ax = plt.subplots(1)
        for i, pat in enumerate(self.paterns[idxOfImage]):
            if patidx == None or i == patidx:
                a1, a2, a3, a4 = patToArw(pat)
                ax.add_patch(a1)
                ax.add_patch(a2)
                ax.add_patch(a3)
                ax.add_patch(a4)
            
        img = Utils.rotate_image(self.getData(idxOfImage), -Utils.rtod(self.avgAng(idxOfImage)))
        ax.imshow(img, cmap="Greys", vmin = np.mean(img)*0.5, vmax = np.mean(img)/0.5)
        
    def checkImgAlignement(self, idx_of_image):

        img1 = self.getData(0)

        img3 = self.correctedImg(idx_of_image) + img1 
        fig,ax = plt.subplots(1)
        ax.imshow(img3, cmap="Greys", vmin = np.mean(img3)*0.5, vmax = np.mean(img3)/0.5)
        
    def show(self, idx, ang = 0, drift = np.zeros(2)):
        img = self.getData(idx)
        Utils.rotate_image(img, -Utils.rtod(ang))
        
        img = np.roll(img, -int(drift[0]+0.5), axis=1)
        img = np.roll(img, -int(drift[1]+0.5), axis=0)
        
        
        fig,ax = plt.subplots(1)
        ax.imshow(img, cmap="Greys", vmin = np.mean(img)*0.5, vmax = np.mean(img)/0.5)
        
        
        
    def add(self, idximg1, idximg2):
        img = (self.getData(idximg1) + self.getData(idximg2))/2
        fig,ax = plt.subplots(1)
        ax.imshow(img, cmap="Greys", vmin = np.mean(img)*0.5, vmax = np.mean(img)/0.5)
        
    def showStarsOfTwoImg(self, idxOfSecond):
        img1 = self.getData()
        img2 = self.getData(idxOfSecond)
        s1 = self.starsPosition[0]
        s2 = self.starsPosition[idxOfSecond]
        s2 = s2 - self.avgDrif(idxOfSecond)
        
        fig,ax = plt.subplots(1)
        
        for i in range(s1.shape[0]):
            c = Circle(s1[i], radius = 25, fill=False)
            ax.add_patch(c)
            
        for i in range(s2.shape[0]):
            c = Circle(s2[i], radius = 25, fill=False, color="r")
            ax.add_patch(c)
        
        ax.imshow((img1+img2)*0.5, cmap="Greys", vmin = np.mean((img1+img2)*0.5)*0.5, vmax = np.mean((img1+img2)*0.5)/0.5)
        
        
        
    def div(self, idx_of_image):
        img = self.getData(idx_of_image)
        dr = self.avgDrif(idx_of_image)
        img1 = self.getData(0)
        img2 = np.roll(img, -int(dr[0]+0.5), axis = 1)
        img2 = np.roll(img2, -int(dr[1]+0.5), axis = 0)
        
        img3 = (img1 / img2)
        fig,ax = plt.subplots(1)
        ax.imshow(img3, cmap="Greys", vmin = np.mean(img3)*0.5, vmax = np.mean(img3)/0.5)
        return img3

    def imshowstar(self, idx = 0):
        
        stars = self.starsPosition[idx]
        img = Utils.rotate_image(self.getData(idx), -Utils.rtod(self.avgAng(idx)))
        fig,ax = plt.subplots(1)
        for i in stars:
            c = Circle(i, radius = 10, fill=False)
            ax.add_patch(c)
        ax.imshow(img, cmap="Greys", vmin = np.mean(img[img>0])*0.5, vmax = np.mean(img[img>0])/0.5)
        
    def imshowstarrot(self, idx = 0):
        stars = self.starsPosition[idx]
        img = self.getData(idx)
        img = Utils.rotate_image(img, -Utils.rtod(self.avgAng(idx)))
        fig,ax = plt.subplots(1)
        for i in stars:
            c = Circle(i, radius = 25, fill=False)
            ax.add_patch(c)
        ax.imshow(img, cmap="Greys", vmin = np.mean(img)*0.5, vmax = np.mean(img)/0.5)
        
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
        
    def getAstPositionAtImg(self, idx):
        return np.asarray(self.asteroidsPositionFirstImage) - (self.getImg(idx).getTime(self.key, self.format) - self.getImg().getTime(self.key,  self.format)).to_value("sec")*np.asarray(self.asteroidsSpeed)
    
    def nofa(self):
        return np.asarray(self.asteroidsPositionFirstImage).shape[0]