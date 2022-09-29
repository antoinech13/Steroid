# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:29:07 2022

@author: antoine
"""

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arrow
from matplotlib.backend_bases import MouseButton
import cv2

class Utils:
    def patToArw(pat):
        s1 = pat.t1.s1
        s2 = pat.t1.s2
        s3 = pat.t1.s3
        s4 = pat.t2.s3
        
        d2 = s2 - s1
        d3 = s3 - s1
        d4 = s4 - s1
        c = [random.random(), random.random(), random.random()]
        a2 = Arrow(s1[0], s1[1], d2[0], d2[1], color=c)
        a3 = Arrow(s1[0], s1[1], d3[0], d3[1], color=c)
        a4 = Arrow(s1[0], s1[1], d4[0], d4[1], color=c)
        
        return a2, a3, a4

    def imshowpatern(img, paterns):
        fig,ax = plt.subplots(1)
        for i in paterns:
            a1, a2, a3 = Utils.patToArw(i)
            ax.add_patch(a1)
            ax.add_patch(a2)
            ax.add_patch(a3)
            
        ax.imshow(img, cmap="Greys", vmin = np.mean(img)*0.5, vmax = np.mean(img)/0.5)
        
    def polarRot(starArray, shape, angle):
        
        if np.isnan(angle):
            print('angle is Nan')
            return
        
        center = ((np.asarray(shape) / 2) + 0.5).astype(int)
       
        stars = starArray - center 
        rot = np.asarray([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])
        return np.dot(stars, rot) + center

    def dist(star1, star2 = [0,0]):
        return np.sum((np.asarray(star2) - np.asarray(star1))**2)**0.5
    
    def insert(arr, idx, val):
        arr = np.insert(arr, idx, val)
        return np.delete(arr, -1)
    
    def isPresent(lis, val):
        for i in lis:
            if i == val:
                return True
        return False
    
    def correctedImg(img, drift = np.zeros((2)), angle = 0):
        
        img = Utils.rotate_image(img, angle)
        img = np.roll(img, -int(drift[0]+0.5), axis=1)
        img = np.roll(img, -int(drift[1]+0.5), axis=0)
        
        return img
    
    def imshowstar(img, stars, angle = None, drift = np.zeros((2)), isImageToMove = True):
        fig,ax = plt.subplots(1)
        
        if angle != None and angle != 0.0:
            if isImageToMove:
                img = Utils.rotate_image(img, Utils.rtod(angle))
            else:
                stars = Utils.polarRot(stars, img.shape, angle)
                
        if (drift != np.zeros((2))).all():
            if isImageToMove:
                img = np.roll(img, -int(drift[0]+0.5), axis=1)
                img = np.roll(img, -int(drift[1]+0.5), axis=0)
            else:
                stars -= drift
                
        for i in stars:
            c = Circle(i, radius = 20, fill=False)
            ax.add_patch(c)
            
        ax.imshow(img, cmap="Greys", vmin = np.nanmean(img)*0.5, vmax = np.nanmean(img)/0.5)

    def getPosFromImage(event, posList, isClose):
        if event.button is MouseButton.LEFT:
           posList.append(np.asarray([event.xdata, event.ydata]))
           print("Position x: ", posList[-1][0], ' y: ', posList[-1][1])
           
        if event.button is MouseButton.RIGHT:
            isClose[0] = True
            plt.close()
    
    def add(img1, img2):
        fig,ax = plt.subplots(1)
        ax.imshow((img1 + img2)/2, cmap="Greys", vmin = np.mean((img1 + img2)/2)*0.5, vmax = np.mean((img1 + img2)/2)/0.5)
        plt.show()
            
    def centred(xy, img, eps):
        
        for i in range(len(xy)):
            x = int(xy[i][0] + 0.5)
            y = int(xy[i][1] + 0.5)
            eps = int(eps + 0.5)

            cent = img[y - eps : y + eps, x - eps : x + eps]

            x = np.linspace(0, 2*eps - 1, 2*eps)            
            idxArr = np.zeros((2, cent.shape[0], cent.shape[1]))
            idxArr[0,:], idxArr[1,:] = np.meshgrid(x,x) 
            cent = cent - np.median(cent)
            cent[cent<0] = 0
            idxArr = (np.sum(np.sum(idxArr * cent**2, axis = 1), axis = 1) / np.sum(cent**2) + 0.5).astype(int)
     
            # Utils.imshow(cent)
            ##pos = np.flip(np.asarray(np.where(cent == np.max(cent)))[:,0].reshape((2)))
            ##xy[i] = np.asarray(xy[i]) + (pos - eps)
            
            xy[i] = np.asarray(xy[i]) + idxArr - eps
            
        return xy

    def imshow(img):
        
        plt.figure()
        plt.imshow(img, cmap='Greys', vmin = np.median(img[img>0])*0.5, vmax = np.median(img[img>0]) / 0.5)
        
    def imshowap(img, r, pos):
      
        img = img[int((pos[1] - r)+0.5): int((pos[1] + r)+0.5), int((pos[0] - r) + 0.5) : int((pos[0] + r) + 0.5)]
        c = Circle(np.asarray(img.shape)/2, radius = 50, fill=False)
        fig,ax = plt.subplots(1)
        ax.add_patch(c)
        ax.imshow(img, cmap='Greys', vmin = np.median(img)*0.5, vmax = np.median(img) / 0.5)
      
        return np.sum(img)
    
    
    def binn(by, file):
        
        pointer = 0
        
        newSize = int(file.shape[0] / by)
        limit  = newSize * by
        
        if len(file.shape) > 1:
            print(file.shape, newSize)
            newData = np.zeros((newSize, file.shape[1]), dtype=file.dtype)
        else:
            newData = np.zeros((newSize), dtype=file.dtype)
        
        while pointer < limit:
            
            for i in range(by):
                if len(file.shape) > 1:
                    newData[int(pointer / by), :] = file[i + pointer, :] + newData[int(pointer / by), :]
                else :
                    newData[int(pointer / by)] = file[i + pointer] + newData[int(pointer / by)]
                    
            pointer += by
    
        return newData / by    
    
    def mask(arr, mask, isAst):
        
        arr = np.asarray(arr)
        newArr = []
    
        for numOfAst in range(mask.shape[0]):
            if not isAst:
            
                element = []
                for nbOfObjects in range(arr.shape[0]):
                    element.append(arr[nbOfObjects][mask[numOfAst] == True])
                newArr.append(np.asarray(element))
            
            else:
    
                newArr.append(arr[numOfAst][mask[numOfAst] == True])
                
        return newArr
    
    def rescalMultiBis(asts, stars):
        mmax = np.nanmax(asts)
        mmin = np.nanmin(asts)
        offset = []
        for i in range(stars.shape[1]):
            mmax2 = np.nanmax(stars[:,i])
            mmin2 = np.nanmin(stars[:,i])
                
            offset.append(np.min([mmax - mmin2, mmin - mmax2]) - 0.01)
        
        return np.asarray(offset)
    
    def rescalMulti(asts, stars):
        std = np.nanstd(asts, axis = 0)
        med = np.nanmedian(asts, axis = 0)
        offset = []
        for i in range(stars.shape[1]):
            std2 = np.nanstd(stars[:,i])
            med2 = np.nanmedian(stars[:,i])
            
            offset.append(np.min(np.asarray([(med + 2*std) - (med2 - 2*std2), (med - 2*std) - (med2 + 2*std2)]), axis = 0))
        
        return np.asarray(offset)
    
    def rotate_image(image, angle):
        image_center = tuple(np.array(image.shape[1::-1]) / 2)
        rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
        result = cv2.warpAffine(image.astype(float), rot_mat, image.shape[1::-1], flags=cv2.INTER_LINEAR)
        return result
    
    
    def rtod(r):
        return r*180 / np.pi
    
    def mad(tab):
        """Median std"""
        med = np.median(tab)
        return np.sqrt(np.median((tab - med)**2))
    
    def hist(tab):
        hist = np.zeros((2**16))
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                hist[tab[i,j]] += 1
        
        return hist
    
    def histCum(tab):
        hist = Utils.hist(tab)
        hc = np.zeros((2**16))
        for i in range(hist.shape[0]):
           if i == 0: 
               hc[0] = hist[0]
           else:
               hc[i] = hc[i-1] + hist[i]
        
        return hc
    
    def findStars(image, tresh = None):
        image = image.astype(np.uint16)
       
        if tresh == None:
            tresh =  1.5*np.median(image)
        
        (thresh, im_bw) = cv2.threshold(image, tresh, 2**8-1, cv2.THRESH_BINARY)
        
        im_bw = cv2.medianBlur(im_bw, 5)
     
        contours, hierarchy = cv2.findContours(np.uint8(im_bw),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        center = []
        for c in contours:
            try:
               M = cv2.moments(c)
               cX = int(M["m10"] / M["m00"])
               cY = int(M["m01"] / M["m00"])
               center.append(np.array([cX, cY]))
        
            except:
                pass
        return np.asarray(center)
    
    def linearBtwTwoPoints(p1, p2):
        '''
        compute linear coeficients a and b for array of points.
        Parameters
        ----------
        p1 : 2d arrays 
            DESCRIPTION. list of points 
        p2 : 2d arrays
            DESCRIPTION. list of points

        Returns
        -------
        list of coef a and b evaluated between p1 and p2 

        '''
        
        dif = p2 - p1
        
        a = dif[:,1] / dif[:,0]
        b = -a*p1[:,0] + p1[:,1]
        
        return a, b
    
    def dotProdWithSharedStart(start, end1, end2):
        """
        Compute the dot product of the vectors pointing from start to end1 and end2

        Parameters
            start: starting point of both vectors
            end1: end point of first vector
            end2: end point of second vector

        Returns
            dot product of (end1 - start) and (end2 - start)
        """
        return (end1[0] - start[0]) * (end2[0] - start[0]) + (end1[1] - start[1]) * (end2[1] - start[1])
    
    def slidingAvg(arr, win, axis = 0):
        
        if len(arr.shape) == 1:
            arr = arr.reshape(arr.shape[0],1)
        
     
        if axis == 0: 
            avg = np.zeros((arr.shape[0] - win + 1, arr.shape[1] ))
            orientation = 1
        elif axis == 1:
             orientation = 0
             avg = np.zeros((arr.shape[1]  - win + 1, arr.shape[0] ))
        elif axis > 1:
            print('can only treat 1d or 2d arr')
            return
            
                
        for i in range(orientation):
            if axis == 0:
                avg[:, i] = np.convolve(arr[:, i], np.ones(win), 'valid') / win
            elif axis == 1:
                avg[:, i] = np.convolve(arr[i, :], np.ones(win), 'valid') / win
            
        return avg
        
    def stojd(time):
        return time / 86400
        