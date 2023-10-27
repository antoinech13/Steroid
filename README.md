# Steroid

Steroid is a python package dedicated to help users to developpe their own automatic or semi-automatic procedure in way to correct frame and/or to performe photometry.

  ## Summary

  

  ## Code Structur

  ## Classes Description

  in this section, we will details all classes of the code and we will details some methods that could be usefull for users

  ### Corrector

  **Description:**

  Astronomical images from the same sequence are rarely aligned with each other. It is common to observe, at best, a drift in both the x and y directions between each image, and at worst, a field rotation. This misalignment can have several origins, but as the  main cause, we can note the type of telescope mount (equatorial or azimuthal), as well as the quality of the mechanics, the presence or absence of guiding, the alignment, a meridian flip, etc. The "Corrector" class is therefore aimed at estimating the drift and the rotation angle between each image in the same sequence. It also provides several functionalities to the user, allowing either simple position correction or direct image correction. In the case of photometric studies, it is preferable not to correct the images directly. Indeed, due to the discrete nature of an image, a rotation of it will introduce undesirable artifacts in the image. In the case of an amateur who simply wants to do astrophotography, it is possible to directly correct the image and apply interpolation to obtain an image without visible artifacts.

  **Constructor:**
  
  ***Corrector(seqManager, flatSeq = None, biasSeq = None, darkSeq = None, exposurKey = None):***
   
  The constructor of the class *Corrector* take, as input: 
  
  -  (mandotory) a sequence of images (see the data structure SeqManager)
  -  (optional) a list of path (list of string) for the flat sequence, a list of path (list of string) for the bias sequence, a list of path (list of string) for the dark sequence and a string which correspond to the fit header key of the exposure (usually exposure is store in fits header under the key EXPOSURE or EXPTIME)


  **Methods:**
  
  
   ***getImgShape(idx = 0, idx_HDU = 0):*** 
 
  -  description: return the shape of an image of the sequence
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (tuple)

  ***getImgCenter(idx_img = 0, idx_HDU = 0):***

  -  description: return the coordinate of the center of an image of the sequence
  -  input: (INT) idx_img of image in the sequence, (INT) idx_HDU in the image
  -  return: (tuple)

  ***getImg(idx = 0):***

  -  description: return an object of type Fit (see data structure Fit)
  -  input: (INT) idx of image in the sequence
  -  return: (Fit)

  ***getData(idx = 0, HDU = 0)***

  -  desciption: return the data of the raw image.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array)

  ***getReducedData(idx = 0, HDU = 0)***

  -  desciption: return the data of the reduced image.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array)
  
  ***getHeader(idx = 0, HDU = 0):***
  
  -  desciption: return the header of the image at idx in the sequence.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (STRING)

  ***histogram(idx = 0 , idx_HDU = 0):***
  
  -  desciption: return the histogram of the image at idx in the sequence and at the HDU idx.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) histogram, (numpy.array) bin edgesarray (see numpy.histogram)

  ***getStarsListOfImg(idx):***

  -  desciption: return the position of stars detected at the iamge idx.
  -  return: (numpy.array) stars position

  ***computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  description: compute the drift and the angle of rotation for each images of the sequence and store them in two list of lengh = to the sequence lengh
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw
    
  ***medDrif(idx):***

  - description: drift is estimated between all stars detected. this function return the median value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (array) 2d arry of drift in both axis

  ***avgDrif(idx):***

  - description: drift is estimated between all stars detected. this function return the average value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (array) 2d arry of drift in both axis

  
   ***medAng(idx):***

  - description: angle is estimated between all stars detected. this function return the median value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (FLOAT) angle of rotation

  ***avgAng(idx):***

  - description: angle is estimated between all stars detected. this function return the average value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (FLOAT) angle of rotation

  ***correctStarsFromRot(arrayToCorrect, idx, coefMultAngle = -1)***

  -   description: according to a given array of positions, this function will correct each position according to the drif and angle of the image idx. The coefMultAngle take 1 or -1 and only give the direction of rotation. (different value frome 1 or -1 will influence the angle of rotation)
  -   input: (2d array) array of position to correct, (INT) idx of image for which to correct, (INT) coeficient multiply to the angle
  -   return: (2d array) new position of objects located at positions arrayToCorrect, according to the image idx drift and angle.

  ***correctedImg(idx = 0, HDU_idx = 0):***

  -  desciption: return the image corrected
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) the corrected image

  ***getSuperImg(idx_ims = None, HDU_idx = 0):***

  -  desciption: return the average combination of all images of the sequence (after correction)
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) the combined image

  ***rejectBadData():***

  - description: reject all data where drift and/or angle was not found

 ***imshowstar(idx = 0):***

 -  description: method to display image at idx and show objects detected
 -  input: (INT) idx of the image in the sequence
 
***checkPatterns(idxOfImage = 0, patidx = None):***

-  description: method to display image idx and, if patterns idx (patidx) set to None, will show all patterns. If patterns idx set to a value, will only show the pattern selected
-  input: (INT) idxOfImage image idex in the sequence. (INT) patidx index of the patterns of all patterns of the image. If set to None, will show all patterns


### Detector

**Description:**

This class is dedicated to detect moving object. it's internaly stor a list of moving objects position and and other list of their speed along x and y axis. with inital poistions and speed, it's easy to determine the position of moving objects on each frames.


**Constructor:**

***Detector(imageSeq, flatSeq = None, biasSeq = None, darkSeq = None):***

 The constructor of the class *Detector* take, as input: 
  
  -  (mandotory) a list of path (list of string) of the main image sequence
  -  (optional) a list of path (list of string) for the flat sequence, a list of path (list of string) for the bias sequence and a list of path (list of string) for the dark sequence.

**Methods:**

***computeImagesDrift(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  description: call the function computeImagesCorrection from *Corrector* than reject all bad data (without drift or/and angle values detected)
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw

***findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2):***

  -  description: will find slow moving object based on method comparing present object or not from one of the first and one of the last frame of the sequence. To help to not overdetect to much, this algorythm is helped by a convolutional neural network based on AlexNet. This method will feed a list of moving object position on the initial frame and a list a object speed along x and y axis.
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw. (INT) epsilon which correspond to the tolerence

***fasterAst():***

  - description: return the index of the faster moving object
  - return: (INT) the index of the faster asteroid in the list

***slowestAst():***

  - description: return the index of the slowest moving object
  - return: (INT) the index of the slowest asteroid in the list

***getAstPositionAtImg(idx):***

  - description: return the position of moving objects on the image at the idx
  - input: (INT) idex of the image where to get positions
  - return: (numpy.array) array of moving objects position

***nofa():***

  - description: return the number of moving object detected
  - return: (INT) number of moving object detected

***astSpeed(idx = 0):***

  - description: return the speed of moving objects on the image at the idx
  - input: (INT) idex of the image where to get positions
  - return: (numpy.array) array of speed on x and y axis of moving objects

                                                           
