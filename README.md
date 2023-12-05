# Steroid

Steroid is a Python package designed to assist users in creating their own automated or semi-automated procedures for correcting frames and/or conducting photometry. It is tailored to facilitate the development of efficient workflows in these domains.

  ## Table of contents
  1. [Code Structure](#code-structure)
  2. [Classes Description](#classes-description)
     1. [Corrector](#corrector)
        1. [Description](#corrector-description)
        2. [Constructor](#corrector-constructor)
        3. [Methods](#corrector-methods)
     2. [Detector](#detector)
        1. [Description](#detector-description)
        2. [Constructor](#detector-constructor)
        3. [Methods](#detector-methods)
  3. [Data structure](#datastruct)
     1. [Triangle](#datastruct-triangle)
        1. [Description](#datastruct-triangle-description)
        2. [Constructor](#datastruct-triangle-constructor)
        3. [Methods](#datastruct-triangle-methods)
     2. [Pattern](#datastruct-pattern)
        1. [Description](#datastruct-pattern-description)
        2. [Constructor](#datastruct-pattern-constructor)
        3. [Methods](#datastruct-pattern-methods)
     3. [SeqManager](#datastruct-seqmanager)
        1. [Description](#datastruct-seqmanager-description)
        2. [Constructor](#datastruct-seqmanager-constructor)
        3. [Methods](#datastruct-seqmanager-methods)
     4. [Appertures](#datastruct-appertures)
        1. [Description](#datastruct-appertures-description)
        2. [Constructor](#datastruct-appertures-constructor)
        3. [Methods](#datastruct-appertures-methods)
     5. [Fit](#datastruct-fit)
        1. [Description](#datastruct-fit-description)
        2. [Constructor](#datastruct-fit-constructor)
        3. [Methods](#datastruct-fit-methods)
  4. [Photometry](#photometry)
     1. [Description](#photometry-description)
     2. [Methods](#photometry-methods)
     3. [How to use](#photometry-howtouse)
 
##Code Structure <a name="code-structure"></a>

##Classes Description <a name="classes-description"></a>


  In this section, we will elaborate on all the classes within the code, offering insights into their functionalities. Additionally, we will provide details on certain methods that may prove useful for users.
  
  <details>

  <summary id="corrector"> 
    
  ### Corrector <a name="corrector"></a>
  
  </summary>
  

  **Description:** <a name="corrector-description"></a>

  
  Astronomical images from the same sequence are seldom perfectly aligned with each other. It is common to observe, at the very least, a drift in both the x and y directions between each image, and at worst, a field rotation. This misalignment can stem from       various sources, with primary factors including the type of telescope mount (equatorial or azimuthal), mechanical quality, presence or absence of guiding, alignment issues, meridian flips, and more.

  The "Corrector" class is specifically crafted to estimate the drift and rotation angle between each image in a given sequence. It offers several functionalities to the user, allowing for either a straightforward position correction or direct image correction.   In the context of photometric studies, it is advisable not to correct the images directly. This is because, due to the discrete nature of an image, rotation can introduce undesirable artifacts. For amateur astronomers engaged in astrophotography, there is an    option to directly correct the image and apply interpolation to produce an image without visible artifacts.

  **Constructor:** <a name="corrector-constructor"></a>
  
  ***Corrector(seqManager, flatSeq = None, biasSeq = None, darkSeq = None, exposurKey = None):***
   
  The constructor of the Corrector class takes the following inputs:

  -  (mandatory) A sequence of images (refer to the data structure SeqManager).
  -  (optional) A list of paths (string list) for the flat sequence, a list of paths (string list) for the bias sequence, a list of paths (string list) for the dark sequence, and a string corresponding to the FITS header key for the exposure (usually, exposure       is stored in the FITS header under the key EXPOSURE or EXPTIME).


  **Methods:** <a name="corrector-methods"></a>
  
  
   ***getImgShape(idx = 0, idx_HDU = 0):*** 
 
  -  Description: Returns the shape of an image in the sequence.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (tuple)

  ***getImgCenter(idx_img = 0, idx_HDU = 0):***

  -  Description: Returns the coordinates of the center of an image in the sequence.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (tuple)

  ***getImg(idx = 0):***

  -  Description: Returns an object of type Fit (refer to the data structure Fit).
  -  Input: (INT) Index of the image in the sequence.
  -  Return: (Fit)

  ***getData(idx = 0, HDU = 0)***

  -  Description: Returns the data of the raw image.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (numpy.array)

  ***getReducedData(idx = 0, HDU = 0)***
  
  -  Description: Returns the data of the reduced image.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (numpy.array)
  
  ***getHeader(idx = 0, HDU = 0):***
  
  -  Description: Returns the header of the image at the specified index in the sequence.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (STRING)

  ***histogram(idx = 0 , idx_HDU = 0):***
  
  -  Description: Returns the histogram of the image at the specified index in the sequence and at the specified HDU index.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (numpy.array) Histogram, (numpy.array) Bin edges (refer to numpy.histogram)

  ***getStarsListOfImg(idx):***

  -  Description: Returns the positions of stars detected in the image at the specified index.
  -  Return: (numpy.array) Star positions

  ***computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  Description: Computes the drift and the angle of rotation for each image in the sequence and stores them in two lists with lengths equal to the sequence length.
  -  Input: (FLOAT) An offset that can be added to adjust the threshold value. (BOOLEAN) Indicates whether the threshold should be estimated on reduced images or raw images.
    
  ***medDrif(idx):***

  -  Description: Drift is estimated between all detected stars. This function returns the median value of the image at the specified index. **NEEDS TO EXECUTE ***computeImagesCorrection*** FIRST**.
  -  Input: (INT) Index of the image in the sequence.
  -  Output: (array) 2D array of drift in both axes.


  ***avgDrif(idx):***

  -  Description: Drift is estimated between all detected stars. This function returns the average value of the image at the specified index. **NEEDS TO EXECUTE ***computeImagesCorrection*** FIRST**.
  -  Input: (INT) Index of the image in the sequence.
  -  Output: (array) 2D array of drift in both axes.

  
   ***medAng(idx):***

  -  Description: Angle is estimated between all detected stars. This function returns the median value of the image at the specified index. **NEEDS TO EXECUTE ***computeImagesCorrection*** FIRST**.
  -  Input: (INT) Index of the image in the sequence.
  -  Output: (FLOAT) Angle of rotation.

  ***avgAng(idx):***

  -  Description: Angle is estimated between all detected stars. This function returns the average value of the image at the specified index. **NEEDS TO EXECUTE ***computeImagesCorrection*** FIRST**.
  -  Input: (INT) Index of the image in the sequence.
  -  Output: (FLOAT) Angle of rotation.

  ***correctStarsFromRot(arrayToCorrect, idx, coefMultAngle = -1)***

  -  Description: According to a given array of positions, this function corrects each position based on the drift and angle of the image at the specified index. The coefficient coefMultAngle takes values of 1 or -1, determining the       direction of rotation. Different values from 1 or -1 will influence the angle of rotation.
  -  Input: (2D array) Array of positions to correct, (INT) Index of the image to correct for, (INT) Coefficient to multiply to the angle.
  -  Return: (2D array) New positions of objects located at positions in arrayToCorrect, adjusted according to the drift and angle of the image at the specified index.

  ***correctedImg(idx = 0, HDU_idx = 0):***

  -  Description: Returns the corrected image.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (numpy.array) The corrected image.

  ***getSuperImg(idx_ims = None, HDU_idx = 0):***

  -  Description: Returns the average combination of all images in the sequence after correction.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the HDU in the image.
  -  Return: (numpy.array) The combined image.

  ***rejectBadData():***

  -  Description: Rejects all data where the drift and/or angle were not found.

 ***imshowstar(idx = 0):***

  -  Description: Method to display the image at the specified index and show detected objects.
  -  Input: (INT) Index of the image in the sequence.
 
 ***checkPatterns(idxOfImage = 0, patidx = None):***

  -  Description: Method to display the image at the specified index. If the pattern index (patidx) is set to None, it will show all patterns. If set to a specific value, it will only display the selected pattern.
  -  Input: (INT) Index of the image in the sequence, (INT) Index of the pattern. If set to None, it will show all patterns.


 </details id="detector">

  <details>

  <summary> 
    
  ### Detector <a name="detector"></a>
  
  </summary>
 


**Description:** <a name="detector-description"></a>

This class is designed for detecting moving objects. It internally maintains a list of positions for these objects and another list for their speeds along both the x and y axes. With the provision of initial positions and speeds, the task of determining the positions of moving objects in each frame becomes straightforward.


**Constructor:** <a name="detector-constructor"></a>

***Detector(imageSeq, flatSeq = None, biasSeq = None, darkSeq = None):***

  The constructor of the Detector class takes the following inputs:

  -  (mandatory) A list of paths (string list) for the main image sequence.
  -  (optional) A list of paths (string list) for the flat sequence, a list of paths (string list) for the bias sequence, and a list of paths (string list) for the dark sequence.

**Methods:** <a name="detector-methods"></a>

***computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  Description: Calls the function computeImagesCorrection from the Corrector class and subsequently rejects all data with missing drift and/or angle values.
  -  Input: (FLOAT) An offset that can be added to adjust the threshold value. (BOOLEAN) Indicates whether the threshold should be estimated on reduced images or raw images.

***findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2):***

  -  Description: Finds slow-moving objects based on a method that compares the presence of an object from one of the first frames to one of the last frames in the sequence. To prevent excessive detection, this algorithm is assisted by a convolutional neural network based on AlexNet. The method store a list of moving object positions on the initial frame and a list of object speeds along the x and y axes.
  -  Input: (FLOAT) An offset that can be added to adjust the threshold value. (BOOLEAN) Indicates whether the threshold should be estimated on reduced images or raw images. (INT) Epsilon, which corresponds to the tolerance.

***fasterAst():***

  -  Description: Returns the index of the fastest moving object.
  -  Return: (INT) The index of the fastest asteroid in the list.

***slowestAst():***

  -  Description: Returns the index of the slowest moving object.
  -  Return: (INT) The index of the slowest asteroid in the list.

***getAstPositionAtImg(idx):***

  -  Description: Returns the position of moving objects in the image at the specified index.
  -  Input: (INT) Index of the image from which to retrieve positions.
  -  Return: (numpy.array) Array of moving objects positions.

***nofa():***

  -  Description: Returns the number of detected moving objects.
  -  Return: (INT) Number of detected moving objects.

***astSpeed(idx = 0):***

  -  Description: Returns the speed of moving objects in the image at the specified index.
  -  Input: (INT) Index of the image from which to retrieve speeds.
  -  Return: (numpy.array) Array of speeds on the x and y axes of moving objects.


</details>

## Data structure description: <a name="datastruct"></a>

This section is dedicated to discussing some classes stored in the file "data_structures.py"

<details>

  <summary> 
    
  ### Triangle<a name="datastruct-triangle"></a>
  
  </summary>



**Description:** <a name="datastruct-triangle-description"></a>

This class stores information about three stars and represents a triangle. It overloads various operations, including addition, subtraction, division, multiplication, and comparison operations. Additionally, it implements the \_\_str\_\_ method for string representation.

**Constructor:**<a name="datastruct-triangle-constructor"></a>

***Triangle(s1, s2, s3, eps = 2):***

S1, s2, and s3 are represented as numpy arrays. The eps parameter is a tolerance used in the \_\_eq\_\_ method to determine whether two triangles are considered equal or not.


**Methods:**<a name="datastruct-triangle-methods"></a>

***d1():***

-  Description: Returns the Euclidean distance between s1 and s2.
-  Return: (FLOAT) Distance between s1 and s2.

***d2():***

-  Description: Returns the Euclidean distance between s1 and s3.
-  Return: (FLOAT) Distance between s1 and s3.

***d3():***

-  Description: Returns the Euclidean distance between s2 and s3.
-  Return: (FLOAT) Distance between s2 and s3.

***buildVect():***

-  Description: Builds tree vectors v1, v2, and v3 between (s1, s2), (s1, s3), and (s2, s3).
-  Return: (numpy.array, numpy.array, numpy.array) Three vectors v1, v2, and v3.

***getRotationAngle(other):***

-  Description: Computes the angle between the triangle and another one. Caution! This method does not check if both triangles are the same.
-  Input: (Triangle) Another triangle.
-  Return: (FLOAT) The angle of rotation between both triangles.


***computeDistance(other):***

-  Description: Returns the mean distance between the triangle and another one.
-  Input: (Triangle) Another triangle.
-  Return: (numpy.array) Mean distance in x and y of both triangles.

***correctRot(angle, center):***

-  Description: Rotates the positions of s1, s2, and s3 by an angle around a specified center of rotation.
-  Input: (FLOAT) Angle of rotation, (numpy.array) Position of the center of rotation.

</details>

<details>

  <summary> 
    
  ### Pattern: <a name="datastruct-pattern"></a>
  
  </summary>



**Desciption:**<a name="datastruct-pattern-description"></a>

This class stores triangles as a pattern. It overloads addition, subtraction, multiplication, division, comparison operations, and the \_\_str\_\_ method.

**Constructor**<a name="datastruct-pattern-constructor"></a>

***Pattern(t1, t2, t3, t4, t5):***


t1, t2, t3, t4, and t5 are objects of the Triangle class (refer to the data structure class *Triangle*).

**Methods:**<a name="datastruct-pattern-methods"></a>

***computeDistance(other):***

-  Description: Computes the mean distance between two patterns.
-  Input: (Pattern) Another pattern to compute the distance.
-  Return: (numpy.array) Mean distance in x and y between the two patterns.

***computeAngle(other):***

-  Description: Computes the angle of rotation between two patterns.
-  Input: (Pattern) Another pattern.
-  Return: (FLOAT) Angle of rotation between the two patterns.


***correctRot(angle, center):***

-  Description: Rotates t1, t2, t3, t4, and t5 by a specified angle around a center of rotation.
-  Input: (FLOAT) Angle of rotation, (numpy.array) Position of the center of rotation.

</details>

<details>

  <summary> 
    
  ### SeqManager:<a name="datastruct-seqmanager"></a>
  
  </summary>



**Description:**<a name="datastruct-seqmanager-description"></a>

This class stores a list of image paths from the same sequence.

**Constructor:**<a name="datastruct-seqmanager-constructor"></a>

***SeqManager(seq):***

The seq attribute represents a list of paths to raw images (STRING) in this class.

**Methods:**<a name="datastruct-seqmanager-methods"></a>

***getPath(idx):***

-  Description: Returns the path of the image at the specified index.
-  Input: (INT) Index of the image of interest in the sequence.
-  Return: (STRING) The path of the image.

***getFileName(idx):***

-  Description: Returns the name of the image at the specified index.
-  Input: (INT) Index of the image of interest in the sequence.
-  Return: (STRING) The name of the image.

***getImg(idx = 0):***

-  Description: Returns an object of type Fit for the image at the specified index.
-  Input: (INT) Index of the image of interest.
-  Return: (Fit) A data structure of type Fit.

***getHDU(idx = 0, HDU = 0):***

-  Description: Returns the HDU (Header Data Unit) of the image at the specified index.
-  Input: (INT) Index of the image, (INT) HDU index.
-  Return: (astropy.io.fits.hdu.image.PrimaryHDU) HDU of the image at the specified index.

***getInfo(idx = 0):***

-  Description: Displays information about the image at the specified index.
-  Input: (INT) Index of the image of interest.

***getHeader(idx = 0, HDU = 0):***

-  Description: Returns the header at the specified HDU index for the image at the given index.
-  Input: (INT) Index of the image, (INT) Index of the HDU of the image at the specified index.
-  Return: (astropy.io.fits.header.Header) Header of the image at the specified HDU index.

***getExpo(idx, key, HDU = 0):***

-  Description: Returns the exposure from the header of the image at the specified index and HDU. The exposure is determined according to the provided key.
-  Input: (INT) Index of the image of interest, (STRING) Key in the header corresponding to the exposure, (INT) HDU index.
-  Return: (FLOAT) Exposure.

***getData(idx = 0, idx_HDU = 0):***

-  Description: Returns the image at the specified index and HDU as an array.
-  Input: (INT) Index of the image of interest, (INT) Index of the HDU.
-  Return: (numpy.array) The image.

***getCenter(idx_img = 0, idx_HDU = 0):***

-  Description: Returns the center coordinates of an image at the specified index and HDU.
-  Input: (INT) Image index, (INT) Image HDU.
-  Return: (numpy.array) Coordinates of the center of the image at the specified index and HDU.

***getImgShape(idx = 0, idx_HDU = 0):***

-  Description: Returns the shape of the image at the specified index and HDU.
-  Input: (INT) Index of the image of interest, (INT) HDU index.
-  Return: (TUPLE) Image shape.

***getTime(key, forma, idx = 0, HDU = 0):***

-  Description: Gets the time of the image at the specified index and HDU from the header using the provided key and format. If the time is stored as Julian Day in the header (e.g., JD=2458780), set the key and format to JD. For more formats, refer to Time.FORMATS from astropy.time.
-  Input: (STRING) Key of the time in the header, (STRING) Format of the time in the header (refer to Time.FORMATS from astropy.time), (INT) Index of the image, (INT) Index of the HDU.
-  Return: (astropy.time.core.Time) Time of the image.

***pop(idx = -1):***

-  Description: Deletes an image at the specified index from the sequence. By default, if idx is set to -1, the last image is deleted.
-  Input: (INT) Index of the image to delete. (Default: -1 to delete the last image.)

***histogram(idx = 0, idx_HDU = 0):***

-  Description: Returns the histogram of the image at the specified index and HDU.
-  Input: (INT) Index of the image of interest, (INT) Index of the HDU.
-  Return: (numpy.array) Histogram, (numpy.array) Bin edges (refer to numpy.histogram).

</details>

<details>

  <summary> 
    
  ### Appertures:<a name="datastruct-appertures"></a>
  
  </summary>

**Description:** <a name="datastruct-appertures-description"></a>

This data structure is dedicated to managing apertures. It takes as input a 2D numpy array of aperture positions with aperture sizes and can handle photometry.

**Constructor:** <a name="datastruct-appertures-constructor"></a>

***Appertures(positions, idxOfStars = None, r = 3, ri = 6, re = 8):*** 

- Positions: 2D numpy array of positions of apertures for all objects. The first rows should represent targets, and the last rows should represent reference stars for differential photometry if needed.
  - `idxOfStars`: (INT) Index of the row in positions where the positions of reference stars' apertures are stored.
  - `r`: (FLOAT) Inner radius of the apertures.
  - `ri`: (FLOAT) Radius of the dead area of the apertures.
  - `re`: (FLOAT) Radius of the background aperture.

**Methods:** <a name="datastruct-appertures-methods"></a>

***photom(img, key, forma, center = False, exposure = None):***

- Description: Perform the photometry and allow users to center the time at the midpoint of exposure if the time in the header is set at the beginning of exposure.
- Input:
  - `img`: (FIT) FIT object of the image used for photometry
  - `key`: (STRING) Keyword of the time in the header
  - `forma`: (STRING) Format of the time in the header
  - `center`: (BOOLEAN) Set to true to center the time in case if the time in header was taken at the beginning of exposure
  - `exposure`: (FLOAT) Exposure time
- Output: (astropy.table.table.QTable) Resume of the photometry
  
</details>


<details>

  <summary>

  ### Fit: <a name="datastruct-fit"></a>
    
  </summary>


  **Description:** <a name="datastruct-fit-description"></a>

  This structure is dedicated to managing FIT images, and numerous methods are implemented to handle various operations on images.

  **Constructor:** <a name="datastruct-fit-constructor"></a>

  ***Fit(path, dark = 0, flat = 1, bias = 0, darkExp = None, exposurKey = None):***

  -  `path`: (STRING) Path of the image in the user's system
  -  `dark`: (NUMPY.ARRAY) Master dark
  -  `flat`: (NUMPY.ARRAY) Master flat
  -  `bias`: (NUMPY.ARRAY) Master bias
  -  `darkExp`: (FLOAT) Exposure of dark images
  -  `exposureKey`: (STRING) The key in the header where exposure is stored

  **Methods** <a name="datastruct-fit-methods"></a>

  ***getHDU(i = 0):***

  -  Description: Get HDU of the image.
  -  Input: (INT) Index of the HDU to get.
  -  Return: (astropy.io.fits.hdu.image.PrimaryHDU)
  
  ***getInfo():***

  -  Description: Print information of the image.

  ***getHeader(HDU = 0):***

  -  Description: Get the header of the HDU.
  -  Input: (INT) HDU index.
  -  Return: (astropy.io.fits.header.Header)

  ***getExposure(self, key, HDU = 0):***

  -  Description: Get the exposure.
  -  Input: (STRING) Key in the header corresponding to the exposure, (INT) Index of the HDU of interest.
  -  Output: (FLOAT)

  ***getTime(key, forma, HDU = 0):***

  -   Description: Get the time from the header.
  -   Input: (STRING) Key in the header corresponding to the time, (STRING) Format of the time stored in the header, (INT) Index of the HDU of interest.
  -   Output: (FLOAT)

  ***getData(idx_HDU = 0):***

  -  Description: Get the image data as a matrix.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (NUMPY.ARRAY) Matrix of the image.

  ***getReducedData(HDU = 0):***

  -  Description: Get the reduced image data as a matrix.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (NUMPY.ARRAY) Matrix of the reduced image.

  ***getCenter(idx_HDU = 0):***

  -  Description: Return the center of the image.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (TUPLE)

  ***getShape(idx_HDU = 0):***

  -  Description: Get the shape of the image.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (numpy.ndarray)

  ***getTresh(reduced = False, display = False):***

  -  Description: Method to automatically determine the best threshold value to binarize the image.
  -  Input: (BOOLEAN) If set to true, will evaluate threshold on reduced frame. (BOOLEAN) If set to true, will plot information to help debug.
  -  Output: (FLOAT) Threshold value.

  ***findStars(tresh = None, onReduced = False):***

  -  Description: Find all objects (not only stars) present on frames.
  -  Input: (FLOAT) Threshold value. If set to None, will be set to 1.5 times the median. (BOOLEAN) If set to true, will find objects on reduced frame.
  -  Output: (NUMPY.ARRAY) x, y coordinates of object centers in the frame.

  ***histogram(idx_HDU = 0):***

  -  Description: Compute histogram of the HDU of interest.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (TUPLE(NUMPY.ARRAY, NUMPY.ARRAY)) The first array corresponds to the histogram, and the second to the bin_edges (see [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)).

  ***reducedHistogram(idx_HDU = 0):***

  -  Description: Compute histogram of the HDU of interest on the reduced frame.
  -  Input: (INT) Index of the HDU of interest.
  -  Output: (TUPLE(NUMPY.ARRAY, NUMPY.ARRAY)) The first array corresponds to the histogram, and the second to the bin_edges (see [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)).

</details>


## Photometry:<a name="photometry"></a>


### Description:<a name="photometry-description"></a>

In the context of our studies on asteroids from the main belt, we have utilized the various functions of Steroid to develop our own photometric software capable of automatic or semi-automatic processing. With this new tool, we have significantly enhanced our ability to generate light curves. This software is integrated into Steroid in the `photometry.py` file. In this section, we will demonstrate how to use it.

<details>

  <summary> 
    
  ### Methods:<a name="photometry-methods"></a>
  
  </summary>


**Constructor:**<a name="photometry-methods-constructor"></a>

***Photometry(detector = None)***

*Photometry* takes only one optional parameter of type *Detector*. Why optional? Because *Photometry* also includes functions to save photometry but also functions to load. If users want to rework on some light curves already processed, they don't need to redo all the work. *Photometry* can reload previous light curves. In this kind of situation, the user doesn't need any *Detector* as the photometry was already done. They just need to build an empty *Photometry* object and use the method ***readCsv(path)***.


**Methods:** <a name="photometry-methods-methods"></a>

***start(nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000)***

- Description: Launch the photometry according to some input parameters.
- Input: 
  - (INT) Number of reference stars (only in case of an automatic procedure).
  - (BOOLEAN) Center or not appertures of the center of brightness.
  - (FLOAT) Maximum value that automatically selected reference stars should not overstep.
  - (FLOAT) The threshold to detect stars in the context of stars passages.

***plotDif(refS = 0, ast = -1, yRange = None, binning = 1, resc = True, forma = 'jd', xtick = None, inMag = True, rmExtremPoint = False, cStd = 2, deg = 4, displayRmFit = False, starPassage = False, markerSize = 100, lineWidths = 5)*** 

- Description: Perform a plot of differential photometry.
- Input:
  - (INT) `refS` is the index of the star selected as reference.
  - (INT) `ast` is, in the case of multiple asteroids, the index of the asteroid that we want to plot. If set to -1, all asteroids will be plotted.
  - (list) `yRange` range of the y-axis.
  - (INT) `binning`. Use to bin the light curve. Automatically chosen if set to -1.
  - (BOOLEAN) `resc`. Rescale stars' light curves close to the asteroid's light curves.
  - (STRING) `forma`. Format of the time. Refer to `Time.FORMATS` from `astropy.time`.
  - (array) `xticks`. New x ticks.
  - (BOOLEAN) `inMag`. If True, the y-axis displays in magnitude. If False, the y-axis displays in instrumental flux.
  - (BOOLEAN) `rmExtremPoint`. If True, will remove extreme points. To remove extreme points, the algorithm will fit a polynomial, then normalize asteroid's light curves with the polynomial. Each point out of [median - C x Std, median + C x Std] are removed.
  - (FLOAT) `cStd`. This corresponds to C.
  - (INT) `deg`. Degree of the polynomial.
  - (BOOLEAN) `displayRmFit`. If True, display more plots to monitor `rmExtremPoint`.
  - (BOOLEAN) `starPassage`. If True, will remove star's passages.
  - (INT) `markerSize`. Corresponds to the size of the marker.
  - (INT) `lineWidths`. Corresponds to the thickness of the marker.

***toDat(path, filename, binning = 1, forma = 'mjd', refS = -1, deg = 4, cStd = 2, displayRmFit = False)***


markdown
Copy code
-  Description: Write files with extension .1, .2, .3, and .4. For each of them, the first column is the time. For others columns, .1 corresponds to data in instrumental flux, .2 corresponds to data in magnitude, .3 corresponds to differential photometry, and .4 corresponds to differential photometry with averaged reference stars.
-  Input: 
   - `path`: (STRING) Path where to save those files.
   - `filename`: (STRING) Name to give to files.
   - `Other` parameters are the same as ***plotDif***.

***log(path, name = "log.txt")***
  
-  Description: Write a log file with information on data rejected, star passages data, FWHM detected on each frame, etc.
-  Input: 
   - `path`: (STRING) Path to save the log file.
   - `name`: (STRING) Name given to the log file. Don't forget the extension.


***toGif(path)***

-  Description: Write a .gif image of all frames with apertures.
-  Input: (STRING) Path + file name with extension (e.g., r"C:/.../myGif.gif").

***toCsv(path)***

-  Description: Write a CSV file summarizing the photometry. It can be used as a backup with the method ***readCsv*** (see below).
-  Input: (STRING) Path + file name with extension (e.g., r"C:/.../myCsv.csv").
  
***readCsv(path)***   

-  Description: Load a CSV file produced with the method ***toCsv***. Can be used to rework plots.
-  Input: (STRING) Path + file name with extension (e.g., r"C:/.../myCsv.csv").

***showAp(idx)***

-  Description: Display the image at the specified index and show aperture positions.
-  Input: (INT) `idx`: Index of the image in the image sequence.
  
***checkBox(ofs)***

-  Description: Display one of the first images of the sequence and show asteroid positions at the beginning and at the end of the sequence, along with box vertices and all objects detected inside the boxes (except asteroids).
-  Input: (FLOAT) `ofs`: Offset added to the threshold to detect objects. To debug star passages, it should be identical to the starPassageOfs parameter of the *start* method.

***log(path, name = "log.txt)***

- description: save a log.txt file containing informations on images not took into account because of bad detection, FWHM computed on each images and also images where star passages were detected.
- input: (STRING) `path` of the directory (exemple: "/my/dir/"), (STRING) `name` of the file. by default it is set to "log.txt"

</details>


<details>

  <summary> 
    
  ### How to use ? <a name="photometry-howtouse"></a>
  
  </summary>




**First step:**

The first thing to use photometry is to import it:

    from photometry import Photometry

*Photometry* object constructor take, as an optional parameter, an object *Detector*. 

If the photometry was not yet done, then the users need to provide an object *Detector*. It is therefore important to include it:

    from detector import Detector


**Second step:** 

The next step is to build an object *Detector* (constructor is describ in the [Detector](#detector) section). To do this, we will use glob

     import glob

an exemple of a piece of code that can be use to build an object *Detector*:

  ~~~

#-----------set up all list of path for the raw data and bias, dark and flat data--------------

    path = r"C:\...\directory_of_your_data/"
   
    seq = glob.glob(path + "*target_repetable_name_pattern*.f*t*") #return a list of path of all file which contain target_repetable_name_pattern in their name
 
    dark = glob.glob(path + "*dark_repetable_name_pattern*.f*t*") #return a list of path of all file which contain dark_repetable_name_pattern in their name
    flat = glob.glob(path + "*flat_repetable_name_pattern*.f*t*") #return a list of path of all file which contain flat_repetable_name_pattern in their name
    bias = glob.glob(path + "*bias_repetable_name_pattern*.f*t*") #return a list of path of all file which contain bias_repetable_name_pattern in their name

#-----------------------check if bias, dark and flat data was found-----------------------------

    if len(dark) == 0:        # Check if dark data was found. If not, set up dark variable to the optional default value of Detector attribute darkSeq 
        dark = None
        print('DARK EMPTY')
    if len(bias) == 0:        # Check if bias data was found. If not, set up dark variable to the optional default value of Detector attribute biasSeq 
        bias = None
        print('BIAS EMPTY')
    if len(flat) == 0:        # Check if flat data was found. If not, set up dark variable to the optional default value of Detector attribute flatSeq 
        flat = None
        print('FLAT EMPTY')

#-------------------construct Detector object----------------------------------------------------

    d = Detector(seq, flatSeq = flat, biasSeq = bias , darkSeq = dark)
  ~~~

**Third step:**

The next step is to launch *Detector* methods to correct images and to detect asteroids:

    d.computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)
    d.findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2)

***computeImagesDrift*** could me internally called in ***findAsteroid*** but we choose to let it like this to give more flexibility to users. Indeed, In this process, ***computeImagesDrift*** will take more time has it has to detect stars from all images of the sequence. Also it's the proccesse the less sensitive to the detection treshold. Indeed, it's only need 5 common stars on each frames to be able to correct images. According to this, users can earn time of execution setting up a high value of **offsetTreshStarsDetection**. For the same reason, **treshOnReduced** can be set to False.  

On the other hand, ***findAsteroid*** is really sensityve to the data quality and to the treshold. More close to the optimal value the treshold will be and better the algorythm will perform. Therefor **offsetTreshStarsDetection** should be small for small adjustement and **treshOnReduced** should be set to True.

**Fourth step**

The next step is to build an object *Photometry* and to launch the photometric processe.

    phot = Photometry(d)
    phot.start(nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000)

**nbOfStars** is mandatory for now but it's only used in the case that you choose automatic reference stars selection. It's the number of reference stars that the code will search. 

**center** is a boolean and is a paramter to allow the code to center appertures on the "center of intensity" (in reference to the center of masse equation where the masse was changed by the intensity of pixels) or not. To be clear, apperture position will all time be set according to the initial positions, to the image drift and angle of rotation and for moving object, to the speed. But, if center is set to true, after placing all appertures, the code will simply perform a centring. In case of starpassages, the apperture will probably stay fixe on the stars the time that the astroid will pass in front but after, it will come back centred on the asteroid when the star will leave the apperture feild. Moreover, with the algorythm set up to delete stars passages, all the time where the apperture will stay focuse on the star will not be present on the final lightcurve.

**maxVal** correspond to the maximum value that reference stars automaticly select should not overstep.

**starPassageOfs** have the same function as **treshOnReduced** for objects and methodes dedicated to detect objects. This one is dedicated to stars Passages detection. It's, by default, set to 15000, which detect only bright stars but can be set much lower to detect fainter stars.

**Fifth step:**

As the photometry is done, the next step is to save results. To save plots, when a plot is done using ***plotDif***, a image of it can be save using the button on the window of the plot.
To save the data, the best function is probably ***ToCsv*** which will save a csv file containing all infomations about the photometry, plus, information on stars passages if detected.
The CSV file could be use as a backup using the function ***readCsv*** which will load in a **Photometry** the photometry previously done. 
The function ***log*** will build a txt file containing additional informations.


    phot.toCsv("/path/of/your/directory/myCsv.csv") #save CSV file

in an other script (or in the same), to reload a CSV file:

    phot = Photometry()  # this line is only needed if you didn't build any Photometry object before (in a new script for exemple)
                         # no need detector object as here we are just loading previous photometry. 
                         
    phot.readCsv("/path/of/your/directory/myCsv.csv")



**Conclusion**

The final code should look like this:

 ~~~

from photometry import Photometry
from detector import Detector

import glob

#-----------Set up all list of path for the raw data and bias, dark and flat data--------------

path = r"C:\...\directory_of_your_data/"
   
seq = glob.glob(path + "*target_repetable_name_pattern*.f*t*") 

dark = glob.glob(path + "*dark_repetable_name_pattern*.f*t*") 
flat = glob.glob(path + "*flat_repetable_name_pattern*.f*t*")
bias = glob.glob(path + "*bias_repetable_name_pattern*.f*t*") 


if len(dark) == 0:      
   dark = None
   print('DARK EMPTY')
if len(bias) == 0:        
   bias = None
   print('BIAS EMPTY')
if len(flat) == 0:        
   flat = None
   print('FLAT EMPTY')

#-------------------Construct Detector object----------------------------------------------------

d = Detector(seq, flatSeq = flat, biasSeq = bias , darkSeq = dark)

#----------------Estimate images correction and perform asteroids detection----------------------

d.computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)
d.findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2)

#---------------Construct Photometry object and launch the photometry----------------------------

phot = Photometry(d)
phot.start(nbOfStars = 3, center = True, maxVal = 30000, starPassageOfs = 15000)

phot.plotDif(refS = 0, ast = -1, yRange = None, binning = 1, resc = True, forma = 'jd', xtick = None, inMag = True, rmExtremPoint = False, cStd = 2, deg = 4, displayRmFit = False, starPassage = False, markerSize = 100, lineWidths = 5)
phot.toCsv("/path/of/your/directory/myCsv.csv")
 ~~~

if loading is need in another file:

    from photometry import Photometry

    phot = Photometry()  # this line is only needed if you didn't build any Photometry object before (in a new script for exemple)
                         # no need detector object as here we are just loading previous photometry. 
                         
    phot.readCsv("/path/of/your/directory/myCsv.csv")

**Caution:** in case where users use semi-automatic procedures, the selection on images are done with the left click and when selection is finish press the right click.


</details>
