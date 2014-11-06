# File: makeFakeGalaxy.py

import numpy as np
import galsim
import pyfits as fits


def makeGalaxy(flux, gal, psfImage, galType='sersic', drawMethod='auto',
               trunc=10.0):
    """
    Function called by task to make galaxy images

    INPUTS:
    flux = calibrated total flux
    gal = galaxy parameters in record (np.recarray)
    psfImage = np.ndarray of psfImage
    galType = type of galaxy we want to make, this has to agree with what's in the
    record array, options now are: 'sersic' (single sersic),
    'dsersic' (double sersic), and
    'real' (for making galaxies from real HST images)
    All the necessary keywords need to be in the fits catalog,
    including maybe drawMethod and trunc...
    """

    if galType is 'sersic':
        return galSimFakeSersic(flux, gal, psfImage=psfImage, trunc=trunc,
                                drawMethod=drawMethod, returnObj=False)

    if galType is 'dsersic':
        (comp1, comp2) = parseDoubleSersic(flux, gal)
        return galSimFakeDoubleSersic(comp1, comp2, psfImage=psfImage,
                                      trunc=trunc, drawMethod=drawMethod,
                                      returnObj=False)

    if galType is 'real':
        # TODO: For real galaxies, we need to decide which to use: index in the
        # catalog or Object ID.  Now, I just use Index
        (real_galaxy_catalog, index) = parseRealGalaxy(gal)
        if index >= 0:
            index = index
            random = False
        else:
            index = None
            random = True
        return galSimRealGalaxy(flux, real_galaxy_catalog, index=index,
                                psfImage=psfImage, random=random,
                                returnObj=False, drawMethod=drawMethod)

def parseRealGalaxy(gal):

        # If an "index" column is presented, then use the index
        # If not, turn random = True
        try:
            index = gal['index']
        except KeyError:
            index = -1

        # Get the catalog name and the directory
        try:
            cat_name = gal['cat_name']
        except KeyError:
            raise KeyError('Can not find the name of the catlog')
        try:
            cat_dir = gal['cat_dir']
        except KeyError:
            cat_dir = None

        real_galaxy_catalog = galsim.RealGalaxyCatalog(cat_name, dir=cat_dir)

        return real_galaxy_catalog, index


def parseDoubleSersic(tflux, gal):
    """
    Parse the input total flux [tflux] and parameter record array [gal] into two
    parameter records for each component [comp1, comp2]
    """

    # Check if this is a real 2-Sersic record
    try:
        frac1 = float(gal['b2t'])
    except KeyError:
        raise KeyError("No b2t parameter is found in the record!!")
    # Make sure the flux fraction of the first component is reasonable
    if (frac1 <= 0) or (frac1 >=1):
        raise Exception("b2t should be > 0 and <1 !!")
    flux1, flux2 = (tflux * frac1), (tflux * (1.0 - frac1))

    # Check, then read in other parameters
    galID = gal["ID"]
    # Effective radius
    try:
        reff1, reff2 = float(gal["reff_pix1"]), float(gal["reff_pix2"])
    except KeyError:
        raise KeyError("reff_pix1 or reff_pix2 is found in the record!!")
    # Sersic index
    try:
        nser1, nser2 = float(gal["sersic_n1"]), float(gal["sersic_n2"])
    except KeyError:
        raise KeyError("sersic_n1 or sersic_n2 is found in the record!!")
    # Axis ratio
    try:
        ba1, ba2 = float(gal["b_a1"]), float(gal["b_a2"])
    except KeyError:
        raise KeyError("b_a1 or b_a2 is found in the record!!")
    # Position angle
    try:
        pa1, pa2 = float(gal["theta1"]), float(gal["theta2"])
    except KeyError:
        raise KeyError("theta1 or theta2 is found in the record!!")


    comp1 = np.array((galID, flux1, nser1, reff1, ba1, pa1),
                     dtype=[('ID','int'), ('mag','float'), ('sersic_n','float'),
                            ('reff_pix','float'), ('b_a','float'),
                            ('theta','float')])
    comp2 = np.array((galID, flux2, nser2, reff2, ba2, pa2),
                     dtype=[('ID','int'), ('mag','float'), ('sersic_n','float'),
                            ('reff_pix','float'), ('b_a','float'),
                            ('theta','float')])

    return comp1, comp2


def arrayToGSObj(imgArr, scale=1.0, norm=False):
    # TODO : Check the scale here
    # According to the GalSim Doxygen
    # If provided, use this as the pixel scale for the Image; this will override
    # the pixel scale stored by the provided Image, in any. If scale is None,
    # then take the provided image's pixel scale. [default: None]
    """
    Convert an input 2-D array into a GalSim Image object
    """
    if norm:
        return galsim.InterpolatedImage(galsim.image.Image(imgArr),
                                         scale=scale, normalization="flux")
    else:
        return galsim.InterpolatedImage(galsim.image.Image(imgArr),
                                         scale=scale)


def galSimDrawImage(galObj, size=0, method="auto", addPoisson=False):
    """
    "Draw" a GalSim Object into an GalSim Image using certain method, and with
    certain size
    """
    # TODO : Think about the scale here:
    # By default scale=None
    # According to GalSim Doxygen :
    # If provided, use this as the pixel scale for the image. If scale is None
    # and image != None, then take the provided image's pixel scale. If scale is
    # None and image == None, then use the Nyquist scale. If scale <= 0
    # (regardless of image), then use the Nyquist scale.

    # Generate an "Image" object for the model
    if size > 0:
        imgTemp = galsim.image.Image(size, size)
        galImg = galObj.drawImage(imgTemp, method=method)
    else:
        galImg = galObj.drawImage(method=method)

    # Just an option for test
    if addPoisson:
        galImg.addNoise(galsim.PoissonNoise())

    # Return the Numpy array version of the image
    return galImg.array


def galSimConvolve(galObj, psfObj, size=0, method="auto", returnObj=False):
    """
    Just do convolution using GalSim
    Make sure the inputs are both GalSim GSObj
    The galaxy model should be the first one, and the PSF object is the second
    one; Returns a imgArr or GSObj
    """
    outObj = galsim.Convolve([galObj, psfObj])

    if returnObj:
        return outObj
    else:
        outArr = galSimDrawImage(galObj, size=size, method=method)
        return outArr


def galSimAdd(galObjList, size=0, method="auto", returnArr=False):
    """
    Just add a list of GSObjs together using GalSim
    Make sure all elements in the input list are GSObjs
    """
    if len(galObjList) < 2:
        raise Exception("Should be more than one GSObjs to add !")

    outObj = galsim.Add(galObjList)

    if returnArr:
        outArr = galSimDrawImage(outObj, size=size, method=method)
        return outArr
    else:
        return outObj


def plotFakeGalaxy(galObj, galID=None, suffix=None, size=0, addPoisson=False):

    """
    Generate a PNG image of the model
    By default, the image will be named as 'fake_galaxy.png'
    """

    import matplotlib.pyplot as plt

    if galID is None:
        outPNG = 'fake_galaxy'
    else:
        outPNG = 'fake_galaxy_%i' % galID
    if suffix is not None:
        outPNG = outPNG + '_' + suffix.strip() + '.png'

    plt.figure(1, figsize=(8,8))

    # Use "fft" just to be fast
    plt.imshow(np.arcsinh(galSimDrawImage(galObj, size=size, method="fft",
                                         addPoisson=addPoisson)))
    plt.savefig(outPNG)


def galSimFakeSersic(flux, gal, psfImage=None, scaleRad=False, returnObj=True,
                     expAll=False, devAll=False, plotFake=False, trunc=0,
                     drawMethod="auto", addPoisson=False):
    """
    Make a fake single Sersic galaxy using the galSim.Sersic function

    Inputs: total flux of the galaxy, and a record array that stores the
    necessary parameters [reffPix, nSersic, axisRatio, posAng]

    Output: a 2-D image array of the galaxy model  OR
            a GalSim object of the model

    Options:
        psfImage:     PSF image for convolution
        trunc:        Flux of Sersic models will truncate at trunc * reffPix
                      radius; trunc=0 means no truncation
        drawMethod:   The method for drawImage: ['auto', 'fft', 'real_space']
        addPoisson:   Add Poisson noise
        plotFake:     Generate a PNG figure of the model
        expAll:       Input model will be seen as nSersic=1
        devAll:       Input model will be seen as nSersic=4
        returnObj:    If TRUE, will return the GSObj, instead of the image array
    """

    # Convert the numpy.float32 into normal float format
    nSersic   = float(gal["sersic_n"])
    reffPix   = float(gal["reff_pix"])
    axisRatio = float(gal["b_a"])
    posAng    = float(gal["theta"])

    # Truncate the flux at trunc x reffPix
    if trunc > 0:
        trunc = int(trunc * reffPix)

    # Make sure Sersic index is not too large
    if nSersic > 6.0:
        raise ValueError("Sersic index is too large! Should be <= 6.0")
    # Check the axisRatio value
    if axisRatio <= 0.24:
        raise ValueError("Axis Ratio is too small! Should be >= 0.24")

    # Make the Sersic model based on flux, re, and Sersic index
    if nSersic == 1.0 or expAll:
        if scaleRad:
            serObj = galsim.Exponential(scale_radius=reffPix)
        else:
            serObj = galsim.Exponential(half_light_radius=reffPix)
        if expAll:
            print " * This model is treated as a n=1 Exponential disk : %d" % (gal["ID"])
    elif nSersic == 4.0 or devAll:
        serObj = galsim.DeVaucouleurs(half_light_radius=reffPix, trunc=trunc)
        if devAll:
            print " * This model is treated as a n=4 De Vaucouleurs model: %d" % (gal["ID"])
    else:
        serObj = galsim.Sersic(nSersic, half_light_radius=reffPix, trunc=trunc)

    # If necessary, apply the Axis Ratio (q=b/a) using the Shear method
    if axisRatio < 1.0:
        serObj = serObj.shear(q=axisRatio, beta=0.0*galsim.degrees)

    # If necessary, apply the Position Angle (theta) using the Rotate method
    if posAng != 0.0 or posAng != 180.0:
        serObj = serObj.rotate(posAng*galsim.degrees)

    # Convolve the Sersic model using the provided PSF image
    if psfImage is not None:
        # Convert the PSF Image Array into a GalSim Object
        # Norm=True by default
        psfObj  = arrayToGSObj(psfImage, norm=True)
        serFinal = galsim.Convolve([serObj, psfObj])
    else:
        serFinal = serObj

    # Pass the flux to the object
    serFinal = serFinal.withFlux(float(flux))

    # Make a PNG figure of the fake galaxy to check if everything is Ok
    # TODO: For test, should be removed later
    if plotFake:
        plotFakeGalaxy(serFinal, galID=gal['ID'])

    # Now, by default, the function will just return the GSObj
    if returnObj:
        return serFinal
    else:
        return galSimDrawImage(serFinal, method=drawMethod,
                               addPoisson=addPoisson)


def galSimFakeDoubleSersic(comp1, comp2, psfImage=None, trunc=0, returnObj=True,
                           devExp=False, plotFake=False, drawMethod='auto',
                           addPoisson=False):
    """
    Make a fake double Sersic galaxy using the galSim.Sersic function

    Inputs: total flux of the galaxy, and a record array that stores the
    necessary parameters [reffPix, nSersic, axisRatio, posAng]

    Output: a 2-D image array of the galaxy model  OR
            a GalSim object of the model

    Options:
        psfImage:     PSF image for convolution
        trunc:        Flux of Sersic models will truncate at trunc * reffPix
                      radius; trunc=0 means no truncation
        drawMethod:   The method for drawImage: ['auto', 'fft', 'real_space']
        addPoisson:   Add Poisson noise
        plotFake:     Generate a PNG figure of the model
        devexp:       The first component will be seen as a nSersic=4 bulge;
                      And, the second one will be seen as a nSersic=1 disk
        returnObj:    If TRUE, will return the GSObj, instead of the image array
    """

    # Get the flux of both components
    flux1 = float(comp1['mag'])
    flux2 = float(comp2['mag'])
    #tflux = flux1 + flux2

    # If devExp = True : Treat the first component as an n=4 DeVaucouleurs bulge
    #                    and, the second component as an n=1 Exponential disk
    if devExp:
        serModel1 = galSimFakeSersic(flux1, comp1, returnObj=True, devAll=True,
                                     trunc=trunc)
        serModel2 = galSimFakeSersic(flux2, comp2, returnObj=True, expAll=True,
                                     trunc=trunc)
    else:
        serModel1 = galSimFakeSersic(flux1, comp1, returnObj=True, trunc=trunc)
        serModel2 = galSimFakeSersic(flux2, comp2, returnObj=True, trunc=trunc)

    # Combine these two components
    doubleSersic = galSimAdd([serModel1, serModel2])

    # Convolve the Sersic model using the provided PSF image
    if psfImage is not None:
        # Convert the PSF Image Array into a GalSim Object
        # Norm=True by default
        psfObj   = arrayToGSObj(psfImage, norm=True)
        dserFinal = galsim.Convolve([doubleSersic, psfObj])
    else:
        dserFinal = doubleSersic

    # Make a PNG figure of the fake galaxy to check if everything is Ok
    # TODO: For test, should be removed later
    if plotFake:
        if devExp:
            plotFakeGalaxy(dserFinal, galID=comp1['ID'], suffix='devexp')
        else:
            plotFakeGalaxy(dserFinal, galID=comp1['ID'], suffix='double')

    # Now, by default, the function will just return the GSObj
    if returnObj:
        return dserFinal
    else:
        return galSimDrawImage(dserFinal, method=drawMethod,
                               addPoisson=addPoisson)


def galSimRealGalaxy(flux, real_galaxy_catalog, index=None, psfImage=None,
                     random=False, returnObj=True, plotFake=False,
                     drawMethod='auto', addPoisson=False):

    """
    Real galaxy
    """

    if index is None:
        random = True
    realObj = galsim.RealGalaxy(real_galaxy_catalog, index=index, random=random)
    index = realObj.index

    # Pass the flux to the object
    realObj = realObj.withFlux(flux)

    # Convolve the Sersic model using the provided PSF image
    if psfImage is not None:
        # Convert the PSF Image Array into a GalSim Object
        # Norm=True by default
        psfObj   = arrayToGSObj(psfImage, norm=True)
        realFinal = galsim.Convolve([realObj, psfObj])
    else:
        realFinal = realFinal

    # Make a PNG figure of the fake galaxy to check if everything is Ok
    # TODO: For test, should be removed later
    if plotFake:
        plotFakeGalaxy(realFinal, galID=index, suffix='realga')

    # Now, by default, the function will just return the GSObj
    if returnObj:
        return realFinal
    else:
        return galSimDrawImage(realFinal, method=drawMethod,
                               addPoisson=addPoisson)


def testMakeFake(galList, asciiTab=False, single=True, double=True, real=True):

    # Make a fake Gaussian PSF
    psfGaussian = galsim.Gaussian(fwhm=2.0)
    psfImage    = psfGaussian.drawImage().array

    # Test SingleSersic
    if single:
        if asciiTab:
            galData = np.loadtxt(galList, dtype=[('ID','int'),
                                                 ('mag','float'),
                                                 ('sersic_n','float'),
                                                 ('reff_pix','float'),
                                                 ('b_a','float'),
                                                 ('theta','float')])
        else:
            galData = fits.open(galList)[1].data

        for igal, gal in enumerate(galData):

            flux = 10.0 ** ((27.0 - gal['mag']) / 2.5)

            print '\n---------------------------------'
            print " Input Flux : ", flux
            print " Input Parameters : ", gal["sersic_n"], gal["reff_pix"]
            print "                    ", gal["b_a"], gal["theta"]

            galArray = galSimFakeSersic(flux, gal, psfImage=psfImage,
                                        plotFake=True, returnObj=False,
                                        trunc=12.0, drawMethod="fft")

            print " Output Flux : ", np.sum(galArray)
            print " Shape of the Output Array : ", galArray.shape
            print '---------------------------------'

    # Test DoubleSersic
    if double:
        if asciiTab:
            raise Exception("For now, only FITS input is allowed !!")
        else:
            galData = fits.open(galList)[1].data

        for igal, gal in enumerate(galData):

            flux = 10.0 ** ((27.0 - gal['mag']) / 2.5)

            print '\n---------------------------------'
            print " Input Flux : ", flux

            (comp1, comp2) = parseDoubleSersic(flux, gal)

            # TODO: Get error when the axis ratio is small: 0.2?
            # RuntimeError: Solve error: Too many iterations in bracketLowerWithLimit()
            # It seems like that GalSim has some issues with highly elliptical
            # objects. Although different, see this one:
            # https://github.com/GalSim-developers/GalSim/issues/384
            # It seems that b/a = 0.25 is fine, so right now, just change the
            # lower limit of b/a to 0.25

            print " Flux for Component 1 : ", comp1['mag']
            print " Flux for Component 2 : ", comp2['mag']
            print " Comp 1 Parameters : %5.2f  %8.2f" % (comp1["sersic_n"],
                                                         comp1["reff_pix"])
            print "                     %5.2f  %8.2f" % (comp1["b_a"],
                                                         comp1["theta"])
            print " Comp 2 Parameters : %5.2f  %8.2f" % (comp2["sersic_n"],
                                                         comp2["reff_pix"])
            print "                     %5.2f  %8.2f" % (comp2["b_a"],
                                                         comp2["theta"])

            doubleArray = galSimFakeDoubleSersic(comp1, comp2,
                                                 psfImage=psfImage,
                                                 trunc=12, returnObj=False,
                                                 devExp=True, plotFake=True,
                                                 drawMethod='auto')

            print " Output Flux : ", np.sum(doubleArray)
            print " Shape of the Output Array : ", doubleArray.shape
            print '---------------------------------'

    # Test RealGalaxy
    if real:

        # Make a special PSF for real galaxy
        psfReal      = galsim.Gaussian(fwhm=0.2)
        psfRealImage = psfReal.drawImage().array
        # TODO : Scale seems to be a problem, should check again !!

        if asciiTab:
            raise Exception("For now, only FITS input is allowed !!")
        else:
            galData = fits.open(galList)[1].data

        for igal, gal in enumerate(galData):

            (real_galaxy_catalog, index) = parseRealGalaxy(gal)
            if index >= 0:
                index = index
                random = False
            else:
                index = None
                random = True

            flux = 10.0 ** ((27.0 - gal['mag']) / 2.5)

            realArray = galSimRealGalaxy(flux, real_galaxy_catalog, index=index,
                                         psfImage=psfRealImage, random=random,
                                         plotFake=True, returnObj=False,
                                         drawMethod='auto')

            print '\n---------------------------------'
            print " Input Flux : ", flux

            print " Output Flux : ", np.sum(realArray)
            print " Shape of the Output Array : ", realArray.shape
            print '---------------------------------'
