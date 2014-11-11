import galsim
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os.path
import subprocess
from randSerCat import serGen

def getFitsFile(flux, nser, reff, ba, pa, zp=26.0, scale=1.0, method='auto',
                psfConv=False, addNoise=False, fwhm=2.0):

    if psfConv:
        psfObj = galsim.Gaussian(fwhm=fwhm)
        psfImg = psfObj.drawImage(scale=scale)
        galsim.fits.write(psfImg, file_name='psf.fits')

    trunc = int(10.0 * reff)

    serTrun = galsim.Sersic(nser, half_light_radius=reff, trunc=trunc)
    serTrun = serTrun.withFlux(flux)
    serTrun = serTrun.shear(q=ba, beta=0.0*galsim.degrees)
    serTrun = serTrun.rotate(pa*galsim.degrees)
    if psfConv:
        serTrun = galsim.Convolve([serTrun, psfObj])
    serImg  = serTrun.drawImage(method=method, scale=scale)
    if addNoise:
        serImg.addNoise(galsim.PoissonNoise())
        serImg.addNoise(galsim.CCDNoise())
    galsim.fits.write(serImg, file_name="temp.fits")

    xs, ys = serImg.array.shape
    if (xs % 2 ==0):
        xc = (xs + 1) / 2.0
    else:
        xc = xs / 2.0
    if (ys % 2 ==0):
        yc = (ys + 1) / 2.0
    else:
        yc = ys / 2.0

    return xs, ys, xc, yc


def writeInput(xs, ys, xc, yc, zp=26.0, scale=1.0):

    f = open('temp.in', 'w')
    f.write('===============================================================================\n')
    f.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
    f.write('A) temp.fits        # Input data image (FITS file)\n')
    f.write('B) out.fits         # Output data image block\n')
    f.write('C) none             # Sigma image name (made from data if blank or "none") \n')
    f.write('D) psf.fits         # Input PSF image and (optional) diffusion kernel\n')
    f.write('E) 1                # PSF fine sampling factor relative to data \n')
    f.write('F) none             # Bad pixel mask (FITS image or ASCII coord list)\n')
    f.write('G) none             # File with parameter constraints (ASCII file) \n')
    f.write('H) 1 %d 1 %d # Image region to fit\n' % (xs, ys))
    f.write('I) %d  %d           # Size of the convolution box\n' % (xs, ys))
    f.write('J) %7.3f            # Magnitude photometric zeropoint\n' % zp)
    f.write('K) %7.3f %7.3f      # Plate scale [arcsec per pixel]\n' % (scale, scale))
    f.write('O) regular             # Display type (regular, curses, both)\n')
    f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n')
    f.write('\n')
    f.write('# \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
    f.write('# -----------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('# Object number: 1\n')
    f.write(' 0) sersic                 #  object type\n')
    f.write(' 1) %5.1f %5.1f  1 1       #  position x, y\n' % (xc, yc))
    f.write(' 3) 20.000      1          #  Integrated magnitude	\n')
    f.write(' 4)  4.000      1          #  R_e (half-light radius)   [pix]\n')
    f.write(' 5)  1.0000     1          #  Sersic index n (de Vaucouleurs n=4) \n')
    f.write(' 9)  0.8        1          #  axis ratio (b/a)  \n')
    f.write('10)  0.0        1          #  position angle (PA) [deg: Up=0, Left=90]\n')
    f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
    f.write('\n')
    f.write('================================================================================\n')
    f.close()

    return None

def runGalfit(galfit='/home/hs/code/galfit/galfit'):

    p = subprocess.Popen([galfit + ' temp.in'], shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    #for line in p.stdout.readlines():
    #    print line

    retval = p.wait()

    if not os.path.exists('out.fits'):
        raise Exception("The GALFIT run failed ! Check !")

    return None

def getResults():

    head = fits.open('out.fits')[2].header

    mag  = float(head['1_MAG'].split()[0])
    re   = float(head['1_RE'].split()[0])
    nser = float(head['1_N'].split()[0])
    ba   = float(head['1_AR'].split()[0])

    return mag, re, nser, ba


if __name__ == '__main__':

    zp = 26.0

    models = serGen(20, minMag=17.0, maxMag=24.0, nSer=1.0, ba=1.0, pa=0.0)

    nModel = models.shape[0]

    print "# mag_in mag_out re_in re_out n_in n_out ba_in ba_out  "

    for i in range(nModel):

        mag = models['mag'][i]
        flux = 10.0 ** ((zp - mag) / 2.50)
        reff = models['reff_pix'][i]
        nser = models['sersic_n'][i]
        ba   = models['b_a'][i]
        pa   = models['theta'][i]

        (xs, ys, xc, yc) = getFitsFile(flux, nser, reff, ba, pa)

        writeInput(xs, ys, xc, yc)

        runGalfit()

        (mag_out, reff_out, nser_out, ba_out) = getResults()

        print "%7.3f %7.3f  %7.3f %7.3f  %7.3f %7.3f  %7.3f %7.3f" % (mag, mag_out,
                                                reff, reff_out, nser, nser_out,
                                                ba, ba_out)

