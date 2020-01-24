import sys
import os
import os.path
import warnings
warnings.filterwarnings("ignore")
import scipy
import pyds9
import numpy as np
import pandas as pd
import subprocess as sp
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import LogNorm
from tkinter import Tk, Label, Button, Entry, IntVar, END, W, E, N
from tkinter import messagebox as tkMessageBox

home_dir = os.getcwd()

def GetSize(x, y, R, theta, ell, ncol, nrow):
    '''this subroutine get the maximun
    and minimim pixels for Kron and sky ellipse'''
    # k Check
    q = (1 - ell)
    bim = q * R
    theta = theta * (np.pi / 180)  # rads!!

    # getting size
    xmin = x - np.sqrt((R**2)   * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)
    xmax = x + np.sqrt((R**2)   * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)
    ymin = y - np.sqrt((R**2)   * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)
    ymax = y + np.sqrt((R**2)   * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        if isinstance(xmin,np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol
    if mask.any():
        if isinstance(xmax,np.ndarray):
            xmax[mask] = ncol
        else:
            xmax = ncol

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin,np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax,np.ndarray):
            ymax[mask] = nrow
        else:
            ymax = nrow

    return (xmin, xmax, ymin, ymax)

def GetPng(Image):
    "Converts FITS file into a PNG image with axis coordinates, inverted colormap, log/zmax style"
    filename = get_pkg_data_filename(Image)

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    counts = hdu.data.ravel()

    n, bins, patches = plt.hist(counts, bins=512, range=(min(counts), max(counts)), color='black')
    n_max = np.argsort(n)[::-1]
    n_idx = n_max[0:20]
    lim_min = min(bins[n_idx])
    lim_max = max(bins[n_idx])

    bri = 33 # brightness, source: docs.opencv.org/3.4/d3/dc1/tutorial_basic_linear_transform.html
    con = 0.98 # contrast, > 0

    plt.subplot(projection=wcs)
    plt.imshow(con*hdu.data+bri, cmap='gray_r', norm=LogNorm(lim_min, lim_max))

    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.grid(color='black', ls='solid', alpha=0.1)
    plt.title(Image)
    plt.tight_layout()
    plt.savefig('%s.png' % (Image), bbox_inches='tight', dpi=300)

def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"

    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)

    hdu = fits.open(Image)
    dat = hdu[0].data[ylo - 1:yhi, xlo - 1:xhi]
    hdu[0].data = dat
    hdu.writeto(Imageout, clobber=True)
    hdu.close()

def GetAxis(Image):
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow

######################################### getting sizes ##################################

def main():

    if len(sys.argv[1:]) != 2:
        print ('Missing arguments')
        print ("Usage:\n %s [Image] [Catalog]" % (sys.argv[0]))
        print ("Example:\n %s Image.fits catalog.cat" % (sys.argv[0]))
        sys.exit()

    print("Sorting and getting sizes for objects")

    image   = sys.argv[1]
    catalog = sys.argv[2]

    #GetAxis(image):
    hdu = fits.open(image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()

    KronScale = 1
    SexSort = "new_catalog.cat"

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        catalog, delimiter="", unpack=True)

    n   =   n.astype(int)
    flg = flg.astype(int)

    Rkron = KronScale * ai * kr

    Bim = (1 - e) * Rkron
    Area = np.pi * Rkron * Bim*(-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)
    f_out = open(SexSort, "w")
    index = mg.argsort()

    for i in index:
        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(n[i],
        alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i],
        theta[i], bkgd[i], idx[i], flg[i],
        np.int(np.round(sxmin[i])), np.int(np.round(sxmax[i])),
        np.int(np.round(symin[i])), np.int(np.round(symax[i])))
        f_out.write(line)

    f_out.close()

######################################### image creation ##################################

    print("Generating object images")
    catalog2 = SexSort

    #reads table with astropy's ascii function and later turns it to a pandas dataframe
    cat2 = ascii.read(catalog2)
    df = cat2.to_pandas()
    #renames columns
    df.columns = ['Num','RA','Dec','XPos','YPos','Mag','Kron','FluxRad','IsoArea',
                  'AIm', 'E','Theta','Background','Class','Flag','XMin','XMax','YMin','YMax']

    lim_mag  = 17
    lim_flag = 0.6
    #Filters table by magnitude (keeps objects with mag <= lim_lag)
    df = df[(df.Mag <= lim_mag)]
    #Filters table by flag (keeps objects with flag <= 0.6)
    #(0 galaxy, 1 = star, object is galaxy if < lim_flag)
    df = df[(df.Flag <= lim_flag)]

    #Index Reset
    df = df.reset_index(drop=True)

    # Catalog size limit for TESTING PURPOSES#
    df = df[0:30]

    #renames columns for later use
    [Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, E, Theta,
    Background, Class, Flag, XMin, XMax, YMin, YMax] = [df.Num, df.RA, df.Dec,
    df.XPos, df.YPos, df.Mag, df.Kron, df.FluxRad, df.IsoArea, df.AIm, df.E,
    df.Theta, df.Background, df.Class, df.Flag,df.XMin, df.XMax, df.YMin, df.YMax]

    def fitslist(Num):
        "Makes a list of .fits file names using the object number"
        temp = list()
        for index in range(0, len(Num)):
            x = str(Num[index])
            temp.append('Image-' + x + '.fits')
        return temp

    # runs previous function and saves list
    fitnames = fitslist(Num)

    # Creates home directory if it doesn't exist
    if not os.path.isdir(home_dir):
        os.makedirs(home_dir)
        print("Home directory %s was created." %home_dir)

    # Makes a .fits file for every object
    for item in fitnames:
        filePath = f"{home_dir}%s" % (item)
        open(filePath, 'w')

    # Makes a text file with every object's .fits file name
    #with open(f'{home_dir}lista_fits.txt', 'w') as f:
    #    for item in fitnames:
    #        f.write("%s\n" % item)

    # make a table with coordinates of every image
    coords_table = 'list_coords.txt'
    df_xy = pd.DataFrame({'Name': fitnames,
                       'X_coord': XPos,
                       'Y_coord': YPos,
                            'RA': RA,
                           'Dec': Dec})

    # rounds numbers for aesthetic purposes
    #df_xy = df_xy.round({'X_coord': 3,'Y_coord': 3,'RA': 6,'Dec': 6})

    # saves table df_xy in a file
    df_xy.to_csv(path_or_buf=coords_table, sep='\t', index=False, header=True, mode='w')
    #############

    Angle = Theta - 90
    AR = 1 - E
    KronScale = 1
    RKron = KronScale * AIm * Kron
    Sky = Background

    Num = Num.astype(int)
    Flag = Flag.astype(int)

    #####AGRANDAR IMAGEN (marcos) USANDO:

    XSize = XMax - XMin
    YSize = YMax - YMin

    # enlarge fit area

    FitBox = 6

    XSize = FitBox * XSize
    YSize = FitBox * YSize

    ##  60 pixels is the minimum area to fit (arbitrary number):
    min_area = 60

    masksize = XSize < min_area

    if masksize.any():
        XSize[masksize] = min_area

    masksize = YSize < min_area

    if masksize.any():
        YSize[masksize] = min_area

    #  Calculate the (x,y) position of the current object relative to
    #  the tile in which it lives.

    XFit = XPos
    YFit = YPos

    #  Calculate fitting box needed to plug into galfit header:

    XLo = XFit - XSize / 2
    XLo = XLo.astype(int)

    maskxy = XLo <= 0
    if maskxy.any():
        XLo[maskxy] = 1

    XHi = XFit + XSize / 2
    XHi = XHi.astype(int)

    NCol, NRow = GetAxis(image)

    maskxy = XHi > NCol  #This does not affect the code at all
    if maskxy.any():
        XHi[maskxy] = NCol

    YLo = YFit - YSize / 2
    YLo = YLo.astype(int)

    maskxy = YLo <= 0
    if maskxy.any():
        YLo[maskxy] = 1

    YHi = YFit + YSize / 2
    YHi = YHi.astype(int)

    maskxy = YHi > NRow  # same as above but for y axis
    if maskxy.any():
        YHi[maskxy] = NRow

    ##### Creates FITS and PNG images with objects cropped from original image

    for i in range(0, len(Num)):
        GetFits(image, fitnames[i], XLo[i], XHi[i], YLo[i], YHi[i])

    for i in fitnames:
        GetPng(i)

if __name__ == '__main__':
    main()

##################################### GUI #####################################

print("Ready to classify")

image  = sys.argv[1]

# table file with coordinates
coords = 'list_coords.txt'
# final table file name with merged data
file_out='list_final.txt'

# folder and file extensions to use
input_folder = home_dir + '/' #current location (os.getcwd())
file_extensions = ['.fits'] # must be fits

# list of all files with the chosen file extension in input folder
file_names = [fn for fn in sorted(os.listdir(input_folder))
              if any(fn.endswith(ext) for ext in file_extensions)]
file_names.sort()
# removes main image from list
file_names.remove(image)

# creates a list of the files, path included
paths = [input_folder + x for x in file_names]

d = pyds9.DS9()
d.set("cmap invert")
d.set("zoom to fit")
d.set("scale log")
d.set("scale zmax")
d.set("file %s" % (paths[0]))

class GUI:
    def __init__(self, master):
        self.master = master

        self.index = 0
        self.paths = paths
        self.n_paths = len(paths)

        self.table  = pd.DataFrame({'Name': file_names})
        self.obj_id = []

        master.title("Classify Galaxies with Python and DS9")

        self.label2 = Label(master, text='Select a class:')
        self.E  = Button(master, text='E', command=self.E)
        self.S0 = Button(master, text='S0',command=self.S0)
        self.S0B= Button(master, text='S0B',command=self.S0B)
        self.S  = Button(master, text='S', command=self.S)
        self.SB = Button(master, text='SB',command=self.SB)
        self.label = Label(master, text=" ")
        self.q  = Button(master, text='Quit', command=self.ask_quit)

        self.label2.grid(row=1, columnspan=5, sticky=N)
        self.E.grid(row=2,  column=1)
        self.S0.grid(row=2, column=2)
        self.S0B.grid(row=2,column=3)
        self.S.grid(row=2,  column=4)
        self.SB.grid(row=2, column=5)
        self.label.grid(row=3, columnspan=5)
        self.q.grid(row=4, column=1)

    def show_next_image(self):
        '''displays next image, or saves data when classification is over'''

        if len(self.obj_id) < self.n_paths: # opens next image
            self.index += 1
            d = pyds9.DS9()
            d.set("cmap invert")
            d.set("zoom to fit")
            d.set("scale log")
            d.set("scale zmax")
            d.set("file %s" % (paths[self.index]))

        if len(self.obj_id) == self.n_paths: # saves data into file
            #printing lengths to verify code working as intended
            #print('# of classes:', len(self.obj_id))
            #print('# of images: ', self.n_paths)

            # loads coordinates dataframe
            dfcoords = ascii.read(coords, header_start=0)
            dfcoords = dfcoords.to_pandas()

            # creates class dataframe
            df = pd.DataFrame({'Name': file_names,'Class': self.obj_id})

            # merge coordinates and class dataframes, then sorts new df by name
            df_m = dfcoords.merge(df, on='Name', how='outer') # outer: use union of keys from both frames
            df_m.sort_values(['Name'], ascending=True, inplace=True)

            # saves newly merged dataframe
            df_m.to_csv(path_or_buf=file_out, sep='\t', index=False, header=True, mode='w')

            print('Classification successful.')
            print('Data has been saved to %s.' %(file_out))
            tkMessageBox.showinfo('Classification successful', 'Data has been saved to %s.' %(file_out))

        if len(self.obj_id) > (self.n_paths):
            print('Everything is already classified.')
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def E(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('E')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S0(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S0')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S0B(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S0B')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def SB(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('SB')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def ask_quit(self):
        if tkMessageBox.askokcancel("Quit", "Are you sure you want to quit?\nData will be lost if you're not finished."):
            root.destroy()

root = Tk()
gui = GUI(root)
root.mainloop()
