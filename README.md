# Classify Galaxies with Python and DS9
A Python GUI tool to classify galaxies with the help of DS9.

This script uses a FITS image and its SExtractor catalog to create cropped images of each galaxy, in order to display them in DS9. DS9 allows the user analyze each galaxy in better detail and classify it correctly. Tkinter is used to have a GUI that displays each possible classification, with an option to quit when requested. 

The classification data will only be saved in a new file once all images are classified. More classification options can be added easily to the GUI, but may be cumbersome to do so. Some variables (listed below) may be changed within the code to adjust to different needs, e.g, image/catalog destination, file names.

The images and tables created are not deleted after the tool is closed.

## Usage
python CGPD.py [image.fits] [catalog]

## SExtractor catalog requirements
In the same order, without headers:

\#   1 NUMBER -                 Running object number

\#   2 ALPHA_J2000 -            Right ascension of barycenter (J2000)                      [deg]

\#   3 DELTA_J2000 -            Declination of barycenter (J2000)                          [deg]

\#   4 X_IMAGE -                Object position along x                                    [pixel]

\#   5 Y_IMAGE -                Object position along y                                    [pixel]

\#   6 MAG_APER -               Fixed aperture magnitude vector                            [mag]

\#   7 KRON_RADIUS -            Kron apertures in units of A or B

\#   8 FLUX_RADIUS -            Fraction-of-light radii                                    [pixel]

\#   9 ISOAREA_IMAGE -          Isophotal area above Analysis threshold                    [pixel**2]

\#  10 A_IMAGE -                Profile RMS along major axis                               [pixel]

\#  11 ELLIPTICITY -            1 - B_IMAGE/A_IMAGE

\#  12 THETA_IMAGE -            Position angle (CCW/x)                                     [deg]

\#  13 BACKGROUND -             Background at centroid position                            [count]

\#  14 CLASS_STAR -             S/G classifier output

\#  15 FLAGS -                  Extraction flags

## Some variables

lim_mag: Magnitude limit for galaxies in catalog (default = 17).

lim_flag: SExtractor flag limit (star = 1, galaxy = 0, default = 0.6).

catalog size limit: The dataframe axis may be changed for testing purposes.

home_dir: Option to change image/catalog directory.

min_area: Minimum area in pixels for galaxy images (default = 60).

coords_table: Name of dataframe that includes image name, pixel and RA/Dec coordinates.

file_out: Name of final table with classifications for each galaxy included.

d.set's: DS9 commands used to help analyze each galaxy.

## Some functions

GetSize: Gets galaxy box dimensions, used to create a new FITS image of cropped galaxy.

GetPng: Creates a PNG image of the galaxy, displaying it in DS9 fashion with FK5 coordinates. Using it makes the process slower and can be disabled.

GetFits: Creates cropped image of galaxy using GetSize outputs.

GetAxis: Gets axis of image's header.

