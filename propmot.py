import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
import pdb
import numpy as np
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
import sys
from scipy import stats
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
import zuds
from iminuit import Minuit
from iminuit.util import describe, make_func_code
from iminuit.cost import LeastSquares


zuds.init_db()


#Once testing is done, use sys.argv[0] here

file_list = os.listdir('/home/conor/.data/newdata/')



#Read in ZTF science image FITS files and sort them

files = []
coords = []

for file in file_list:
	if 'sciimg.cat' in file:
		coords.append('/home/conor/.data/newdata/'+file)
	elif 'sciimg.fits' in file:
		files.append('/home/conor/.data/newdata/' + file)
coords.sort()
files.sort()

#Read in files until one with a decent seeing and mag limit is found

i = 0
for image in files:
	f = fits.open(image)[0]
	if f.header['MAGLIM'] > 20.5 and f.header['SEEING'] <= 2:
		base_image = image
		base_coord = coords[i]
		magb = f.header['MAGZP']
		date_base = f.header['OBSJD']
		
		break
	i+=1
	

#Read in base image to pull Gaia data to calibrate proper motions later


sci = zuds.ScienceImage.from_file(base_image)

gaia = sci.gaia_dr2_calibrators()


gpmra = gaia['pmra']
gpmdec = gaia['pmdec']
gpmra_err = gaia['pmra_error']
gpmdec_err = gaia['pmdec_error']
gra = gaia['ra']
gdec = gaia['dec']

#Put data into a form that SkyCoord will accept

gra2 = []
gdec2 = []

for i in range(len(gra)):
	gra2.append(gra[i])
	gdec2.append(gdec[i])
	
#Read in base cat file and read in base coords

f = fits.open(base_coord)

ra_base = f[2].data['XWIN_WORLD']
dec_base = f[2].data['YWIN_WORLD']
ra_basserr = np.sqrt(f[2].data['ERRX2WIN_IMAGE']) * 1.012 * 1000
dec_baseerr = np.sqrt(f[2].data['ERRY2WIN_IMAGE']) * 1.012 * 1000



files.remove(base_image)
coords.remove(base_coord)

ra = []
dec = []
ra_err = []
dec_err = []
mag = []
t_frame = []
sn = []
dates = []
dates.append(date_base)



#Save all dates in one list that will be accessed later to assign each location a date

for file in files:

	fit = fits.open(file)

	dates.append(fit[0].header['OBSJD'])

	
count = 0
for reg in coords:
	f = fits.open(reg)


	for i in range(len(f[2].data['XWIN_WORLD'])):
		#Sometimes when using stacked data, there are very large reported magnitudes
		if f[2].data['MAG_AUTO'][i] + magb < 50:
			ra.append(f[2].data['XWIN_WORLD'][i])
			dec.append(f[2].data['YWIN_WORLD'][i])
			mag.append(f[2].data['MAG_AUTO'][i] + magb)
			#Err units should be in mas while ra and dec should be in degrees, which will be changed to mas later
			ra_err.append(np.sqrt(f[2].data['ERRX2WIN_IMAGE'][i])*1.012*1000)
			dec_err.append(np.sqrt(f[2].data['ERRY2WIN_IMAGE'][i])*1.012*1000)
			sn.append(f[2].data['FLUX_AUTO'][i]/f[2].data['FLUXERR_AUTO'][i])
	

	#This is used later to assign dates	
	if count == 0:
		t_frame.append(0)
		t_frame.append(t_frame[-1] + len(f[2].data['XWIN_WORLD'])-1)
	else:
		t_frame.append(t_frame[-1] + len(f[2].data['XWIN_WORLD']))


	count+=1




stars = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_base*u.degree,dec=dec_base*u.degree),0.000277778*u.degree)



#Match stars in the base cat file to the other ZTF cat files

x = []
time = []
ra_matched = []
dec_matched = []
ra_bmatched = []
dec_bmatched = []
mag_matched = []
raerr_matched = []
raerr_bmatched = []
decerr_matched = []
sn_matched = []
decerr_bmatched = []
for j in range(len(sep2d.arcsecond)):
	x.append(sep2d.arcsecond[j])
	ra_matched.append(ra[idxself[j]])
	ra_bmatched.append(ra_base[idxs[j]])
	dec_matched.append(dec[idxself[j]])
	dec_bmatched.append(dec_base[idxs[j]])
	mag_matched.append(mag[idxself[j]])
	raerr_matched.append(ra_err[idxself[j]])
	raerr_bmatched.append(ra_basserr[idxs[j]])
	decerr_matched.append(dec_err[idxself[j]])
	sn_matched.append(sn[idxself[j]])
	

#Determine date of star location by cross matching its position in ra with the time frame list

for k in range(len(idxself)):
	for l in range(len(t_frame)):
		if idxself[k] <= t_frame[l] and idxself[k] > t_frame[l-1]:
			#Units should be in yrs
			time.append((dates[l-1]-date_base)/365)
		



def line(x, a, b):  #Simple straight line model with explicit parameters
    return a + b * x



chunksra = []
chunksdec = []
chunkspmra = []
chunkspmdec = []
chunkspmraerr = []
chunkspmdecerr = []
chunkspmr = []
chunkspmr2 = []
chunksmag = []
chunksx = []
chunkssn = []
count = 0

badsn = []
badmag = []

index = 0

#Code to group stars tht are within one arcsecond with each other in chunks

for i in range(len(x)):
	ra_chunk = []
	dec_chunk = []
	t_chunk = []
	mag_chunk =[]
	x_chunk = []
	raerr_chunk =[]
	decerr_chunk = []
	ra_real = []
	dec_real = []
	sn_chunk = []
	#If locations are within one arcsecond, continue
	if i!= len(x)-1 and i!=0 and np.abs(ra_matched[i] - ra_matched[i-1]) < 0.000277778:
		continue
	elif i == 0:
		continue
	else:
		for l in range(index,i):
			ra_real.append(ra_matched[l])
			dec_real.append(dec_matched[l])
			#Units in ra and dec chunks should be in mas now
			ra_chunk.append((ra_matched[l]-ra_bmatched[l])*3600*1000)
			dec_chunk.append((dec_matched[l]-dec_bmatched[l])*3600*1000)
			sn_chunk.append(sn_matched[l])
			x_chunk.append(x[l])
			mag_chunk.append(mag_matched[l])
			t_chunk.append(time[l])
			raerr_chunk.append(raerr_matched[l])
			decerr_chunk.append(decerr_matched[l])
		
		index = i
		#Want at least 10% of files to be represented
		if len(t_chunk) < 10:
			continue


		new_chunk = []
		newt_chunk1 = []
		newt_chunk = []
		new_chunk1 = []
		new_chunk2 = []
		new_chunk3 = []
		new_chunk4 = []
		new_mag = []
		for i in range(len(x_chunk)):
			if np.abs(x_chunk[i] - np.median(x_chunk)) < 2 * 1.428 * stats.median_absolute_deviation(x_chunk - np.median(x_chunk)):
				new_chunk.append(x_chunk[i])
				newt_chunk.append(t_chunk[i])
				new_chunk1.append(ra_chunk[i])
				new_chunk2.append(dec_chunk[i])
				new_mag.append(mag_chunk[i])
				new_chunk3.append(raerr_chunk[i])
				new_chunk4.append(decerr_chunk[i])



		#Run Minuit on outlier-free ra and dec chunks
		
		least_squares = LeastSquares(newt_chunk, new_chunk1, new_chunk3,line)


		m = Minuit(least_squares, a=0,b=0)

		
		m.migrad()
		m.hesse()

		least_squares1 = LeastSquares(newt_chunk, new_chunk2, new_chunk4,line)

		m1 = Minuit(least_squares1, a=0,b=0)

		
		m1.migrad()
		m1.hesse()

		

		#Run linear regression with scipy.stats
		slope, intercept, r_value, p_value, std_err = stats.linregress(newt_chunk,new_chunk)
		slope1, intercept1, r_value1, p_value1,std_err1 = stats.linregress(newt_chunk,new_chunk1)
		slope2, intercept2, r_value2, p_value2,std_err2 = stats.linregress(newt_chunk,new_chunk2)

		stars = SkyCoord(ra=[np.median(ra_real)]*u.degree,dec=[np.median(dec_real)]*u.degree)
		idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=gra2*u.degree,dec=gdec2*u.degree),0.000277778*u.degree)
		
		#Check if stars are blended together
		if len(idxs) > 1:
			continue


		chunksra.append(np.median(ra_real))
		chunksdec.append(np.median(dec_real))
		chunkspmra.append(m.values['b'])
		chunkspmdec.append(m1.values['b'])
		chunkspmraerr.append(m.errors['b'])
		chunkspmdecerr.append(m1.errors['b'])
		chunksmag.append(np.median(new_mag))
		chunkssn.append(np.median(sn_chunk))
			




gmra_matched = []
pmra_matched = []
gmraerr_matched = []
pmraerr_matched = []
sn_matched = []
mag_matched = []

#Cross check with Gaia proper motions

stars = SkyCoord(ra=chunksra*u.degree,dec=chunksdec*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=gra2*u.degree,dec=gdec2*u.degree),0.000277778*u.degree)
for i in range(len(idxs)):
	pmra_matched.append(chunkspmra[idxself[i]])
	gmra_matched.append(gpmra[idxs[i]])
	gmraerr_matched.append(gpmra_err[idxs[i]])
	pmraerr_matched.append(chunkspmraerr[idxself[i]])
	sn_matched.append(chunkssn[idxself[i]])
	mag_matched.append(chunksmag[idxself[i]])

model = []

for i in range(len(pmra_matched)):
	model.append((pmra_matched[i]-gmra_matched[i])/pmraerr_matched[i])

count = 0

for i in range(len(model)):
	if np.abs(model[i]) < 2:
		count +=1
		
print(count/len(idxs))




#plt.errorbar(gmra_matched,model,yerr = pmraerr_matched,xerr = gmraerr_matched,fmt = 'o',ecolor = 'b',capthick=1,capsize = 2)
plt.scatter(gmra_matched,model)
#ax.add_artist(plt.scatter(mag2,x2,color = 'r'))
#plt.plot(gmra_matched,gmra_matched)
plt.xlabel("Gaia Proper Motion in RA (mas/yr)")
plt.ylabel("Simga Offset from True Values (Using ZTF Errors)")
plt.title('Gaia vs ZTF Proper Motions')
plt.show()

