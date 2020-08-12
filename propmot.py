import glob
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
from itertools import zip_longest
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

file_list = os.listdir('/home/conor/.data/newdata/')


files = []
coords = []

for file in file_list:
	if 'sciimg.cat' in file:
		coords.append('/home/conor/.data/newdata/'+file)
	elif 'sciimg.fits' in file:
		files.append('/home/conor/.data/newdata/' + file)
coords.sort()
files.sort()



matches = []
dates = []

check = '/home/conor/.data/newdata/ztf_20200715280428_000762_zg_c01_o_q1_sciimg.fits'

lgo = fits.open(check)

magb = lgo[0].header['MAGZP']

sci = zuds.ScienceImage.from_file(check)

yep = sci.gaia_dr2_calibrators()


gpmra = yep['pmra']
gpmdec = yep['pmdec']
gpmra_err = yep['pmra_error']
gpmdec_err = yep['pmdec_error']
gra = yep['ra']
gdec = yep['dec']


gra2 = []
gdec2 = []

for i in range(len(gra)):
	gra2.append(gra[i])
	gdec2.append(gdec[i])


i = 0
for image in files:
	f = fits.open(image)[0]
	if f.header['MAGLIM'] > 20.5 and f.header['SEEING'] <= 2:
		base_image = image
		base_coord = coords[i]
		break
	i+=1
	

f = fits.open(base_coord)


ra_base = f[2].data['XWIN_WORLD']
dec_base = f[2].data['YWIN_WORLD']
ra_basserr = f[2].data['ERRAWIN_IMAGE'] * 1.012
dec_baseerr = f[2].data['ERRBWIN_IMAGE'] * 1.012


fit = fits.open(base_image)
date_base = fit[0].header['OBSJD']

dates.append(date_base)


'''sci = zuds.ScienceImage.from_file('/home/conor/.data/newdata/ztf_20200715280428_000762_zg_c01_o_q1_sciimg.fits')

ra_g = []
dec_g = []
pmra_g = []
pmdec_g = []

columns = sci.gaia_dr2_calibrators()
rag = columns['ra']
decg = columns['dec']
pmrag = columns['pmra']
pmdecg = columns['pmdec']


for i in range(len(rag)):
    ra_g.append(rag[i])
    dec_g.append(decg[i])
    pmra_g.append(pmrag[i])
    pmdec_g.append(pmdecg[i])'''


files.remove(base_image)
coords.remove(base_coord)

ra = []
dec = []
ra_err = []
dec_err = []
mag = []
t_frame = []
sn = []

for file in files:

	fit = fits.open(file)

	dates.append(fit[0].header['OBSJD'])

count = 0
for reg in coords:
	f = fits.open(reg)


	for i in range(len(f[2].data['XWIN_WORLD'])):
		if f[2].data['MAG_AUTO'][i] + magb < 50:
			ra.append(f[2].data['XWIN_WORLD'][i])
			dec.append(f[2].data['YWIN_WORLD'][i])
			mag.append(f[2].data['MAG_AUTO'][i] + magb)
			ra_err.append(np.sqrt(f[2].data['ERRX2WIN_IMAGE'][i])*1.012*1000)
			dec_err.append(np.sqrt(f[2].data['ERRY2WIN_IMAGE'][i])*1.012*1000)
			sn.append(f[2].data['FLUX_AUTO'][i]/f[2].data['FLUXERR_AUTO'][i])
	

		
	if count ==0:
		t_frame.append(0)
		t_frame.append(t_frame[-1] + len(f[2].data['XWIN_WORLD'])-1)
	else:
		t_frame.append(t_frame[-1] + len(f[2].data['XWIN_WORLD']))


	count+=1

x = []



stars = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_base*u.degree,dec=dec_base*u.degree),0.000277778*u.degree)
	

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
	



for k in range(len(idxself)):
	for l in range(len(t_frame)):
		if idxself[k] <= t_frame[l] and idxself[k] > t_frame[l-1]:
			time.append((dates[l-1]-date_base)/365)
		
pmra = []
pmdec = []
index = 0


def line(x, a, b):  # simple straight line model with explicit parameters
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
	if i!=0 and np.abs(ra_matched[i] - ra_matched[i-1]) < 0.000277778:
		continue
	elif i == 0:
		continue
	else:
		for l in range(index,i):
			ra_real.append(ra_matched[l])
			dec_real.append(dec_matched[l])
			ra_chunk.append((ra_matched[l]-ra_bmatched[l])*3600*1000)
			dec_chunk.append((dec_matched[l]-dec_bmatched[l])*3600*1000)
			sn_chunk.append(sn_matched[l])
			x_chunk.append(x[l])
			mag_chunk.append(mag_matched[l])
			t_chunk.append(time[l])
			raerr_chunk.append(raerr_matched[l])
			decerr_chunk.append(decerr_matched[l])
		
		index = i
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





		least_squares = LeastSquares(newt_chunk, new_chunk1, new_chunk3,line)


		m = Minuit(least_squares, a=0,b=0)

		
		m.migrad()
		m.hesse()

		least_squares1 = LeastSquares(newt_chunk, new_chunk2, new_chunk4,line)

		m1 = Minuit(least_squares1, a=0,b=0)

		
		m1.migrad()
		m1.hesse()

		

		
		'''yes = []
		yes1 = []

		for i in range(len(newt_chunk)):
			yes.append(newt_chunk[i]*m.values['b']+m.values['a'])
			yes1.append(newt_chunk[i]*m1.values['b']+m1.values['a'])
		if np.median(new_mag) >16 and np.median(new_mag) < 19:
			plt.errorbar(new_chunk1,new_chunk2,xerr= new_chunk4, yerr = new_chunk3,fmt = 'o',ecolor = 'b',capthick=1,capsize = 2)
			plt.scatter(yes,yes1,color = 'r')
			plt.title('Postions of Predicted (Red) and Measured (Blue)')
			plt.ylabel('Declination offset from base image (arcsec)')
			plt.xlabel('Right ascension offset from base image (arcsec)')
			plt.show()'''

		slope, intercept, r_value, p_value, std_err = stats.linregress(newt_chunk,new_chunk)
		slope1, intercept1, r_value1, p_value1,std_err1 = stats.linregress(newt_chunk,new_chunk1)
		slope2, intercept2, r_value2, p_value2,std_err2 = stats.linregress(newt_chunk,new_chunk2)

		stars = SkyCoord(ra=[np.median(ra_real)]*u.degree,dec=[np.median(dec_real)]*u.degree)
		idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=gra2*u.degree,dec=gdec2*u.degree),0.000277778*u.degree)

		if len(idxs) > 1:
			continue


		'''if m.values['b'] > 1000:

			y = []

			for i in range(len(newt_chunk)):
				y.append(newt_chunk[i]*m.values['b'] + m.values['a'])



			plt.errorbar(newt_chunk,new_chunk1,yerr = new_chunk3,ls = 'none')
			plt.plot(newt_chunk,y)
			plt.show()'''


		chunksra.append(np.median(ra_real))
		chunksdec.append(np.median(dec_real))
		chunkspmra.append(m.values['b'])
		chunkspmdec.append(m1.values['b'])
		chunkspmraerr.append(m1.errors['b'])
		chunkspmdecerr.append(m.errors['b'])
		chunksmag.append(np.median(new_mag))
		chunkssn.append(np.median(sn_chunk))

		
		'''if slope1 * 365 * 3600000> 100:
			
			print(slope1*365*3600000)


			print(ra_chunk)
			print(dec_chunk)

			print(slope1 * 365 * 3600000)
			print(slope2 * 365 * 3600000)
			print(std_err1* 365 * 3600000)
			print(std_err2* 365 * 3600000)'''
			




gmra_matched = []
pmra_matched = []
gmraerr_matched = []
pmraerr_matched = []
sn_matched = []
mag_matched = []



ras = []


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
model1 = []
gmra_matched1 = []
for i in range(len(model)):
	if np.abs(model[i]) < 2:
		count +=1
	if np.abs(sn_matched[i]) < 15:
		model1.append(model[i])
		gmra_matched1.append(gmra_matched[i])
		
print(count/len(idxs))




#plt.errorbar(gmra_matched,model,yerr = pmraerr_matched,xerr = gmraerr_matched,fmt = 'o',ecolor = 'b',capthick=1,capsize = 2)
plt.scatter(gmra_matched1,model1)
#ax.add_artist(plt.scatter(mag2,x2,color = 'r'))
#plt.plot(gmra_matched,gmra_matched)
plt.xlabel("Gaia Proper Motion in RA (mas/yr)")
plt.ylabel("Simga Offset from True Values")
plt.title('Gaia vs ZTF Proper Motions')
plt.show()

