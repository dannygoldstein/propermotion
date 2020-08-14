import sys
import zuds
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.optimize import fmin_l_bfgs_b
from scipy.stats import linregress

PIX_SCALE = 1.012 * u.arcsec

# initialize the database
zuds.init_db()

# get the list of science images to read from the command line
image_list = sys.argv[1]

# load in image paths
image_paths = np.genfromtxt(image_list, encoding='ascii', dtype=None)

# load in images, filter out bad images, and sort by date (ascending)
images = list(
    filter(
        lambda img: img.header['MAGLIM'] > 20.5 and img.header['SEEING'] < 2.5,
        sorted(
            [zuds.ScienceImage.from_file(p) for p in image_paths],
            key=lambda img: img.mjd
        )
    )
)

# make catalogs
for i in images:
    i.mask_image = zuds.MaskImage.from_file(i.local_path.replace('sci', 'msk'))
    i.catalog = zuds.PipelineFITSCatalog.from_image(i)

# define the base image as the chronologically first image, pop it off the stack
base_image = images[0]
images = images[1:]

# get proper motions from gaia
gaia = base_image.gaia_dr2_calibrators()

# keep only the objects in gaia that are in a reasonable magnitude range
gaia = gaia[(gaia['phot_rp_mean_mag'] > 15) & (gaia['phot_rp_mean_mag'] < 21)]

# match the base catalog to gaia
base_coords = SkyCoord(ra=base_image.catalog.data['XWIN_WORLD'],
                       dec=base_image.catalog.data['YWIN_WORLD'],
                       unit='deg')
gaia_coords = SkyCoord(ra=gaia['ra'], dec=gaia['dec'], unit='deg')
idx_gaia, idx_base, _, _ = base_coords.search_around_sky(gaia_coords, 1 * u.arcsec)

# prune the base and gaia catalogs, keeping only the matches in both catalogs
base_catalog = base_image.catalog.data[idx_base]
gaia_catalog = gaia[idx_gaia]
base_coords = base_coords[idx_base]

# match each catalog to the base catalog, storing matches in a dictionary of
# the following form:
#       base_star         : [(   matched_image,          matched_row         )]
#  base_catalog_row_index : [(image_list_index_other, other_catalog_row_index)]


# iterate over the remaining images
match_dictionary = {ind: [] for ind in range(len(base_coords))}
for img_ind, other in enumerate(images):
    match_coords = SkyCoord(ra=other.catalog.data['XWIN_WORLD'],
                            dec=other.catalog.data['YWIN_WORLD'],
                            unit='deg')
    idx_other, idx_base, _, _ = base_coords.search_around_sky(
        match_coords,
        1 * u.arcsec
    )

    for cat_idxother, cat_idxbase in zip(idx_other, idx_base):
        match_dictionary[cat_idxbase].append((img_ind, cat_idxother))

# fit a line to the motion of each star
for base_index in match_dictionary:

    # initialize the output plot
    fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(4, 8))

    # fit RA and Dec separately, then combine
    for key, a in zip(['ra', 'dec'], ax):
        sex_key = 'XWIN_WORLD' if key == 'ra' else 'YWIN_WORLD'
        err_key = 'ERRX2WIN_IMAGE' if key == 'ra' else 'ERRY2WIN_IMAGE'

        # get the ztf position and error along this axis for all matched stars
        y_data_deg = np.asarray([
            images[img_idx].catalog.data[row_idx][sex_key]
            for img_idx, row_idx in match_dictionary[base_index]
        ])

        var_pix2 = np.asarray([
            images[img_idx].catalog.data[row_idx][err_key]
            for img_idx, row_idx in match_dictionary[base_index]
        ])

        # get the mjd of the image of each matched star
        x = np.asarray([images[img_idx].mjd for img_idx, _ in
                        match_dictionary[base_index]])

        # convert from variance to sigma
        sigma_deg = (np.sqrt(var_pix2) * PIX_SCALE).to('deg').value

        def objective_function(parameters):
            """The chi-squared of the model against the data. Minimize this."""
            slope, intercept = parameters
            y_mod_deg = slope * x + intercept
            chi = (y_mod_deg - y_data_deg) / sigma_deg
            return np.sum(chi * chi)

        # minimize the function
        guess = linregress(x, y_data_deg)
        minres = fmin_l_bfgs_b(objective_function, guess[:2], approx_grad=True)
        slope_fit, intercept_fit = minres[0]

        # plot the result
        a.plot(x, slope_fit * x + intercept_fit)
        a.errorbar(x, y_data_deg, yerr=sigma_deg, marker='.', linestyle='none')
        a.set_ylabel(f'{key}_ZTF (deg)')
        a.minorticks_on()
        a.tick_params(which='both', direction='in', bottom=True, top=True,
                      left=True, right=True)

        if key == 'dec':
            a.set_xlabel('MJD [days]')

    fig.tight_layout()
    fig.savefig(f'{base_index:05d}.fit.pdf')
