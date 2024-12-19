# # Setup libraries
import os
import numpy as np
from astropy.table import Table, join, vstack
from tqdm import tqdm


# # Cleaning and calibration functions
#
# These functions follow the cuts and calibration defined in Mandelbaum 2017. The code for them was taken from [CLMM's](https://github.com/LSSTDESC/CLMM) `Example4_Fit_Halo_mass_to_HSC_data.ipynb`. All credit belongs to the original authors.
#
# Observation: the cut based on `photoz_risk_best` is not done for P(z) pdf catalog.
def make_cuts(catalog_in, is_pdf):
    # We consider some cuts in Mandelbaum et al. 2018 (HSC SSP Y1 shear catalog).
    select = catalog_in["detect_is_primary"] == True
    select &= catalog_in["icmodel_flux_flags"] == False
    select &= catalog_in["iclassification_extendedness"] > 0.5
    select &= catalog_in["icmodel_mag_err"] <= 2.5 / np.log(10.0) / 10.0
    select &= (
        catalog_in["ishape_hsm_regauss_e1"] ** 2
        + catalog_in["ishape_hsm_regauss_e2"] ** 2
        < 4.0
    )
    select &= catalog_in["icmodel_mag"] <= 24.5
    select &= catalog_in["iblendedness_abs_flux"] < (10 ** (-0.375))
    select &= (
        catalog_in["ishape_hsm_regauss_resolution"] >= 0.3
    )  # similar to extendedness
    select &= catalog_in["ishape_hsm_regauss_sigma"] <= 0.4
    # Note "zbest" minimizes the risk of the photo-z point estimate being far away from the true value.
    # Details: https://hsc-release.mtk.nao.ac.jp/doc/wp-content/uploads/2017/02/pdr1_photoz_release_note.pdf
    if not is_pdf:
        select &= catalog_in["photoz_risk_best"] < 0.5

    catalog_out = catalog_in[select]

    return catalog_out


# Reference: Mandelbaum et al. 2018 "The first-year shear catalog of the Subaru Hyper Suprime-Cam Subaru Strategic Program Survey".
# Section A.3.2: "per-object galaxy shear estimate".
def apply_shear_calibration(catalog_in):
    e1_0 = catalog_in["ishape_hsm_regauss_e1"]
    e2_0 = catalog_in["ishape_hsm_regauss_e2"]
    e_rms = catalog_in["ishape_hsm_regauss_derived_rms_e"]
    m = catalog_in["ishape_hsm_regauss_derived_shear_bias_m"]
    c1 = catalog_in["ishape_hsm_regauss_derived_shear_bias_c1"]
    c2 = catalog_in["ishape_hsm_regauss_derived_shear_bias_c2"]
    weight = catalog_in["ishape_hsm_regauss_derived_shape_weight"]

    R = 1.0 - np.sum(weight * e_rms**2.0) / np.sum(weight)
    m_mean = np.sum(weight * m) / np.sum(weight)

    g1 = (e1_0 / (2.0 * R) - c1) / (1.0 + m_mean)
    g2 = (e2_0 / (2.0 * R) - c2) / (1.0 + m_mean)

    return g1, g2


for root, dirs, files in os.walk("clusters"):
    for file in files:
        if file.endswith("_raw_shear_catalog.fits"):
            raw_catalog = Table.read(os.path.join(root, file))

            catalog = make_cuts(raw_catalog, is_pdf=False)
            catalog_pdf = make_cuts(raw_catalog, is_pdf=True)
            g1, g2 = apply_shear_calibration(catalog)
            g1_pdf, g2_pdf = apply_shear_calibration(catalog_pdf)
            catalog["e1"] = g1
            catalog["e2"] = g2
            catalog_pdf["e1"] = g1_pdf
            catalog_pdf["e2"] = g2_pdf

            catalog.write(
                os.path.join(
                    root, file.replace("_raw_shear_catalog", "_shear_catalog")
                ),
                overwrite=True,
            )
            catalog_pdf.write(
                os.path.join(
                    root, file.replace("_raw_shear_catalog", "_shear_catalog_pdf")
                ),
                overwrite=True,
            )
