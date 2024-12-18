# s16a mass and center fit with NumCosmo

Most of these steps follow [HSC's own guide](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/). To get this data you must register an account on [HSC's website](https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new).

## Getting s16a catalog from HSC

Visit HSC's [CAS search](https://hsc-release.mtk.nao.ac.jp/datasearch/) and use the following SQL query to download the s16a shear catalog:

```
select
 b.*, 
 c.ira, c.idec, c.tract,
 a.ishape_hsm_regauss_e1, a.ishape_hsm_regauss_e2, 
 a.ishape_hsm_regauss_resolution, a.ishape_hsm_regauss_sigma, 
 d.photoz_best, d.photoz_risk_best, d.photoz_std_best,
 e.icmodel_mag, e.icmodel_mag_err, 
 e.detect_is_primary, 
 e.iclassification_extendedness, 
 e.icmodel_flux_flags, 
 e.icmodel_flux, e.icmodel_flux_err, 
 c.iblendedness_abs_flux
from
 s16a_wide.meas2 a
 inner join s16a_wide.weaklensing_hsm_regauss b using (object_id)
 inner join s16a_wide.meas c using (object_id)
 inner join s16a_wide.photoz_frankenz d using (object_id)
 inner join s16a_wide.forced e using (object_id)
```

You must place the catalog on the root of this repository and name it `s16a_shear_catalog.fits`.

See [Mandelbaum 2018](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract) for more information on the catalog.


## Getting the s16a frankenz P(z) catalog

To get the s16a photo-z pdf data, we need to download the fits files for each of the tracts and fields from their [file system](https://hsc-release.mtk.nao.ac.jp/archive/filetree/s16a-shape-catalog/Sirius/). We also need to download the redshift bins definitions from there. To download the wanted `frankenz` P(z), we run the notebook `scrape_pz.ipynb`. Remember to add your credentials in the notebook.


## Matching the shear and cluster catalogs

We are focusing on the cluster analyzed in [Hamana 2020](https://arxiv.org/abs/2004.00170) on this work. The file `hamana_clusters.fits` contains their information. The complete CAMIRA s16a cluster catalog can be downloaded from [Oguri's repository](https://github.com/oguri/cluster_catalogs).

To match the clusters to the shear catalog, run `match_catalogs.ipynb`.


## Calibrating the catalogs

To calibrate the catalogs, we follow [Mandelbaum 2018](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract)  and run `calibrate_catalogs.ipynb`.
This notebook saves two catalogs: one for usage with photo-z best and another with photo-z pdf.


## Creating NumCosmo experiments

Finally, we create the NumCosmo experiment files for each of the clusters running `generate_experiments.ipynb`.


## Fitting cluster mass and center

There are two ways to fit the clusters. The simplest option is to run the numcosmo CLI app. We may also load the experiment files in python and then run the fit. This option provides more customizilability. We provide the notebook `cluster_fitting.ipynb` with examples on how to use both options.
