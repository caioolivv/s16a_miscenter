# # Setup libraries
import os
import sys
import numpy as np
from astropy.table import Table, join, vstack
from tqdm import tqdm
from scipy.optimize import fsolve

from numcosmo_py import nc
from numcosmo_py import ncm

__name__ = "NcContext"

ncm.cfg_init()
ncm.cfg_set_log_handler(lambda msg: sys.stdout.write(msg) and sys.stdout.flush())

# # Load the cluster and shear catalogs
#
# You may choose different cluster catalogs, but be sure that they have the columns: `wl_name`,
# `ra`, and `dec` (you may have to adjust the case of these column names).
cluster_catalog = Table.read("hamana_clusters.fits")
shear_catalog = Table.read("s16a_shear_catalog.fits")


# # Create individual shear catalogs for each cluster
#
# We create a directory `clusters` and inside of it create a new `cluster_name` directory for each
# cluster. Each directory then holds the `fits` file with the raw individual shear catalog. These
# are created by selecting all galaxies within a square of 0.4 degrees centered on the cluster
# coordinates. Then we select the pz catalogs for the present tracts and join the shear and pz
# catalogs using the `object_id` column.
if not os.path.exists("clusters"):
    os.makedirs("clusters")

for cluster in tqdm(cluster_catalog):
    ra = cluster["ra"]
    dec = cluster["dec"]
    z = cluster["z"]

    cosmo = nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.omega_x2omega_k()

    cosmo["H0"] = 69.7
    cosmo["Omegab"] = 0.0464
    cosmo["Omegac"] = 0.235
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.00

    prim = nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745
    prim["n_SA"] = 0.9660

    reion = nc.HIReionCamb.new()

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    dist = nc.Distance.new(6.0)
    halo_position = nc.HaloPosition.new(dist)

    halo_position.prepare(cosmo)

    halo_position["ra"] = ra
    halo_position["dec"] = dec
    halo_position["z"] = z

    half_box_side = fsolve(
        lambda sep: halo_position.projected_radius_from_ra_dec(cosmo, ra, dec + sep / 2)
        - 5.0,
        0.5,
    )[0]

    cluster_shear_catalog = shear_catalog[
        (np.abs(shear_catalog["ira"] - ra) < half_box_side / np.cos(np.radians(dec)))
        & (np.abs(shear_catalog["idec"] - dec) < half_box_side)
    ]

    pz_catalogs = []

    for tract in np.unique(cluster_shear_catalog["tract"]):
        for root, dirs, files in os.walk("pz"):
            for f in files:
                if str(tract) in f:
                    _ = Table.read(root + "/" + f)
                    pz_catalogs.append(_)

    pz_catalog = vstack(pz_catalogs)

    cluster_shear_catalog = join(pz_catalog, cluster_shear_catalog, keys="object_id")

    if not os.path.exists(f"clusters/{cluster['wl_name']}"):
        os.makedirs(f"clusters/{cluster['wl_name']}")

    cluster_shear_catalog.write(
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_raw_shear_catalog.fits",
        overwrite=True,
    )
