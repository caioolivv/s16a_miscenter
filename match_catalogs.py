# # Setup libraries
import os
import numpy as np
from astropy.table import Table, join, vstack
from tqdm import tqdm


# # Load the cluster and shear catalogs
#
# You may choose different cluster catalogs, but be sure that they have the columns: `name`, `ra`,
# and `dec` (you may have to adjust the case of these column names).
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

    cluster_shear_catalog = shear_catalog[
        (np.abs(shear_catalog["ira"] - ra) < 0.2 / np.cos(np.radians(dec)))
        & (np.abs(shear_catalog["idec"] - dec) < 0.2)
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

    if not os.path.exists(f"clusters/{cluster['name']}"):
        os.makedirs(f"clusters/{cluster['name']}")

    cluster_shear_catalog.write(
        f"clusters/{cluster['name']}/{cluster['name']}_raw_shear_catalog.fits",
        overwrite=True,
    )
