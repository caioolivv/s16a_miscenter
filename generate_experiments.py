# ## Setup libraries
import sys

from numcosmo_py import nc
from numcosmo_py import ncm
import numpy as np
from astropy.table import Table
from tqdm import tqdm
from scipy.optimize import fsolve

__name__ = "NcContext"

ncm.cfg_init()
ncm.cfg_set_log_handler(lambda msg: sys.stdout.write(msg) and sys.stdout.flush())

cluster_catalog = Table.read("hamana_clusters.fits")
pz_bins = Table.read("pz/s16a_pz_pdf_bins.fits")

for cluster in tqdm(cluster_catalog):
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

    cosmo.param_set_desc("H0", {"fit": False})
    cosmo.param_set_desc("Omegac", {"fit": False})
    cosmo.param_set_desc("Omegab", {"fit": False})
    cosmo.param_set_desc("w", {"fit": False})
    cosmo.param_set_desc("Omegak", {"fit": False})
    prim.param_set_desc("ln10e10ASA", {"fit": False})
    prim.param_set_desc("n_SA", {"fit": False})
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    min_radius = 0.3 * cosmo.h()
    max_radius = 3.0 * cosmo.h()

    dist = nc.Distance.new(6.0)
    halo_mass_summary = nc.HaloCMParam.new(nc.HaloMassSummaryMassDef.CRITICAL, 200.0)
    density_profile = nc.HaloDensityProfileNFW.new(halo_mass_summary)
    surface_mass_density = nc.WLSurfaceMassDensity.new(dist)
    halo_position = nc.HaloPosition.new(dist)

    surface_mass_density.prepare(cosmo)
    halo_position.prepare(cosmo)

    halo_mass_summary.param_set_desc("log10MDelta", {"fit": True})
    halo_position.param_set_desc("ra", {"fit": True})
    halo_position.param_set_desc("dec", {"fit": True})

    ra = cluster["ra"]
    dec = cluster["dec"]
    z = cluster["z"]
    log10mass = cluster["log10M200_wmap"]
    c = cluster["c200_wmap"]

    halo_position["ra"] = ra
    halo_position["dec"] = dec
    halo_position["z"] = z

    sep_Mpc = halo_position.projected_radius_from_ra_dec(
        cosmo, cluster["ra"], cluster["dec"] + cluster["sep"]
    )

    half_box_side = abs(
        fsolve(
            lambda sep: halo_position.projected_radius_from_ra_dec(cosmo, ra, dec + sep)
            - (max_radius + sep_Mpc / 2),
            0.6,
        )[0]
    )

    ra_min = ra - half_box_side / np.cos(np.radians(dec))
    ra_max = ra + half_box_side / np.cos(np.radians(dec))
    dec_min = dec - half_box_side
    dec_max = dec + half_box_side

    if not isinstance(c, np.float64):
        c = 4.0

    halo_mass_summary["log10MDelta"] = log10mass
    halo_mass_summary["cDelta"] = c

    side = abs(
        fsolve(
            lambda sep: halo_position.projected_radius_from_ra_dec(cosmo, ra, dec + sep)
            - 0.4,
            0.6,
        )[0]
    )

    halo_position.param_set_desc(
        "ra",
        {
            "lower-bound": float(ra - side / np.cos(np.radians(dec))),
            "upper-bound": float(ra + side / np.cos(np.radians(dec))),
        },
    )
    halo_position.param_set_desc(
        "dec",
        {
            "lower-bound": float(dec - side),
            "upper-bound": float(dec + side),
        },
    )

    galaxy_position = nc.GalaxySDPositionFlat.new(ra_min, ra_max, dec_min, dec_max)
    galaxy_true_redshift = nc.GalaxySDTrueRedshiftLSSTSRD.new()
    galaxy_redshift_obs = nc.GalaxySDObsRedshiftGauss.new(galaxy_true_redshift)
    galaxy_shape = nc.GalaxySDShapeGaussHSC.new()

    z_data = nc.GalaxySDObsRedshiftData.new(galaxy_redshift_obs)
    p_data = nc.GalaxySDPositionData.new(galaxy_position, z_data)
    s_data = nc.GalaxySDShapeData.new(galaxy_shape, p_data)

    shear_catalog = Table.read(
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_shear_catalog.fits"
    )

    cut_shear_catalog_dict = {col: [] for col in s_data.required_columns()}

    for i in range(len(shear_catalog)):
        radius = halo_position.projected_radius_from_ra_dec(
            cosmo, shear_catalog["ira"][i], shear_catalog["idec"][i]
        )
        sep_Mpc = halo_position.projected_radius_from_ra_dec(
            cosmo, cluster["ra"], cluster["dec"] + cluster["sep"]
        )

        if radius > max_radius or radius < min_radius:
            continue

        cut_shear_catalog_dict["ra"].append(shear_catalog["ira"][i])
        cut_shear_catalog_dict["dec"].append(shear_catalog["idec"][i])
        cut_shear_catalog_dict["epsilon_int_1"].append(0.0)
        cut_shear_catalog_dict["epsilon_int_2"].append(0.0)
        cut_shear_catalog_dict["epsilon_obs_1"].append(shear_catalog["e1"][i])
        cut_shear_catalog_dict["epsilon_obs_2"].append(shear_catalog["e2"][i])
        cut_shear_catalog_dict["sigma_int"].append(
            shear_catalog["ishape_hsm_regauss_derived_rms_e"][i]
        )
        cut_shear_catalog_dict["sigma_obs"].append(
            shear_catalog["ishape_hsm_regauss_derived_sigma_e"][i]
        )
        cut_shear_catalog_dict["z"].append(shear_catalog["photoz_best"][i])
        cut_shear_catalog_dict["zp"].append(shear_catalog["photoz_best"][i])
        cut_shear_catalog_dict["sigma_z"].append(0.05)  # TODO: Get the real value

    wl_obs = nc.GalaxyWLObs.new(
        nc.GalaxyWLObsCoord.CELESTIAL,
        len(cut_shear_catalog_dict["ra"]),
        s_data.required_columns(),
    )

    for i in range(len(cut_shear_catalog_dict["ra"])):
        for key, value in cut_shear_catalog_dict.items():
            wl_obs.set(key, i, value[i])

    data_cluster = nc.DataClusterWL.new()
    data_cluster.set_obs(wl_obs)
    data_cluster.set_init(True)

    mset = ncm.MSet.new_array(
        [
            cosmo,
            density_profile,
            surface_mass_density,
            halo_position,
            galaxy_position,
            galaxy_redshift_obs,
            galaxy_shape,
        ]
    )
    dataset = ncm.Dataset.new_array([data_cluster])
    likelihood = ncm.Likelihood.new(dataset)

    mset.prepare_fparam_map()

    experiment = ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    ser = ncm.Serialize.new(ncm.SerializeOpt.CLEAN_DUP)

    ser.to_binfile(
        dataset,
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_experiment.dataset.gvar",
    )
    ser.dict_str_to_yaml_file(
        experiment,
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_experiment.yaml",
    )


for cluster in tqdm(cluster_catalog):
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

    cosmo.param_set_desc("H0", {"fit": False})
    cosmo.param_set_desc("Omegac", {"fit": False})
    cosmo.param_set_desc("Omegab", {"fit": False})
    cosmo.param_set_desc("w", {"fit": False})
    cosmo.param_set_desc("Omegak", {"fit": False})
    prim.param_set_desc("ln10e10ASA", {"fit": False})
    prim.param_set_desc("n_SA", {"fit": False})
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    min_radius = 0.3 * cosmo.h()
    max_radius = 3.0 * cosmo.h()

    dist = nc.Distance.new(6.0)
    halo_mass_summary = nc.HaloCMParam.new(nc.HaloMassSummaryMassDef.CRITICAL, 200.0)
    density_profile = nc.HaloDensityProfileNFW.new(halo_mass_summary)
    surface_mass_density = nc.WLSurfaceMassDensity.new(dist)
    halo_position = nc.HaloPosition.new(dist)

    surface_mass_density.prepare(cosmo)
    halo_position.prepare(cosmo)

    halo_mass_summary.param_set_desc("log10MDelta", {"fit": True})
    halo_position.param_set_desc("ra", {"fit": True})
    halo_position.param_set_desc("dec", {"fit": True})

    ra = cluster["ra"]
    dec = cluster["dec"]
    z = cluster["z"]
    log10mass = cluster["log10M200_wmap"]
    c = cluster["c200_wmap"]

    halo_position["ra"] = ra
    halo_position["dec"] = dec
    halo_position["z"] = z

    sep_Mpc = halo_position.projected_radius_from_ra_dec(
        cosmo, cluster["ra"], cluster["dec"] + cluster["sep"]
    )

    half_box_side = abs(
        fsolve(
            lambda sep: halo_position.projected_radius_from_ra_dec(cosmo, ra, dec + sep)
            - (max_radius + sep_Mpc / 2),
            0.6,
        )[0]
    )

    ra_min = ra - half_box_side / np.cos(np.radians(dec))
    ra_max = ra + half_box_side / np.cos(np.radians(dec))
    dec_min = dec - half_box_side
    dec_max = dec + half_box_side

    if not isinstance(c, np.float64):
        c = 4.0

    halo_mass_summary["log10MDelta"] = log10mass
    halo_mass_summary["cDelta"] = c

    halo_position["ra"] = ra
    halo_position["dec"] = dec
    halo_position["z"] = z

    side = abs(
        fsolve(
            lambda sep: halo_position.projected_radius_from_ra_dec(cosmo, ra, dec + sep)
            - 0.4,
            0.6,
        )[0]
    )

    halo_position.param_set_desc(
        "ra",
        {
            "lower-bound": float(ra - side / np.cos(np.radians(dec))),
            "upper-bound": float(ra + side / np.cos(np.radians(dec))),
        },
    )
    halo_position.param_set_desc(
        "dec",
        {
            "lower-bound": float(dec - side),
            "upper-bound": float(dec + side),
        },
    )

    galaxy_position = nc.GalaxySDPositionFlat.new(ra_min, ra_max, dec_min, dec_max)
    galaxy_redshift_obs = nc.GalaxySDObsRedshiftPz.new()
    galaxy_shape = nc.GalaxySDShapeGaussHSC.new()

    z_data = nc.GalaxySDObsRedshiftData.new(galaxy_redshift_obs)
    p_data = nc.GalaxySDPositionData.new(galaxy_position, z_data)
    s_data = nc.GalaxySDShapeData.new(galaxy_shape, p_data)

    shear_catalog = Table.read(
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_shear_catalog_pdf.fits"
    )

    cut_shear_catalog_dict = {col: [] for col in s_data.required_columns()}
    cut_shear_catalog_dict["pz"] = []

    for i in range(len(shear_catalog)):
        radius = halo_position.projected_radius_from_ra_dec(
            cosmo, shear_catalog["ira"][i], shear_catalog["idec"][i]
        )
        sep_Mpc = halo_position.projected_radius_from_ra_dec(
            cosmo, cluster["ra"], cluster["dec"] + cluster["sep"]
        )

        if radius > max_radius or radius < min_radius:
            continue

        cut_shear_catalog_dict["ra"].append(shear_catalog["ira"][i])
        cut_shear_catalog_dict["dec"].append(shear_catalog["idec"][i])
        cut_shear_catalog_dict["epsilon_int_1"].append(0.0)
        cut_shear_catalog_dict["epsilon_int_2"].append(0.0)
        cut_shear_catalog_dict["epsilon_obs_1"].append(shear_catalog["e1"][i])
        cut_shear_catalog_dict["epsilon_obs_2"].append(shear_catalog["e2"][i])
        cut_shear_catalog_dict["sigma_int"].append(
            shear_catalog["ishape_hsm_regauss_derived_rms_e"][i]
        )
        cut_shear_catalog_dict["sigma_obs"].append(
            shear_catalog["ishape_hsm_regauss_derived_sigma_e"][i]
        )
        cut_shear_catalog_dict["z"].append(shear_catalog["photoz_best"][i])

        lower_index = 0
        upper_index = len(shear_catalog["P(z)"][i]) - 1

        for j in range(len(shear_catalog["P(z)"][i])):
            if shear_catalog["P(z)"][i][j] / max(shear_catalog["P(z)"][i]) > 1e-3:
                lower_index = j
                break

        for j in range(len(shear_catalog["P(z)"][i]) - 1, -1, -1):
            if shear_catalog["P(z)"][i][j] / max(shear_catalog["P(z)"][i]) > 1e-3:
                upper_index = j
                break

        xv = ncm.Vector.new(upper_index - lower_index + 1)
        yv = ncm.Vector.new(upper_index - lower_index + 1)

        for j in range(0, upper_index - lower_index + 1):
            xv.set(j, np.array(pz_bins["BINS"])[j + lower_index])
            yv.set(j, shear_catalog["P(z)"][i][j + lower_index])

        pz_spline = ncm.SplineCubicNotaknot.new_full(xv, yv, True)

        pz_spline.prepare()

        norm = pz_spline.eval_integ(xv.get(0), xv.get(xv.len() - 1))

        for j in range(0, yv.len()):
            yv.set(j, yv.get(j) / norm)

        pz_spline.prepare()

        cut_shear_catalog_dict["pz"].append(pz_spline)

    wl_obs = nc.GalaxyWLObs.new(
        nc.GalaxyWLObsCoord.CELESTIAL,
        len(cut_shear_catalog_dict["ra"]),
        s_data.required_columns(),
    )

    for i in range(len(cut_shear_catalog_dict["ra"])):
        for key, value in cut_shear_catalog_dict.items():
            if key == "pz":
                continue
            wl_obs.set(key, i, value[i])

        wl_obs.set_pz(i, cut_shear_catalog_dict["pz"][i])

    data_cluster = nc.DataClusterWL.new()
    data_cluster.set_obs(wl_obs)
    data_cluster.set_init(True)

    mset = ncm.MSet.new_array(
        [
            cosmo,
            density_profile,
            surface_mass_density,
            halo_position,
            galaxy_position,
            galaxy_redshift_obs,
            galaxy_shape,
        ]
    )
    dataset = ncm.Dataset.new_array([data_cluster])
    likelihood = ncm.Likelihood.new(dataset)

    mset.prepare_fparam_map()

    experiment = ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    ser = ncm.Serialize.new(ncm.SerializeOpt.CLEAN_DUP)

    ser.to_binfile(
        dataset,
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_experiment_pdf.dataset.gvar",
    )
    ser.dict_str_to_yaml_file(
        experiment,
        f"clusters/{cluster['wl_name']}/{cluster['wl_name']}_experiment_pdf.yaml",
    )
