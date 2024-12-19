# # Setup libraries
import os
import requests
from requests.auth import HTTPBasicAuth
from bs4 import BeautifulSoup
from tqdm import tqdm


# # Setup options
#
# Adjust `photoz_method` to choose your preferred photo-z estimate. Also, be sure to add your
# credentials to `auth` before running this notebook.

# Choose your photo-z method
PHOTO_Z_METHOD = "frankenz"
# Insert your HSC user and password
auth = HTTPBasicAuth("user", "pass")

PZ_URL = "https://hsc-release.mtk.nao.ac.jp/archive/filetree/s16a-shape-catalog/"


# # Downloading the data
#
# Now we create the `pz` directory and download the photo-z pdf bins definitions:
pz_bins_request = requests.get(PZ_URL, auth=auth, timeout=10)
pz_bins_soup = BeautifulSoup(pz_bins_request.text, "html.parser")

if not os.path.exists("pz"):
    os.makedirs("pz")

for file in pz_bins_soup.find_all("a"):
    if PHOTO_Z_METHOD in file.get("href"):
        pz_bins_definition_request = requests.get(
            f"{PZ_URL}{file.get('href')}", stream=True, auth=auth, timeout=10
        )

        with open(f"pz/s16a_{PHOTO_Z_METHOD}_bins.fits", "wb") as f:
            for chunk in pz_bins_definition_request.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)


# Finally, we iterate through the HSC fields and download the tract data for the chosen photo-z
# method. Each `fits` object is downloaded in its own field directory inside `pz`.
pz_fields_request = requests.get(PZ_URL + "Sirius/", auth=auth, timeout=10)
pz_fields_soup = BeautifulSoup(pz_fields_request.content, "html.parser")

for field in tqdm(pz_fields_soup.find_all("a")):
    if field["href"] == "../":
        continue

    pz_tracts_request = requests.get(
        PZ_URL + "Sirius/" + field["href"], auth=auth, timeout=10
    )
    pz_tracts_soup = BeautifulSoup(pz_tracts_request.content, "html.parser")

    if not os.path.exists(f"pz/{field['href']}"):
        os.makedirs(f"pz/{field['href']}")

    for file in tqdm(pz_tracts_soup.find_all("a")):
        if "frankenz" not in file["href"]:
            continue

        pz_pdf_request = requests.get(
            PZ_URL + "Sirius/" + field["href"] + file["href"],
            stream=True,
            auth=auth,
            timeout=10,
        )

        with open(f"pz/{field['href']}{file['href']}", "wb") as f:
            for chunk in pz_pdf_request.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
