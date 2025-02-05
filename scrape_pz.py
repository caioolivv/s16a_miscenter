# # Setup libraries
import os
import subprocess
import requests
from requests.auth import HTTPBasicAuth
from bs4 import BeautifulSoup
from tqdm import tqdm

PHOTO_Z_METHOD = "mlz"
auth = HTTPBasicAuth("user", "pass")
PZ_URL = "https://hsc-release.mtk.nao.ac.jp/archive/filetree/s16a-shape-catalog/"

pz_bins_request = requests.get(PZ_URL, auth=auth, timeout=10)
pz_bins_soup = BeautifulSoup(pz_bins_request.text, "html.parser")

if not os.path.exists("pz"):
    os.makedirs("pz")

for file in pz_bins_soup.find_all("a"):
    if PHOTO_Z_METHOD in file.get("href"):
        pz_bins_definition_request = requests.get(
            f"{PZ_URL}{file.get('href')}", stream=True, auth=auth, timeout=10
        )

        with open("pz/s16a_pz_pdf_bins.fits", "wb") as f:
            for chunk in pz_bins_definition_request.iter_content(
                chunk_size=10 * 1024 * 1024
            ):
                if chunk:
                    f.write(chunk)

pz_fields_request = requests.get(PZ_URL + "Sirius/", auth=auth, timeout=10)
pz_fields_soup = BeautifulSoup(pz_fields_request.content, "html.parser")

URIS_FILE = "pz/uris.txt"

with open(URIS_FILE, "w", encoding="UTF8") as f:
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
            if PHOTO_Z_METHOD not in file["href"]:
                continue

            fits_url = PZ_URL + "Sirius/" + field["href"] + file["href"]

            file_name = f"pz/{field['href']}{file['href'].split('_')[0]}_pz_pdf.fits"

            f.write(fits_url + "\n")
            f.write(f"  out={file_name}\n")

subprocess.run(
    [
        "aria2c",
        f"--http-user={auth.username}",
        f"--http-passwd={auth.password}",
        "--max-connection-per-server=16",
        "--split=16",
        "--max-tries=100",
        "--retry-wait=10",
        f"--input-file={URIS_FILE}",
    ],
    check=True,
)

subprocess.run(["rm", URIS_FILE], check=True)
