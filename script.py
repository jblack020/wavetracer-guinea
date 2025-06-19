import earthaccess
import pathlib

tiles = [
    "N12W016", "N12W015", "N12W014", "N12W013", "N12W012", "N12W011", "N12W010", "N12W009", "N12W008",
    "N11W016", "N11W015", "N11W014", "N11W013", "N11W012", "N11W011", "N11W010", "N11W009", "N11W008",
    "N10W016", "N10W015", "N10W014", "N10W013", "N10W012", "N10W011", "N10W010", "N10W009", "N10W008",
    "N09W015", "N09W014", "N09W013", "N09W012", "N09W011", "N09W010", "N09W009", "N09W008",
    "N08W014", "N08W013", "N08W012", "N08W011", "N08W010", "N08W009", "N08W008",
    "N07W014", "N07W013", "N07W012", "N07W011", "N07W010", "N07W009",
    "N06W012", "N06W011", "N06W010", "N06W009", "N06W008",
    "N05W011", "N05W010", "N05W009", "N05W008",
    "N04W010", "N04W009", "N04W008",
]

# ------------------------------------------------------------------
# 1. Authenticate once (looks in ~/.netrc or prompts interactively)
# ------------------------------------------------------------------
# <-- interactive prompt if needed :contentReference[oaicite:0]{index=0}
earthaccess.login()

# ------------------------------------------------------------------
# 2. Build the HTTPS links you want
#    (LP-DAAC’s on-prem “Data-Pool” keeps SRTMGL3 v003 under this date folder)
# ------------------------------------------------------------------
base = "https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL3.003/2000.02.11"
urls = [f"{base}/{t}.SRTMGL3.hgt.zip" for t in tiles]

# ------------------------------------------------------------------
# 3. Grab the files
#    - pass `provider="LPDAAC"` because we’re giving raw URLs
# ------------------------------------------------------------------
# downloads in parallel  :contentReference[oaicite:1]{index=1}
earthaccess.download(urls, "./srtm_tiles", provider="LPDAAC")
