import requests, json, sys
# find latest tectonic release windows asset
r = requests.get("https://api.github.com/repos/tectonic-typesetting/tectonic/releases", timeout=40)
rels = r.json()
asset_url=None; tag=None
for rel in rels:
    for a in rel.get("assets",[]):
        n=a["name"].lower()
        if "windows" in n and ("x86_64" in n or "amd64" in n) and n.endswith(".zip"):
            asset_url=a["browser_download_url"]; tag=rel["tag_name"]; break
    if asset_url: break
print("tag:", tag); print("asset:", asset_url)
if asset_url:
    with open("data/external/tectonic.zip","wb") as f:
        f.write(requests.get(asset_url, timeout=180).content)
    import zipfile, os
    with zipfile.ZipFile("data/external/tectonic.zip") as z:
        z.extractall("data/external/tectonic_bin")
    print("extracted:", os.listdir("data/external/tectonic_bin"))
