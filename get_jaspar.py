import requests, pathlib
s=requests.Session(); s.headers.update({"User-Agent":"POLARIS/1.0"})
urls=[
 "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
 "https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
]
dest=pathlib.Path("data/external/JASPAR2024_CORE_vertebrates_nr.txt")
ok=False
for u in urls:
    try:
        r=s.get(u,timeout=120)
        if r.status_code==200 and len(r.text)>1000:
            dest.write_text(r.text,encoding="utf-8"); ok=True
            print("downloaded",u,"->",len(r.text),"bytes"); break
        print("status",r.status_code,u[:60])
    except Exception as e: print("ERR",str(e)[:80])
if ok:
    txt=dest.read_text(encoding="utf-8")
    print("first 18 lines:")
    print("\n".join(txt.splitlines()[:18]))
    print("total motif headers (>):", sum(1 for l in txt.splitlines() if l.startswith(">")))
