#!/usr/bin/env python3
"""NCBI GEO GSE97358 source inspection. Raw data are not committed or uploaded."""
import gzip, hashlib, json, pathlib, time, urllib.request
root = pathlib.Path(__file__).parent
inp = root / "raw"
out = root / "results"
inp.mkdir(exist_ok=True)
out.mkdir(exist_ok=True)
files = ["GSE97358_fibroblasts_samples_info.txt.gz", "GSE97358_168_fibroblast_samples_gene_counts.txt.gz"]
base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97358/suppl/"
alt = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97358&format=file&file="
def acquire(name):
    target = inp / name
    if target.is_file() and target.stat().st_size > 200: return target
    for url in [base + name, alt + name]:
        for attempt in range(3):
            try:
                req = urllib.request.Request(url,headers={"User-Agent":"GSE97358-donor-paired-NB-scientific-reproduction/1.0"})
                with urllib.request.urlopen(req,timeout=120) as response, target.open("wb") as fh:
                    while True:
                        chunk = response.read(1 << 20)
                        if not chunk: break
                        fh.write(chunk)
                with gzip.open(target,"rb") as fh: fh.read(100)
                return target
            except Exception as e:
                print("SOURCE_RETRY",name,type(e).__name__,str(e)[:120],flush=True)
                target.unlink(missing_ok=True)
                time.sleep(2+attempt)
    raise RuntimeError("Official GEO file download unavailable: " + name)
def details(path):
    with gzip.open(path,"rt",encoding="utf-8-sig",errors="replace") as fh:
        lines = [fh.readline().rstrip("\n\r") for _ in range(5)]
    delimiter = max(("\t",",",";"), key=lambda x:lines[0].count(x))
    columns = lines[0].split(delimiter)
    return {"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "compressed_bytes":path.stat().st_size,
            "delimiter":repr(delimiter),
            "header_fields":len(columns),
            "header_first_12":columns[:12],
            "first_data_row_fields":len(lines[1].split(delimiter)),
            "first_data_row_first_8":lines[1].split(delimiter)[:8] if "gene_counts" in path.name else [],
            "second_data_row_first_4":lines[2].split(delimiter)[:4] if "gene_counts" in path.name else [],
            "rows_3_to_5_first_cell":[line.split(delimiter)[0] for line in lines[2:5]]}
summary = {"source":"GSE97358 official GEO supplementary",
           "status":"layout-and-hashes-only; no statistical result",
           "files":{p.name:details(p) for p in map(acquire,files)}}
(out/"source_profile.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
