#!/usr/bin/env python3
"""Fetch or reuse the seven original GEO gzip inputs for Project 30 ED3R1.

Only exact, SHA-256 verified public GEO sources are admitted. No third-party
binary is redistributed in the manuscript/reproducibility ZIP.

Usage:
  python download_sources.py --out source_data
  python download_sources.py --out source_data --from-dir /path/to/existing/source_files
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, time, urllib.request, gzip
from pathlib import Path
from replay import INPUTS, EXPECTED_SHA


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def validate(path:Path,name:str):
    if not path.is_file() or path.stat().st_size<100: raise ValueError(f'Missing or implausibly small source file: {name}')
    with path.open('rb') as f:
        if f.read(2)!=b'\x1f\x8b':raise ValueError('Not a gzip file: '+name)
    with gzip.open(path,'rb') as f:
        for _ in iter(lambda:f.read(1<<20), b''):pass
    digest=sha(path)
    if digest!=EXPECTED_SHA[name]:raise ValueError('SHA-256 mismatch: '+name)
    return digest


def run(out:Path,from_dir:Path|None):
    out.mkdir(parents=True,exist_ok=True)
    results=[]
    for accession,names in INPUTS.items():
        for name in names:
            dest=out/name
            url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/suppl/{name}'
            if dest.exists():
                try:
                    digest=validate(dest,name)
                    status='CACHED_SHA_VERIFIED'
                    results.append(dict(accession=accession,file=name,status=status,bytes=dest.stat().st_size,sha256=digest,official_url=url))
                    print(accession, name, status, flush=True)
                    continue
                except Exception as exc:
                    raise RuntimeError(f'Existing file failed validation: {dest}: {exc}')
            tmp=dest.with_suffix(dest.suffix+'.part')
            try:
                if from_dir is not None:
                    src=from_dir/name
                    if not src.is_file():raise FileNotFoundError(src)
                    shutil.copyfile(src,tmp)
                    status='LOCAL_COPY_SHA_VERIFIED'
                else:
                    req=urllib.request.Request(url,headers={'User-Agent':'Project30 GEO academic verification (public data)'})
                    with urllib.request.urlopen(req,timeout=90) as source, tmp.open('wb') as fd:
                        size=0
                        while True:
                            data=source.read(1<<20)
                            if not data:break
                            size+=len(data)
                            if size>150*1024*1024:raise ValueError('Source unexpectedly large')
                            fd.write(data)
                    status='DOWNLOADED_SHA_VERIFIED'
                digest=validate(tmp,name)
                tmp.replace(dest)
                results.append(dict(accession=accession,file=name,status=status,bytes=dest.stat().st_size,sha256=digest,official_url=url))
                print(accession, name, status,flush=True)
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
    ledger=out/'GEO_SOURCE_ADMISSION.json'
    ledger.write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
    return len(results)

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--out',type=Path,default=Path('source_data'))
    a.add_argument('--from-dir',type=Path,default=None)
    args=a.parse_args()
    n=run(args.out,args.from_dir)
    print('VERIFIED',n,'of',sum(map(len,INPUTS.values())),'official source files; no source included in public code release.')
