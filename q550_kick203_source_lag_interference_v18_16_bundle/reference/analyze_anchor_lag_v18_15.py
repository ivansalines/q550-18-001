#!/usr/bin/env python3
"""Q550 v18.15 independent anchor-aware age partition and opposing-lag check.

Reads *unchanged* v18.15 original sendback and authentic nested v18.7
checkpoints. No Q550 ODE integration and no physical-metric claims.
All numerical comparisons use the initial co-moving Euclidean 10D metric.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator

NEST = ['inputs/v18_13_original_sendback.zip',
        'inputs/v18_12_original_sendback.zip',
        'inputs/v18_11_original_sendback.zip',
        'inputs/v18_10_original_sendback.zip',
        'inputs/v18_9_original_sendback.zip',
        'inputs/v18_7_original_sendback.zip']
WIDTHS = (0.5, 0.25, 0.125, 0.0625)


def bytes_member(archive, suffix):
    names = [n for n in archive.namelist() if n.endswith(suffix)]
    if len(names) != 1:
        raise RuntimeError(f'Could not uniquely locate {suffix}: {len(names)} matches')
    return archive.read(names[0])


def verify_manifest(z, base=''):
    name = base + 'SHA256SUMS_OUTPUT.json'
    if name not in z.namelist():
        raise RuntimeError(f'Missing checksums: {name}')
    expected = json.loads(z.read(name))['files']
    for path, sha in expected.items():
        if hashlib.sha256(z.read(base+path)).hexdigest() != sha:
            raise RuntimeError(f'SHA256 mismatch: {base+path}')
    if z.testzip() is not None:
        raise RuntimeError('ZIP CRC mismatch')
    return len(expected)


def load_original(zip_path):
    srcbytes = zip_path.read_bytes()
    with zipfile.ZipFile(io.BytesIO(srcbytes)) as z:
        first = z.namelist()[0].split('/')[0] + '/'
        manifest_count = verify_manifest(z, first)
        summary = json.loads(z.read(first+'results/summary.json'))
        nested = z.read(first+'inputs/v18_14_original_sendback.zip')
    for path in NEST:
        with zipfile.ZipFile(io.BytesIO(nested)) as z:
            # All old sendback archives are original root-relative archives.
            verify_manifest(z)
            nested = z.read(path)
    with zipfile.ZipFile(io.BytesIO(nested)) as z:
        verify_manifest(z)
        data={}
        for tol in ('standard','tight'):
            with np.load(io.BytesIO(z.read(f'checkpoint/{tol}_sign_crossing.npz')),allow_pickle=False) as npz:
                data[tol]={
                    'xi':npz['xi_nodes'].copy(),
                    'A':npz['integral_source_frame_11D'][:,:10].copy()
                }
    return data, summary, hashlib.sha256(srcbytes).hexdigest(), manifest_count


def project(e, x):
    unit = e/np.linalg.norm(e)
    return x - np.sum(x*unit,axis=-1)[...,None]*unit


def calculate(A, xi, z, width, initial):
    interp=PchipInterpolator(xi,A,axis=0)
    e=A[80]
    stop=z if initial=='full' else z-1.0
    ages=np.arange(0,stop,width)
    ages=np.r_[ages, stop]
    slabs=project(e,interp(z-ages[:-1])-interp(z-ages[1:]))
    w=project(e,interp(z))
    rho=np.linalg.norm(w)
    axial=slabs @ (w/rho)
    magnitudes=np.linalg.norm(slabs,axis=1)
    neg=-np.minimum(axial,0).sum()
    denom=np.abs(axial).sum()
    vc=1-rho/magnitudes.sum()
    lc=1-rho/denom
    if np.linalg.norm(slabs.sum(axis=0)-w)>2e-9:
        raise RuntimeError('Projected source closure failed')
    gram=slabs@slabs.T
    gram[np.tril_indices(len(gram))]=np.inf
    a,b=np.unravel_index(np.argmin(gram),gram.shape)
    centers=(ages[:-1]+ages[1:])/2
    qpre=project(e,interp(1)-interp(0))
    if np.linalg.norm(qpre)>2e-9:
        raise RuntimeError('Projected pre-anchor history is not zero')
    return {
        'xi':z,'width':width,'subset':initial, 'number_slabs':len(slabs),
        'negative_axial_mass':float(neg),
        'longitudinal_cancellation':float(lc),
        'vector_cancellation':float(vc),
        'strongest_opposing_age_a':float(centers[a]),
        'strongest_opposing_age_b':float(centers[b]),
        'strongest_pair_age_separation':float(abs(centers[b]-centers[a])),
        'strongest_pair_dot':float(gram[a,b]),
        'normalized_pair_dot_per_width_sq':float(gram[a,b]/(width*width)),
        'pre_anchor_net_norm':float(np.linalg.norm(qpre)),
        'closure_residual':float(np.linalg.norm(slabs.sum(axis=0)-w)),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    data,summary,sha,nverified=load_original(args.input)
    rows=[]
    for tol in ('standard','tight'):
        d=data[tol];xi=d['xi'];A=d['A']
        if xi.shape != (641,) or A.shape != (641,10) or abs(xi[80]-1)>1e-10:
            raise RuntimeError('Expected original Q550 v18.7 source grid')
        for j,z in enumerate(summary['events_xi'][tol],1):
            for width in WIDTHS:
                for subset in ('full','post_anchor'):
                    rows.append({'tolerance':tol,'event':j,
                                 **calculate(A,xi,float(z),width,'full' if subset=='full' else 'post_anchor')})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print(json.dumps({'status':'completed','input_sha256':sha,'source_manifest_count':nverified,
                      'output':str(args.output),'rows':len(rows)},indent=2))

if __name__=='__main__':main()
