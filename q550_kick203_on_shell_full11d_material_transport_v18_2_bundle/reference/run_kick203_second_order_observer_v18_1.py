#!/usr/bin/env python3
"""Q550 kick 203 v18.1: second-order one-sided B observers and convergence.

A is frozen to the original v17.9 hquarter frames; all M transport comes from
archived v17.9 full11D arrays. Only section B's tangent difference changes.
The first-order branch is compared to the archived v18.0 result. Comparisons
across different B step sizes therefore hold A, M, histories and restart fixed.
This is an off-shell Shapley diagnostic, not physical flip certification.
"""
from __future__ import annotations
import argparse
import io
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'reference'), str(ROOT/'scripts')]
import history_bridge_v16 as old
import q550_state_impulse_response_v13 as v13
from frame_geometry import shapley_three_factors

VERSION = '2026-10-16-Q550-kick203-second-order-one-sided-observer-v18.1'
STENCILS = ('hhalf', 'hquarter', 'heighth', 'hsixteenth')


def read_npz_from_zip(path: Path, basename: str) -> dict:
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith('/' + basename)]
        if len(names) != 1:
            raise RuntimeError(f'{path}: expected exactly one {basename}; found {names}')
        with np.load(io.BytesIO(z.read(names[0])), allow_pickle=False) as d:
            return {k: np.asarray(d[k]) for k in d.files}


def unit(v):
    v = np.asarray(v, float)
    r = float(np.linalg.norm(v))
    return v / r if r > 1e-25 else np.zeros_like(v)


def angle_deg(a, b):
    aa, bb = np.linalg.norm(a), np.linalg.norm(b)
    if min(aa, bb) < 1e-25:
        return float('nan')
    return float(np.degrees(np.arccos(np.clip((a @ b) / (aa * bb), -1.0, 1.0))))


def tangent_differences(history, s, h):
    """Return first/second-order one-sided, and original centered Gram frames.

    The 2h sample for second-order forward/backward is ALWAYS on the same side
    of s as its associated difference. Selecting backward before the restart
    and forward after avoids sampling tangents across the B restart.
    """
    t0 = history.tangent(float(s))
    tp1 = history.tangent(float(s + h))
    tp2 = history.tangent(float(s + 2 * h))
    tm1 = history.tangent(float(s - h))
    tm2 = history.tangent(float(s - 2 * h))
    derivatives = {
        'centered': (tp1 - tm1) / (2.0 * h),
        'forward1': (tp1 - t0) / h,
        'backward1': (t0 - tm1) / h,
        'forward2': (-3.0 * t0 + 4.0 * tp1 - tp2) / (2.0 * h),
        'backward2': (3.0 * t0 - 4.0 * tm1 + tm2) / (2.0 * h),
    }
    return {name: old.gram_frame(t0, d) for name, d in derivatives.items()}


def obtain_crossing_and_heights(base_root, kick):
    events = pd.read_csv(base_root / v13.DEPENDENCIES['events'][0]).set_index('event_n')
    a0 = float(events.loc[390, 'support_center'])
    a1 = float(events.loc[433, 'support_center']) - a0
    b0 = float(events.loc[472, 'support_center'])
    b1 = float(events.loc[513, 'support_center']) - b0
    if min(a1, b1) <= 0:
        raise RuntimeError('invalid section map slopes')
    shared = json.loads((base_root / v13.DEPENDENCIES['shared'][0]).read_text())
    hhalf = 0.125 * float(shared['p_line']) / (int(shared['arc_samples']) - 1)
    heights = {name: hhalf / (2 ** j) for j, name in enumerate(STENCILS)}
    crossing = (a0 + a1 * kick - b0) / b1
    return a0, a1, b0, b1, crossing, heights


def compute(args):
    args.out.mkdir(parents=True, exist_ok=True)
    # Never mix previous addendum arrays even if their filenames look alike.
    archived_manifest = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    for file_path in (args.v179_sendback, args.v180_sendback, args.singular_vectors):
        entry = archived_manifest['archived_input_sha256'].get(file_path.name)
        if entry is None or v13.hash_file(file_path) != entry['sha256']:
            raise RuntimeError(f'Archived input SHA256 mismatch: {file_path}')
    original_hash = v13.hash_file(args.source_bundle_zip)
    if original_hash != v13.SOURCE_SHA:
        raise RuntimeError(f'Original source hash mismatch: {original_hash}')
    provenance = v13.validate_inputs(args.base_root, args.source_bundle_zip)
    archive = read_npz_from_zip(args.v179_sendback, 'full_frames_cocycle_arrays.npz')
    orig_shap = read_npz_from_zip(args.v179_sendback, 'oriented_shapley_arrays.npz')
    prev = read_npz_from_zip(args.v180_sendback, 'observer_shapley_arrays.npz')
    pa = v13.load_module(args.base_root / v13.DEPENDENCIES['audit'][0], 'v181_audit')
    support = v13.load_module(args.base_root / v13.DEPENDENCIES['support'][0], 'v181_support')
    tube = v13.load_module(args.base_root / v13.DEPENDENCIES['tube'][0], 'v181_tube')
    params = pa.ModelParameters(gamma_rho=0.3)
    k0 = (55.0 / 6.0) * (10.0 / 550.0)
    a0, a1, b0, b1, crossing, heights = obtain_crossing_and_heights(args.base_root, args.kick)

    key0 = f'standard__kick{args.kick:g}__eps0__sign+0'
    xi = np.asarray(archive[key0 + '__xi'], float)
    if not np.all(np.diff(xi) > 0):
        raise RuntimeError('non-monotone original xi grid')
    prior_xi = np.asarray(prev['xi'], float)
    idx = np.where((abs(xi - crossing) <= args.micro_radius) |
                   np.isclose(xi, 211.0, atol=1e-9, rtol=0) |
                   np.isclose(xi, 211.25, atol=1e-9, rtol=0) |
                   np.isclose(xi, 212.10, atol=1e-9, rtol=0))[0]
    selected_xi = xi[idx]
    if len(idx) < 20:
        raise RuntimeError(f'insufficient selected samples: {len(idx)}')
    if len(prior_xi) != len(selected_xi) or not np.allclose(prior_xi, selected_xi, atol=1e-11, rtol=0):
        raise RuntimeError('v18.0 and v18.1 xi selections differ: refuse regression')

    sA = a0 + a1 * xi
    sB = b0 + b1 * xi
    kick_support_s = float(a0 + a1 * args.kick)
    stop = float(max(a0 + a1 * xi.max(), b0 + b1 * xi.max()) + max(.1, 2 * max(heights.values())))

    def rhs(s, y):
        return np.asarray(support.support_rhs_float(np.asarray(y, float), kappa0=k0, params=params, audit=pa), float)

    def jac(s, y):
        return np.asarray(tube.float_support_jacobian(np.asarray(y, float), kappa0=k0, params=params, audit=pa, base=support), float)

    orbit = pa.load_orbit(550., args.base_root / v13.DEPENDENCIES['orbit'][0])
    tail = pa.collect_uniform_tail(orbit, tail_periods=70, steps_per_period=160,
                                   samples_per_rotation=720, initial_x_perturbation=1e-6,
                                   params=params, progress_every=1000000)
    archived_spline = pa.build_state_spline(tail)
    seed = max(float(tail.s[0]), float(tail.s[-1]) - .5)
    if seed >= min(sB):
        raise RuntimeError('baseline seed outside original support overlap')
    direction = old.columns_for_vectors(args.singular_vectors, args.kick, [1])[1]

    rows = []
    regress = []
    output_arrays = {'xi': selected_xi, 'crossing_xi': np.array([crossing]),
                     'B_support_heights': np.array([heights[x] for x in STENCILS])}
    max_closure = 0.
    for tol in args.tolerance_sets:
        rt, at = {'standard': (1e-12, 1e-14), 'tight': (3e-13, 3e-15)}[tol]
        baseline = old.dense_ivp(rhs, seed, stop, archived_spline(seed), rt, at)
        y0 = np.asarray(baseline.sol(kick_support_s), float).reshape(11)
        norm_y0 = max(1.0, float(np.linalg.norm(y0)))
        for eps in args.eps:
            variants = {}
            for sign in (+1, -1):
                history_ivp = old.dense_ivp(rhs, kick_support_s, stop,
                                            y0 + sign * eps * norm_y0 * direction, rt, at)
                history = old.History(baseline, history_ivp, kick_support_s, rhs, jac)
                tag = f'{tol}__kick{args.kick:g}__eps{eps:g}__sign{sign:+d}'
                if not np.array_equal(xi, archive[tag + '__xi']):
                    raise RuntimeError(f'xi grid mismatch: {tag}')
                all_B = {}
                for bname, h in heights.items():
                    frames = [tangent_differences(history, float(sB[i]), h) for i in idx]
                    all_B[bname] = {name: np.stack([f[name] for f in frames])
                                    for name in frames[0]}
                    if bname in ('hhalf', 'hquarter'):
                        saved = archive[tag + '__' + bname + '__BB'][idx]
                        err = float(np.max(np.linalg.norm(all_B[bname]['centered'] - saved, axis=(1, 2))))
                        regress.append(dict(tolerance=tol, eps=eps, sign=sign,
                                            comparison='legacy_B_centered_' + bname, max_error=err))
                        if err > args.frame_gate:
                            raise RuntimeError(f'v17.9 central frame regression failed: {tag}, {bname}, {err}')
                variants[sign] = {
                    'A': archive[tag + '__hquarter__BA'][idx],
                    'M': archive[tag + '__M'][idx],
                    'B': all_B,
                }

            key = f'{tol}__kick{args.kick:g}__eps{eps:g}__A_hquarter'
            primary = f'{tol}__kick{args.kick:g}__eps{eps:g}__hquarter'
            # One frozen 11D reference axis from archived, centered section A.
            cent_archive = orig_shap[primary + '__shapley_ABC'][idx]
            anchor_id = np.where((selected_xi > crossing + .025) &
                                 (selected_xi < crossing + .035))[0]
            if len(anchor_id) == 0:
                raise RuntimeError('frozen A anchor missing')
            frozen_A = unit(cent_archive[int(anchor_id[0]), 0])
            output_arrays[key + '__frozen_A_11D'] = frozen_A

            for bname, h in heights.items():
                pairs = {}
                for mode in ('forward1', 'backward1', 'non_straddling1',
                             'forward2', 'backward2', 'non_straddling2'):
                    shares = np.full((len(idx), 3, 11), np.nan)
                    deltas = np.full((len(idx), 11), np.nan)
                    for j, i in enumerate(idx):
                        off = float(xi[i] - crossing)
                        if mode.startswith('non_straddling') and abs(off) < 1e-11:
                            continue
                        if mode.startswith('non_straddling'):
                            side = 'backward' if off < 0 else 'forward'
                            real_mode = side + mode[-1]
                            # Guard that full O(h) or O(2h) support never crosses restart.
                            sb = float(sB[i])
                            max_sample = sb - (h * (2 if mode[-1] == '2' else 1)) if off < 0 else sb + (h * (2 if mode[-1] == '2' else 1))
                            if (sb < kick_support_s and max_sample >= kick_support_s) or (sb > kick_support_s and max_sample <= kick_support_s):
                                raise RuntimeError('non-straddling stencil spans the restart')
                        else:
                            real_mode = mode
                        neg, pos = variants[-1], variants[+1]
                        minus = (neg['A'][j], neg['B'][bname][real_mode][j], neg['M'][j])
                        plus = (pos['A'][j], pos['B'][bname][real_mode][j], pos['M'][j])
                        share, delta, error = shapley_three_factors(minus, plus)
                        if error > args.closure_gate:
                            raise RuntimeError(f'closure failed: {key}/{bname}/{mode}, {error}')
                        shares[j] = share
                        deltas[j] = delta
                        max_closure = max(max_closure, error)
                        A, B, M = share
                        S = A + B
                        along = unit(A)
                        rows.append(dict(tolerance=tol, eps_rel=eps, A_stencil='hquarter',
                                         B_stencil=bname, B_h_support_s=float(h),
                                         B_h_xi=float(h / b1), observer=mode,
                                         xi=float(xi[i]), relative_xi=off,
                                         A_norm=float(np.linalg.norm(A)),
                                         B_norm=float(np.linalg.norm(B)),
                                         S_norm=float(np.linalg.norm(S)),
                                         M_norm=float(np.linalg.norm(M)),
                                         delta_norm=float(np.linalg.norm(delta)),
                                         angle_AB_deg=angle_deg(A, B),
                                         angle_SM_deg=angle_deg(S, M),
                                         AB_cancel=float(np.linalg.norm(S) / max(np.linalg.norm(A) + np.linalg.norm(B), 1e-300)),
                                         S_along_A=float(S @ along),
                                         M_along_A=float(M @ along),
                                         delta_along_A=float(delta @ along),
                                         delta_along_frozen_A=float(delta @ frozen_A),
                                         delta_transverse_frozen_A=float(np.linalg.norm(delta - (delta @ frozen_A) * frozen_A)),
                                         closure=float(error)))
                    pairs[mode] = deltas
                    output_arrays[key + '__B_' + bname + '__' + mode + '__dlocal_delta'] = deltas
                    output_arrays[key + '__B_' + bname + '__' + mode + '__shapley_ABC'] = shares
                if bname == 'hquarter':
                    # Exact regression to original v18.0 for both one-sided first-order differences.
                    for mode in ('forward1', 'backward1', 'non_straddling1'):
                        v18_mode = mode.removesuffix('1') if mode != 'non_straddling1' else 'non_straddling'
                        saved_key = primary + '__' + v18_mode + '__dlocal_delta'
                        if saved_key not in prev:
                            raise RuntimeError(f'v18.0 regression input lacks {saved_key}')
                        err = float(np.nanmax(np.linalg.norm(pairs[mode] - prev[saved_key], axis=1)))
                        regress.append(dict(tolerance=tol, eps=eps, sign='pair',
                                            comparison='v18.0_' + v18_mode + '_B_hquarter', max_error=err))
                        if err > args.shapley_gate:
                            raise RuntimeError(f'v18.0 regression failed: {saved_key} {err}')
                print('computed', tol, f'eps={eps:g}', f'A=hquarter B={bname}',
                      f'points={len(idx)}', flush=True)

    data = pd.DataFrame(rows).sort_values(['tolerance', 'eps_rel', 'B_stencil', 'observer', 'xi'])
    data.to_csv(args.out / 'observer_orders.csv', index=False, float_format='%.16g')
    pd.DataFrame(regress).to_csv(args.out / 'regression_checks.csv', index=False, float_format='%.16g')
    np.savez_compressed(args.out / 'second_order_full11D_arrays.npz', **output_arrays)

    # Convergence diagnostics: for each point and fixed (tol,eps), compare
    # successively finer one-sided approximations of the exact same 11D delta.
    conv = []
    for (tol, eps), _ in data.groupby(['tolerance', 'eps_rel']):
        prefix = f'{tol}__kick{args.kick:g}__eps{eps:g}__A_hquarter'
        for mode in ('non_straddling1', 'non_straddling2'):
            dv = [output_arrays[prefix + '__B_' + b + '__' + mode + '__dlocal_delta']
                  for b in STENCILS]
            for j, x in enumerate(selected_xi):
                ref = dv[3][j]
                if not np.isfinite(ref).all():
                    continue
                diffs = [float(np.linalg.norm(d[j] - ref)) for d in dv[:3]]
                stepdiffs = [float(np.linalg.norm(dv[k][j] - dv[k+1][j])) for k in range(3)]
                refnorm = float(np.linalg.norm(ref))
                conv.append(dict(tolerance=tol, eps_rel=eps, observer=mode,
                                 xi=float(x), relative_xi=float(x-crossing),
                                 reference_delta_norm=refnorm,
                                 diff_hhalf_to_h16=diffs[0], diff_hquarter_to_h16=diffs[1],
                                 diff_h8_to_h16=diffs[2],
                                 normalized_diff_h8_h16=diffs[2]/max(refnorm, 1e-300),
                                 step_difference_hhalf_hquarter=stepdiffs[0],
                                 step_difference_hquarter_h8=stepdiffs[1],
                                 step_difference_h8_h16=stepdiffs[2],
                                 consecutive_ratio_coarse=stepdiffs[0]/max(stepdiffs[1],1e-300),
                                 consecutive_ratio_fine=stepdiffs[1]/max(stepdiffs[2],1e-300),
                                 angle_hhalf_h16_deg=angle_deg(dv[0][j],ref),
                                 angle_hquarter_h16_deg=angle_deg(dv[1][j],ref),
                                 angle_h8_h16_deg=angle_deg(dv[2][j],ref),
                                 signed_projection_h16=float(ref @ output_arrays[prefix+'__frozen_A_11D'])))
    cdf = pd.DataFrame(conv)
    cdf.to_csv(args.out/'residual_convergence.csv', index=False, float_format='%.16g')

    # Direct first-order versus second-order difference at IDENTICAL B steps.
    # This makes their comparison independent of the h16 reference choice.
    cross_order=[]
    for tol in args.tolerance_sets:
      for eps in args.eps:
        prefix=f'{tol}__kick{args.kick:g}__eps{eps:g}__A_hquarter'
        for b in STENCILS:
          a=output_arrays[prefix+'__B_'+b+'__non_straddling1__dlocal_delta']
          d=output_arrays[prefix+'__B_'+b+'__non_straddling2__dlocal_delta']
          axis=output_arrays[prefix+'__frozen_A_11D']
          for j,x in enumerate(selected_xi):
            if not (np.isfinite(a[j]).all() and np.isfinite(d[j]).all()):continue
            cross_order.append(dict(tolerance=tol,eps_rel=eps,B_stencil=b,
                xi=float(x),relative_xi=float(x-crossing),
                delta1_norm=float(np.linalg.norm(a[j])),
                delta2_norm=float(np.linalg.norm(d[j])),
                delta1_vs_2_absolute=float(np.linalg.norm(a[j]-d[j])),
                delta1_vs_2_relative_to_second=float(np.linalg.norm(a[j]-d[j])/max(np.linalg.norm(d[j]),1e-300)),
                delta1_vs_2_angle_deg=angle_deg(a[j],d[j]),
                delta1_frozen_signed=float(a[j]@axis),
                delta2_frozen_signed=float(d[j]@axis)))
    pd.DataFrame(cross_order).to_csv(args.out/'order1_vs_order2.csv',index=False,float_format='%.16g')

    primary = data[(data.tolerance == 'standard') & (data.eps_rel == 1e-8)]
    for suffix, ylabel, field in [
        ('angle', 'angle(A,B) [deg]', 'angle_AB_deg'),
        ('norm', '||Delta||', 'delta_norm'),
        ('signed', 'Delta along frozen A', 'delta_along_frozen_A'),
    ]:
        fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), sharey=(suffix == 'angle'))
        for ax, side, mode in zip(axs, ['before crossing', 'after crossing'], ['backward2','forward2']):
            ss = primary[(primary.relative_xi < 0 if side=='before crossing' else primary.relative_xi > 0)]
            for stencil in STENCILS:
                seg = ss[(ss.B_stencil == stencil) & (ss.observer == mode)]
                ax.plot(seg.relative_xi, seg[field], label=stencil, linewidth=1.2)
            ax.axvline(0, color='gray', linestyle=':', linewidth=.8)
            ax.grid(alpha=.18)
            ax.set_xlim(-args.micro_radius if side=='before crossing' else 0,
                        0 if side=='before crossing' else args.micro_radius)
            ax.set_title(side)
            ax.set_xlabel('xi minus B crossing')
        axs[0].set_ylabel(ylabel)
        axs[1].legend(fontsize='small')
        fig.suptitle(f'Q550 kick203 — second-order B observer / fixed A=hquarter / {suffix}')
        fig.tight_layout()
        fig.savefig(args.out / f'order2_{suffix}_versus_B_step.png', dpi=180)
        plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, mode in zip(axs, ('non_straddling1','non_straddling2')):
        seg = cdf[(cdf.tolerance=='standard') & (cdf.eps_rel==1e-8) & (cdf.observer == mode)]
        for col,label in [('angle_hquarter_h16_deg','hquarter versus h16'),
                          ('angle_h8_h16_deg','h8 versus h16')]:
            ax.plot(seg.relative_xi, seg[col], label=label, linewidth=1.2)
        ax.axvline(0, color='gray', linestyle=':', linewidth=.8)
        ax.set_title('first order' if mode.endswith('1') else 'second order')
        ax.set_xlabel('xi minus B crossing')
        ax.set_ylabel('angle between 11D residuals [deg]')
        ax.grid(alpha=.18)
        ax.set_xlim(-args.micro_radius,args.micro_radius)
        ax.legend(fontsize='small')
    fig.tight_layout()
    fig.savefig(args.out / 'direction_convergence_orders_1_vs_2.png', dpi=180)
    plt.close(fig)

    summaries = []
    for mode in ('non_straddling1','non_straddling2'):
        seg = cdf[(cdf.tolerance=='standard')&(cdf.eps_rel==1e-8)&(cdf.observer==mode)]
        near=seg[(abs(seg.relative_xi) <= .01) & (abs(seg.relative_xi) >= 1e-8)]
        summaries.append(dict(observer=mode,n_near=int(len(near)),
                              max_angle_h8_vs_h16=float(near.angle_h8_h16_deg.max()),
                              median_angle_h8_vs_h16=float(near.angle_h8_h16_deg.median()),
                              median_relative_diff_h8_h16=float(near.normalized_diff_h8_h16.median()),
                              median_fine_difference_ratio=float(near.consecutive_ratio_fine.median())))
    summary = dict(version=VERSION,crossing_xi=float(crossing),
                   sample_count_total=int(len(xi)),sample_count_zoom=int(len(idx)),
                   frozen_A_stencil='hquarter',
                   B_support_heights=heights,
                   B_halfwidth_in_xi={s:float(h/b1) for s,h in heights.items()},
                   source_sha256=original_hash, original_source_provenance=provenance,
                   max_legacy_frame_reproduction=float(pd.DataFrame(regress).max_error.max()),
                   max_shapley_closure=float(max_closure),
                   cross_order_csv="order1_vs_order2.csv",
                   convergence_near_event=summaries,
                   reading=('No two physical histories are recomputed differently across observers. '
                            'Frozen section A and full11D M preserve v17.9 genealogy; only B tangent '
                            'differences change. Odd/even first/second-order numerical convergence of '
                            'the off-shell Shapley residual is tested. Neither the observer-induced '
                            '90-degree band nor a sign change in a moving frame proves physical torsion.'))
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame(summaries).to_csv(args.out/'near_event_convergence_summary.csv',index=False)
    (args.out/'REPORT.md').write_text(
        '# Q550 kick203 v18.1 — one-sided second-order observer convergence\n\n'
        f'B crossing: `xi_B = {crossing:.15f}`. Same frozen A=hquarter and full11D M in all comparisons.\n\n'
        'B observer derivatives (t is the 11D unit tangent):\n\n'
        '- Forward first: `(t(s+h)-t(s))/h`\n'
        '- Backward first: `(t(s)-t(s-h))/h`\n'
        '- Forward second: `(-3t(s)+4t(s+h)-t(s+2h))/(2h)`\n'
        '- Backward second: `(3t(s)-4t(s-h)+t(s-2h))/(2h)`\n\n'
        'Only backward stencils are used before the event and only forward stencils after the event.\n'
        'Four B steps: historical hhalf and hquarter, and extensions heighth and hsixteenth.\n\n'
        'Tests: original central frames reproduce v17.9, and the first-order hquarter residuals '\
        'reproduce v18.0. Shapley is algebraically checked in full11D.\n\n'
        f'Max regression error: {summary["max_legacy_frame_reproduction"]:.4g}.\n'
        f'Max Shapley closure: {summary["max_shapley_closure"]:.4g}.\n\n'
        'No physical torsion/flip certification: off-shell allocation and observer sensitivity remain.\n'
        'See second_order_full11D_arrays.npz, observer_orders.csv, residual_convergence.csv and figures.\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='original_source_provenance'},indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-root', type=Path, default=Path('/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle'))
    parser.add_argument('--source-bundle-zip', type=Path, default=Path('/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip'))
    parser.add_argument('--v179-sendback', type=Path, default=ROOT/'inputs/q550_kick203_event_micro_v17_9_sendback.zip')
    parser.add_argument('--v180-sendback', type=Path, default=ROOT/'inputs/q550_kick203_one_sided_observer_v18_0_sendback.zip')
    parser.add_argument('--singular-vectors', type=Path, default=ROOT/'inputs/singular_direction_vectors.csv')
    parser.add_argument('--out', type=Path, default=ROOT/'results/kick203_v18_1')
    parser.add_argument('--kick', type=float, default=203.)
    parser.add_argument('--micro-radius', type=float, default=.04)
    parser.add_argument('--eps', type=float, nargs='+', default=[1e-8,3e-8])
    parser.add_argument('--tolerance-sets', nargs='+', choices=('standard','tight'), default=['standard','tight'])
    parser.add_argument('--frame-gate', type=float, default=1e-8)
    parser.add_argument('--shapley-gate', type=float, default=1e-8)
    parser.add_argument('--closure-gate', type=float, default=1e-11)
    compute(parser.parse_args())
