from pathlib import Path
import json, time
import numpy as np
import pandas as pd
import requests
import freesasa
from Bio.PDB import PDBParser, is_aa

# Isolated audit rerun for the only identified constant error in the historical
# protein workflow: GLN/GLU maximum accessible surface areas were transposed.
# This script holds the AlphaFold v6 coordinates, FreeSASA settings, exposure
# threshold and residue-class definitions fixed and recomputes both the old and
# corrected surface classifications from the same residue-level absolute SASA.

PROBE_RADIUS_A = 1.4
SURFACE_RASA_THRESHOLD = 0.20
PATCH_CA_DISTANCE_A = 10.0
PATCH_AREA_PER_RESIDUE_A2 = 12.5

MAX_ASA_BASE = {
    "ALA":129.0, "ARG":274.0, "ASN":195.0, "ASP":193.0, "CYS":167.0,
    "GLY":104.0, "HIS":224.0, "ILE":197.0, "LEU":201.0, "LYS":236.0,
    "MET":224.0, "PHE":240.0, "PRO":159.0, "SER":155.0, "THR":172.0,
    "TRP":285.0, "TYR":263.0, "VAL":174.0,
}
MAX_ASA_OLD = dict(MAX_ASA_BASE, GLN=223.0, GLU=225.0)
MAX_ASA_CORRECTED = dict(MAX_ASA_BASE, GLN=225.0, GLU=223.0)

STANDARD_AA_3 = set(MAX_ASA_CORRECTED)
POSITIVE = {"ARG","LYS","HIS"}
NEGATIVE = {"ASP","GLU"}
HYDROPHOBIC = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","TYR","PRO"}
AROMATIC = {"PHE","TRP","TYR","HIS"}
POLAR_OR_CHARGED = {"ARG","LYS","HIS","ASP","GLU","ASN","GLN","SER","THR","TYR","CYS"}
HBOND_DONOR = {"ARG","LYS","HIS","ASN","GLN","SER","THR","TYR","TRP","CYS"}
HBOND_ACCEPTOR = {"ASP","GLU","ASN","GLN","HIS","SER","THR","TYR","CYS"}

OUT = Path("audit_rerun/output")
PDBDIR = OUT / "alphafold_v6_pdbs"
OUT.mkdir(parents=True, exist_ok=True)
PDBDIR.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "HOF-protein-GlnGlu-correction-audit/1.0"})


def normalize_bio_resid(residue):
    _, resseq, icode = residue.id
    icode = str(icode).strip()
    return f"{resseq}{icode}" if icode else str(resseq)


def normalize_fs_resid(x):
    return str(x).strip().replace(" ", "")


def get_pdb(acc):
    path = PDBDIR / f"AF-{acc}-F1-model_v6.pdb"
    if path.exists() and path.stat().st_size > 1000:
        return path
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"
    r = session.get(url, timeout=90)
    if r.status_code != 200 or len(r.content) < 1000:
        # Fallback to API only if the frozen v6 URL is unavailable.
        api = session.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=90)
        api.raise_for_status()
        recs = api.json()
        rec = recs[0] if isinstance(recs, list) else recs
        url = rec.get("pdbUrl")
        if not url:
            raise RuntimeError(f"No AlphaFold PDB URL for {acc}")
        r = session.get(url, timeout=90)
        r.raise_for_status()
    path.write_bytes(r.content)
    return path


def parse_residue_map(pdb_path):
    structure = PDBParser(QUIET=True).get_structure("prot", str(pdb_path))
    residue_map = {}
    for model in structure:
        for chain in model:
            chain_id = str(chain.id).strip()
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                resname = residue.get_resname().upper().strip()
                if resname not in STANDARD_AA_3:
                    continue
                resid = normalize_bio_resid(residue)
                ca = residue["CA"].coord.astype(float) if "CA" in residue else None
                residue_map[(chain_id, resid)] = {"resname": resname, "ca_coord": ca}
        break
    return residue_map


def absolute_residue_sasa(pdb_path, residue_map):
    params = freesasa.Parameters({
        "algorithm": freesasa.LeeRichards,
        "probe-radius": PROBE_RADIUS_A,
        "n-slices": 20,
    })
    fs_structure = freesasa.Structure(str(pdb_path))
    fs_result = freesasa.calc(fs_structure, params)
    residue_areas = fs_result.residueAreas()
    records = []
    for chain_id, residues in residue_areas.items():
        chain_id = str(chain_id).strip()
        for resid_key, ra in residues.items():
            resid = normalize_fs_resid(resid_key)
            meta = residue_map.get((chain_id, resid))
            if meta is None:
                cands = [v for (ch, rid), v in residue_map.items() if rid == resid]
                if len(cands) == 1:
                    meta = cands[0]
            if meta is None or meta["resname"] not in STANDARD_AA_3:
                continue
            records.append({
                "chain": chain_id,
                "resid": resid,
                "resname": meta["resname"],
                "sasa": float(ra.total),
                "ca_coord": meta["ca_coord"],
            })
    if not records:
        raise RuntimeError(f"No mapped FreeSASA residue areas for {pdb_path.name}")
    return float(fs_result.totalArea()), records


def connected_components(coords, cutoff):
    if not coords:
        return []
    coords = np.asarray(coords, dtype=float)
    d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    adj = d <= cutoff
    visited = np.zeros(len(coords), dtype=bool)
    comps = []
    for i in range(len(coords)):
        if visited[i]:
            continue
        stack = [i]
        visited[i] = True
        comp = [i]
        while stack:
            u = stack.pop()
            for v in np.where(adj[u])[0]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(int(v))
                    comp.append(int(v))
        comps.append(comp)
    return comps


def patch_stats(surface, residue_set, prefix):
    coords = [r["ca_coord"] for r in surface if r["resname"] in residue_set and r["ca_coord"] is not None]
    comps = connected_components(coords, PATCH_CA_DISTANCE_A)
    sizes = sorted((len(c) for c in comps), reverse=True)
    mx = sizes[0] if sizes else 0
    return {
        f"{prefix}_patch_count": int(len(comps)),
        f"max_{prefix}_patch_size_res": int(mx),
        f"max_{prefix}_patch_area_A2": float(mx * PATCH_AREA_PER_RESIDUE_A2),
    }


def classify(records, max_asa):
    surface = []
    for r in records:
        rasa = r["sasa"] / max_asa[r["resname"]]
        if np.isfinite(rasa) and rasa >= SURFACE_RASA_THRESHOLD:
            surface.append(r)
    if not surface:
        raise RuntimeError("No surface residues")
    n = len(surface)
    total = len(records)
    def frac(s):
        return sum(r["resname"] in s for r in surface) / n
    def density(s):
        return 100.0 * sum(r["resname"] in s for r in surface) / n
    out = {
        "surface_residue_count": int(n),
        "surface_residue_fraction": float(n / total),
        "surface_pos_residue_frac": float(frac(POSITIVE)),
        "surface_neg_residue_frac": float(frac(NEGATIVE)),
        "surface_hydrophobic_frac": float(frac(HYDROPHOBIC)),
        "surface_aromatic_frac": float(frac(AROMATIC)),
        "surface_polar_frac": float(frac(POLAR_OR_CHARGED)),
        "surface_hbond_donor_density": float(density(HBOND_DONOR)),
        "surface_hbond_acceptor_density": float(density(HBOND_ACCEPTOR)),
    }
    out["surface_net_charge_index"] = out["surface_pos_residue_frac"] - out["surface_neg_residue_frac"]
    out.update(patch_stats(surface, POSITIVE, "pos"))
    out.update(patch_stats(surface, NEGATIVE, "neg"))
    out.update(patch_stats(surface, HYDROPHOBIC, "hyd"))
    return out, surface


hist = pd.read_csv("protein_descriptors_FINAL_100_UI_MIN.csv")
assert len(hist) == 100 and hist["uniprot_id"].nunique() == 100

ranking_fields = [
    "surface_pos_residue_frac", "surface_neg_residue_frac",
    "surface_hydrophobic_frac", "surface_aromatic_frac",
    "surface_hbond_donor_density", "surface_hbond_acceptor_density",
]
compare_fields = ranking_fields + ["surface_residue_count", "surface_residue_fraction"]

corrected_rows = []
delta_rows = []
flip_rows = []
repro_rows = []
failures = []

for i, h in hist.iterrows():
    acc = h["uniprot_id"]
    print(f"[{i+1:03d}/100] {acc}", flush=True)
    try:
        pdb_path = get_pdb(acc)
        residue_map = parse_residue_map(pdb_path)
        total_sasa, records = absolute_residue_sasa(pdb_path, residue_map)
        old, old_surface = classify(records, MAX_ASA_OLD)
        new, new_surface = classify(records, MAX_ASA_CORRECTED)
        row = {"protein_name": h["protein_name"], "uniprot_id": acc, "sasa_total_A2": total_sasa,
               "modeled_residue_count_freesasa": len(records), **new}
        corrected_rows.append(row)

        oldkeys = {(r["chain"], r["resid"], r["resname"]) for r in old_surface}
        newkeys = {(r["chain"], r["resid"], r["resname"]) for r in new_surface}
        for key in sorted(oldkeys ^ newkeys):
            rec = next(r for r in records if (r["chain"], r["resid"], r["resname"]) == key)
            flip_rows.append({
                "uniprot_id": acc, "protein_name": h["protein_name"],
                "chain": key[0], "resid": key[1], "resname": key[2],
                "absolute_sasa_A2": rec["sasa"],
                "old_rasa": rec["sasa"] / MAX_ASA_OLD[key[2]],
                "corrected_rasa": rec["sasa"] / MAX_ASA_CORRECTED[key[2]],
                "old_surface": key in oldkeys, "corrected_surface": key in newkeys,
            })

        d = {"protein_name": h["protein_name"], "uniprot_id": acc}
        changed = False
        ranking_changed = False
        for f in sorted(set(new) | set(old)):
            ov, nv = old[f], new[f]
            d[f"old_{f}"] = ov
            d[f"corrected_{f}"] = nv
            if isinstance(ov, (int, np.integer)) and isinstance(nv, (int, np.integer)):
                diff = int(nv) - int(ov)
                is_changed = diff != 0
            else:
                diff = float(nv) - float(ov)
                is_changed = abs(diff) > 1e-12
            d[f"delta_{f}"] = diff
            changed = changed or is_changed
            if f in ranking_fields:
                ranking_changed = ranking_changed or is_changed
        d["any_surface_descriptor_changed"] = changed
        d["any_ranking_input_changed"] = ranking_changed
        delta_rows.append(d)

        rr = {"protein_name": h["protein_name"], "uniprot_id": acc}
        maxerr = 0.0
        for f in compare_fields:
            hv = float(h[f])
            ov = float(old[f])
            err = abs(hv - ov)
            rr[f"abs_error_{f}"] = err
            maxerr = max(maxerr, err)
        rr["max_abs_error"] = maxerr
        rr["historical_reproduced"] = maxerr < 1e-9
        repro_rows.append(rr)
    except Exception as e:
        failures.append({"uniprot_id": acc, "protein_name": h["protein_name"], "error": repr(e)})
        print("FAILED", acc, repr(e), flush=True)
    time.sleep(0.05)

corr = pd.DataFrame(corrected_rows)
delta = pd.DataFrame(delta_rows)
flips = pd.DataFrame(flip_rows)
repro = pd.DataFrame(repro_rows)
fail = pd.DataFrame(failures)

corr.to_csv(OUT / "protein_surface_descriptors_corrected.csv", index=False)
delta.to_csv(OUT / "protein_surface_descriptor_delta.csv", index=False)
flips.to_csv(OUT / "gln_glu_surface_residue_flips.csv", index=False)
repro.to_csv(OUT / "historical_surface_reproduction_check.csv", index=False)
fail.to_csv(OUT / "rerun_failures.csv", index=False)

summary = {
    "requested_proteins": 100,
    "successful_proteins": int(len(corr)),
    "failed_proteins": int(len(fail)),
    "historical_reproduced_at_1e-9": int(repro["historical_reproduced"].sum()) if len(repro) else 0,
    "proteins_with_any_surface_descriptor_change": int(delta["any_surface_descriptor_changed"].sum()) if len(delta) else 0,
    "proteins_with_ranking_input_change": int(delta["any_ranking_input_changed"].sum()) if len(delta) else 0,
    "surface_residue_flips": int(len(flips)),
    "old_constants": {"GLN": 223.0, "GLU": 225.0},
    "corrected_constants": {"GLN": 225.0, "GLU": 223.0},
    "probe_radius_A": PROBE_RADIUS_A,
    "surface_rasa_threshold": SURFACE_RASA_THRESHOLD,
    "patch_CA_distance_A": PATCH_CA_DISTANCE_A,
    "patch_area_per_residue_A2": PATCH_AREA_PER_RESIDUE_A2,
    "coordinate_source": "frozen AlphaFold model_v6 URLs; API fallback only if v6 URL unavailable",
}
(OUT / "protein_surface_correction_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

if len(fail):
    raise SystemExit(f"Rerun incomplete: {len(fail)} protein(s) failed")
