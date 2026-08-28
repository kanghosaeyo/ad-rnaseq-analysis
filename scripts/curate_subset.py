"""
Parse GSE121212 sample headers into a clean metadata table,
and subset the count matrix to the homogeneous AD/CTRL cohort.

Excludes:
  - All PSO_* columns (psoriasis cohort, not needed for this analysis)
  - AD_032 through AD_037 (the acute/chronic sub-study cohort — a
    methodologically distinct protocol from the main lesional/non-lesional
    design; see project notes for why these are excluded)

Keeps:
  - AD_004 through AD_031: paired lesional / non-lesional samples (21 patients)
  - CTRL_*_healthy: all 38 healthy control samples
"""

from pathlib import Path
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

RAW_COUNTS_PATH = REPO_ROOT / "data" / "raw" / "GSE121212_readcount.txt"
OUT_METADATA_PATH = REPO_ROOT / "data" / "processed" / "metadata.csv"
OUT_COUNTS_PATH = REPO_ROOT / "data" / "processed" / "counts_subset.csv"

# Patients belonging to the separate acute/chronic sub-study,
# excluded to avoid mixing two different biopsy-timing protocols.
CHRONIC_SUBSTUDY_PATIENTS = {"032", "033", "034", "035", "036", "037"}


def parse_sample_name(sample_name: str) -> dict:
    """
    Parse a column name like 'AD_004_lesional' or 'CTRL_002_healthy'
    into its group, patient_id, and condition.
    """
    parts = sample_name.split("_")
    group = parts[0]  # AD, CTRL, or PSO
    patient_number = parts[1]
    patient_id = f"{group}_{patient_number}"
    condition = "_".join(parts[2:])  # handles 'non-lesional', 'chronic_lesion', etc.
    return {
        "sample_id": sample_name,
        "group": group,
        "patient_number": patient_number,
        "patient_id": patient_id,
        "condition": condition,
    }


def build_metadata(columns: list[str]) -> pd.DataFrame:
    records = [parse_sample_name(c) for c in columns]
    return pd.DataFrame(records)


def apply_inclusion_filter(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only:
      - group == 'CTRL' (all healthy controls), or
      - group == 'AD' AND patient_id not in the chronic sub-study set
        AND condition in {'lesional', 'non-lesional'}
    """
    is_ctrl = meta["group"] == "CTRL"
    is_clean_ad = (
        (meta["group"] == "AD")
        & (~meta["patient_number"].isin(CHRONIC_SUBSTUDY_PATIENTS))
        & (meta["condition"].isin(["lesional", "non-lesional"]))
    )
    return meta[is_ctrl | is_clean_ad].reset_index(drop=True)


def main():
    if not RAW_COUNTS_PATH.exists():
        raise FileNotFoundError(
            f"Raw counts file not found at {RAW_COUNTS_PATH}. "
            f"See data/raw/README.md for download instructions."
        )

    counts = pd.read_csv(RAW_COUNTS_PATH, sep="\t", index_col=0)
    print(f"Loaded raw matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")

    meta = build_metadata(counts.columns.tolist())
    meta_clean = apply_inclusion_filter(meta)

    print("\nSample counts by group/condition (clean cohort):")
    print(meta_clean.groupby(["group", "condition"]).size())

    # Sanity check: every AD patient should have exactly 2 samples (paired).
    ad_patient_counts = (
        meta_clean[meta_clean["group"] == "AD"].groupby("patient_id").size()
    )
    unpaired = ad_patient_counts[ad_patient_counts != 2]
    if len(unpaired) > 0:
        print(f"\nWARNING: unpaired AD patients (not exactly 2 samples): "
              f"{unpaired.to_dict()}")
        print("These will be kept but won't contribute to the paired design "
              "term cleanly — review before running DESeq2.")

    counts_clean = counts[meta_clean["sample_id"].tolist()]

    OUT_METADATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    meta_clean.to_csv(OUT_METADATA_PATH, index=False)
    counts_clean.to_csv(OUT_COUNTS_PATH)

    print(f"\nSaved metadata -> {OUT_METADATA_PATH} ({meta_clean.shape[0]} samples)")
    print(f"Saved count matrix -> {OUT_COUNTS_PATH} "
          f"({counts_clean.shape[0]} genes x {counts_clean.shape[1]} samples)")


if __name__ == "__main__":
    main()