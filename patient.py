# John Hancock kyg7rx
import csv
import re
import matplotlib.pyplot as plt
from scipy import stats

# ---------- Data model ----------
class Patient:
    all_patients = []  # collect every created Patient

    def __init__(self, donor_id, apoe, age_onset, age_diag, sex):
        self.donor_id = donor_id
        self.apoe = apoe
        self.age_onset = age_onset
        self.age_diag = age_diag
        self.sex = sex
        Patient.all_patients.append(self)

    def __repr__(self):
        return (f"Patient(ID={self.donor_id}, APOE={self.apoe}, "
                f"Onset={self.age_onset}, Diagnosis={self.age_diag}, Sex={self.sex})")


# ---------- Helpers ----------
def _norm_header(s: str) -> str:
    """Normalize a header for resilient matching (lowercase, alnum only)."""
    return re.sub(r'[^0-9a-z]+', '', (s or '').lower())

def _resolve(headers, candidates_norm_keys):
    """
    Resolve an actual header name from a list of candidate normalized keys.
    Returns the matched real header or None.
    """
    norm_map = {_norm_header(h): h for h in headers if h is not None}
    for key in candidates_norm_keys:
        if key in norm_map:
            return norm_map[key]
    return None

def _clean_missing(val):
    """Map common missing markers to None; strip whitespace."""
    v = (val or "").strip()
    return None if v.lower() in {"", "na", "n/a", "nan", "none", "null", "."} else v

def _clean_number_like(val):
    """
    Parse numbers: '81.0' -> '81' for display; but we return as float for analysis.
    Returns float or None.
    """
    v = _clean_missing(val)
    if v is None:
        return None
    try:
        return float(v)
    except:
        return None


# ---------- Loader (skips rows with NO onset-of-symptoms value) ----------
def load_patients_from_csv(filepath="NO DATE GENOTYPE Metadata (1).csv"):
    """
    Create Patient objects from the CSV, skipping rows with NO onset-of-symptoms value.
    """
    # Clear any prior loaded patients if you re-run
    Patient.all_patients.clear()

    with open(filepath, newline='', encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames or []

        # Resolve columns robustly (tolerate spacing/case variations)
        donor_col = _resolve(headers, ["donorid", "donor_id", "id", "donoridnumber"])
        apoe_col  = _resolve(headers, ["apoegenotype", "apoe", "apoe_genotype", "apoegenotype", "apoe(genotype)"])
        onset_col = _resolve(headers, ["ageofonsetcognitivesymptoms", "ageofonset", "onsetage", "ageatcognitivesymptomonset"])
        diag_col  = _resolve(headers, ["ageofdementiadiagnosis", "ageatdiagnosis", "diagnosisage"])
        sex_col   = _resolve(headers, ["sex", "gender"])

        for row in reader:
            donor_id  = _clean_missing(row.get(donor_col, "") if donor_col else None)
            apoe_raw  = _clean_missing(row.get(apoe_col,  "") if apoe_col  else None)
            age_onset = _clean_number_like(row.get(onset_col, "") if onset_col else None)
            age_diag  = _clean_number_like(row.get(diag_col,  "") if diag_col  else None)
            sex       = _clean_missing(row.get(sex_col,   "") if sex_col   else None)

            # Skip entries with NO onset-of-symptoms value
            if age_onset is None:
                continue

            # Normalize APOE string formatting (e.g., strip spaces)
            apoe = apoe_raw.replace(" ", "") if apoe_raw else None

            Patient(donor_id, apoe, age_onset, age_diag, sex)


# ---------- Grouping, Printing, T-test, and Plot ----------
def apoe_group_label(apoe: str) -> str:
    """
    Per instructions:
      - Group '3/4' and '4/4' as 'E4 Present'
      - All other genotypes -> 'E4 Not Present'
    """
    if not apoe:
        return "E4 Not Present"
    apoe_norm = apoe.strip()
    return "E4 Present" if apoe_norm in {"3/4", "4/4"} else "E4 Not Present"

def print_transformed_data():
    """
    Print the transformed data: Donor ID, E4-group label, Age of Onset.
    """
    print("DonorID\tGroup\t\tAgeOnset")
    for p in Patient.all_patients:
        group = apoe_group_label(p.apoe)
        # Align group label for readability
        pad = "\t" if group == "E4 Present" else "\t"
        print(f"{p.donor_id or ''}\t{group}{pad}{int(p.age_onset) if p.age_onset.is_integer() else p.age_onset}")

def analyze_and_plot():
    """
    Create two groups (E4 Present vs E4 Not Present), run Student's t-test,
    and make a 2-bar plot of mean Age of Onset for each group.
    """
    e4_present = []
    e4_not_present = []

    for p in Patient.all_patients:
        group = apoe_group_label(p.apoe)
        onset = p.age_onset
        if onset is None:
            continue
        if group == "E4 Present":
            e4_present.append(onset)
        else:
            e4_not_present.append(onset)

    # Basic safety checks
    if len(e4_present) == 0 or len(e4_not_present) == 0:
        print("Not enough data to compare groups.")
        print(f"E4 Present count: {len(e4_present)}, E4 Not Present count: {len(e4_not_present)}")
        return

    # Means
    mean_present = sum(e4_present) / len(e4_present)
    mean_not = sum(e4_not_present) / len(e4_not_present)

    # --- Bar plot with two bars (no specific colors per instruction constraints) ---
    labels = ["E4 Present", "E4 Not Present"]
    means = [mean_present, mean_not]

    plt.figure()
    plt.bar(labels, means)
    plt.ylabel("Mean Age of Onset (years)")
    plt.title("APOE (Grouped) vs Age of Onset")
    plt.tight_layout()
    plt.show()

    # --- Student's t-test (Welch's version: robust to unequal variances) ---
    tval, pval = stats.ttest_ind(e4_present, e4_not_present, equal_var=False)

    print("\n--- Statistical Comparison (Student's t-test, Welch) ---")
    print(f"Group sizes: E4 Present n={len(e4_present)}, E4 Not Present n={len(e4_not_present)}")
    print(f"Means: E4 Present = {mean_present:.2f}, E4 Not Present = {mean_not:.2f}")
    print(f"t-value = {tval:.4f}")
    print(f"p-value = {pval:.6f}")
    if pval < 0.05:
        print("Result: Statistically significant difference at α = 0.05")
    else:
        print("Result: Not statistically significant at α = 0.05")


# ---------- Run everything ----------
if __name__ == "__main__":
    # 1) Load patients from the provided CSV, skipping rows with no onset age
    load_patients_from_csv("NO DATE GENOTYPE Metadata (1).csv")

    # 2) Print transformed data (APOE collapsed into two groups; print age of onset)
    print_transformed_data()

    # 3) Make two-bar plot and run t-test
    analyze_and_plot()
