# John Hancock kyg7rx
import csv
import re
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

# ---------- Data model ----------
class Patient:
    all_patients = []  # collect every created Patient

    def __init__(self, donor_id, apoe, age_onset, age_diag, sex, education):
        self.donor_id = donor_id
        self.apoe = apoe
        self.age_onset = age_onset
        self.age_diag = age_diag
        self.sex = sex
        self.years_education = education
        Patient.all_patients.append(self)

    def __repr__(self):
        return (f"Patient(ID={self.donor_id}, APOE={self.apoe}, "
                f"Onset={self.age_onset}, Diagnosis={self.age_diag}, "
                f"Sex={self.sex}, Education={self.years_education})")


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
    Parse numbers: '81.0' -> 81.0 (float). Returns float or None.
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
        apoe_col  = _resolve(headers, ["apoegenotype", "apoe", "apoe_genotype", "apoe(genotype)"])
        onset_col = _resolve(headers, ["ageofonsetcognitivesymptoms", "ageofonset", "onsetage", "ageatcognitivesymptomonset"])
        diag_col  = _resolve(headers, ["ageofdementiadiagnosis", "ageatdiagnosis", "diagnosisage"])
        sex_col   = _resolve(headers, ["sex", "gender"])
        education_col = _resolve(headers, ["education", "yearsofeducation", "yoe", "education(years)"])

        for row in reader:
            donor_id  = _clean_missing(row.get(donor_col, "") if donor_col else None)
            apoe_raw  = _clean_missing(row.get(apoe_col,  "") if apoe_col  else None)
            age_onset = _clean_number_like(row.get(onset_col, "") if onset_col else None)
            age_diag  = _clean_number_like(row.get(diag_col,  "") if diag_col  else None)
            sex       = _clean_missing(row.get(sex_col,   "") if sex_col   else None)
            education = _clean_number_like(row.get(education_col, "") if education_col else None)

            # Skip entries with NO onset-of-symptoms value
            if age_onset is None:
                continue

            # Normalize APOE string formatting (e.g., strip spaces)
            apoe = apoe_raw.replace(" ", "") if apoe_raw else None

            Patient(donor_id, apoe, age_onset, age_diag, sex, education)


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
        age_disp = int(p.age_onset) if (p.age_onset is not None and float(p.age_onset).is_integer()) else p.age_onset
        print(f"{p.donor_id or ''}\t{group}{pad}{age_disp}")

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

# Standard deviation (instead of SEM)
    sd_present = np.std(e4_present, ddof=1)
    sd_not = np.std(e4_not_present, ddof=1)

# --- Bar plot with two bars (no specific colors) ---
    labels = ["E4 Present", "E4 Not Present"]
    means = [mean_present, mean_not]
    errors = [sd_present, sd_not]  # use SD

    plt.figure()
    plt.bar(labels, means, yerr=errors, capsize=5)
    plt.ylabel("Mean Age of Onset (years)")
    plt.title("APOE (Grouped) vs Age of Onset")
    plt.tight_layout()
    plt.show()

    # --- Student's t-test (Welch) ---
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


# ---------- Scatter + Linear Regression ----------
def plot_education_vs_onset():
    """
    Scatter of Years of Education vs. Age of Onset with a linear regression line.
    Prints slope, intercept, R^2, p-value, and std error.
    """
    x = []  # years of education
    y = []  # age of onset

    for p in Patient.all_patients:
        if p.years_education is not None and p.age_onset is not None:
            x.append(p.years_education)
            y.append(p.age_onset)

    if len(x) == 0:
        print("No data available for education vs. onset.")
        return

    # Need at least 2 points for a regression
    if len(x) < 2:
        print("Not enough points for regression (need at least 2).")
        return

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    # --- Linear regression ---
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    # Use a smooth line across the span of x for nicer display
    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = slope * x_line + intercept
    r2 = r_value ** 2

    # --- Plot scatter + regression line ---
    plt.figure()
    plt.scatter(x, y, label="Data")
    plt.plot(x_line, y_line, label=f"Fit: y = {slope:.2f}x + {intercept:.2f}")
    plt.xlabel("Years of Education")
    plt.ylabel("Age of Onset (years)")
    plt.title("Education vs. Age of Onset (with Linear Regression)")
    plt.legend()
    # Put R^2 on the plot
    plt.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=plt.gca().transAxes,
             ha="left", va="top")
    plt.tight_layout()
    plt.show()

    # --- Report stats ---
    print("\n--- Linear Regression: Education vs. Age of Onset ---")
    print(f"Slope = {slope:.4f}")
    print(f"Intercept = {intercept:.4f}")
    print(f"R-squared = {r2:.4f}")


# ---------- Run everything ----------
if __name__ == "__main__":
    # 1) Load patients from the provided CSV, skipping rows with no onset age
    load_patients_from_csv("NO DATE GENOTYPE Metadata (1).csv")

    # 2) Print transformed data (APOE collapsed into two groups; print age of onset)
    print_transformed_data()

    # 3) Make two-bar plot and run t-test
    analyze_and_plot()

    # 4) Scatter + regression with R^2
    plot_education_vs_onset()
