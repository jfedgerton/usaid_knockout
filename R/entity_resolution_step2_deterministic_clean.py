###############################################################################
# Entity Resolution Step 2: Deterministic regex cleaning
#
# Input:  Data/entity_resolution/raw_names_summary.csv
# Output: Data/entity_resolution/cleaned_names.csv
#         Data/entity_resolution/deterministic_matches.csv
#
# Applies a cleaning pipeline to produce cleaned_name from raw_name:
#   1. Unicode NFKD normalization + strip accents + lowercase
#   2. Strip punctuation
#   3. Collapse whitespace
#   4. Remove legal suffixes
#   5. Expand aid-sector abbreviations (both directions)
#   6. Remove leading "the"
#   7. Final trim
#
# Then groups by cleaned_name to find deterministic merges.
###############################################################################

import csv
import re
import unicodedata
import os
from collections import defaultdict

print("=== Entity Resolution Step 2: Deterministic cleaning ===")

# --------------------------------------------------------------------------
# 0. Paths
# --------------------------------------------------------------------------

INPUT_PATH = "Data/entity_resolution/raw_names_summary.csv"
OUTPUT_CLEANED = "Data/entity_resolution/cleaned_names.csv"
OUTPUT_MATCHES = "Data/entity_resolution/deterministic_matches.csv"

# --------------------------------------------------------------------------
# 1. Abbreviation expansion map
#    Key = lowercase abbreviation or canonical full form
#    Value = canonical cleaned form (what both should map to)
#
#    Strategy: if the lowered name exactly matches an abbreviation key,
#    replace it with the canonical form. The canonical form itself will
#    also pass through cleaning and match.
# --------------------------------------------------------------------------

# Map: abbreviation (lowercase) -> canonical full name (lowercase, no accents)
ABBREV_TO_CANONICAL = {
    "usaid": "united states agency for international development",
    "dfid": "department for international development",
    "fcdo": "foreign commonwealth and development office",
    "giz": "deutsche gesellschaft fur internationale zusammenarbeit",
    "jica": "japan international cooperation agency",
    "afd": "agence francaise de developpement",
    "sida": "swedish international development cooperation agency",
    "cida": "canadian international development agency",
    "norad": "norwegian agency for development cooperation",
    "msf": "medecins sans frontieres",
    "icrc": "international committee of the red cross",
    "ifrc": "international federation of red cross and red crescent societies",
    "who": "world health organization",
    "unicef": "united nations childrens fund",
    "undp": "united nations development programme",
    "unhcr": "united nations high commissioner for refugees",
    "wfp": "world food programme",
    "fao": "food and agriculture organization",
    "iom": "international organization for migration",
    "ocha": "office for the coordination of humanitarian affairs",
    "unfpa": "united nations population fund",
    "bmz": "federal ministry for economic cooperation and development germany",
}

# Additional full-name aliases that should map to the same canonical name
# These handle cases where the OECD CRS uses a specific long-form name
# that differs from the abbreviation expansion target
FULL_NAME_ALIASES = {
    # DFID variants
    "department for international development united kingdom":
        "department for international development",
    # GIZ variants (the CRS uses the English name)
    "german corporation for international cooperation":
        "deutsche gesellschaft fur internationale zusammenarbeit",
    # CIDA / Global Affairs Canada (successor org)
    "global affairs canada":
        "canadian international development agency",
    # FAO full UN name
    "food and agriculture organization of the united nations":
        "food and agriculture organization",
    # World Vision variants
    "world vision international": "world vision",
    "world vision inc": "world vision",
    # FCDO is successor to DFID — keep them separate since they are
    # genuinely different orgs at different time periods.
    # But flag: "foreign commonwealth and development office" stays as-is.
}

# --------------------------------------------------------------------------
# 2. Legal suffixes to strip (whole-word at end of string)
# --------------------------------------------------------------------------

LEGAL_SUFFIXES = [
    "inc", "ltd", "llc", "corp", "plc", "gmbh", "sarl",
    "bv", "nv", "srl", "pty", "pvt", "limited",
    "incorporated", "corporation",
]
# NOTE: Removed "foundation", "association", "sa", "ag", "ev", "co" because:
#   - "foundation" and "association" are often part of org proper names
#     (e.g., "International Development Association", "United Nations Foundation")
#   - "sa", "ag", "co" are too short and cause false removals
#     (e.g., "Energie du Mali SA" is fine but "co" clips valid words)

# Build regex: match any of these as a whole word at the end of the string
LEGAL_SUFFIX_RE = re.compile(
    r"\b(" + "|".join(re.escape(s) for s in LEGAL_SUFFIXES) + r")\s*$"
)

# --------------------------------------------------------------------------
# 3. Cleaning functions
# --------------------------------------------------------------------------

def normalize_unicode(s):
    """NFKD decomposition, strip combining marks (accents), lowercase."""
    nfkd = unicodedata.normalize("NFKD", s)
    stripped = "".join(c for c in nfkd if not unicodedata.combining(c))
    return stripped.lower()


def strip_punctuation(s):
    """Remove commas, periods, hyphens, parentheses, quotes, slashes, ampersands.
    Replace with space to avoid merging words."""
    s = re.sub(r"[,.\-\(\)\"\'/\\&;:!?#@\[\]\{\}]", " ", s)
    return s


def collapse_whitespace(s):
    """Multiple spaces -> single space. Strip leading/trailing."""
    return re.sub(r"\s+", " ", s).strip()


def remove_legal_suffixes(s):
    """Strip trailing legal suffixes (whole-word, case-insensitive)."""
    # Apply repeatedly in case there are stacked suffixes (e.g., "Corp Inc")
    for _ in range(3):
        new_s = LEGAL_SUFFIX_RE.sub("", s).strip()
        if new_s == s:
            break
        s = new_s
    return s


def expand_abbreviation(s):
    """If the cleaned name exactly matches a known abbreviation, expand it.
    Only exact matches to avoid false positives (e.g., 'SIDA' = Spanish AIDS)."""
    s_stripped = s.strip()
    if s_stripped in ABBREV_TO_CANONICAL:
        return ABBREV_TO_CANONICAL[s_stripped]
    return s


def apply_full_name_aliases(s):
    """Map known full-name variants to canonical forms."""
    s_stripped = s.strip()
    if s_stripped in FULL_NAME_ALIASES:
        return FULL_NAME_ALIASES[s_stripped]
    return s


def remove_leading_the(s):
    """Remove 'the' at start of name."""
    return re.sub(r"^the\s+", "", s)


def clean_name(raw_name):
    """Full cleaning pipeline."""
    s = normalize_unicode(raw_name)           # 1. Unicode + lowercase
    s = strip_punctuation(s)                   # 2. Strip punctuation
    s = collapse_whitespace(s)                 # 3. Collapse whitespace
    s = remove_legal_suffixes(s)               # 4. Remove legal suffixes
    s = collapse_whitespace(s)                 # 4b. Clean up after suffix removal
    s = expand_abbreviation(s)                 # 5a. Expand abbreviations
    s = apply_full_name_aliases(s)             # 5b. Map full-name aliases
    s = remove_leading_the(s)                  # 6. Remove leading "the"
    s = collapse_whitespace(s)                 # 7. Final trim
    return s

# --------------------------------------------------------------------------
# 4. Read input
# --------------------------------------------------------------------------

if not os.path.exists(INPUT_PATH):
    print(f"ERROR: {INPUT_PATH} not found. Run Step 1 first.")
    exit(1)

raw_names = []
with open(INPUT_PATH, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        raw_names.append(row["raw_name"])

print(f"Read {len(raw_names)} unique raw names")

# --------------------------------------------------------------------------
# 5. Apply cleaning pipeline
# --------------------------------------------------------------------------

cleaned_pairs = []  # list of (raw_name, cleaned_name)

for raw in raw_names:
    cleaned = clean_name(raw)
    cleaned_pairs.append((raw, cleaned))

# --------------------------------------------------------------------------
# 6. Write cleaned_names.csv
# --------------------------------------------------------------------------

with open(OUTPUT_CLEANED, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["raw_name", "cleaned_name"])
    for raw, cleaned in cleaned_pairs:
        writer.writerow([raw, cleaned])

print(f"Wrote {len(cleaned_pairs)} rows to {OUTPUT_CLEANED}")

# --------------------------------------------------------------------------
# 7. Group by cleaned_name to find deterministic merges
# --------------------------------------------------------------------------

cleaned_to_raw = defaultdict(list)
for raw, cleaned in cleaned_pairs:
    cleaned_to_raw[cleaned].append(raw)

# Only report groups where multiple raw names map to same cleaned name
merges = []
for cleaned, raws in sorted(cleaned_to_raw.items()):
    merges.append({
        "cleaned_name": cleaned,
        "n_raw_variants": len(raws),
        "raw_variants": "|".join(sorted(raws)),
    })

with open(OUTPUT_MATCHES, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["cleaned_name", "n_raw_variants", "raw_variants"])
    writer.writeheader()
    for row in merges:
        writer.writerow(row)

# --------------------------------------------------------------------------
# 8. Report statistics
# --------------------------------------------------------------------------

n_total = len(raw_names)
n_cleaned_unique = len(cleaned_to_raw)
n_merged = sum(1 for raws in cleaned_to_raw.values() if len(raws) > 1)
n_raw_merged = sum(len(raws) for raws in cleaned_to_raw.values() if len(raws) > 1)

print(f"\n=== Deterministic Cleaning Summary ===")
print(f"Input raw names:          {n_total}")
print(f"Unique cleaned names:     {n_cleaned_unique}")
print(f"Names reduced by:         {n_total - n_cleaned_unique}")
print(f"Merged groups (2+ raws):  {n_merged}")
print(f"Raw names in merged groups: {n_raw_merged}")

print(f"\nTop 20 merged groups (most variants):")
multi = [(cleaned, raws) for cleaned, raws in cleaned_to_raw.items() if len(raws) > 1]
multi.sort(key=lambda x: -len(x[1]))
for cleaned, raws in multi[:20]:
    print(f"  [{len(raws)} variants] {cleaned}")
    for r in sorted(raws):
        print(f"    <- {r}")

# Show names that changed (not just case)
print(f"\nSample of names changed by cleaning (first 30):")
changed_count = 0
for raw, cleaned in cleaned_pairs:
    raw_lower = raw.lower().strip()
    if raw_lower != cleaned and changed_count < 30:
        print(f"  {raw}")
        print(f"    -> {cleaned}")
        changed_count += 1

print(f"\n=== Step 2 complete ===")
print(f"Output files:")
print(f"  {OUTPUT_CLEANED}")
print(f"  {OUTPUT_MATCHES}")
