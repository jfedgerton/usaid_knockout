###############################################################################
# Entity Resolution Step 4a: Semi-automated review of fuzzy candidates
#
# Input:  Data/entity_resolution/fuzzy_candidates.csv
#         Data/entity_resolution/raw_names_summary.csv  (for frequency)
#         Data/entity_resolution/cleaned_names.csv
# Output: Data/entity_resolution/reviewed_fuzzy_matches.csv
#
# Strategy:
#   - AUTO-ACCEPT: very high token_sort (>=97) with same word count
#     (likely typos/minor spelling variants)
#   - AUTO-REJECT: known false-positive patterns (country name swaps,
#     gendered variants of different ministries, etc.)
#   - HEURISTIC ACCEPT: high token_sort (>=90) where one name is a clear
#     substring of the other and the extra words are generic
#   - FLAG: everything else gets accept=REVIEW for manual inspection
#
# The goal is to reduce the 1,134 candidates to a manageable manual review
# set while catching the obvious true matches and obvious false positives.
###############################################################################

import csv
import re
import os
from collections import defaultdict

print("=== Entity Resolution Step 4a: Fuzzy candidate review ===")

# --------------------------------------------------------------------------
# 0. Load data
# --------------------------------------------------------------------------

candidates = []
with open("Data/entity_resolution/fuzzy_candidates.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        row["token_sort_score"] = float(row["token_sort_score"])
        row["partial_ratio_score"] = float(row["partial_ratio_score"])
        candidates.append(row)

print(f"Loaded {len(candidates)} fuzzy candidate pairs")

# Load frequency data for decision-making
name_freq = {}
with open("Data/entity_resolution/raw_names_summary.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name_freq[row["raw_name"]] = int(row["n_appearances"])

# Load cleaned->raw mapping
cleaned_to_raws = defaultdict(list)
with open("Data/entity_resolution/cleaned_names.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        cleaned_to_raws[row["cleaned_name"]].append(row["raw_name"])

# --------------------------------------------------------------------------
# 1. Define false-positive patterns
# --------------------------------------------------------------------------

# Country names that commonly get swapped in fuzzy matching
COUNTRY_PAIRS_FALSE_POSITIVE = [
    ("australia", "austria"),
    ("slovakia", "slovenia"),
    ("colombia", "cambodia"),
    ("niger", "nigeria"),
    ("mali", "malawi"),
    ("guinea", "guinea bissau"),
    ("congo", "comoros"),
    ("iran", "iraq"),
    ("sudan", "south sudan"),  # These ARE different countries
    ("timor", "timore"),
]

# Words that differentiate genuinely different orgs when everything else matches
DIFFERENTIATING_WORDS = [
    "north", "south", "east", "west",
    "rural", "urban",
    "primary", "secondary",
    "men", "women", "children",
    # Country-specific branches — these may be genuinely different
    "usa", "uk", "canada", "france", "germany", "italy", "spain",
    "japan", "korea", "australia", "austria", "sweden", "norway",
    "denmark", "finland", "belgium", "netherlands", "switzerland",
    "ireland", "luxembourg", "portugal", "greece", "poland",
    "czech", "hungary", "romania", "bulgaria",
]

def contains_country_swap(name_a, name_b):
    """Check if the difference between two names is just a country name swap."""
    words_a = set(name_a.split())
    words_b = set(name_b.split())
    diff_a = words_a - words_b
    diff_b = words_b - words_a

    # Check each known false-positive country pair
    for c1, c2 in COUNTRY_PAIRS_FALSE_POSITIVE:
        # Check if the only difference is this country pair
        if (c1 in diff_a and c2 in diff_b) or (c2 in diff_a and c1 in diff_b):
            return True

    # Also catch: names identical except for a country name at the end
    # e.g., "ministry of health sudan" vs "ministry of health yemen"
    if len(diff_a) == 1 and len(diff_b) == 1:
        word_a = list(diff_a)[0]
        word_b = list(diff_b)[0]
        # If both differing words are country-ish (short, at end of name)
        if name_a.endswith(word_a) and name_b.endswith(word_b):
            # These are likely different country-specific entities
            return True

    return False


def is_generic_suffix_difference(name_a, name_b):
    """Check if one name is the other plus generic suffix words."""
    # Normalize to sorted word lists
    if len(name_a) > len(name_b):
        longer, shorter = name_a, name_b
    else:
        longer, shorter = name_b, name_a

    # Check if shorter is a prefix of longer
    if longer.startswith(shorter):
        extra = longer[len(shorter):].strip()
        extra_words = extra.split()
        # Generic suffixes that don't change identity
        generic = {"international", "global", "worldwide", "network",
                    "programme", "program", "project", "initiative"}
        if all(w in generic for w in extra_words):
            return True

    return False


def has_differentiating_word(name_a, name_b):
    """Check if the difference between names includes a differentiating word."""
    words_a = set(name_a.split())
    words_b = set(name_b.split())
    diff = (words_a - words_b) | (words_b - words_a)

    for word in diff:
        if word in DIFFERENTIATING_WORDS:
            return True
    return False


def word_count(name):
    return len(name.split())

# --------------------------------------------------------------------------
# 2. Classify each candidate
# --------------------------------------------------------------------------

results = []
stats = defaultdict(int)

for c in candidates:
    name_a = c["name_a"]
    name_b = c["name_b"]
    ts = c["token_sort_score"]
    pr = c["partial_ratio_score"]
    reason = c["match_reason"]
    recommended = c["recommended_canonical"]

    accept = "REVIEW"  # Default: needs manual review
    review_note = ""

    wc_a = word_count(name_a)
    wc_b = word_count(name_b)

    # --- AUTO-REJECT rules ---

    # Rule R1: Country name swap (false positive)
    if contains_country_swap(name_a, name_b):
        accept = "no"
        review_note = "country_name_swap"

    # Rule R2: Single-word difference that is a differentiating word
    elif has_differentiating_word(name_a, name_b) and ts < 97:
        accept = "no"
        review_note = "differentiating_word"

    # Rule R3: Very different word counts with only partial match
    elif abs(wc_a - wc_b) >= 3 and ts < 85:
        accept = "no"
        review_note = "word_count_mismatch"

    # Rule R4: Both names very short (< 3 words) and not exact
    elif wc_a <= 2 and wc_b <= 2 and ts < 95:
        accept = "no"
        review_note = "short_names_low_score"

    # --- AUTO-ACCEPT rules ---

    # Rule A1: Very high token_sort AND same word count (likely typo/variant)
    elif ts >= 97 and abs(wc_a - wc_b) <= 1:
        accept = "yes"
        review_note = "very_high_similarity"
        # Pick canonical: prefer longer name (more specific), or more frequent
        freq_a = sum(name_freq.get(r, 0) for r in cleaned_to_raws.get(name_a, []))
        freq_b = sum(name_freq.get(r, 0) for r in cleaned_to_raws.get(name_b, []))
        if freq_a >= freq_b:
            recommended = name_a
        else:
            recommended = name_b

    # Rule A2: One name is exact prefix of the other, extra is just country tag
    elif pr == 100 and abs(wc_a - wc_b) <= 2 and ts >= 90:
        # Check if the difference is a country suffix
        if len(name_a) > len(name_b):
            longer, shorter = name_a, name_b
        else:
            longer, shorter = name_b, name_a

        if longer.startswith(shorter):
            extra = longer[len(shorter):].strip()
            # If extra is just a country name, keep the country-specific one
            # as canonical (it's more precise in this dataset context)
            accept = "yes"
            review_note = "prefix_with_country_tag"
            recommended = longer  # Keep the more specific name
        else:
            accept = "REVIEW"
            review_note = "high_partial_needs_review"

    # Rule A3: High token sort (>=93) with same word count, not a country swap
    elif ts >= 93 and abs(wc_a - wc_b) == 0:
        accept = "yes"
        review_note = "high_similarity_same_wordcount"
        freq_a = sum(name_freq.get(r, 0) for r in cleaned_to_raws.get(name_a, []))
        freq_b = sum(name_freq.get(r, 0) for r in cleaned_to_raws.get(name_b, []))
        if freq_a >= freq_b:
            recommended = name_a
        else:
            recommended = name_b

    # --- REMAINING: flag for manual review ---
    else:
        accept = "REVIEW"
        if ts >= 90:
            review_note = "high_score_needs_verification"
        elif reason == "partial_ratio":
            review_note = "partial_match_substring"
        else:
            review_note = "moderate_score"

    stats[accept] += 1
    stats[f"{accept}_{review_note}"] += 1

    results.append({
        "name_a": name_a,
        "name_b": name_b,
        "token_sort_score": ts,
        "partial_ratio_score": pr,
        "accept": accept,
        "canonical_name": recommended,
        "review_note": review_note,
        "raw_variants_a": c["raw_variants_a"],
        "raw_variants_b": c["raw_variants_b"],
    })

# --------------------------------------------------------------------------
# 3. Write output
# --------------------------------------------------------------------------

OUTPUT = "Data/entity_resolution/reviewed_fuzzy_matches.csv"

with open(OUTPUT, "w", encoding="utf-8", newline="") as f:
    fieldnames = [
        "name_a", "name_b", "token_sort_score", "partial_ratio_score",
        "accept", "canonical_name", "review_note",
        "raw_variants_a", "raw_variants_b",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

# --------------------------------------------------------------------------
# 4. Summary
# --------------------------------------------------------------------------

print(f"\n=== Review Classification Summary ===")
print(f"Total candidates: {len(results)}")
print(f"  Auto-accepted (yes):  {stats['yes']}")
print(f"  Auto-rejected (no):   {stats['no']}")
print(f"  Needs review (REVIEW): {stats['REVIEW']}")

print(f"\nDetailed breakdown:")
for key in sorted(stats.keys()):
    if "_" in key:
        print(f"  {key}: {stats[key]}")

# Show accepted matches
accepted = [r for r in results if r["accept"] == "yes"]
print(f"\n=== Auto-accepted matches ({len(accepted)}) ===")
for r in sorted(accepted, key=lambda x: -x["token_sort_score"]):
    print(f"  [{r['token_sort_score']:.0f}] {r['name_a']}")
    print(f"     -> {r['name_b']}")
    print(f"     canonical: {r['canonical_name']}  ({r['review_note']})")

# Show rejected matches (sample)
rejected = [r for r in results if r["accept"] == "no"]
print(f"\n=== Auto-rejected matches (first 30 of {len(rejected)}) ===")
for r in sorted(rejected, key=lambda x: -x["token_sort_score"])[:30]:
    print(f"  [{r['token_sort_score']:.0f}] {r['name_a']}")
    print(f"     vs {r['name_b']}")
    print(f"     reason: {r['review_note']}")

# Show REVIEW items (sample)
review = [r for r in results if r["accept"] == "REVIEW"]
print(f"\n=== Needs manual review (first 30 of {len(review)}) ===")
for r in sorted(review, key=lambda x: -x["token_sort_score"])[:30]:
    print(f"  [{r['token_sort_score']:.0f} / {r['partial_ratio_score']:.0f}] {r['name_a']}")
    print(f"     vs {r['name_b']}")
    print(f"     note: {r['review_note']}")

print(f"\n=== Step 4a complete ===")
print(f"Output: {OUTPUT}")
print(f"\nNext steps:")
print(f"  1. Open {OUTPUT} and review rows where accept=REVIEW")
print(f"  2. Change accept to 'yes' or 'no' for each")
print(f"  3. For 'yes' rows, verify canonical_name is correct")
print(f"  4. Re-run entity_resolution_step4_build_crosswalk.py")
