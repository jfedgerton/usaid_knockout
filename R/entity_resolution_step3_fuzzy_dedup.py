###############################################################################
# Entity Resolution Step 3: Fuzzy deduplication using rapidfuzz
#
# Input:  Data/entity_resolution/cleaned_names.csv
#         Data/entity_resolution/raw_names_summary.csv  (for frequency info)
# Output: Data/entity_resolution/fuzzy_candidates.csv
#
# Pairwise fuzzy matching on cleaned names:
#   1. Block by first word to avoid O(n^2)
#   2. token_sort_ratio for word-reordering robustness
#   3. partial_ratio for substring containment
#   4. Minimum length filter (>=8 chars) to avoid short-name false positives
#   5. Output candidate pairs for manual review — NO auto-merging
###############################################################################

import csv
import os
from collections import defaultdict
from rapidfuzz import fuzz

print("=== Entity Resolution Step 3: Fuzzy deduplication ===")

# --------------------------------------------------------------------------
# 0. Config
# --------------------------------------------------------------------------

INPUT_CLEANED = "Data/entity_resolution/cleaned_names.csv"
INPUT_SUMMARY = "Data/entity_resolution/raw_names_summary.csv"
OUTPUT_FUZZY = "Data/entity_resolution/fuzzy_candidates.csv"

MIN_TOKEN_SORT = 85       # Minimum token_sort_ratio to flag as candidate
MIN_PARTIAL = 90           # Minimum partial_ratio to flag as candidate (substring)
MIN_LENGTH = 8             # Only fuzzy-match names with >= 8 chars
PARTIAL_MIN_LENGTH = 12    # partial_ratio only for names >= 12 chars (avoid short substring matches)

# --------------------------------------------------------------------------
# 1. Read cleaned names (deduplicated)
# --------------------------------------------------------------------------

if not os.path.exists(INPUT_CLEANED):
    print(f"ERROR: {INPUT_CLEANED} not found. Run Step 2 first.")
    exit(1)

raw_to_cleaned = {}
cleaned_to_raws = defaultdict(list)

with open(INPUT_CLEANED, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        raw = row["raw_name"]
        cleaned = row["cleaned_name"]
        raw_to_cleaned[raw] = cleaned
        cleaned_to_raws[cleaned].append(raw)

# Unique cleaned names
unique_cleaned = sorted(set(raw_to_cleaned.values()))
print(f"Read {len(raw_to_cleaned)} raw->cleaned mappings")
print(f"Unique cleaned names: {len(unique_cleaned)}")

# --------------------------------------------------------------------------
# 2. Read frequency info for choosing canonical names
# --------------------------------------------------------------------------

name_frequency = {}
if os.path.exists(INPUT_SUMMARY):
    with open(INPUT_SUMMARY, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name_frequency[row["raw_name"]] = int(row["n_appearances"])

def get_frequency(cleaned_name):
    """Get total frequency across all raw variants of a cleaned name."""
    raws = cleaned_to_raws.get(cleaned_name, [])
    return sum(name_frequency.get(r, 0) for r in raws)

# --------------------------------------------------------------------------
# 3. Build blocks by first word
# --------------------------------------------------------------------------

blocks = defaultdict(list)
for name in unique_cleaned:
    first_word = name.split()[0] if name.split() else name
    blocks[first_word].append(name)

n_blocks = len(blocks)
block_sizes = [len(v) for v in blocks.values()]
print(f"Blocking by first word: {n_blocks} blocks")
print(f"  Block size range: {min(block_sizes)}-{max(block_sizes)}")
print(f"  Blocks with >1 name: {sum(1 for s in block_sizes if s > 1)}")

# Also build a second blocking scheme by first 3 characters
# to catch cases where the first word differs but names are similar
blocks_3char = defaultdict(list)
for name in unique_cleaned:
    if len(name) >= 3:
        blocks_3char[name[:3]].append(name)

# --------------------------------------------------------------------------
# 4. Pairwise fuzzy matching within blocks
# --------------------------------------------------------------------------

candidates = []
seen_pairs = set()  # Avoid duplicate pairs from overlapping blocks

def process_block(names_in_block, block_label=""):
    """Compare all pairs within a block."""
    n = len(names_in_block)
    if n < 2:
        return

    for i in range(n):
        name_a = names_in_block[i]
        len_a = len(name_a)

        for j in range(i + 1, n):
            name_b = names_in_block[j]
            len_b = len(name_b)

            # Skip if both names are too short
            if len_a < MIN_LENGTH and len_b < MIN_LENGTH:
                continue

            # Skip if already seen (from another block)
            pair_key = (min(name_a, name_b), max(name_a, name_b))
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)

            # Compute scores
            token_sort = fuzz.token_sort_ratio(name_a, name_b)
            partial = fuzz.partial_ratio(name_a, name_b)

            # Flag as candidate?
            is_token_match = token_sort >= MIN_TOKEN_SORT
            is_partial_match = (
                partial >= MIN_PARTIAL
                and len_a >= PARTIAL_MIN_LENGTH
                and len_b >= PARTIAL_MIN_LENGTH
            )

            if is_token_match or is_partial_match:
                # Choose recommended canonical: prefer more frequent, then shorter
                freq_a = get_frequency(name_a)
                freq_b = get_frequency(name_b)
                if freq_a > freq_b:
                    recommended = name_a
                elif freq_b > freq_a:
                    recommended = name_b
                else:
                    recommended = name_a if len_a <= len_b else name_b

                candidates.append({
                    "name_a": name_a,
                    "name_b": name_b,
                    "token_sort_score": round(token_sort, 1),
                    "partial_ratio_score": round(partial, 1),
                    "recommended_canonical": recommended,
                    "match_reason": "token_sort" if is_token_match else "partial_ratio",
                    "raw_variants_a": "|".join(sorted(cleaned_to_raws.get(name_a, []))),
                    "raw_variants_b": "|".join(sorted(cleaned_to_raws.get(name_b, []))),
                })


# Process first-word blocks
print("\nProcessing first-word blocks...")
processed = 0
for block_key, names_in_block in blocks.items():
    process_block(names_in_block, block_key)
    processed += 1
    if processed % 500 == 0:
        print(f"  Processed {processed}/{n_blocks} blocks, {len(candidates)} candidates so far")

print(f"After first-word blocks: {len(candidates)} candidates")

# Process 3-char blocks (catches cross-first-word matches)
print("Processing 3-char blocks...")
n_3char = len(blocks_3char)
processed = 0
for block_key, names_in_block in blocks_3char.items():
    process_block(names_in_block, block_key)
    processed += 1
    if processed % 500 == 0:
        print(f"  Processed {processed}/{n_3char} 3-char blocks, {len(candidates)} candidates so far")

print(f"After 3-char blocks: {len(candidates)} candidates")

# --------------------------------------------------------------------------
# 5. Cross-block matching for known problem patterns
#    Short canonical org names that may appear as prefixes
# --------------------------------------------------------------------------

print("Running cross-block checks for known org patterns...")

# Orgs that commonly appear as prefix of longer names
PREFIX_PATTERNS = [
    "world vision", "care", "mercy corps", "save the children",
    "red cross", "oxfam", "plan international", "action against hunger",
    "medecins sans frontieres", "doctors without borders",
    "international rescue committee", "catholic relief services",
    "world food programme", "world health organization",
]

for prefix in PREFIX_PATTERNS:
    matching = [n for n in unique_cleaned if n.startswith(prefix) and n != prefix]
    if prefix in unique_cleaned and matching:
        for other in matching:
            pair_key = (min(prefix, other), max(prefix, other))
            if pair_key not in seen_pairs:
                seen_pairs.add(pair_key)
                token_sort = fuzz.token_sort_ratio(prefix, other)
                partial = fuzz.partial_ratio(prefix, other)
                freq_prefix = get_frequency(prefix)
                freq_other = get_frequency(other)
                recommended = prefix if freq_prefix >= freq_other else other
                candidates.append({
                    "name_a": prefix,
                    "name_b": other,
                    "token_sort_score": round(token_sort, 1),
                    "partial_ratio_score": round(partial, 1),
                    "recommended_canonical": recommended,
                    "match_reason": "prefix_check",
                    "raw_variants_a": "|".join(sorted(cleaned_to_raws.get(prefix, []))),
                    "raw_variants_b": "|".join(sorted(cleaned_to_raws.get(other, []))),
                })

print(f"Total candidates after all passes: {len(candidates)}")

# --------------------------------------------------------------------------
# 6. Write output
# --------------------------------------------------------------------------

candidates.sort(key=lambda x: -x["token_sort_score"])

with open(OUTPUT_FUZZY, "w", encoding="utf-8", newline="") as f:
    fieldnames = [
        "name_a", "name_b", "token_sort_score", "partial_ratio_score",
        "recommended_canonical", "match_reason", "raw_variants_a", "raw_variants_b",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in candidates:
        writer.writerow(row)

print(f"\nWrote {len(candidates)} candidate pairs to {OUTPUT_FUZZY}")

# --------------------------------------------------------------------------
# 7. Summary statistics
# --------------------------------------------------------------------------

print(f"\n=== Fuzzy Matching Summary ===")
print(f"Total candidate pairs: {len(candidates)}")

if candidates:
    by_reason = defaultdict(int)
    for c in candidates:
        by_reason[c["match_reason"]] += 1
    for reason, count in sorted(by_reason.items()):
        print(f"  {reason}: {count} pairs")

    scores = [c["token_sort_score"] for c in candidates]
    print(f"\nToken sort score distribution:")
    for threshold in [95, 90, 85, 80, 75]:
        n = sum(1 for s in scores if s >= threshold)
        print(f"  >= {threshold}: {n} pairs")

    print(f"\nTop 30 candidate pairs (highest token_sort_score):")
    for c in candidates[:30]:
        print(f"  [{c['token_sort_score']:.0f} / {c['partial_ratio_score']:.0f}] "
              f"{c['name_a']}")
        print(f"    vs {c['name_b']}")
        print(f"    reason: {c['match_reason']}, "
              f"recommended: {c['recommended_canonical']}")

print(f"\n=== Step 3 complete ===")
print(f"Output: {OUTPUT_FUZZY}")
print(f"\n*** IMPORTANT: Review fuzzy_candidates.csv manually before proceeding ***")
print(f"*** Some high-scoring pairs WILL be false positives ***")
