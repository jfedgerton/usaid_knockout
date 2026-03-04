###############################################################################
# Entity Resolution Step 4: Build initial crosswalk from deterministic matches
#
# Input:  Data/entity_resolution/cleaned_names.csv
#         Data/entity_resolution/deterministic_matches.csv
#         Data/entity_resolution/fuzzy_candidates.csv (for reference)
# Output: Data/entity_resolution/name_crosswalk.csv
#
# This script builds the crosswalk in two parts:
#   1. All raw names get their cleaned_name as canonical (match_type = "exact"
#      if 1:1, "deterministic" if many:1)
#   2. A reviewed_fuzzy_matches.csv file (you create manually) can be read
#      to add fuzzy merges
#
# After running this, you'll manually review fuzzy_candidates.csv, create
# reviewed_fuzzy_matches.csv, and re-run this script.
###############################################################################

import csv
import os
from collections import defaultdict

print("=== Entity Resolution Step 4: Building crosswalk ===")

# --------------------------------------------------------------------------
# 1. Read deterministic cleaning results
# --------------------------------------------------------------------------

cleaned_names = {}  # raw_name -> cleaned_name
with open("Data/entity_resolution/cleaned_names.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        cleaned_names[row["raw_name"]] = row["cleaned_name"]

print(f"Loaded {len(cleaned_names)} raw->cleaned mappings")

# Group cleaned names by their cleaned form
cleaned_groups = defaultdict(list)
for raw, cleaned in cleaned_names.items():
    cleaned_groups[cleaned].append(raw)

# --------------------------------------------------------------------------
# 2. Build base crosswalk from deterministic matches
# --------------------------------------------------------------------------

crosswalk = []

for cleaned, raws in sorted(cleaned_groups.items()):
    match_type = "deterministic" if len(raws) > 1 else "exact"
    for raw in raws:
        crosswalk.append({
            "raw_name": raw,
            "canonical_name": cleaned,
            "match_type": match_type,
        })

print(f"Base crosswalk: {len(crosswalk)} entries")
print(f"  exact:          {sum(1 for c in crosswalk if c['match_type'] == 'exact')}")
print(f"  deterministic:  {sum(1 for c in crosswalk if c['match_type'] == 'deterministic')}")

# --------------------------------------------------------------------------
# 3. Apply fuzzy matches if reviewed file exists
# --------------------------------------------------------------------------

FUZZY_REVIEWED = "Data/entity_resolution/reviewed_fuzzy_matches.csv"

if os.path.exists(FUZZY_REVIEWED):
    print(f"\nFound {FUZZY_REVIEWED} — applying fuzzy merges")

    # Expected columns: name_a, name_b, accept (yes/no), canonical_name
    # If accept=yes, both name_a and name_b map to canonical_name
    fuzzy_merges = {}  # old_cleaned -> new_canonical

    with open(FUZZY_REVIEWED, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("accept", "").lower().strip() in ("yes", "y", "true", "1"):
                canonical = row["canonical_name"].strip()
                name_a = row["name_a"].strip()
                name_b = row["name_b"].strip()

                # Map both to the canonical
                if name_a != canonical:
                    fuzzy_merges[name_a] = canonical
                if name_b != canonical:
                    fuzzy_merges[name_b] = canonical

    print(f"  {len(fuzzy_merges)} fuzzy merge mappings loaded")

    # Update crosswalk: if a raw_name's cleaned_name is in fuzzy_merges,
    # update its canonical to the fuzzy target
    n_updated = 0
    for entry in crosswalk:
        if entry["canonical_name"] in fuzzy_merges:
            entry["canonical_name"] = fuzzy_merges[entry["canonical_name"]]
            entry["match_type"] = "fuzzy"
            n_updated += 1

    print(f"  Updated {n_updated} crosswalk entries via fuzzy merges")
else:
    print(f"\nNo {FUZZY_REVIEWED} found — skipping fuzzy merges")
    print(f"To add fuzzy matches:")
    print(f"  1. Review Data/entity_resolution/fuzzy_candidates.csv")
    print(f"  2. Create {FUZZY_REVIEWED} with columns:")
    print(f"     name_a, name_b, accept, canonical_name")
    print(f"  3. Re-run this script")

# --------------------------------------------------------------------------
# 4. Apply manual overrides if file exists
# --------------------------------------------------------------------------

MANUAL_OVERRIDES = "Data/entity_resolution/manual_overrides.csv"

if os.path.exists(MANUAL_OVERRIDES):
    print(f"\nFound {MANUAL_OVERRIDES} — applying manual merges")

    # Expected columns: raw_name, canonical_name
    manual_map = {}
    with open(MANUAL_OVERRIDES, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            manual_map[row["raw_name"].strip()] = row["canonical_name"].strip()

    n_manual = 0
    for entry in crosswalk:
        if entry["raw_name"] in manual_map:
            entry["canonical_name"] = manual_map[entry["raw_name"]]
            entry["match_type"] = "manual"
            n_manual += 1

    # Add any manual mappings for raw names not already in crosswalk
    existing_raws = {e["raw_name"] for e in crosswalk}
    for raw, canonical in manual_map.items():
        if raw not in existing_raws:
            crosswalk.append({
                "raw_name": raw,
                "canonical_name": canonical,
                "match_type": "manual",
            })
            n_manual += 1

    print(f"  Applied {n_manual} manual overrides")
else:
    print(f"\nNo {MANUAL_OVERRIDES} found — skipping manual overrides")

# --------------------------------------------------------------------------
# 5. Write final crosswalk
# --------------------------------------------------------------------------

OUTPUT = "Data/entity_resolution/name_crosswalk.csv"

crosswalk.sort(key=lambda x: x["raw_name"])

with open(OUTPUT, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["raw_name", "canonical_name", "match_type"])
    writer.writeheader()
    for row in crosswalk:
        writer.writerow(row)

# --------------------------------------------------------------------------
# 6. Summary
# --------------------------------------------------------------------------

n_total = len(crosswalk)
n_canonical = len(set(e["canonical_name"] for e in crosswalk))
n_changed = sum(1 for e in crosswalk if e["raw_name"] != e["canonical_name"])

match_types = defaultdict(int)
for e in crosswalk:
    match_types[e["match_type"]] += 1

print(f"\n=== Crosswalk Summary ===")
print(f"Total raw names:       {n_total}")
print(f"Unique canonical:      {n_canonical}")
print(f"Names that changed:    {n_changed}")
print(f"Match type breakdown:")
for mt, count in sorted(match_types.items()):
    print(f"  {mt}: {count}")

print(f"\nOutput: {OUTPUT}")
print(f"\n=== Step 4 complete ===")
