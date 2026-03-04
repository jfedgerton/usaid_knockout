###############################################################################
# Step 4b: Fix problematic auto-accepts in reviewed_fuzzy_matches.csv
#
# Issues:
# 1. "community of canarias spain" vs "community of cantabria spain" — different
#    Spanish autonomous communities, incorrectly auto-accepted at score 94
# 2. prefix_with_country_tag matches — too risky to auto-accept because
#    "UNDP" and "UNDP Mali" could be different nodes in the same network.
#    Move all to REVIEW for manual inspection.
###############################################################################

import csv

INPUT = "Data/entity_resolution/reviewed_fuzzy_matches.csv"
OUTPUT = INPUT  # overwrite

rows = []
with open(INPUT, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

n_fixed = 0

for row in rows:
    # Fix 1: canarias vs cantabria
    if ("canarias" in row["name_a"] and "cantabria" in row["name_b"]) or \
       ("cantabria" in row["name_a"] and "canarias" in row["name_b"]):
        if row["accept"] == "yes":
            row["accept"] = "no"
            row["review_note"] = "fixed_different_regions"
            n_fixed += 1

    # Fix 2: Move prefix_with_country_tag from yes to REVIEW
    if row["review_note"] == "prefix_with_country_tag" and row["accept"] == "yes":
        row["accept"] = "REVIEW"
        row["review_note"] = "prefix_with_country_tag_needs_review"
        n_fixed += 1

with open(OUTPUT, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

# Stats
accepted = sum(1 for r in rows if r["accept"] == "yes")
rejected = sum(1 for r in rows if r["accept"] == "no")
review = sum(1 for r in rows if r["accept"] == "REVIEW")

print(f"Fixed {n_fixed} entries")
print(f"Final counts: yes={accepted}, no={rejected}, REVIEW={review}")
print(f"Auto-accepted matches:")
for r in rows:
    if r["accept"] == "yes":
        print(f"  [{r['token_sort_score']}] {r['name_a']}  ->  {r['name_b']}")
