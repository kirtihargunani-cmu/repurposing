# One-time seed of the 79o answer cache from the stored top-10 / top-3 literature tables (user 2026-11-08: no re-pull).
import json, re, urllib.parse
import pandas as pd
R = "results/lincs/"; LC = "data/interim/lincs/"
pm, ct = {}, {}
tf = LC + "literature_fda_all_titles.json"; titles = json.load(open(tf))
for f in ["literature_v2_consensus", "literature_fda_top_v2"]:
    l = pd.read_csv(R + f + ".tsv", sep="\t")
    for r in l.itertuples():
        if isinstance(r.pubmed_query, str) and r.pubmed_count == r.pubmed_count:
            ids = []
            for j in (1, 2, 3):
                u = getattr(r, f"paper{j}")
                if isinstance(u, str):
                    pid = u.rstrip("/").split("/")[-1]; ids.append(pid)
                    tt = getattr(r, f"paper{j}_title")
                    if isinstance(tt, str) and pid not in titles:
                        titles[pid] = tt
            pm[r.pubmed_query] = (int(r.pubmed_count), ids)
        if isinstance(r.trials_search_url, str) and r.trial_count == r.trial_count:
            qs = urllib.parse.parse_qs(urllib.parse.urlparse(r.trials_search_url).query)
            st = []
            for j in (1, 2, 3):
                u, info = getattr(r, f"trial{j}"), getattr(r, f"trial{j}_info")
                if isinstance(u, str) and isinstance(info, str):
                    head, title = info.split(": ", 1); phase, status = head.split(", ", 1)
                    st.append((u.rstrip("/").split("/")[-1], phase, status, title))
            ct[(qs["intr"][0], qs["cond"][0])] = (int(r.trial_count), st)
with open(LC + "literature_query_pubmed.jsonl", "a") as fh:
    for q, (n, ids) in pm.items():
        fh.write(json.dumps(dict(q=q, n=n, ids=ids)) + "\n")
with open(LC + "literature_query_trials.jsonl", "a") as fh:
    for (d, c), (n, st) in ct.items():
        fh.write(json.dumps(dict(drug=d, cond=c, n=n, st=st)) + "\n")
json.dump(titles, open(tf, "w"))
print("seeded pubmed", len(pm), "trials", len(ct), "titles", len(titles))
