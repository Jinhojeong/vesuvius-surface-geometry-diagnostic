"""Score flummoxjr's harness against exact truth, exactly as he asked:
per d12 bin, (a) fraction of step-3-REJECTED pairs whose two points genuinely
lie on different sheets (turn_id), (b) false-accept rate on ACCEPTED sites
(pairs the guard kept whose points are actually the same sheet in truth)."""
import csv
from collections import defaultdict
BINS=[(2,3,"2-3"),(3,4,"3-4"),(4,5,"4-5"),(5,6.01,"5-6.01")]
def b(d):
    for lo,hi,name in BINS:
        if lo<=d<hi: return name
    return None

rej=defaultdict(lambda:[0,0]); acc=defaultdict(lambda:[0,0])
unresolved_r=unresolved_a=0
for row in csv.DictReader(open("/mnt/vesuvius/fp_out/facing_pairs_cc_rejected.csv")):
    if row["site_type"]!="tight": continue
    bb=b(float(row["d12"]))
    if bb is None: continue
    g=row["gt_diff_sheet"]
    if g=="" or g is None: unresolved_r+=1; continue
    rej[bb][0]+=1; rej[bb][1]+=int(g=="1")
for row in csv.DictReader(open("/mnt/vesuvius/fp_out/facing_pairs_sites.csv")):
    if row["site_type"]!="tight": continue
    bb=b(float(row["d12"]))
    if bb is None: continue
    g=row["gt_diff_sheet"]
    if g=="" or g is None: unresolved_a+=1; continue
    acc[bb][0]+=1; acc[bb][1]+=int(g=="0")
print("unresolved gt: rejected %d, accepted %d"%(unresolved_r,unresolved_a))
print("\n%-8s %26s %30s"%("d12 bin","rejected: n, truly-diff %","accepted: n, false-accept %"))
for lo,hi,name in BINS:
    rn,rd=rej[name]; an,af=acc[name]
    print("%-8s %14d  %8s %19d  %10s"%(name,rn,
        "%.1f%%"%(100*rd/rn) if rn else "-",an,"%.2f%%"%(100*af/an) if an else "-"))
tr=sum(v[0] for v in rej.values()); td=sum(v[1] for v in rej.values())
ta=sum(v[0] for v in acc.values()); tf=sum(v[1] for v in acc.values())
print("%-8s %14d  %8s %19d  %10s"%("all",tr,"%.1f%%"%(100*td/tr),ta,"%.2f%%"%(100*tf/ta) if ta else "-"))
print("\nper-case rejected truly-diff, to see spread:")
pc=defaultdict(lambda:[0,0])
for row in csv.DictReader(open("/mnt/vesuvius/fp_out/facing_pairs_cc_rejected.csv")):
    if row["site_type"]!="tight" or row["gt_diff_sheet"]=="": continue
    pc[row["case"]][0]+=1; pc[row["case"]][1]+=int(row["gt_diff_sheet"]=="1")
for c in sorted(pc): print("  %-18s %6d rejects, %.1f%% truly different"%(c,pc[c][0],100*pc[c][1]/pc[c][0]))
