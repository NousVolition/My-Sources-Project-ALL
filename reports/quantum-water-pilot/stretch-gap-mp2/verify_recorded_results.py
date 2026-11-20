"""Check the published records and arithmetic; do not run MP2 or Lean."""
import argparse
import ast
import hashlib
import json
import math
import re
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
read=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b,tol=1e-10):
    assert abs(a-b)<=tol,(a,b,tol)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protected-root",type=Path,default=ROOT.parents[2])
    args=ap.parse_args()
    manifest=read(DATA/"addition-checksums.json")["files"]
    for name,digest in manifest.items():
        assert sha(DATA/name)==digest,name
    protocol=read(DATA/"protocol.json")
    prereg=read(DATA/"preregistration-receipt.json")
    run=read(DATA/"execution-receipt.json")
    start=read(DATA/"run-start.json")
    result=read(DATA/"comparison-results.json")
    ph=sha(DATA/"protocol.json")
    assert ph==prereg["protocol_sha256"]==run["protocol_sha256"]==result["protocol_sha256"]
    assert datetime.fromisoformat(protocol["registered_before_any_MP2_energy_or_optimization_utc"])<datetime.fromisoformat(start["started_at_utc"])
    assert datetime.fromisoformat(prereg["checksum_corrected_before_calculations_utc"])<datetime.fromisoformat(start["started_at_utc"])
    assert sha(DATA/"run_comparison.py")==run["runner_sha256"]==start["runner_sha256"]
    assert sha(DATA/"analyze_comparison.py")==result["analyzer_sha256"]
    assert sha(DATA/"implementation-receipt.json")==start["implementation_receipt_sha256"]
    assert sha(DATA/"preregistration-receipt.json")==start["preregistration_receipt_sha256"]
    assert run["status"]=="completed"
    records={}
    for row in run["evaluations"]:
        assert sha(DATA/row["file"])==row["sha256"]
        r=read(DATA/row["file"])
        assert r["protocol_sha256"]==ph and r["RHF_converged"]
        records[row["id"]]=r
    assert len(records)==run["energy_gradient_evaluation_count"]==91
    for h in read(DATA/"hessians.json"):
        step=h["step_bohr"]
        H=h["raw_hessian_Eh_per_bohr2"]
        for refs in h["evaluation_references"]:
            plus,minus=records[refs["plus"]],records[refs["minus"]]
            pg=sum(plus["gradient_Eh_per_bohr"],[])
            mg=sum(minus["gradient_Eh_per_bohr"],[])
            for row in range(9):
                close((pg[row]-mg[row])/(2*step),H[row][refs["axis"]],1e-12)
        for i in range(9):
            for j in range(9):
                close((H[i][j]+H[j][i])/2,h["symmetrized_hessian_Eh_per_bohr2"][i][j],1e-12)
    primary=result["primary_results"]["H2O"]
    gap=primary["antisymmetric_stretch_cm1"]-primary["symmetric_stretch_cm1"]
    close(gap,primary["gap_cm1"])
    test=protocol["primary_gap_test"]
    assert test["reference_gap_cm1"]==111 and test["absolute_gap_error_limit_cm1"]==10
    assert test["pass_interval_cm1_inclusive"]==[101,121]
    assert len(result["numerical_checks"])==36
    for c in result["numerical_checks"]:
        assert c["passes"]==(c["value"]<=c["limit"]) and c["passes"]
    assert result["all_numerical_checks_pass"] and result["H2O_numerical_validity_passes"]
    error=abs(gap-111)
    close(error,result["H2O_absolute_gap_error_cm1"])
    assert error>10 and result["primary_H2O_outcome"]=="FAIL"
    follow=read(DATA/"response-check/verification-results.json")
    fp=read(DATA/"response-check/protocol.json")
    assert sha(DATA/"comparison-results.json")==fp["parent_results_sha256"]
    assert ph==fp["parent_protocol_sha256"]
    assert sha(DATA/"response-check/protocol.json")==follow["protocol_sha256"]
    assert sha(DATA/"response-check/run_response_check.py")==follow["runner_sha256"]
    assert datetime.fromisoformat(fp["written_before_followup_calculations_utc"])<datetime.fromisoformat(follow["started_at_utc"])
    follow_records={}
    for row in follow["evaluations"]:
        assert sha(DATA/"response-check"/row["file"])==row["sha256"]
        r=read(DATA/"response-check"/row["file"])
        assert r["RHF_converged"] and r["response"]["gmres_info"]==0
        assert r["response"]["true_preconditioned_residual_norm"]<=fp["validation_limits"]["maximum_true_preconditioned_response_residual"]
        follow_records[r["tag"]]=r
    assert len(follow_records)==follow["evaluation_count"]==19
    fh=read(DATA/"response-check/hessian.json")
    for axis in range(9):
        pg=sum(follow_records[f"axis-{axis}-plus"]["gradient_Eh_per_bohr"],[])
        mg=sum(follow_records[f"axis-{axis}-minus"]["gradient_Eh_per_bohr"],[])
        for row in range(9):
            close((pg[row]-mg[row])/(2*fh["step_bohr"]),fh["raw_hessian_Eh_per_bohr2"][row][axis],1e-12)
    assert len(follow["numerical_checks"])==6 and follow["all_validation_checks_pass"]
    for c in follow["numerical_checks"]:
        assert c["passes"]==(c["value"]<=c["limit"]) and c["passes"]
    close(abs(follow["results"]["H2O"]["gap_cm1"]-gap),follow["results"]["H2O"]["gap_change_cm1"])
    assert not (101<=follow["results"]["H2O"]["gap_cm1"]<=121)
    preserved=0
    for name,digest in protocol["preservation"]["protected_published_files_sha256"].items():
        if name=="README.md":
            continue
        assert sha(args.protected_root/name)==digest,name
        preserved+=1
    assert preserved==15
    path_pattern=re.compile(rb"(?i)\b[A-Z]:[\\/]|file:[/]{2,3}(?=[A-Za-z/])|/(?:Users|home|mnt|tmp|private/var)/")
    for f in ROOT.rglob("*"):
        if not f.is_file() or "__pycache__" in f.parts:
            continue
        assert not path_pattern.search(f.read_bytes()),f.relative_to(ROOT).as_posix()
        if f.suffix==".py":
            ast.parse(f.read_text(encoding="utf-8-sig"))
        if f.suffix==".json":
            read(f)
        if f.suffix==".md":
            for link in re.findall(r"\]\(([^)]+)\)",f.read_text(encoding="utf-8-sig")):
                if "://" not in link and not link.startswith("#"):
                    assert (f.parent/link.split("#")[0]).exists(),(f.name,link)
    print(json.dumps({"recorded_package_verification":"PASS","MP2_H2O_gap_test":"FAIL",
                      "H2O_gap_cm1":gap,"pass_interval_cm1":[101,121],
                      "primary_numerical_checks_passed":36,"response_followup_checks_passed":6,
                      "archived_files_verified":len(manifest),"gradient_records_verified":110,
                      "published_xTB_study_files_unchanged":preserved,
                      "new_MP2_run":False,"Lean_run":False,"local_folder_path_matches":0}))
if __name__=="__main__":
    main()
