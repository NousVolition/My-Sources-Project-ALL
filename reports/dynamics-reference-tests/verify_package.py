"""Validate deliverable integrity and internal report links without network access."""
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parent


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets=[]
    def handle_starttag(self,tag,attrs):
        for name,value in attrs:
            if name in ("href","src") and value and not value.startswith(("https:","http:","#")):
                self.targets.append(value)


def main():
    manifest=json.loads((ROOT/"manifest.json").read_text())
    failures=[]
    for relative,expected in manifest.items():
        path=ROOT/relative
        if not path.is_file() or path.stat().st_size!=expected["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest()!=expected["sha256"]:
            failures.append(relative)
    internal=[];missing=[];outside=[]
    root_resolved=ROOT.resolve()
    for report in ROOT.rglob("*.html"):
        parser=Links()
        parser.feed(report.read_text(encoding="utf-8"))
        for link in parser.targets:
            path=(report.parent/unquote(urlsplit(link).path)).resolve()
            reference=f"{report.relative_to(ROOT).as_posix()}: {link}"
            if not path.is_relative_to(root_resolved):
                outside.append(reference)
                continue
            internal.append(reference)
            if not path.is_file(): missing.append(reference)
    result=dict(files_checked=len(manifest),checksum_failures=failures,internal_links_checked=len(internal),
                missing_internal_targets=missing,optional_links_outside_package=outside,passed=not failures and not missing,
                external_package_reference="../water-biology-study/report.html is optional and not bundled.")
    (ROOT/"package-verification.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    if not result["passed"]: sys.exit(1)


if __name__=="__main__": main()
