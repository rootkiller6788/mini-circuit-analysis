# build_all.py - Generate all files for mini-transient-analysis
import os

BASE = r"F:\nano-everything\mini-electronic-info\1. mini-circuit-analysis\mini-transient-analysis"

def w(rel, content):
    path = os.path.join(BASE, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    n = content.count("\n")
    print(f"  {rel}: {n} lines")

