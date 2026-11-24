"""Resolve release inputs and experiment outputs independently of cwd."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
OPD_CACHE = REPO_ROOT / "cache/opd/decomp_cache.json"
RUNTIME_CACHE = REPO_ROOT / "cache/runtime/time_cache.json"
OPS_COEFFICIENTS = REPO_ROOT / "cache/ops"
OPS_RESULTS = REPO_ROOT / "test/OPS"
OPS_TEMP = OPS_RESULTS / ".tmp"


def add_ere_module_path():
    """Prefer the extension built by this release; also allow PYTHONPATH."""
    build = REPO_ROOT / "build"
    if build.is_dir() and str(build) not in sys.path:
        sys.path.insert(0, str(build))
