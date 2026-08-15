"""Step 5: blinding.

A secret constant is drawn once and stored in ``blind/secret_offset.json``,
which is gitignored.  Every residual passed through :func:`apply` is shifted by
it.  The analyst fixes all analysis choices while looking only at blinded
numbers, commits the frozen configuration, and only then unblinds.

Two properties make this a real blind rather than a ritual:

* The offset is drawn from a distribution wide compared with any plausible
  signal, so "the blinded answer looks like zero" carries no information.
* A SHA-256 commitment to the offset is written to ``blind/commitment.json``,
  which IS tracked by git.  The commitment is created at the same moment as the
  offset, so the offset cannot be quietly redrawn after seeing a result -- the
  committed hash would no longer match.

:func:`apply` deliberately never returns or logs the offset value.
"""

from __future__ import annotations

import hashlib
import json
import secrets
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np

from . import config as cfg

SECRET_PATH = cfg.BLIND_DIR / "secret_offset.json"
COMMIT_PATH = cfg.BLIND_DIR / "commitment.json"

# Width of the blinding distribution, magnitudes.  Chosen to be large compared
# with the systematic floor we expect (~1e-3 mag) and with any signal we could
# plausibly claim, so a blinded result near zero is uninformative.
BLIND_SCALE_MAG = 0.05


@dataclass
class BlindState:
    exists: bool
    committed: bool
    created_utc: str | None
    commitment_sha256: str | None
    unblinded: bool


def _digest(value: float, salt: str) -> str:
    return hashlib.sha256(f"{salt}:{value!r}".encode()).hexdigest()


def create(force: bool = False) -> BlindState:
    """Draw the secret offset and write the public commitment. Idempotent."""
    if SECRET_PATH.exists() and not force:
        return status()

    salt = secrets.token_hex(16)
    # Uniform, not Gaussian: bounded support means the unblinded correction can
    # never be so large that it looks like a mistake.
    offset = float(np.random.default_rng(
        secrets.randbits(63)).uniform(-BLIND_SCALE_MAG, BLIND_SCALE_MAG))

    created = datetime.now(timezone.utc).isoformat()
    SECRET_PATH.write_text(json.dumps(
        {"offset_mag": offset, "salt": salt, "created_utc": created}, indent=2))
    COMMIT_PATH.write_text(json.dumps({
        "commitment_sha256": _digest(offset, salt),
        "salt": salt,
        "scale_mag": BLIND_SCALE_MAG,
        "distribution": "uniform(-scale, +scale)",
        "created_utc": created,
        "unblinded": False,
        "note": ("SHA-256 of 'salt:repr(offset)'. Verify after unblinding with "
                 "pipeline.blind.verify()."),
    }, indent=2))
    return status()


def _load_secret() -> dict:
    if not SECRET_PATH.exists():
        raise FileNotFoundError(
            f"{SECRET_PATH} not found -- run pipeline.blind.create() first")
    return json.loads(SECRET_PATH.read_text())


def apply(resid: np.ndarray) -> np.ndarray:
    """Return blinded residuals. Never returns or logs the offset itself."""
    return np.asarray(resid, dtype=float) + _load_secret()["offset_mag"]


def is_unblinded() -> bool:
    if not COMMIT_PATH.exists():
        return False
    return bool(json.loads(COMMIT_PATH.read_text()).get("unblinded", False))


def status() -> BlindState:
    commit = json.loads(COMMIT_PATH.read_text()) if COMMIT_PATH.exists() else {}
    return BlindState(
        exists=SECRET_PATH.exists(),
        committed=COMMIT_PATH.exists(),
        created_utc=commit.get("created_utc"),
        commitment_sha256=commit.get("commitment_sha256"),
        unblinded=bool(commit.get("unblinded", False)),
    )


def unblind(confirm: str) -> float:
    """Reveal the offset. Requires the literal confirmation string.

    Records the unblinding time in the tracked commitment file, so the moment
    of unblinding is part of the git history.
    """
    if confirm != "I have frozen the analysis and committed it":
        raise ValueError("refusing to unblind without the confirmation string")
    secret = _load_secret()
    commit = json.loads(COMMIT_PATH.read_text())
    if _digest(secret["offset_mag"], secret["salt"]) != commit["commitment_sha256"]:
        raise RuntimeError(
            "COMMITMENT MISMATCH: the secret offset does not hash to the "
            "committed value. The blind is broken and the result is void.")
    commit["unblinded"] = True
    commit["unblinded_utc"] = datetime.now(timezone.utc).isoformat()
    commit["revealed_offset_mag"] = secret["offset_mag"]
    COMMIT_PATH.write_text(json.dumps(commit, indent=2))
    return float(secret["offset_mag"])


def verify() -> bool:
    secret = _load_secret()
    commit = json.loads(COMMIT_PATH.read_text())
    return _digest(secret["offset_mag"], secret["salt"]) == commit["commitment_sha256"]
