#!/usr/bin/env python
"""Count the wide-distance (1250 pc) sample on the same probe partitions.

If the wide sample is only a few times the 500 pc sample, it is cheaper to pull
the wide superset once and subset client-side than to run two full downloads.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import logging
import numpy as np
from pipeline import config as cfg
from pipeline.adql import count_query
from pipeline.tap import run_adql

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

PROBE = [10, 95, 180]
for dmax in (cfg.CUTS.distance_max_pc, cfg.CUTS.distance_max_pc_wide):
    ns = []
    for p in PROBE:
        n = int(run_adql(count_query(p, distance_max_pc=dmax),
                         name=f"count_d{int(dmax)}_p{p}")["n"].iloc[0])
        ns.append(n)
    print(f"d < {dmax:6.0f} pc : per-partition {ns}  mean {np.mean(ns):8.0f} "
          f"-> total ~{np.mean(ns)*cfg.N_PARTITIONS:,.0f}")
