#!/usr/bin/env python
"""Repair UTF-8 text that was round-tripped through a latin-1 decode.

Windows PowerShell's Get-Content/Set-Content pipeline reads UTF-8 bytes as the
system codepage and writes them back as UTF-8, turning every non-ASCII
character into mojibake ("Omega" -> "Î©").  The transformation is exactly
    utf8_bytes -> decode(latin-1) -> encode(utf-8)
so it inverts cleanly with
    text.encode('latin-1').decode('utf-8')

Never edit these files from PowerShell again; use WSL sed or a UTF-8-aware
editor.
"""
from __future__ import annotations

import sys
from pathlib import Path


MARKERS = ("Ã", "Î", "â€", "Â")

# Windows PowerShell uses the ANSI codepage (cp1252 on most systems), NOT
# latin-1. They differ in 0x80-0x9F, which is exactly where the bytes for
# em-dash, curly quotes and the Greek block land -- so latin-1 fails to invert.
CODECS = ("cp1252", "latin-1", "cp1250", "cp1254")


def repair(path: Path) -> bool:
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    text = raw.decode("utf-8", errors="strict")
    if not any(m in text for m in MARKERS):
        print(f"{path.name}: already clean")
        return False

    best, best_codec, best_score = None, None, 10 ** 9
    for codec in CODECS:
        try:
            cand = text.encode(codec, errors="strict").decode("utf-8", errors="strict")
        except (UnicodeEncodeError, UnicodeDecodeError):
            continue
        score = sum(cand.count(m) for m in MARKERS)
        if score < best_score:
            best, best_codec, best_score = cand, codec, score

    if best is None:
        print(f"{path.name}: NOT reversible by any codec -- left untouched")
        return False
    before = sum(text.count(m) for m in MARKERS)
    path.write_text(best, encoding="utf-8")
    print(f"{path.name}: repaired via {best_codec} "
          f"({before} markers -> {best_score})")
    return True


if __name__ == "__main__":
    targets = sys.argv[1:] or ["PAPER.md"]
    for t in targets:
        repair(Path(t))
