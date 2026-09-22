#!/usr/bin/env python3

import argparse
import csv
import re

import pandas as pd


RANKS = {
    "d": "domain",
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
TOKEN = re.compile(r'(?:^|,)([dkpcofgs]):([^,(]+?)(?:\(([0-9.]+)\))?(?=,|$)')


def parse_taxonomy(value: str) -> dict[str, object]:
    parsed: dict[str, object] = {}
    for rank_code, label, confidence in TOKEN.findall(value or ""):
        rank = RANKS[rank_code]
        parsed[rank] = label.strip().strip('"')
        parsed[f"{rank}_confidence"] = float(confidence) if confidence else None
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert USEARCH SINTAX output to a rank table")
    parser.add_argument("--sintax", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    with open(args.sintax, encoding="utf-8", newline="") as stream:
        for fields in csv.reader(stream, delimiter="\t"):
            if not fields:
                continue
            record: dict[str, object] = {"ASV": fields[0]}
            record.update(parse_taxonomy(fields[1] if len(fields) > 1 else ""))
            rows.append(record)

    ordered = ["ASV"]
    for rank in RANKS.values():
        ordered.extend([rank, f"{rank}_confidence"])
    pd.DataFrame(rows).reindex(columns=ordered).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
