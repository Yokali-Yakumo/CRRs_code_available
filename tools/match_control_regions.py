#!/usr/bin/env python3
"""
Generate genome-wide control regions matched to each input BED interval by
length and GC content.

Features
--------
- Generates N-fold controls for every input interval.
- Matches each control to its source interval, rather than only matching the
  global length/GC distributions.
- Supports symmetric or asymmetric GC and length tolerances.
- Excludes the input BED itself and one or more blacklist BED files.
- Rejects sequences with excessive non-ACGT bases.
- Avoids overlap among generated controls by default.
- Parallel candidate search using multiple worker processes.

Dependency
----------
pysam
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import math
import multiprocessing as mp
import os
import random
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pysam


@dataclass(frozen=True)
class BedRecord:
    chrom: str
    start: int
    end: int
    name: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Candidate:
    target_index: int
    replicate_index: int
    chrom: str
    start: int
    end: int
    gc: float


# Worker globals. Each process opens its own FASTA handle.
_W_FASTA: Optional[pysam.FastaFile] = None
_W_CHROMS: Sequence[str] = ()
_W_CHROM_LENGTHS: Dict[str, int] = {}
_W_EXCLUSION: Dict[str, Tuple[List[int], List[int]]] = {}
_W_GC_LOWER = 0.02
_W_GC_UPPER = 0.02
_W_LEN_LOWER = 0.10
_W_LEN_UPPER = 0.10
_W_MAX_ATTEMPTS = 10000
_W_MAX_N_FRAC = 0.0
_W_SAME_CHROM = False
_W_BASE_SEED = 1


def open_text(path: str):
    """Open plain text or gzip-compressed text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def read_bed(path: str, require_name: bool = False) -> List[BedRecord]:
    records: List[BedRecord] = []
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: BED line has fewer than 3 tab-delimited columns")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: BED start/end must be integers") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid interval {fields[0]}:{start}-{end}")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"region_{len(records) + 1}"
            if require_name and len(fields) < 4:
                raise ValueError(f"{path}:{line_no}: BED4 name column is required")
            records.append(BedRecord(fields[0], start, end, name))
    if not records:
        raise ValueError(f"No valid BED intervals found in {path}")
    return records


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return []
    merged: List[List[int]] = [[sorted_intervals[0][0], sorted_intervals[0][1]]]
    for start, end in sorted_intervals[1:]:
        last = merged[-1]
        if start <= last[1]:
            if end > last[1]:
                last[1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def build_static_index(records: Iterable[BedRecord]) -> Dict[str, Tuple[List[int], List[int]]]:
    grouped: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for rec in records:
        grouped[rec.chrom].append((rec.start, rec.end))

    result: Dict[str, Tuple[List[int], List[int]]] = {}
    for chrom, intervals in grouped.items():
        merged = merge_intervals(intervals)
        result[chrom] = ([x[0] for x in merged], [x[1] for x in merged])
    return result


def static_overlaps(index: Dict[str, Tuple[List[int], List[int]]], chrom: str, start: int, end: int) -> bool:
    data = index.get(chrom)
    if data is None:
        return False
    starts, ends = data
    i = bisect.bisect_left(starts, end) - 1
    return i >= 0 and ends[i] > start


class DynamicIntervalIndex:
    """Simple per-chromosome sorted interval index for accepted controls."""

    def __init__(self) -> None:
        self._starts: Dict[str, List[int]] = defaultdict(list)
        self._ends: Dict[str, List[int]] = defaultdict(list)

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts = self._starts.get(chrom)
        if not starts:
            return False
        ends = self._ends[chrom]
        i = bisect.bisect_left(starts, end) - 1
        return i >= 0 and ends[i] > start

    def add(self, chrom: str, start: int, end: int) -> None:
        starts = self._starts[chrom]
        ends = self._ends[chrom]
        pos = bisect.bisect_left(starts, start)
        starts.insert(pos, start)
        ends.insert(pos, end)


def gc_content(sequence: str) -> Tuple[Optional[float], float]:
    seq = sequence.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    valid = a + c + g + t
    if valid == 0:
        return None, 1.0
    non_acgt_fraction = 1.0 - valid / len(seq)
    return (g + c) / valid, non_acgt_fraction


def choose_weighted_chromosome(rng: random.Random, chroms: Sequence[str], lengths: Dict[str, int], region_len: int) -> Optional[str]:
    eligible: List[Tuple[str, int]] = []
    total = 0
    for chrom in chroms:
        n_positions = lengths[chrom] - region_len + 1
        if n_positions > 0:
            eligible.append((chrom, n_positions))
            total += n_positions
    if total <= 0:
        return None
    pick = rng.randrange(total)
    cumulative = 0
    for chrom, weight in eligible:
        cumulative += weight
        if pick < cumulative:
            return chrom
    return eligible[-1][0]


def init_worker(
    fasta_path: str,
    chroms: Sequence[str],
    chrom_lengths: Dict[str, int],
    exclusion_index: Dict[str, Tuple[List[int], List[int]]],
    gc_lower: float,
    gc_upper: float,
    len_lower: float,
    len_upper: float,
    max_attempts: int,
    max_n_frac: float,
    same_chrom: bool,
    base_seed: int,
) -> None:
    global _W_FASTA, _W_CHROMS, _W_CHROM_LENGTHS, _W_EXCLUSION
    global _W_GC_LOWER, _W_GC_UPPER, _W_LEN_LOWER, _W_LEN_UPPER
    global _W_MAX_ATTEMPTS, _W_MAX_N_FRAC, _W_SAME_CHROM, _W_BASE_SEED

    _W_FASTA = pysam.FastaFile(fasta_path)
    _W_CHROMS = tuple(chroms)
    _W_CHROM_LENGTHS = chrom_lengths
    _W_EXCLUSION = exclusion_index
    _W_GC_LOWER = gc_lower
    _W_GC_UPPER = gc_upper
    _W_LEN_LOWER = len_lower
    _W_LEN_UPPER = len_upper
    _W_MAX_ATTEMPTS = max_attempts
    _W_MAX_N_FRAC = max_n_frac
    _W_SAME_CHROM = same_chrom
    _W_BASE_SEED = base_seed


def find_one_candidate(task: Tuple[int, int, int, str, int, float]) -> Tuple[int, int, Optional[Candidate]]:
    """Find one candidate for one target/replicate/round."""
    if _W_FASTA is None:
        raise RuntimeError("Worker FASTA was not initialized")

    target_index, replicate_index, round_index, target_chrom, target_len, target_gc = task
    seed = (
        _W_BASE_SEED
        + 1_000_003 * target_index
        + 97_409 * replicate_index
        + 7_919 * round_index
    ) & 0xFFFFFFFFFFFF
    rng = random.Random(seed)

    min_len = max(1, math.floor(target_len * (1.0 - _W_LEN_LOWER)))
    max_len = max(min_len, math.ceil(target_len * (1.0 + _W_LEN_UPPER)))
    min_gc = max(0.0, target_gc - _W_GC_LOWER)
    max_gc = min(1.0, target_gc + _W_GC_UPPER)

    for _ in range(_W_MAX_ATTEMPTS):
        candidate_len = rng.randint(min_len, max_len)

        if _W_SAME_CHROM:
            chrom = target_chrom
            if chrom not in _W_CHROM_LENGTHS or _W_CHROM_LENGTHS[chrom] < candidate_len:
                continue
        else:
            chrom = choose_weighted_chromosome(rng, _W_CHROMS, _W_CHROM_LENGTHS, candidate_len)
            if chrom is None:
                continue

        chrom_len = _W_CHROM_LENGTHS[chrom]
        start = rng.randint(0, chrom_len - candidate_len)
        end = start + candidate_len

        if static_overlaps(_W_EXCLUSION, chrom, start, end):
            continue

        seq = _W_FASTA.fetch(chrom, start, end)
        candidate_gc, non_acgt_fraction = gc_content(seq)
        if candidate_gc is None or non_acgt_fraction > _W_MAX_N_FRAC:
            continue
        if min_gc <= candidate_gc <= max_gc:
            return target_index, replicate_index, Candidate(
                target_index=target_index,
                replicate_index=replicate_index,
                chrom=chrom,
                start=start,
                end=end,
                gc=candidate_gc,
            )

    return target_index, replicate_index, None


def validate_args(args: argparse.Namespace) -> None:
    for label, value in (
        ("--gc-tol", args.gc_tol),
        ("--gc-lower", args.gc_lower),
        ("--gc-upper", args.gc_upper),
        ("--length-tol", args.length_tol),
        ("--length-lower", args.length_lower),
        ("--length-upper", args.length_upper),
        ("--max-n-frac", args.max_n_frac),
    ):
        if value is not None and value < 0:
            raise ValueError(f"{label} must be >= 0")
    if args.fold < 1:
        raise ValueError("--fold must be >= 1")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.max_attempts < 1 or args.max_rounds < 1:
        raise ValueError("--max-attempts and --max-rounds must be >= 1")
    if args.max_n_frac > 1:
        raise ValueError("--max-n-frac must be <= 1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate genome-wide control BED regions matched by length and GC content.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-f", "--fasta", required=True, help="Reference genome FASTA; .fai will be created if absent")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument(
        "-b", "--blacklist", action="append", default=[],
        help="Blacklist BED file to exclude; may be supplied multiple times",
    )
    parser.add_argument("--fold", type=int, default=2, help="Number of controls per input interval")

    parser.add_argument(
        "--gc-tol", type=float, default=0.02,
        help="Symmetric absolute GC tolerance; 0.02 means +/- 2 percentage points",
    )
    parser.add_argument("--gc-lower", type=float, default=None, help="Allowed GC decrease; overrides lower side of --gc-tol")
    parser.add_argument("--gc-upper", type=float, default=None, help="Allowed GC increase; overrides upper side of --gc-tol")

    parser.add_argument(
        "--length-tol", type=float, default=0.10,
        help="Symmetric relative length tolerance; 0.10 means +/- 10%%",
    )
    parser.add_argument("--length-lower", type=float, default=None, help="Allowed relative length decrease; overrides lower side of --length-tol")
    parser.add_argument("--length-upper", type=float, default=None, help="Allowed relative length increase; overrides upper side of --length-tol")

    parser.add_argument("--threads", type=int, default=12, help="Number of parallel worker processes")
    parser.add_argument("--max-attempts", type=int, default=10000, help="Random candidates tested per task in each round")
    parser.add_argument("--max-rounds", type=int, default=20, help="Refill rounds for failed or conflicting controls")
    parser.add_argument("--max-n-frac", type=float, default=0.0, help="Maximum fraction of non-ACGT bases in a control")
    parser.add_argument("--seed", type=int, default=1, help="Random seed")
    parser.add_argument("--same-chrom", action="store_true", help="Require each control to come from the same chromosome as its target")
    parser.add_argument(
        "--chrom-regex", default=None,
        help="Only sample chromosomes whose names match this regular expression",
    )
    parser.add_argument(
        "--allow-control-overlap", action="store_true",
        help="Allow generated controls to overlap one another",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    validate_args(args)

    gc_lower = args.gc_tol if args.gc_lower is None else args.gc_lower
    gc_upper = args.gc_tol if args.gc_upper is None else args.gc_upper
    len_lower = args.length_tol if args.length_lower is None else args.length_lower
    len_upper = args.length_tol if args.length_upper is None else args.length_upper

    if not os.path.exists(args.fasta):
        raise FileNotFoundError(args.fasta)
    if not os.path.exists(args.fasta + ".fai"):
        print(f"[INFO] FASTA index not found; creating {args.fasta}.fai", file=sys.stderr)
        pysam.faidx(args.fasta)

    targets = read_bed(args.input)
    blacklist_records: List[BedRecord] = []
    for blacklist_path in args.blacklist:
        blacklist_records.extend(read_bed(blacklist_path))

    fasta = pysam.FastaFile(args.fasta)
    chrom_lengths = {chrom: length for chrom, length in zip(fasta.references, fasta.lengths)}

    chroms = list(fasta.references)
    if args.chrom_regex:
        pattern = re.compile(args.chrom_regex)
        chroms = [chrom for chrom in chroms if pattern.fullmatch(chrom)]
        if not chroms:
            raise ValueError("--chrom-regex excluded every FASTA sequence")

    target_gc: List[float] = []
    for rec in targets:
        if rec.chrom not in chrom_lengths:
            raise ValueError(f"Input chromosome {rec.chrom!r} is absent from FASTA")
        if rec.end > chrom_lengths[rec.chrom]:
            raise ValueError(
                f"Input interval {rec.chrom}:{rec.start}-{rec.end} exceeds chromosome length {chrom_lengths[rec.chrom]}"
            )
        seq = fasta.fetch(rec.chrom, rec.start, rec.end)
        gc, _ = gc_content(seq)
        if gc is None:
            raise ValueError(f"Input interval {rec.chrom}:{rec.start}-{rec.end} contains no A/C/G/T bases")
        target_gc.append(gc)
    fasta.close()

    if args.same_chrom:
        missing = sorted({rec.chrom for rec in targets if rec.chrom not in chroms})
        if missing:
            raise ValueError(
                "--same-chrom was requested, but these target chromosomes are excluded by --chrom-regex: "
                + ", ".join(missing)
            )

    # Exclude both original target regions and explicit blacklist regions.
    exclusion_index = build_static_index([*targets, *blacklist_records])

    pending = {(target_index, replicate_index) for target_index in range(len(targets)) for replicate_index in range(1, args.fold + 1)}
    accepted: Dict[Tuple[int, int], Candidate] = {}
    accepted_index = DynamicIntervalIndex()

    print(
        f"[INFO] targets={len(targets)}, requested_controls={len(pending)}, threads={args.threads}",
        file=sys.stderr,
    )
    print(
        f"[INFO] GC tolerance=-{gc_lower:.4f}/+{gc_upper:.4f}; "
        f"length tolerance=-{len_lower:.4f}/+{len_upper:.4f}",
        file=sys.stderr,
    )

    initargs = (
        args.fasta,
        tuple(chroms),
        chrom_lengths,
        exclusion_index,
        gc_lower,
        gc_upper,
        len_lower,
        len_upper,
        args.max_attempts,
        args.max_n_frac,
        args.same_chrom,
        args.seed,
    )

    context = mp.get_context("spawn")
    with context.Pool(processes=args.threads, initializer=init_worker, initargs=initargs) as pool:
        for round_index in range(args.max_rounds):
            if not pending:
                break

            tasks = [
                (
                    target_index,
                    replicate_index,
                    round_index,
                    targets[target_index].chrom,
                    targets[target_index].length,
                    target_gc[target_index],
                )
                for target_index, replicate_index in sorted(pending)
            ]
            chunksize = max(1, len(tasks) // max(1, args.threads * 4))
            next_pending = set()

            for target_index, replicate_index, candidate in pool.imap_unordered(
                find_one_candidate, tasks, chunksize=chunksize
            ):
                key = (target_index, replicate_index)
                if candidate is None:
                    next_pending.add(key)
                    continue
                if (
                    not args.allow_control_overlap
                    and accepted_index.overlaps(candidate.chrom, candidate.start, candidate.end)
                ):
                    next_pending.add(key)
                    continue

                accepted[key] = candidate
                if not args.allow_control_overlap:
                    accepted_index.add(candidate.chrom, candidate.start, candidate.end)

            pending = next_pending
            print(
                f"[INFO] round={round_index + 1}, accepted={len(accepted)}, remaining={len(pending)}",
                file=sys.stderr,
            )

    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)

    with open(args.output, "wt", encoding="utf-8") as out:
        out.write(
            "#chrom\tstart\tend\tcontrol_id\tsource_id\tsource_chrom\tsource_start\tsource_end"
            "\tsource_length\tcontrol_length\tsource_gc\tcontrol_gc\tgc_difference\n"
        )
        for key in sorted(accepted):
            target_index, replicate_index = key
            candidate = accepted[key]
            target = targets[target_index]
            control_id = f"{target.name}__control_{replicate_index}"
            out.write(
                f"{candidate.chrom}\t{candidate.start}\t{candidate.end}\t{control_id}\t{target.name}\t"
                f"{target.chrom}\t{target.start}\t{target.end}\t{target.length}\t{candidate.end - candidate.start}\t"
                f"{target_gc[target_index]:.6f}\t{candidate.gc:.6f}\t{candidate.gc - target_gc[target_index]:.6f}\n"
            )

    unmatched_path = args.output + ".unmatched.tsv"
    if pending:
        with open(unmatched_path, "wt", encoding="utf-8") as out:
            out.write("source_id\tsource_chrom\tsource_start\tsource_end\treplicate_index\n")
            for target_index, replicate_index in sorted(pending):
                target = targets[target_index]
                out.write(
                    f"{target.name}\t{target.chrom}\t{target.start}\t{target.end}\t{replicate_index}\n"
                )
        print(
            f"[WARNING] {len(pending)} controls were not found. See {unmatched_path}. "
            "Consider increasing tolerances, --max-attempts, or --max-rounds.",
            file=sys.stderr,
        )
    else:
        if os.path.exists(unmatched_path):
            os.remove(unmatched_path)

    print(f"[INFO] Wrote {len(accepted)} controls to {args.output}", file=sys.stderr)
    return 0 if not pending else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(1)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
