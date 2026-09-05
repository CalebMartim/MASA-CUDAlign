#!/usr/bin/env python3
"""Generate reproducible synthetic DNA as a single-record FASTA file.

Examples:
    python3 generate_dna.py 1000000 --seed 42 -o sequence_a.fasta
    python3 generate_dna.py 1000000 --seed 43 -o sequence_b.fasta
    python3 generate_dna.py 1000000 --seed 42 -o reproduced_a.fasta
    sha256sum sequence_a.fasta reproduced_a.fasta

Requires Python 3.8+; uses only the standard library. Existing output files
are never overwritten. Without -o, FASTA bytes go to stdout. A JSON summary
with the FASTA SHA-256 goes to stderr after successful generation.

Reproducibility contract, algorithm synthetic-dna-v1:
  * Seed is an unsigned 64-bit integer, serialized as 8 big-endian bytes.
  * For counters 0, 1, ..., compute SHA-256 of:
        b"synthetic-dna-v1\\0" + seed_bytes + counter_as_8_big_endian_bytes
  * Read digest bytes in order, each as four 2-bit groups, most significant
    first. Map 00 -> A, 01 -> C, 10 -> G, 11 -> T.
  * Concatenate blocks and take exactly SIZE bases.
  * Emit ASCII header ">synthetic-dna-v1 seed=SEED length=SIZE\\n",
    sequence lines of 80 bases, and a final LF. No timestamps or paths.

Same algorithm version, integer seed, and size guarantee identical FASTA
bytes on conforming Python implementations across machines and operating
systems. No dependency on random.seed(), Python hash(), locale, or native
newline conversion. Preserve this file and the printed checksum with your
experiment. A future change to the generation/format rules requires a new
algorithm version.

Different sizes with the same seed intentionally share a sequence prefix
(their headers and whole-file checksums differ). Different seeds create
separate pseudorandom sequences, not a controlled-similarity alignment pair.
Bases have approximately equal frequencies; exact 25% composition is not
enforced. This is a synthetic benchmark model, not a biological simulator.
Generation streams with bounded memory, so large sequences need not fit RAM.
"""

import argparse
import hashlib
import json
import sys


VERSION = "synthetic-dna-v1"
DOMAIN = b"synthetic-dna-v1\0"
ALPHABET = b"ACGT"
EXPANSION = tuple(
    bytes(ALPHABET[(value >> shift) & 3] for shift in (6, 4, 2, 0))
    for value in range(256)
)


def size_argument(value):
    try:
        number = int(value, 10)
    except ValueError:
        raise argparse.ArgumentTypeError("size must be a positive integer")
    if not 1 <= number <= 2**64:
        raise argparse.ArgumentTypeError("size must be between 1 and 2^64 bases")
    return number


def seed_argument(value):
    try:
        number = int(value, 10)
    except ValueError:
        raise argparse.ArgumentTypeError("seed must be an unsigned 64-bit integer")
    if not 0 <= number < 2**64:
        raise argparse.ArgumentTypeError("seed must be between 0 and 2^64 - 1")
    return number


def fasta_chunks(size, seed):
    """Yield canonical FASTA bytes without constructing the entire sequence."""
    yield (">{} seed={} length={}\n".format(VERSION, seed, size)).encode("ascii")
    prefix = DOMAIN + seed.to_bytes(8, "big")
    pending = b""
    remaining = size
    counter = 0
    while remaining:
        digest = hashlib.sha256(prefix + counter.to_bytes(8, "big")).digest()
        block = b"".join(EXPANSION[value] for value in digest)
        take = min(remaining, len(block))
        pending += block[:take]
        remaining -= take
        counter += 1
        while len(pending) >= 80:
            yield pending[:80] + b"\n"
            pending = pending[80:]
    if pending:
        yield pending + b"\n"


def write_fasta(output, size, seed):
    checksum = hashlib.sha256()
    for chunk in fasta_chunks(size, seed):
        output.write(chunk)
        checksum.update(chunk)
    output.flush()
    return checksum.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Generate deterministic A/C/G/T DNA in FASTA format.",
        epilog="Same version + size + seed = identical bytes. No external packages.",
    )
    parser.add_argument("size", type=size_argument, help="sequence length in bases")
    parser.add_argument("--seed", required=True, type=seed_argument,
                        help="explicit reproducibility seed (0 through 2^64 - 1)")
    parser.add_argument("-o", "--output", help="new output file; default: stdout")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    try:
        if args.output:
            with open(args.output, "xb") as output:
                checksum = write_fasta(output, args.size, args.seed)
        else:
            checksum = write_fasta(sys.stdout.buffer, args.size, args.seed)
    except OSError as error:
        parser.exit(1, "generation failed: {}\n".format(error))
    summary = {"algorithm": VERSION, "length": args.size, "seed": args.seed,
               "fasta_sha256": checksum}
    print(json.dumps(summary, sort_keys=True), file=sys.stderr)


if __name__ == "__main__":
    main()

