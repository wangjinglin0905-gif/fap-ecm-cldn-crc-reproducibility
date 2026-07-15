from __future__ import annotations

import argparse
import hashlib
import math
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("output", type=Path)
    parser.add_argument("--size", type=int)
    parser.add_argument("--segments", type=int, default=8)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--checksum")
    parser.add_argument("--retries", type=int, default=8)
    return parser.parse_args()


def remote_size(url: str) -> int:
    request = urllib.request.Request(
        url,
        headers={"Range": "bytes=0-0", "User-Agent": "FAP-CLDN-reanalysis/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        content_range = response.headers.get("Content-Range")
        if content_range and "/" in content_range:
            return int(content_range.rsplit("/", 1)[1])
        content_length = response.headers.get("Content-Length")
        if content_length:
            return int(content_length)
    raise RuntimeError("Remote file size could not be determined")


def download_segment(
    url: str,
    output: Path,
    start: int,
    end: int,
    retries: int,
) -> Path:
    expected_size = end - start + 1
    if output.exists() and output.stat().st_size == expected_size:
        return output

    for attempt in range(1, retries + 1):
        temporary = output.with_suffix(output.suffix + ".tmp")
        request = urllib.request.Request(
            url,
            headers={
                "Range": f"bytes={start}-{end}",
                "User-Agent": "FAP-CLDN-reanalysis/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                if response.status != 206:
                    raise RuntimeError(
                        f"Server returned HTTP {response.status} for range {start}-{end}"
                    )
                with temporary.open("wb") as handle:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
            if temporary.stat().st_size != expected_size:
                raise RuntimeError(
                    f"Segment size mismatch: {temporary.stat().st_size} != {expected_size}"
                )
            temporary.replace(output)
            return output
        except Exception:
            if temporary.exists():
                temporary.unlink()
            if attempt == retries:
                raise
            time.sleep(min(30, 2**attempt))
    raise RuntimeError("Unreachable retry state")


def checksum(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    size = args.size or remote_size(args.url)
    segment_count = max(1, min(args.segments, size))
    segment_size = math.ceil(size / segment_count)
    part_dir = args.output.parent / f".{args.output.name}.parts"
    part_dir.mkdir(parents=True, exist_ok=True)

    ranges = []
    for index in range(segment_count):
        start = index * segment_size
        end = min(size - 1, start + segment_size - 1)
        if start <= end:
            ranges.append((index, start, end))

    completed = []
    with ThreadPoolExecutor(max_workers=min(args.workers, len(ranges))) as executor:
        futures = {
            executor.submit(
                download_segment,
                args.url,
                part_dir / f"part_{index:03d}",
                start,
                end,
                args.retries,
            ): (index, start, end)
            for index, start, end in ranges
        }
        for future in as_completed(futures):
            index, start, end = futures[future]
            completed.append(future.result())
            print(
                f"segment {index + 1}/{len(ranges)} verified ({start}-{end})",
                flush=True,
            )

    assembled = args.output.with_suffix(args.output.suffix + ".assembling")
    with assembled.open("wb") as destination:
        for index, _, _ in ranges:
            part_path = part_dir / f"part_{index:03d}"
            with part_path.open("rb") as source:
                while True:
                    block = source.read(1024 * 1024)
                    if not block:
                        break
                    destination.write(block)

    if assembled.stat().st_size != size:
        raise RuntimeError(f"Assembled size mismatch: {assembled.stat().st_size} != {size}")

    if args.checksum:
        algorithm, expected = args.checksum.split(":", 1)
        observed = checksum(assembled, algorithm)
        if observed.lower() != expected.lower():
            raise RuntimeError(f"Checksum mismatch: {observed} != {expected}")
        print(f"{algorithm} verified: {observed}", flush=True)

    os.replace(assembled, args.output)
    print(f"complete: {args.output} ({size} bytes)", flush=True)


if __name__ == "__main__":
    main()
