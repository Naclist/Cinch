"""Deterministic triangular pair blocks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator


@dataclass(frozen=True)
class PairBlock:
    block_id: int
    pairs: tuple[tuple[int, int], ...]


def iter_pair_blocks(n_loci: int, max_pairs: int) -> Iterator[PairBlock]:
    """Yield every unordered pair exactly once using bounded pair storage."""

    if n_loci < 0:
        raise ValueError("n_loci must be non-negative")
    if max_pairs < 1:
        raise ValueError("max_pairs must be >= 1")
    block: list[tuple[int, int]] = []
    block_id = 0
    for left in range(n_loci - 1):
        for right in range(left + 1, n_loci):
            block.append((left, right))
            if len(block) == max_pairs:
                yield PairBlock(block_id, tuple(block))
                block_id += 1
                block.clear()
    if block:
        yield PairBlock(block_id, tuple(block))
