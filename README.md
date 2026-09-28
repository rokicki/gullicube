# gullicube

A solver for the N×N×N cube, named for Gulliver, who travelled between the
tiny and the giant: it handles cubes from 2×2×2 to thousands of layers, and
aims for short solutions.

It works for any N ≥ 2. You can give it more time to get a shorter solution,
and it has a fast mode for very large cubes.  Typical performance is
about one move per sticker.

## Build

```sh
make            # builds ./gullicube (C++20; any recent clang or gcc)
```

## Run

Run `gullicube` from this directory: it reads its data files (`*.cls`) from
the same directory the executable is invoked as (not through path
lookup).

```sh
./gullicube                              # full usage
./gullicube -b 256 -r 5000 --verify 10  # 10x10x10, 5000-move scramble, beam 256, independent check
./gullicube -b 64 -o 12 > solution.txt  # random state; solution lines start with a space
./gullicube -i -o 16 < scramble.txt     # read a move scramble from stdin
./gullicube --fast 1024                 # huge cubes: minimum time
./gullicube -2 -b 1024 20               # widths 1, 2, 4, ... until Ctrl-C; keeps the best
```

- **Moves:** one per layer turn (a half turn counts 1). Notation: `3R2` turns
  the third layer from R by a half turn, `U'` is the U face
  counter-clockwise.
- **Seed:** without `--seed` a random state is chosen from system entropy,
  and the seed is printed as `# seed S` so the run can be repeated.
- **Checking:** `--verify` replays scramble + solution on an independent
  sticker simulator. It needs a move scramble (`-r`, `-R` or `-i`).
- **Effort:** the `-b` beam width is the main knob. Wings use 4× that width
  (`--wingfactor`), and the merge gets a time budget of half the pair search
  time (`--mergeshare`).
- **Memory:** about 3 GB (0.6 GB with `--fast`).

## How it solves

1. **Corners** (and, on odd cubes, middle edges and centre centres):
   optimal for just corners, two-phase for corners and middle edges.
2. **Wings**, one orbit per depth, with centres free: a beam search over mined
   wing algorithms, with a finish table of every state one algorithm solves.
   The per-depth solutions are merged to cancel setup turns.
3. **Oblique centre orbits**, solved in pairs (x,y)/(y,x): a beam over pair
   algorithms with exact one-algorithm finishes (finish table), with costs that
   count cancellation between consecutive algorithms.
4. **Diagonal and mid centre orbits**: the same machinery over algorithms that
   move only that orbit.
5. **One merge** of all centre algorithms, reordered within their dependency
   lattice to cancel moves at the seams: a merge beam for small cubes, a
   budgeted greedy merge for large ones.

Every phase has a guaranteed fallback (3-cycle tables), so a solve always
completes.

## Data files

Every algorithm table is closed under the 48 cube symmetries and inversion, so
each file stores one representative per symmetry class. The tables are
expanded when loaded, and the finish tables are rebuilt in parallel at each
start (about 2 s in total).

| file | contents | used for |
|---|---|---|
| `p8s10.cls` | pair algorithms ≤ 8 moves moving ≤ 10 pieces (408k) | pair beam |
| `fast6.cls` | pair algorithms ≤ 6 moves plus single-orbit 8-movers (23k) | pair beam in `--fast` |
| `c12x.cls`, `p9.cls`, `eg8.cls` | pair and single-orbit algorithms (17M) | pair finish table |
| `diag2.cls`, `mid2.cls` | diagonal / mid algorithms ≤ 8 moves (plus 8–12-move ones) | diagonal / mid beams |
| `diag_x10.cls`, `mid_x10.cls` | every diagonal / mid algorithm ≤ 10 moves (49M / 5.3M) | their finish tables |
| `wing_d6c.cls` | wing algorithms (centres free): ≤ 6 moves plus face-turn conjugates (1.2M) | wing beam and finish table |
| `wing3.cls` | brobdicube's wing 3-cycles | older wing solver and parity finisher |
| `pair3.cls`, `diag3.cls`, `mid3.cls` | complete 3-cycle tables | guaranteed fallback |

**How the tables were made:**
- **Pair tables:** twsearch enumerations of 8×8×8 sequences up to 11 moves that
  restore everything except one oblique pair and the diagonals, expanded 96
  ways.
- **Diagonal and mid tables:** those sequences with their two slice depths
  collapsed onto one depth (or onto the middle layer), kept only if they
  change nothing but their own orbit.
- **Wing tables:** an enumeration of face turns plus one slice depth that
  restore the corners and the other wings.
- **Checking:** every algorithm was validated on several cube sizes (even and
  odd, outer and innermost orbits).

## Source

- `src/full`: the solver.
- `src/pair`: the pair-orbit pool format, move groups and the independent
  sticker simulator used by `--verify`.
- `src/brob`: brobdicube's cube representation, move application and corner
  solver (brobdicube still unreleased), unmodified.

## License

This work is dual-licensed under the Mozilla Public License 2.0 and GPL 3.0 (or
any later version). If you use this work, you can choose either (or both) license terms to adhere to.
The license texts are in `LICENSE-MPL.md` and `LICENSE-GPL.md`.

This covers the whole distribution, including the brobdicube sources in `src/brob`
and the data files.

`SPDX-License-Identifier: MPL-2.0 OR GPL-3.0-or-later`
