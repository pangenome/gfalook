# gfalook

A Rust tool for visualizing variation graphs from GFA files.

This is a reimplementation of [odgi viz](https://github.com/pangenome/odgi) that works directly with GFA files, avoiding the need to convert to ODGI format first.

## Features

- Direct GFA file input (no conversion required)
- PNG and SVG output formats
- Memory-efficient processing
- SHA256-based path coloring (matching odgi)
- Binned mode visualization with multiple coloring modes
- Path clustering by similarity
- Edge visualization
- X-axis coordinate display

## Installation

```bash
cargo build --release
```

## Usage

```bash
gfalook -i input.gfa -o output.png [OPTIONS]
```

### Core Options

```
-i, --idx <FILE>           Input GFA file (required)
-o, --out <FILE>           Output PNG/SVG file (required)
-x, --width <N>            Image width in pixels (default: 1500)
-y, --height <N>           Image height in pixels (default: 500)
-a, --path-height <N>      Height per path in pixels
-t, --threads <N>          Number of threads
-v, --verbose              Enable verbose logging
```

### Path Selection

```
-p, --paths-to-display <FILE>  List of paths to display
-I, --ignore-prefix <PREFIX>   Ignore paths with this prefix
-r, --range <PATH:START-END>   Display specific range
```

### Path Appearance

```
-H, --hide-path-names      Hide path names on the left
-n, --no-path-borders      Don't show path borders
-b, --black-path-borders   Draw path borders in black
-R, --pack-paths           Pack paths compactly (2D layout)
-L, --link-path-pieces <F> Draw connector lines between path gaps
```

### Coloring Modes

```
-s, --color-by-prefix <CHAR>    Color by path name prefix
-F, --path-colors <FILE>        Load custom path colors
-m, --color-by-mean-depth       Color by coverage depth (Spectral palette)
-z, --color-by-mean-inversion-rate  Color by strand orientation
-S, --show-strand               Show forward (blue) / reverse (red) strand
-N, --color-by-uncalled-bases   Color by N proportion (black to green)
-J, --highlight-node-ids <FILE> Highlight specific nodes in red
-B, --colorbrewer-palette <SCHEME:N>  Select palette (Spectral, RdBu, etc.)
-G, --no-grey-depth             Use full palette range for low coverage
```

### Gradient Mode

```
-d, --change-darkness      Vary color darkness by position in path
-l, --longest-path         Use longest path for darkness normalization
-u, --white-to-black       White-to-black gradient (with -d)
-A, --alignment-prefix <S> Apply effects only to matching paths
```

### Special Modes

```
-O, --compressed-mode      Single row showing mean coverage (RdBu palette)
-M, --prefix-merges <FILE> Merge paths with same prefix into groups
-k, --cluster-paths        Cluster paths by similarity
    --cluster-threshold    Similarity threshold for clusters
    --cluster-gap <N>      Gap between clusters (default: 10)
```

### X-Axis

```
--x-axis <COORD>           Show x-axis (pangenomic or path name)
--x-ticks <N>              Number of tick marks (default: 10)
--x-axis-absolute          Show absolute coordinates
```

## Examples

Basic visualization:
```bash
gfalook -i graph.gfa -o graph.png -v
```

Depth coloring with RdBu palette:
```bash
gfalook -i graph.gfa -o graph.png -m -B RdBu:11
```

Compressed coverage view:
```bash
gfalook -i graph.gfa -o graph.png -O
```

Pack paths compactly:
```bash
gfalook -i graph.gfa -o graph.png -R
```

Position-based darkness gradient:
```bash
gfalook -i graph.gfa -o graph.png -d -l
```

## Citation

Based on the `viz` command from [ODGI](https://github.com/pangenome/odgi):

> Guarracino A, Heumos S, Naber F, Panber P, Kelleher J, Garrison E. ODGI: understanding pangenome graphs. *Bioinformatics*. 2022;38(13):3319-3326. https://doi.org/10.1093/bioinformatics/btac308

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [odgi](https://github.com/pangenome/odgi) - Optimized Dynamic Genome/Graph Implementation
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
