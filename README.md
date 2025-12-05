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

## Examples

All examples use the `chr6.C4.gfa` test graph.

### Basic visualization

```bash
gfalook -i graph.gfa -o basic.png -x 1000 -y 500
```

![Basic visualization](images/basic.png)

### Depth coloring (`-m`)

Color by coverage depth using Spectral palette:

```bash
gfalook -i graph.gfa -o depth.png -x 1000 -y 500 -m
```

![Depth coloring](images/depth.png)

### Depth coloring with RdBu palette (`-m -B`)

```bash
gfalook -i graph.gfa -o depth_rdbu.png -x 1000 -y 500 -m -B RdBu:11
```

![Depth with RdBu](images/depth_rdbu.png)

### Strand coloring (`-S`)

Show forward (blue) and reverse (red) strand orientation:

```bash
gfalook -i graph.gfa -o strand.png -x 1000 -y 500 -S
```

![Strand coloring](images/strand.png)

### Position-based darkness gradient (`-d -l`)

Color darkness varies by position within each path:

```bash
gfalook -i graph.gfa -o darkness.png -x 1000 -y 500 -d -l
```

![Darkness gradient](images/darkness.png)

### White-to-black gradient (`-d -u`)

```bash
gfalook -i graph.gfa -o white_black.png -x 1000 -y 500 -d -u
```

![White to black](images/white_black.png)

### Link path pieces (`-L`)

Draw connector lines between discontinuous path segments:

```bash
gfalook -i graph.gfa -o links.png -x 1000 -y 500 -L 0.3
```

![Link path pieces](images/links.png)

### Compressed mode (`-O`)

Single row showing mean coverage across all paths:

```bash
gfalook -i graph.gfa -o compressed.png -x 1000 -y 100 -O
```

![Compressed mode](images/compressed.png)

## Options Reference

### Core Options

| Option | Description |
|--------|-------------|
| `-i, --idx <FILE>` | Input GFA file (required) |
| `-o, --out <FILE>` | Output PNG/SVG file (required) |
| `-x, --width <N>` | Image width in pixels (default: 1500) |
| `-y, --height <N>` | Image height in pixels (default: 500) |
| `-a, --path-height <N>` | Height per path in pixels |
| `-t, --threads <N>` | Number of threads |
| `-v, --verbose` | Enable verbose logging |

### Path Selection

| Option | Description |
|--------|-------------|
| `-p, --paths-to-display <FILE>` | List of paths to display |
| `-I, --ignore-prefix <PREFIX>` | Ignore paths with this prefix |
| `-r, --range <PATH:START-END>` | Display specific range |

### Path Appearance

| Option | Description |
|--------|-------------|
| `-H, --hide-path-names` | Hide path names on the left |
| `-n, --no-path-borders` | Don't show path borders |
| `-b, --black-path-borders` | Draw path borders in black |
| `-R, --pack-paths` | Pack paths compactly (2D layout) |
| `-L, --link-path-pieces <F>` | Draw connector lines between path gaps |

### Coloring Modes

| Option | Description |
|--------|-------------|
| `-s, --color-by-prefix <CHAR>` | Color by path name prefix |
| `-F, --path-colors <FILE>` | Load custom path colors |
| `-m, --color-by-mean-depth` | Color by coverage depth |
| `-z, --color-by-mean-inversion-rate` | Color by strand orientation |
| `-S, --show-strand` | Show forward/reverse strand |
| `-N, --color-by-uncalled-bases` | Color by N proportion |
| `-J, --highlight-node-ids <FILE>` | Highlight specific nodes in red |
| `-B, --colorbrewer-palette <SCHEME:N>` | Select palette (Spectral, RdBu, RdYlGn, PiYG, PRGn, RdYlBu, BrBG) |
| `-G, --no-grey-depth` | Use full palette range for low coverage |

### Gradient Mode

| Option | Description |
|--------|-------------|
| `-d, --change-darkness` | Vary color darkness by position |
| `-l, --longest-path` | Use longest path for normalization |
| `-u, --white-to-black` | White-to-black gradient (with -d) |
| `-A, --alignment-prefix <S>` | Apply effects only to matching paths |

### Special Modes

| Option | Description |
|--------|-------------|
| `-O, --compressed-mode` | Single row showing mean coverage |
| `-M, --prefix-merges <FILE>` | Merge paths with same prefix |
| `-k, --cluster-paths` | Cluster paths by similarity |
| `--cluster-threshold <F>` | Similarity threshold for clusters |
| `--cluster-gap <N>` | Gap between clusters (default: 10) |

### X-Axis

| Option | Description |
|--------|-------------|
| `--x-axis <COORD>` | Show x-axis (pangenomic or path name) |
| `--x-ticks <N>` | Number of tick marks (default: 10) |
| `--x-axis-absolute` | Show absolute coordinates |

## Citation

Based on the `viz` command from [ODGI](https://github.com/pangenome/odgi):

> Guarracino A, Heumos S, Naber F, Panber P, Kelleher J, Garrison E. ODGI: understanding pangenome graphs. *Bioinformatics*. 2022;38(13):3319-3326. https://doi.org/10.1093/bioinformatics/btac308

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [odgi](https://github.com/pangenome/odgi) - Optimized Dynamic Genome/Graph Implementation
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
