# gfalook

A Rust tool for visualizing variation graphs from GFA files.

This is a reimplementation of [odgi viz](https://github.com/pangenome/odgi) that works directly with GFA files, avoiding the need to convert to ODGI format first.

## Citation

Based on the `viz` command from [ODGI](https://github.com/pangenome/odgi):

> Guarracino A, Heumos S, Naber F, Panber P, Kelleher J, Garrison E. ODGI: understanding pangenome graphs. *Bioinformatics*. 2022;38(13):3319-3326. https://doi.org/10.1093/bioinformatics/btac308

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

### Path clustering (`-k`)

Automatically cluster and reorder paths by similarity. Colored bars on the left indicate cluster membership. Combined with `-m` for better visibility:

```bash
gfalook -i graph.gfa -o clustered.png -x 1000 -y 500 -k -m
```

![Path clustering](images/clustered.png)

Show only cluster representatives (medoids) with `-K`. Each path label shows the cluster size:

```bash
gfalook -i graph.gfa -o clustered_representatives.png -x 1000 -y 500 -k -K -m
```

![Cluster representatives](images/clustered_representatives.png)

### Dendrogram visualization (`-k -D -m`)

Show hierarchical clustering tree alongside paths with depth coloring:

```bash
gfalook -i graph.gfa -o dendrogram.png -x 1000 -y 500 -k -D -m
```

![Dendrogram](images/dendrogram.png)

### Dendrogram with X-axis (`-k -D -m --x-axis`)

Combine dendrogram, depth coloring, and absolute X-axis coordinates:

```bash
gfalook -i graph.gfa -o dendrogram_xaxis.png -x 1000 -y 500 \
    -k -D -m --x-axis "chm13#chr6:31825251-31908851" --x-axis-absolute
```

![Dendrogram with X-axis](images/dendrogram_xaxis.png)

### UPGMA hierarchical clustering (`--use-upgma`)

Use pure UPGMA hierarchical clustering instead of DBSCAN. This creates clusters by cutting the tree at a height threshold:

```bash
gfalook -i graph.gfa -o upgma.png -x 1000 -y 500 -k --use-upgma -D -m
```

![UPGMA clustering](images/upgma.png)

Control cluster granularity with `--upgma-threshold` (0.0-1.0, lower = more clusters):

```bash
gfalook -i graph.gfa -o upgma_fine.png -x 1000 -y 500 -k --use-upgma --upgma-threshold 0.3 -D -m
```

![UPGMA fine clustering](images/upgma_fine.png)

### X-axis with pangenomic coordinates (`--x-axis pangenomic`)

Display coordinates based on node order in the graph:

```bash
gfalook -i graph.gfa -o xaxis_pangenomic.png -x 1000 -y 500 --x-axis pangenomic --x-ticks 5
```

![X-axis pangenomic](images/xaxis_pangenomic.png)

### X-axis with path reference (`--x-axis <path>`)

Display coordinates based on a reference path (e.g., chm13):

```bash
gfalook -i graph.gfa -o xaxis_chm13.png -x 1000 -y 500 \
    --x-axis "chm13#chr6:31825251-31908851" --x-ticks 5
```

![X-axis chm13](images/xaxis_chm13.png)

### X-axis with absolute coordinates (`--x-axis-absolute`)

Show absolute chromosome coordinates instead of relative positions:

```bash
gfalook -i graph.gfa -o xaxis_chm13_abs.png -x 1000 -y 500 \
    --x-axis "chm13#chr6:31825251-31908851" --x-ticks 5 --x-axis-absolute
```

![X-axis absolute](images/xaxis_chm13_abs.png)

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [odgi](https://github.com/pangenome/odgi) - Optimized Dynamic Genome/Graph Implementation
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
