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
- Path clustering by similarity with dendrogram visualization
- Path annotations with legend (e.g., population labels)
- Edge visualization
- X-axis coordinate display
- Linear reference axis linked to graph columns (LD-plot style)

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

### CpG density coloring (`--color-by-cpg-density`)

Color bins by the density of CG dinucleotides inside each node, from light grey (0 CG per bp) to near-black. `--cpg-density-max` sets the density that maps to near-black (default 0.10 = 10 CG per 100 bp, a typical CpG island). CG dinucleotides that span two nodes are not counted.

```bash
gfalook -i graph.gfa -o cpg.png -x 1000 -y 500 --color-by-cpg-density
```

![CpG density coloring](images/cpg.png)

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

### Minimum cluster size (`--min-cluster-size`)

Merge clusters with fewer than N paths into their nearest cluster (single-linkage distance), starting from the smallest. On this graph, `-k` alone gives 7 clusters (50, 19, 13, 4, 2, 1, 1 paths), `--min-cluster-size 5` gives 4 (53, 19, 13, 5):

```bash
gfalook -i graph.gfa -o clustered_min_size.png -x 1000 -y 500 -k -m --min-cluster-size 5
```

![Minimum cluster size](images/clustered_min_size.png)

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

### Reference axis (`--ref-axis`)

The graph x-axis is not the coordinate system of any one genome, so a plot made on genome coordinates does not line up with it. `--ref-axis <PATH>` draws, under the graph, a linear axis of that path. As in LD plots, lines link each reference position to the graph columns that hold it. From top to bottom the panel shows:

1. a strip with the graph columns covered by each reference range,
2. lines from each range start and end in the graph to the same position on the linear axis,
3. the ranges as boxes on the linear axis, with ticks in reference coordinates.

The linear axis spans the full graph width, from 0 to the path length. A plot on reference coordinates with the same width and x limits lines up with it. `--x-axis-absolute` adds the subpath start from `name:start-end`.

```bash
gfalook -i graph.gfa -o ref_axis.png -x 1000 -y 500 -m \
    --ref-axis "chm13#chr6:31825251-31908851" --x-axis-absolute
```

![Reference axis](images/ref_axis.png)

`--ref-ranges` sets the ranges:

- `parts:N`: N equal parts (default `parts:10`).
- `step:BP`: one range every BP bp.
- `bubbles[:F]`: boundaries at nodes traversed once by at least a fraction F of the displayed paths (default F 0.9). Boundaries closer than `--ref-min-range` bp (default path length / 50) are merged. If no boundary is left, gfalook exits with an error.
- a BED file: path name, start, end, optional name (0-based, half-open, path-relative). The name is drawn under its box when it fits. A row with a different path name or a range beyond the path is an error.

```bash
gfalook -i graph.gfa -o ref_axis_bubbles.png -x 1000 -y 500 -m \
    --ref-axis "chm13#chr6:31825251-31908851" --ref-ranges bubbles
gfalook -i graph.gfa -o genes.png -x 1000 -y 500 -m \
    --ref-axis "chm13#chr6:31825251-31908851" --ref-ranges genes.bed
```

![Reference axis with bubble ranges](images/ref_axis_bubbles.png)

gfalook logs how many reference steps go back in graph order. When the graph order or strand differs from the reference, `--ref-sort <PATH>` orders the graph by it and the lines stop crossing. `--ref-sort` also flips nodes so the reference walks them forward, as `odgi groom -R` does: other nodes take the orientation most paths give them, relative to the reference. This haplotype is on the reverse strand: without `--ref-sort`, 1155 of 1156 steps go back and 1156 are reverse. With it, 0 go back and 0 are reverse.

```bash
gfalook -i graph.gfa -o ref_axis_reverse.png -x 1000 -y 500 -m \
    --ref-axis "HG00438#2#JAHBCA010000042.1:24398231-24449090" --ref-ranges parts:8
gfalook -i graph.gfa -o ref_axis_reverse_sorted.png -x 1000 -y 500 -m \
    --ref-sort "HG00438#2#JAHBCA010000042.1:24398231-24449090" \
    --ref-axis "HG00438#2#JAHBCA010000042.1:24398231-24449090" --ref-ranges parts:8
```

![Reverse-strand reference](images/ref_axis_reverse.png)
![Reverse-strand reference, sorted](images/ref_axis_reverse_sorted.png)

`--ref-sort` does not always help. `--ref-sort` ranks each node at its first visit, so a reference that crosses a repeat twice jumps back and forth in its second copy. For chm13 in this C4 graph, `--ref-sort` raises the steps that go back from 1 to 115. Compare the logged count with and without it.

### Path annotations (`-E`)

Label paths by category (e.g., population, sample type) using an annotation file. The annotation file is a CSV or TSV with two columns: `prefix` and `annotation`. Each prefix is matched against path names (longest match wins), and categories are displayed as a centered legend at the top with colored indicators next to each path.

```bash
gfalook -i graph.gfa -o annotation.png -x 1000 -y 500 \
    -E test/hprc_year1_samples.csv
```

![Path annotations](images/annotation.png)

Example annotation file (`test/hprc_year1_samples.csv`) - includes all 234 HPRC Year 1 samples with their population codes:
```
prefix,annotation
HG00438,CHS
HG00621,CHS
HG01891,ACB
...
CHM13,REF
GRCh38,REF
HG002,NA
HG005,NA
```

### Annotations with clustering (`-E -k -D -m --x-axis`)

Combine path annotations with clustering, dendrogram, depth coloring, and x-axis coordinates:

```bash
gfalook -i graph.gfa -o annotation_clustered.png -x 1000 -y 500 \
    -E test/hprc_year1_samples.csv -k -D -m \
    --x-axis "chm13#chr6:31825251-31908851" --x-axis-absolute
```

![Annotations with clustering](images/annotation_clustered.png)

This example shows all visualization features combined: legend at top, dendrogram on left, cluster bar, annotation bar (with a gap between them), path names, depth coloring, and x-axis with absolute coordinates at the bottom.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [odgi](https://github.com/pangenome/odgi) - Optimized Dynamic Genome/Graph Implementation
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
