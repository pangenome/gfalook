# gfalook

A Rust tool for visualizing variation graphs from GFA files.

This is a reimplementation of [odgi viz](https://github.com/pangenome/odgi) that works directly with GFA files, avoiding the need to convert to ODGI format first.

## Features

- Direct GFA file input (no conversion required)
- Memory-efficient processing
- SHA256-based path coloring (matching odgi)
- Binned mode visualization
- Path name rendering with bitmap font
- Edge visualization
- Multiple coloring modes (-m, -z, -S)

## Installation

```bash
cargo build --release
```

## Usage

```bash
gfalook -i input.gfa -o output.png [OPTIONS]
```

### Options

```
-i, --idx <FILE>           Input GFA file (required)
-o, --out <FILE>           Output PNG file (required)
-x, --width <N>            Image width in pixels (default: 1500)
-y, --height <N>           Image height in pixels (default: 500)
-a, --path-height <N>      Height per path in pixels
-P, --progress             Show progress messages
-H, --hide-path-names      Hide path names on the left
-n, --no-path-borders      Don't show path borders
-b, --black-path-borders   Draw path borders in black
-m, --color-by-mean-depth  Color by coverage depth
-z, --color-by-mean-inversion-rate  Color by strand orientation
-S, --show-strand          Show forward/reverse strand coloring
-s, --color-by-prefix <CHAR>  Color by path name prefix
-F, --path-colors <FILE>   Load custom path colors from file
-p, --paths-to-display <FILE>  List of paths to display
-I, --ignore-prefix <PREFIX>  Ignore paths with this prefix
```

## Examples

Basic visualization:
```bash
gfalook -i graph.gfa -o graph.png -P
```

With depth coloring:
```bash
gfalook -i graph.gfa -o graph.png -m -P
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Related Projects

- [odgi](https://github.com/pangenome/odgi) - Optimized Dynamic Genome/Graph Implementation
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
