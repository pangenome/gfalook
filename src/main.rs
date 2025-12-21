#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use clap::Parser;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "gfalook")]
#[command(about = "Visualize a variation graph in 1D.", long_about = None)]
struct Args {
    // === Input/Output ===
    /// Load the variation graph in GFA format from this FILE.
    #[arg(
        short = 'i',
        long = "idx",
        value_name = "FILE",
        help_heading = "Input/Output"
    )]
    idx: PathBuf,

    /// Write the visualization to this FILE (PNG or SVG based on extension).
    #[arg(
        short = 'o',
        long = "out",
        value_name = "FILE",
        help_heading = "Input/Output"
    )]
    out: PathBuf,

    // === Image Size ===
    /// Set the width in pixels of the output image.
    #[arg(
        short = 'x',
        long = "width",
        value_name = "N",
        default_value_t = 1500,
        help_heading = "Image Size"
    )]
    width: u32,

    /// Set the height in pixels of the output image.
    #[arg(
        short = 'y',
        long = "height",
        value_name = "N",
        default_value_t = 500,
        help_heading = "Image Size"
    )]
    height: u32,

    /// The height in pixels for a path.
    #[arg(
        short = 'a',
        long = "path-height",
        value_name = "N",
        default_value_t = 10,
        help_heading = "Image Size"
    )]
    path_height: u32,

    /// The padding in pixels on the x-axis for a path.
    #[arg(
        short = 'X',
        long = "path-x-padding",
        value_name = "N",
        default_value_t = 0,
        help_heading = "Image Size"
    )]
    path_x_padding: u32,

    // === Clustering ===
    /// Automatically order paths by similarity.
    #[arg(
        short = 'k',
        long = "cluster-paths",
        conflicts_with = "paths_to_display",
        help_heading = "Clustering"
    )]
    cluster_paths: bool,

    /// Similarity threshold for cluster detection (automatic if not specified).
    #[arg(
        long = "cluster-threshold",
        value_name = "F",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    cluster_threshold: Option<f64>,

    /// Use all nodes for clustering instead of only variable nodes.
    #[arg(
        long = "cluster-all-nodes",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    cluster_all_nodes: bool,

    /// Gap in pixels between clusters.
    #[arg(
        long = "cluster-gap",
        value_name = "N",
        requires = "cluster_paths",
        default_value_t = 10,
        help_heading = "Clustering"
    )]
    cluster_gap: u32,

    /// Maximum number of clusters allowed (automatic if not specified).
    #[arg(
        long = "max-clusters",
        value_name = "N",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    max_clusters: Option<usize>,

    /// Show only one representative path (medoid) per cluster.
    #[arg(
        short = 'K',
        long = "cluster-representatives",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    cluster_representatives: bool,

    /// Show dendrogram on the left (hierarchical clustering tree).
    #[arg(
        short = 'D',
        long = "dendrogram",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    dendrogram: bool,

    /// Width of the dendrogram in pixels.
    #[arg(
        long = "dendrogram-width",
        value_name = "PIXELS",
        default_value = "100",
        requires = "dendrogram",
        help_heading = "Clustering"
    )]
    dendrogram_width: u32,

    /// Use pure UPGMA hierarchical clustering instead of DBSCAN.
    /// Clusters are determined by cutting the tree at a height threshold.
    #[arg(
        long = "use-upgma",
        requires = "cluster_paths",
        help_heading = "Clustering"
    )]
    use_upgma: bool,

    /// Height threshold for cutting UPGMA tree (0.0-1.0, default: auto-detect).
    /// Lower values create more clusters, higher values create fewer.
    #[arg(
        long = "upgma-threshold",
        value_name = "THRESHOLD",
        requires = "use_upgma",
        help_heading = "Clustering"
    )]
    upgma_threshold: Option<f64>,

    // === Path Selection ===
    /// List of paths to display in the specified order.
    #[arg(
        short = 'p',
        long = "paths-to-display",
        value_name = "FILE",
        help_heading = "Path Selection"
    )]
    paths_to_display: Option<PathBuf>,

    /// Ignore paths starting with the given PREFIX.
    #[arg(
        short = 'I',
        long = "ignore-prefix",
        value_name = "PREFIX",
        help_heading = "Path Selection"
    )]
    ignore_prefix: Option<String>,

    /// Nucleotide range to visualize: STRING=[PATH:]start-end.
    #[arg(
        short = 'r',
        long = "path-range",
        value_name = "STRING",
        help_heading = "Path Selection"
    )]
    path_range: Option<String>,

    /// Merge paths beginning with prefixes listed in FILE.
    #[arg(
        short = 'M',
        long = "prefix-merges",
        value_name = "FILE",
        help_heading = "Path Selection"
    )]
    prefix_merges: Option<PathBuf>,

    // === Path Appearance ===
    /// Don't show path borders.
    #[arg(
        short = 'n',
        long = "no-path-borders",
        help_heading = "Path Appearance"
    )]
    no_path_borders: bool,

    /// Draw path borders in black (default is white).
    #[arg(
        short = 'b',
        long = "black-path-borders",
        help_heading = "Path Appearance"
    )]
    black_path_borders: bool,

    /// Pack all paths rather than displaying a single path per row.
    #[arg(short = 'R', long = "pack-paths", conflicts_with_all = ["paths_to_display", "compressed_mode", "prefix_merges", "cluster_paths"], help_heading = "Path Appearance")]
    pack_paths: bool,

    /// Show thin links of this relative width to connect path pieces.
    #[arg(
        short = 'L',
        long = "link-path-pieces",
        value_name = "FLOAT",
        help_heading = "Path Appearance"
    )]
    link_path_pieces: Option<f64>,

    // === Path Names ===
    /// Hide the path names on the left of the generated image.
    #[arg(short = 'H', long = "hide-path-names", help_heading = "Path Names")]
    hide_path_names: bool,

    /// Color path names background with the same color as paths.
    #[arg(
        short = 'C',
        long = "color-path-names-background",
        help_heading = "Path Names"
    )]
    color_path_names_background: bool,

    /// Maximum number of characters to display for each path name.
    #[arg(
        short = 'c',
        long = "max-num-of-characters",
        value_name = "N",
        help_heading = "Path Names"
    )]
    max_num_of_characters: Option<usize>,

    // === Coloring ===
    /// Color paths by their names looking at the prefix before the given character.
    #[arg(
        short = 's',
        long = "color-by-prefix",
        value_name = "CHAR",
        help_heading = "Coloring"
    )]
    color_by_prefix: Option<char>,

    /// Read per-path RGB colors from FILE.
    #[arg(
        short = 'F',
        long = "path-colors",
        value_name = "FILE",
        help_heading = "Coloring"
    )]
    path_colors: Option<PathBuf>,

    /// Use red and blue coloring to display forward and reverse alignments.
    #[arg(short = 'S', long = "show-strand", help_heading = "Coloring")]
    show_strand: bool,

    /// Change the color respect to the node strandness (black for forward, red for reverse).
    #[arg(
        short = 'z',
        long = "color-by-mean-inversion-rate",
        help_heading = "Coloring"
    )]
    color_by_mean_inversion_rate: bool,

    /// Change the color with respect to the uncalled bases.
    #[arg(
        short = 'N',
        long = "color-by-uncalled-bases",
        help_heading = "Coloring"
    )]
    color_by_uncalled_bases: bool,

    /// Color nodes listed in FILE in red and all other nodes in grey.
    #[arg(
        short = 'J',
        long = "highlight-node-ids",
        value_name = "FILE",
        help_heading = "Coloring"
    )]
    highlight_node_ids: Option<PathBuf>,

    // === Binned Mode ===
    /// The bin width specifies the size of each bin in the binned mode.
    #[arg(
        short = 'w',
        long = "bin-width",
        value_name = "bp",
        help_heading = "Binned Mode"
    )]
    bin_width: Option<f64>,

    /// Automatically set width so each node/segment gets at least 1 pixel.
    #[arg(long = "show-all-nodes", help_heading = "Binned Mode")]
    show_all_nodes: bool,

    /// Minimum width in pixels for each node (use with --show-all-nodes, default: 1).
    #[arg(
        long = "node-width",
        value_name = "N",
        default_value = "1",
        help_heading = "Binned Mode"
    )]
    node_width: u32,

    /// Change the color with respect to the mean coverage.
    #[arg(
        short = 'm',
        long = "color-by-mean-depth",
        help_heading = "Binned Mode"
    )]
    color_by_mean_depth: bool,

    /// Use the colorbrewer palette specified by SCHEME:N.
    #[arg(
        short = 'B',
        long = "colorbrewer-palette",
        value_name = "SCHEME:N",
        help_heading = "Binned Mode"
    )]
    colorbrewer_palette: Option<String>,

    /// Use the colorbrewer palette for <0.5x and ~1x coverage bins.
    #[arg(short = 'G', long = "no-grey-depth", help_heading = "Binned Mode")]
    no_grey_depth: bool,

    // === Gradient Mode ===
    /// Change the color darkness based on nucleotide position.
    #[arg(short = 'd', long = "change-darkness", help_heading = "Gradient Mode")]
    change_darkness: bool,

    /// Use the longest path length to change the color darkness.
    #[arg(short = 'l', long = "longest-path", help_heading = "Gradient Mode")]
    longest_path: bool,

    /// Change the color darkness from white to black.
    #[arg(short = 'u', long = "white-to-black", help_heading = "Gradient Mode")]
    white_to_black: bool,

    // === Special Modes ===
    /// Compress the view vertically, summarizing path coverage.
    #[arg(short = 'O', long = "compressed-mode", conflicts_with_all = ["cluster_paths", "prefix_merges"], help_heading = "Special Modes")]
    compressed_mode: bool,

    /// Apply alignment related visual motifs to paths which have this name prefix.
    #[arg(
        short = 'A',
        long = "alignment-prefix",
        value_name = "STRING",
        help_heading = "Special Modes"
    )]
    alignment_prefix: Option<String>,

    // === X-Axis ===
    /// Show x-axis with coordinates. Use "pangenomic" for node-order coordinates or a path name for path-based coordinates.
    #[arg(long = "x-axis", value_name = "COORD_SYSTEM", help_heading = "X-Axis")]
    x_axis: Option<String>,

    /// Number of ticks on the x-axis.
    #[arg(
        long = "x-ticks",
        value_name = "N",
        default_value_t = 10,
        help_heading = "X-Axis"
    )]
    x_ticks: u32,

    /// Show absolute coordinates by adding the subpath start position (from name:start-end format). Cannot be used with "pangenomic".
    #[arg(long = "x-axis-absolute", requires = "x_axis", help_heading = "X-Axis")]
    x_axis_absolute: bool,

    // === Annotation ===
    /// Load path annotations from TSV file (columns: prefix, annotation). Prefix matches path names.
    #[arg(
        short = 'E',
        long = "annotation-file",
        value_name = "FILE",
        help_heading = "Annotation"
    )]
    annotation_file: Option<PathBuf>,

    /// Column number for annotation values (1-based). Default: 2 for TSV, 4 for CSV.
    #[arg(
        long = "annotation-column",
        value_name = "N",
        requires = "annotation_file",
        help_heading = "Annotation"
    )]
    annotation_column: Option<usize>,

    /// Width of annotation bar in pixels.
    #[arg(
        long = "annotation-bar-width",
        value_name = "N",
        default_value = "10",
        requires = "annotation_file",
        help_heading = "Annotation"
    )]
    annotation_bar_width: u32,

    /// Height of legend area in pixels.
    #[arg(
        long = "legend-height",
        value_name = "N",
        default_value = "30",
        requires = "annotation_file",
        help_heading = "Annotation"
    )]
    legend_height: u32,

    // === Performance ===
    /// Number of threads to use for parallel operations.
    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        help_heading = "Performance"
    )]
    threads: Option<usize>,

    /// Verbosity level (0 = error, 1 = info, 2 = debug).
    #[arg(
        short = 'v',
        long = "verbose",
        value_name = "N",
        default_value_t = 1,
        help_heading = "Performance"
    )]
    verbose: u8,
}

/// A segment (node) in the graph
#[derive(Debug, Clone)]
struct Segment {
    sequence_len: u64,
    n_count: u64, // Number of uncalled bases (N's) in the sequence
}

/// An edge between two segments
#[derive(Debug, Clone)]
struct Edge {
    from_id: u64,
    from_rev: bool,
    to_id: u64,
    to_rev: bool,
}

/// A step in a path: (segment_id, is_reverse)
#[derive(Debug, Clone)]
struct PathStep {
    segment_id: u64,
    is_reverse: bool,
}

/// A path through the graph
#[derive(Debug, Clone)]
struct GfaPath {
    name: String,
    steps: Vec<PathStep>,
}

/// Minimal graph representation for visualization
struct Graph {
    segments: Vec<Segment>,
    segment_name_to_id: FxHashMap<String, u64>,
    segment_offsets: Vec<u64>,
    total_length: u64,
    paths: Vec<GfaPath>,
    edges: Vec<Edge>,
}

/// Canonical edge key for deduplication
fn edge_key(from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) -> (u64, bool, u64, bool) {
    // Normalize edge direction for deduplication
    if from_id < to_id || (from_id == to_id && !from_rev) {
        (from_id, from_rev, to_id, to_rev)
    } else {
        (to_id, !to_rev, from_id, !from_rev)
    }
}

impl Graph {
    fn new() -> Self {
        Graph {
            segments: Vec::new(),
            segment_name_to_id: FxHashMap::default(),
            segment_offsets: Vec::new(),
            total_length: 0,
            paths: Vec::new(),
            edges: Vec::new(),
        }
    }
}

/// 5x8 bitmap font (matching odgi's font5x8.h)
const FONT_5X8: [[u8; 8]; 128] = {
    let mut font = [[0u8; 8]; 128];
    font[b' ' as usize] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
    font[b'!' as usize] = [0x20, 0x20, 0x20, 0x20, 0x20, 0x00, 0x20, 0x00];
    font[b'"' as usize] = [0x50, 0x50, 0x50, 0x00, 0x00, 0x00, 0x00, 0x00];
    font[b'#' as usize] = [0x50, 0x50, 0xF8, 0x50, 0xF8, 0x50, 0x50, 0x00];
    font[b'$' as usize] = [0x20, 0x78, 0xA0, 0x70, 0x28, 0xF0, 0x20, 0x00];
    font[b'%' as usize] = [0xC0, 0xC8, 0x10, 0x20, 0x40, 0x98, 0x18, 0x00];
    font[b'&' as usize] = [0x40, 0xA0, 0xA0, 0x40, 0xA8, 0x90, 0x68, 0x00];
    font[b'\'' as usize] = [0x20, 0x20, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00];
    font[b'(' as usize] = [0x10, 0x20, 0x40, 0x40, 0x40, 0x20, 0x10, 0x00];
    font[b')' as usize] = [0x40, 0x20, 0x10, 0x10, 0x10, 0x20, 0x40, 0x00];
    font[b'*' as usize] = [0x00, 0x20, 0xA8, 0x70, 0xA8, 0x20, 0x00, 0x00];
    font[b'+' as usize] = [0x00, 0x20, 0x20, 0xF8, 0x20, 0x20, 0x00, 0x00];
    font[b',' as usize] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x20, 0x40];
    font[b'-' as usize] = [0x00, 0x00, 0x00, 0xF8, 0x00, 0x00, 0x00, 0x00];
    font[b'.' as usize] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x20, 0x00];
    font[b'/' as usize] = [0x00, 0x08, 0x10, 0x20, 0x40, 0x80, 0x00, 0x00];
    font[b'0' as usize] = [0x70, 0x88, 0x98, 0xA8, 0xC8, 0x88, 0x70, 0x00];
    font[b'1' as usize] = [0x20, 0x60, 0x20, 0x20, 0x20, 0x20, 0x70, 0x00];
    font[b'2' as usize] = [0x70, 0x88, 0x08, 0x30, 0x40, 0x80, 0xF8, 0x00];
    font[b'3' as usize] = [0xF8, 0x10, 0x20, 0x10, 0x08, 0x88, 0x70, 0x00];
    font[b'4' as usize] = [0x10, 0x30, 0x50, 0x90, 0xF8, 0x10, 0x10, 0x00];
    font[b'5' as usize] = [0xF8, 0x80, 0xF0, 0x08, 0x08, 0x88, 0x70, 0x00];
    font[b'6' as usize] = [0x30, 0x40, 0x80, 0xF0, 0x88, 0x88, 0x70, 0x00];
    font[b'7' as usize] = [0xF8, 0x08, 0x10, 0x20, 0x40, 0x40, 0x40, 0x00];
    font[b'8' as usize] = [0x70, 0x88, 0x88, 0x70, 0x88, 0x88, 0x70, 0x00];
    font[b'9' as usize] = [0x70, 0x88, 0x88, 0x78, 0x08, 0x10, 0x60, 0x00];
    font[b':' as usize] = [0x00, 0x00, 0x20, 0x00, 0x00, 0x20, 0x00, 0x00];
    font[b';' as usize] = [0x00, 0x00, 0x20, 0x00, 0x00, 0x20, 0x20, 0x40];
    font[b'<' as usize] = [0x08, 0x10, 0x20, 0x40, 0x20, 0x10, 0x08, 0x00];
    font[b'=' as usize] = [0x00, 0x00, 0xF8, 0x00, 0xF8, 0x00, 0x00, 0x00];
    font[b'>' as usize] = [0x80, 0x40, 0x20, 0x10, 0x20, 0x40, 0x80, 0x00];
    font[b'?' as usize] = [0x70, 0x88, 0x08, 0x10, 0x20, 0x00, 0x20, 0x00];
    font[b'@' as usize] = [0x70, 0x88, 0xB8, 0xA8, 0xB8, 0x80, 0x70, 0x00];
    font[b'A' as usize] = [0x70, 0x88, 0x88, 0xF8, 0x88, 0x88, 0x88, 0x00];
    font[b'B' as usize] = [0xF0, 0x88, 0x88, 0xF0, 0x88, 0x88, 0xF0, 0x00];
    font[b'C' as usize] = [0x70, 0x88, 0x80, 0x80, 0x80, 0x88, 0x70, 0x00];
    font[b'D' as usize] = [0xE0, 0x90, 0x88, 0x88, 0x88, 0x90, 0xE0, 0x00];
    font[b'E' as usize] = [0xF8, 0x80, 0x80, 0xF0, 0x80, 0x80, 0xF8, 0x00];
    font[b'F' as usize] = [0xF8, 0x80, 0x80, 0xF0, 0x80, 0x80, 0x80, 0x00];
    font[b'G' as usize] = [0x70, 0x88, 0x80, 0xB8, 0x88, 0x88, 0x70, 0x00];
    font[b'H' as usize] = [0x88, 0x88, 0x88, 0xF8, 0x88, 0x88, 0x88, 0x00];
    font[b'I' as usize] = [0x70, 0x20, 0x20, 0x20, 0x20, 0x20, 0x70, 0x00];
    font[b'J' as usize] = [0x38, 0x10, 0x10, 0x10, 0x10, 0x90, 0x60, 0x00];
    font[b'K' as usize] = [0x88, 0x90, 0xA0, 0xC0, 0xA0, 0x90, 0x88, 0x00];
    font[b'L' as usize] = [0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0xF8, 0x00];
    font[b'M' as usize] = [0x88, 0xD8, 0xA8, 0xA8, 0x88, 0x88, 0x88, 0x00];
    font[b'N' as usize] = [0x88, 0xC8, 0xA8, 0x98, 0x88, 0x88, 0x88, 0x00];
    font[b'O' as usize] = [0x70, 0x88, 0x88, 0x88, 0x88, 0x88, 0x70, 0x00];
    font[b'P' as usize] = [0xF0, 0x88, 0x88, 0xF0, 0x80, 0x80, 0x80, 0x00];
    font[b'Q' as usize] = [0x70, 0x88, 0x88, 0x88, 0xA8, 0x90, 0x68, 0x00];
    font[b'R' as usize] = [0xF0, 0x88, 0x88, 0xF0, 0xA0, 0x90, 0x88, 0x00];
    font[b'S' as usize] = [0x70, 0x88, 0x80, 0x70, 0x08, 0x88, 0x70, 0x00];
    font[b'T' as usize] = [0xF8, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x00];
    font[b'U' as usize] = [0x88, 0x88, 0x88, 0x88, 0x88, 0x88, 0x70, 0x00];
    font[b'V' as usize] = [0x88, 0x88, 0x88, 0x88, 0x88, 0x50, 0x20, 0x00];
    font[b'W' as usize] = [0x88, 0x88, 0x88, 0xA8, 0xA8, 0xD8, 0x88, 0x00];
    font[b'X' as usize] = [0x88, 0x88, 0x50, 0x20, 0x50, 0x88, 0x88, 0x00];
    font[b'Y' as usize] = [0x88, 0x88, 0x50, 0x20, 0x20, 0x20, 0x20, 0x00];
    font[b'Z' as usize] = [0xF8, 0x08, 0x10, 0x20, 0x40, 0x80, 0xF8, 0x00];
    font[b'[' as usize] = [0x70, 0x40, 0x40, 0x40, 0x40, 0x40, 0x70, 0x00];
    font[b'\\' as usize] = [0x00, 0x80, 0x40, 0x20, 0x10, 0x08, 0x00, 0x00];
    font[b']' as usize] = [0x70, 0x10, 0x10, 0x10, 0x10, 0x10, 0x70, 0x00];
    font[b'^' as usize] = [0x20, 0x50, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00];
    font[b'_' as usize] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF8, 0x00];
    font[b'`' as usize] = [0x40, 0x20, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00];
    font[b'a' as usize] = [0x00, 0x00, 0x70, 0x08, 0x78, 0x88, 0x78, 0x00];
    font[b'b' as usize] = [0x80, 0x80, 0xB0, 0xC8, 0x88, 0x88, 0xF0, 0x00];
    font[b'c' as usize] = [0x00, 0x00, 0x70, 0x80, 0x80, 0x88, 0x70, 0x00];
    font[b'd' as usize] = [0x08, 0x08, 0x68, 0x98, 0x88, 0x88, 0x78, 0x00];
    font[b'e' as usize] = [0x00, 0x00, 0x70, 0x88, 0xF8, 0x80, 0x70, 0x00];
    font[b'f' as usize] = [0x30, 0x48, 0x40, 0xE0, 0x40, 0x40, 0x40, 0x00];
    font[b'g' as usize] = [0x00, 0x00, 0x78, 0x88, 0x78, 0x08, 0x70, 0x00];
    font[b'h' as usize] = [0x80, 0x80, 0xB0, 0xC8, 0x88, 0x88, 0x88, 0x00];
    font[b'i' as usize] = [0x20, 0x00, 0x60, 0x20, 0x20, 0x20, 0x70, 0x00];
    font[b'j' as usize] = [0x10, 0x00, 0x30, 0x10, 0x10, 0x90, 0x60, 0x00];
    font[b'k' as usize] = [0x80, 0x80, 0x90, 0xA0, 0xC0, 0xA0, 0x90, 0x00];
    font[b'l' as usize] = [0x60, 0x20, 0x20, 0x20, 0x20, 0x20, 0x70, 0x00];
    font[b'm' as usize] = [0x00, 0x00, 0xD0, 0xA8, 0xA8, 0xA8, 0xA8, 0x00];
    font[b'n' as usize] = [0x00, 0x00, 0xB0, 0xC8, 0x88, 0x88, 0x88, 0x00];
    font[b'o' as usize] = [0x00, 0x00, 0x70, 0x88, 0x88, 0x88, 0x70, 0x00];
    font[b'p' as usize] = [0x00, 0x00, 0xF0, 0x88, 0xF0, 0x80, 0x80, 0x00];
    font[b'q' as usize] = [0x00, 0x00, 0x78, 0x88, 0x78, 0x08, 0x08, 0x00];
    font[b'r' as usize] = [0x00, 0x00, 0xB0, 0xC8, 0x80, 0x80, 0x80, 0x00];
    font[b's' as usize] = [0x00, 0x00, 0x70, 0x80, 0x70, 0x08, 0xF0, 0x00];
    font[b't' as usize] = [0x40, 0x40, 0xE0, 0x40, 0x40, 0x48, 0x30, 0x00];
    font[b'u' as usize] = [0x00, 0x00, 0x88, 0x88, 0x88, 0x98, 0x68, 0x00];
    font[b'v' as usize] = [0x00, 0x00, 0x88, 0x88, 0x88, 0x50, 0x20, 0x00];
    font[b'w' as usize] = [0x00, 0x00, 0x88, 0x88, 0xA8, 0xA8, 0x50, 0x00];
    font[b'x' as usize] = [0x00, 0x00, 0x88, 0x50, 0x20, 0x50, 0x88, 0x00];
    font[b'y' as usize] = [0x00, 0x00, 0x88, 0x88, 0x78, 0x08, 0x70, 0x00];
    font[b'z' as usize] = [0x00, 0x00, 0xF8, 0x10, 0x20, 0x40, 0xF8, 0x00];
    font[b'{' as usize] = [0x10, 0x20, 0x20, 0x40, 0x20, 0x20, 0x10, 0x00];
    font[b'|' as usize] = [0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x00];
    font[b'}' as usize] = [0x40, 0x20, 0x20, 0x10, 0x20, 0x20, 0x40, 0x00];
    font[b'~' as usize] = [0x00, 0x00, 0x40, 0xA8, 0x10, 0x00, 0x00, 0x00];
    font
};

const TRAILING_DOTS: [u8; 8] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xA8, 0x00];

/// ColorBrewer Spectral 11-class diverging palette (default for -m)
/// With two grey colors prepended for low coverage (matching odgi)
const COLORBREWER_SPECTRAL_13: [(u8, u8, u8); 13] = [
    (196, 196, 196), // < 0.5x coverage (very low)
    (128, 128, 128), // < 1.5x coverage (low)
    // Spectral 11 (reversed from standard order)
    (158, 1, 66),
    (213, 62, 79),
    (244, 109, 67),
    (253, 174, 97),
    (254, 224, 139),
    (255, 255, 191),
    (230, 245, 152),
    (171, 221, 164),
    (102, 194, 165),
    (50, 136, 189),
    (94, 79, 162),
];

/// ColorBrewer RdBu 11-class diverging palette (used for compressed mode by default)
const COLORBREWER_RDBU_11: [(u8, u8, u8); 11] = [
    (103, 0, 31),
    (178, 24, 43),
    (214, 96, 77),
    (244, 165, 130),
    (253, 219, 199),
    (247, 247, 247),
    (209, 229, 240),
    (146, 197, 222),
    (67, 147, 195),
    (33, 102, 172),
    (5, 48, 97),
];

/// ColorBrewer RdYlGn 11-class diverging palette
const COLORBREWER_RDYLGN_11: [(u8, u8, u8); 11] = [
    (165, 0, 38),
    (215, 48, 39),
    (244, 109, 67),
    (253, 174, 97),
    (254, 224, 139),
    (255, 255, 191),
    (217, 239, 139),
    (166, 217, 106),
    (102, 189, 99),
    (26, 152, 80),
    (0, 104, 55),
];

/// ColorBrewer PiYG 11-class diverging palette
const COLORBREWER_PIYG_11: [(u8, u8, u8); 11] = [
    (142, 1, 82),
    (197, 27, 125),
    (222, 119, 174),
    (241, 182, 218),
    (253, 224, 239),
    (247, 247, 247),
    (230, 245, 208),
    (184, 225, 134),
    (127, 188, 65),
    (77, 146, 33),
    (39, 100, 25),
];

/// ColorBrewer PRGn 11-class diverging palette
const COLORBREWER_PRGN_11: [(u8, u8, u8); 11] = [
    (64, 0, 75),
    (118, 42, 131),
    (153, 112, 171),
    (194, 165, 207),
    (231, 212, 232),
    (247, 247, 247),
    (217, 240, 211),
    (166, 219, 160),
    (90, 174, 97),
    (27, 120, 55),
    (0, 68, 27),
];

/// ColorBrewer RdYlBu 11-class diverging palette
const COLORBREWER_RDYLBU_11: [(u8, u8, u8); 11] = [
    (165, 0, 38),
    (215, 48, 39),
    (244, 109, 67),
    (253, 174, 97),
    (254, 224, 144),
    (255, 255, 191),
    (224, 243, 248),
    (171, 217, 233),
    (116, 173, 209),
    (69, 117, 180),
    (49, 54, 149),
];

/// ColorBrewer BrBG 11-class diverging palette
const COLORBREWER_BRBG_11: [(u8, u8, u8); 11] = [
    (84, 48, 5),
    (140, 81, 10),
    (191, 129, 45),
    (223, 194, 125),
    (246, 232, 195),
    (245, 245, 245),
    (199, 234, 229),
    (128, 205, 193),
    (53, 151, 143),
    (1, 102, 94),
    (0, 60, 48),
];

/// Get a ColorBrewer palette by name. Returns the palette colors or None if not found.
/// Supports: Spectral, RdBu, RdYlGn, PiYG, PRGn, RdYlBu, BrBG (all 11-class)
fn get_colorbrewer_palette(name: &str) -> Option<&'static [(u8, u8, u8)]> {
    match name.to_lowercase().as_str() {
        "spectral" => Some(&COLORBREWER_SPECTRAL_13[2..]), // Skip the 2 grey colors
        "rdbu" => Some(&COLORBREWER_RDBU_11),
        "rdylgn" => Some(&COLORBREWER_RDYLGN_11),
        "piyg" => Some(&COLORBREWER_PIYG_11),
        "prgn" => Some(&COLORBREWER_PRGN_11),
        "rdylbu" => Some(&COLORBREWER_RDYLBU_11),
        "brbg" => Some(&COLORBREWER_BRBG_11),
        _ => None,
    }
}

/// Parse colorbrewer palette argument "SCHEME:N" and return (scheme_name, n)
fn parse_colorbrewer_arg(arg: &str) -> Option<(String, usize)> {
    let parts: Vec<&str> = arg.split(':').collect();
    if parts.len() == 2 {
        if let Ok(n) = parts[1].parse::<usize>() {
            return Some((parts[0].to_string(), n));
        }
    }
    // Just scheme name without :N
    if parts.len() == 1 {
        return Some((parts[0].to_string(), 11));
    }
    None
}

/// ColorBrewer Set1 qualitative palette for cluster indicators
/// 9 distinct colors that are easy to distinguish
const CLUSTER_COLORS: [(u8, u8, u8); 9] = [
    (228, 26, 28),   // red
    (55, 126, 184),  // blue
    (77, 175, 74),   // green
    (152, 78, 163),  // purple
    (255, 127, 0),   // orange
    (255, 255, 51),  // yellow
    (166, 86, 40),   // brown
    (247, 129, 191), // pink
    (153, 153, 153), // grey
];

/// Get color for a cluster ID
fn get_cluster_color(cluster_id: usize) -> (u8, u8, u8) {
    CLUSTER_COLORS[cluster_id % CLUSTER_COLORS.len()]
}

/// ColorBrewer Set2 qualitative palette for annotations (8 pastel colors)
/// Distinct from CLUSTER_COLORS (Set1) to avoid confusion when both are displayed
const ANNOTATION_COLORS: [(u8, u8, u8); 8] = [
    (102, 194, 165), // teal
    (252, 141, 98),  // coral
    (141, 160, 203), // periwinkle
    (231, 138, 195), // pink
    (166, 216, 84),  // lime
    (255, 217, 47),  // yellow
    (229, 196, 148), // tan
    (179, 179, 179), // grey
];

/// ColorBrewer Paired palette for annotations (12 colors for more categories)
const ANNOTATION_COLORS_EXTENDED: [(u8, u8, u8); 12] = [
    (166, 206, 227), // light blue
    (31, 120, 180),  // blue
    (178, 223, 138), // light green
    (51, 160, 44),   // green
    (251, 154, 153), // light red
    (227, 26, 28),   // red
    (253, 191, 111), // light orange
    (255, 127, 0),   // orange
    (202, 178, 214), // light purple
    (106, 61, 154),  // purple
    (255, 255, 153), // light yellow
    (177, 89, 40),   // brown
];

/// Get color for an annotation category
fn get_annotation_color(category_index: usize, total_categories: usize) -> (u8, u8, u8) {
    if total_categories <= 8 {
        ANNOTATION_COLORS[category_index % ANNOTATION_COLORS.len()]
    } else {
        ANNOTATION_COLORS_EXTENDED[category_index % ANNOTATION_COLORS_EXTENDED.len()]
    }
}

/// Parse a GFA file efficiently
fn parse_gfa(path: &PathBuf) -> std::io::Result<Graph> {
    let mut graph = Graph::new();

    info!("Loading GFA file...");

    // First pass: collect segments
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("S\t") {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let name = parts[1].to_string();
                let seq = parts[2];
                let seq_len = seq.len() as u64;
                // Count uncalled bases (N's)
                let n_count = seq.bytes().filter(|&b| b == b'N' || b == b'n').count() as u64;
                let id = graph.segments.len() as u64;
                graph.segment_name_to_id.insert(name, id);
                graph.segments.push(Segment {
                    sequence_len: seq_len,
                    n_count,
                });
            }
        }
    }

    // Calculate segment offsets (linear layout)
    let mut offset = 0u64;
    for seg in &graph.segments {
        graph.segment_offsets.push(offset);
        offset += seg.sequence_len;
    }
    graph.total_length = offset;

    info!(
        "Found {} segments, total length: {} bp",
        graph.segments.len(),
        graph.total_length
    );

    // Use a set to deduplicate edges
    let mut edge_set: std::collections::HashSet<(u64, bool, u64, bool)> =
        std::collections::HashSet::new();

    // Second pass: collect paths and edges (from L-lines)
    let file2 = File::open(path)?;
    let reader2 = BufReader::new(file2);
    for line in reader2.lines() {
        let line = line?;
        if line.starts_with("P\t") {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path_name = parts[1].to_string();
                let segments_str = parts[2];
                let mut steps = Vec::new();

                for seg in segments_str.split(',') {
                    let seg = seg.trim();
                    if seg.is_empty() {
                        continue;
                    }
                    let (name, is_reverse) = if let Some(stripped) = seg.strip_suffix('+') {
                        (stripped, false)
                    } else if let Some(stripped) = seg.strip_suffix('-') {
                        (stripped, true)
                    } else {
                        (seg, false)
                    };
                    if let Some(&id) = graph.segment_name_to_id.get(name) {
                        steps.push(PathStep {
                            segment_id: id,
                            is_reverse,
                        });
                    }
                }

                graph.paths.push(GfaPath {
                    name: path_name,
                    steps,
                });
            }
        } else if line.starts_with("W\t") {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 7 {
                let sample = parts[1];
                let hap = parts[2];
                let seq = parts[3];
                let walk_str = parts[6];

                let path_name = format!("{}#{}#{}", sample, hap, seq);
                let mut steps = Vec::new();

                let mut chars = walk_str.chars().peekable();
                while let Some(c) = chars.next() {
                    if c == '>' || c == '<' {
                        let is_reverse = c == '<';
                        let mut seg_name = String::new();
                        while let Some(&nc) = chars.peek() {
                            if nc == '>' || nc == '<' {
                                break;
                            }
                            seg_name.push(chars.next().unwrap());
                        }
                        if !seg_name.is_empty() {
                            if let Some(&id) = graph.segment_name_to_id.get(&seg_name) {
                                steps.push(PathStep {
                                    segment_id: id,
                                    is_reverse,
                                });
                            }
                        }
                    }
                }

                graph.paths.push(GfaPath {
                    name: path_name,
                    steps,
                });
            }
        } else if line.starts_with("L\t") {
            // Parse edge: L<TAB>from<TAB>from_orient<TAB>to<TAB>to_orient<TAB>overlap
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let from_name = parts[1];
                let from_orient = parts[2];
                let to_name = parts[3];
                let to_orient = parts[4];

                if let (Some(&from_id), Some(&to_id)) = (
                    graph.segment_name_to_id.get(from_name),
                    graph.segment_name_to_id.get(to_name),
                ) {
                    let from_rev = from_orient == "-";
                    let to_rev = to_orient == "-";
                    edge_set.insert(edge_key(from_id, from_rev, to_id, to_rev));
                }
            }
        }
    }

    // Third pass: add edges from consecutive path steps (implicit edges)
    for path in &graph.paths {
        for window in path.steps.windows(2) {
            let from = &window[0];
            let to = &window[1];
            // Edge from end of 'from' to start of 'to'
            // from_rev=true means we're going through from in reverse, so edge starts from beginning
            // to_rev=true means we're entering to in reverse, so edge goes to end
            edge_set.insert(edge_key(
                from.segment_id,
                from.is_reverse,
                to.segment_id,
                to.is_reverse,
            ));
        }
    }

    // Convert edge set to vector
    for (from_id, from_rev, to_id, to_rev) in edge_set {
        graph.edges.push(Edge {
            from_id,
            from_rev,
            to_id,
            to_rev,
        });
    }

    info!(
        "Found {} paths, {} edges",
        graph.paths.len(),
        graph.edges.len()
    );

    Ok(graph)
}

/// Compute SHA256-based path color (matching odgi algorithm exactly)
fn compute_path_color(path_name: &str, color_by_prefix: Option<char>) -> (u8, u8, u8) {
    let hash_input = if let Some(sep) = color_by_prefix {
        path_name.split(sep).next().unwrap_or(path_name)
    } else {
        path_name
    };

    let mut hasher = Sha256::new();
    hasher.update(hash_input.as_bytes());
    let result = hasher.finalize();

    let path_r = result[24];
    let path_g = result[8];
    let path_b = result[16];

    // Normalize to 0-1 range (divide by 255)
    let mut path_r_f = path_r as f32 / 255.0;
    let mut path_g_f = path_g as f32 / 255.0;
    let mut path_b_f = path_b as f32 / 255.0;

    // Normalize by sum
    let sum = path_r_f + path_g_f + path_b_f;
    if sum > 0.0 {
        path_r_f /= sum;
        path_g_f /= sum;
        path_b_f /= sum;
    }

    // Brighten the color (matching odgi's exact formula)
    let max_component = path_r_f.max(path_g_f).max(path_b_f);
    let f = if max_component > 0.0 {
        1.5f32.min(1.0 / max_component)
    } else {
        1.0
    };

    let r_out = (255.0 * (path_r_f * f).min(1.0)).round() as u8;
    let g_out = (255.0 * (path_g_f * f).min(1.0)).round() as u8;
    let b_out = (255.0 * (path_b_f * f).min(1.0)).round() as u8;

    (r_out, g_out, b_out)
}

fn load_path_colors(path: &PathBuf) -> std::io::Result<FxHashMap<String, (u8, u8, u8)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut colors = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let path_name = parts[0].to_string();
            let color_str = parts[1];

            let rgb = if color_str.starts_with('#') {
                let r = u8::from_str_radix(&color_str[1..3], 16).unwrap_or(0);
                let g = u8::from_str_radix(&color_str[3..5], 16).unwrap_or(0);
                let b = u8::from_str_radix(&color_str[5..7], 16).unwrap_or(0);
                (r, g, b)
            } else {
                let rgb_parts: Vec<u8> = color_str
                    .split(',')
                    .filter_map(|s| s.trim().parse().ok())
                    .collect();
                if rgb_parts.len() == 3 {
                    (rgb_parts[0], rgb_parts[1], rgb_parts[2])
                } else {
                    (128, 128, 128)
                }
            };

            colors.insert(path_name, rgb);
        }
    }

    Ok(colors)
}

fn load_paths_to_display(path: &PathBuf) -> std::io::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() {
            paths.push(line.to_string());
        }
    }

    Ok(paths)
}

/// Load node IDs to highlight from a file (one ID per line)
fn load_highlight_node_ids(path: &PathBuf) -> std::io::Result<FxHashSet<u64>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut node_ids = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() {
            if let Ok(id) = line.parse::<u64>() {
                node_ids.insert(id);
            }
        }
    }

    Ok(node_ids)
}

/// Result of path grouping by prefix
struct PathGrouping {
    /// For each original path index, the group index (-1 if not grouped)
    path_to_group: Vec<i64>,
    /// List of valid prefixes (group names)
    prefixes: Vec<String>,
    /// Number of groups
    num_groups: usize,
}

/// Load prefixes and create path groupings
fn load_prefix_merges(path: &PathBuf, paths: &[GfaPath]) -> std::io::Result<PathGrouping> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut prefixes_tmp: Vec<String> = Vec::new();
    let mut seen: FxHashSet<String> = FxHashSet::default();

    // Read all prefixes from file
    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if !line.is_empty() {
            if seen.contains(&line) {
                eprintln!("[gfalook] warning: duplicate prefix found: {}", line);
            } else {
                prefixes_tmp.push(line.clone());
                seen.insert(line);
            }
        }
    }

    let mut path_to_group: Vec<i64> = vec![-1; paths.len()];
    let mut prefixes: Vec<String> = Vec::new();

    // Assign each path to a group based on matching prefix
    for (path_idx, gfa_path) in paths.iter().enumerate() {
        let mut found = false;

        // First search in already validated prefixes
        for (group_idx, prefix) in prefixes.iter().enumerate() {
            if gfa_path.name.starts_with(prefix) {
                path_to_group[path_idx] = group_idx as i64;
                found = true;
                break;
            }
        }

        // If not found, search in all read prefixes
        if !found {
            for prefix in &prefixes_tmp {
                if gfa_path.name.starts_with(prefix) {
                    let group_idx = prefixes.len();
                    prefixes.push(prefix.clone());
                    path_to_group[path_idx] = group_idx as i64;
                    break;
                }
            }
        }
    }

    let num_groups = prefixes.len();
    Ok(PathGrouping {
        path_to_group,
        prefixes,
        num_groups,
    })
}

/// Annotation data loaded from TSV file
struct AnnotationData {
    /// Map from prefix to annotation category
    prefix_to_annotation: FxHashMap<String, String>,
    /// Ordered list of prefixes (sorted by length descending for longest-match-first)
    prefixes: Vec<String>,
    /// Ordered list of unique categories (sorted alphabetically)
    categories: Vec<String>,
    /// Map from category name to assigned color
    category_colors: FxHashMap<String, (u8, u8, u8)>,
}
/// Parse a CSV line handling quoted fields that may contain commas
fn parse_csv_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;

    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c == '"' {
            if in_quotes && i + 1 < chars.len() && chars[i + 1] == '"' {
                // Escaped quote
                current.push('"');
                i += 1;
            } else {
                in_quotes = !in_quotes;
            }
        } else if c == ',' && !in_quotes {
            fields.push(current.trim().to_string());
            current = String::new();
        } else {
            current.push(c);
        }
        i += 1;
    }
    fields.push(current.trim().to_string());
    fields
}

impl AnnotationData {
    /// Find annotation for a path by matching against prefixes (longest match wins)
    fn get_annotation(&self, path_name: &str) -> Option<&String> {
        for prefix in &self.prefixes {
            if path_name.starts_with(prefix) {
                return self.prefix_to_annotation.get(prefix);
            }
        }
        None
    }
}

/// Load path annotations from a TSV/CSV file
/// Expected format: prefix,annotation (first line is header)
/// The prefix column matches path names that start with that prefix
/// Supports both TSV (tab-separated) and CSV (comma-separated) based on file extension
/// annotation_column is 1-based (None = auto: 2 for TSV, 4 for CSV/HPRC format)
fn load_annotations(path: &PathBuf, annotation_column: Option<usize>) -> std::io::Result<AnnotationData> {
    // Read file as bytes and convert lossy to handle non-UTF8 characters
    let bytes = std::fs::read(path)?;
    let content = String::from_utf8_lossy(&bytes);

    // Detect delimiter based on file extension
    let is_csv = path
        .extension()
        .map(|e| e.to_string_lossy().to_lowercase() == "csv")
        .unwrap_or(false);

    // Determine annotation column (0-based internally)
    // Default: column 2 (index 1) for TSV, column 4 (index 3) for CSV
    let ann_col_idx = annotation_column
        .map(|c| c.saturating_sub(1))
        .unwrap_or(if is_csv { 3 } else { 1 });

    let mut prefix_to_annotation: FxHashMap<String, String> = FxHashMap::default();
    let mut categories_set: FxHashSet<String> = FxHashSet::default();
    let mut is_first_line = true;

    for line in content.lines() {
        let line = line.trim();

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        // Skip header line (first non-empty line)
        if is_first_line {
            is_first_line = false;
            continue;
        }

        // Parse fields based on delimiter
        let fields: Vec<String> = if is_csv {
            parse_csv_fields(line)
        } else {
            line.split('\t').map(|s| s.to_string()).collect()
        };

        // Get prefix (column 0) and annotation (specified column)
        if fields.len() > ann_col_idx {
            let prefix = fields[0].clone();
            let annotation = fields[ann_col_idx].clone();

            if !annotation.is_empty() && !prefix.is_empty() {
                categories_set.insert(annotation.clone());
                prefix_to_annotation.insert(prefix, annotation);
            }
        }
    }

    // Sort prefixes by length descending (longest match first)
    let mut prefixes: Vec<String> = prefix_to_annotation.keys().cloned().collect();
    prefixes.sort_by(|a, b| b.len().cmp(&a.len()));

    // Sort categories alphabetically for consistent ordering
    let mut categories: Vec<String> = categories_set.into_iter().collect();
    categories.sort();

    // Assign colors to categories
    let total = categories.len();
    let category_colors: FxHashMap<String, (u8, u8, u8)> = categories
        .iter()
        .enumerate()
        .map(|(i, cat)| (cat.clone(), get_annotation_color(i, total)))
        .collect();

    Ok(AnnotationData {
        prefix_to_annotation,
        prefixes,
        categories,
        category_colors,
    })
}

/// Result of path clustering
struct ClusteringResult {
    ordering: Vec<usize>,
    cluster_ids: Vec<usize>,
    num_clusters: usize,
    representatives: Vec<usize>, // medoid index (into original paths array) per cluster
    cluster_sizes: Vec<usize>,   // member count per cluster
    dendrogram: Option<Dendrogram>, // hierarchical clustering tree
}

/// A node in the dendrogram tree
#[derive(Clone, Debug)]
struct DendrogramNode {
    left: usize,  // index of left child (< n means leaf, >= n means internal node)
    right: usize, // index of right child
    height: f64,  // merge height (distance at which clusters merged)
    #[allow(dead_code)]
    size: usize, // number of leaves in this subtree
}

/// Dendrogram structure for hierarchical clustering visualization
#[derive(Clone, Debug)]
struct Dendrogram {
    nodes: Vec<DendrogramNode>, // internal nodes (n-1 nodes for n leaves)
    leaf_order: Vec<usize>,     // optimal leaf ordering for visualization
    max_height: f64,            // maximum merge height
}

/// Build a dendrogram using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
/// cluster_assignments: DBSCAN cluster IDs for each path (used to constrain merging order)
/// If provided, merging happens within clusters first, then between clusters
fn build_dendrogram(dist_matrix: &[Vec<f64>], cluster_assignments: Option<&[usize]>) -> Dendrogram {
    let n = dist_matrix.len();
    if n == 0 {
        return Dendrogram {
            nodes: Vec::new(),
            leaf_order: Vec::new(),
            max_height: 0.0,
        };
    }
    if n == 1 {
        return Dendrogram {
            nodes: Vec::new(),
            leaf_order: vec![0],
            max_height: 0.0,
        };
    }

    // Working distance matrix (will be modified during clustering)
    let mut dists: Vec<Vec<f64>> = dist_matrix.to_vec();

    // Track which cluster each index belongs to (-1 means merged)
    // Cluster IDs: 0..n are leaves, n..2n-1 are internal nodes
    let mut cluster_id: Vec<isize> = (0..n as isize).collect();
    let mut cluster_sizes: Vec<usize> = vec![1; n];

    // Track DBSCAN cluster for each active node (for constrained merging)
    let dbscan_cluster: Vec<Option<usize>> = if let Some(assignments) = cluster_assignments {
        assignments.iter().map(|&c| Some(c)).collect()
    } else {
        vec![None; n]
    };

    let mut nodes: Vec<DendrogramNode> = Vec::with_capacity(n - 1);
    let mut max_height: f64 = 0.0;

    // Track children (leaf indices) for each cluster ID
    // Key: cluster_id (0..n for leaves, n..2n-1 for internal nodes)
    // Value: list of original leaf indices in this cluster
    let mut children: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        children.insert(i, vec![i]);
    }

    for merge_idx in 0..(n - 1) {
        // Find minimum distance pair
        // If cluster_assignments provided, prefer merging within same DBSCAN cluster first
        let mut min_dist = f64::MAX;
        let mut min_i = 0;
        let mut min_j = 0;
        let mut found_same_cluster = false;

        // First pass: look for merges within same DBSCAN cluster
        if cluster_assignments.is_some() {
            for i in 0..n {
                if cluster_id[i] < 0 {
                    continue;
                }
                for j in (i + 1)..n {
                    if cluster_id[j] < 0 {
                        continue;
                    }
                    // Only consider if both are in the same DBSCAN cluster
                    if dbscan_cluster[i] == dbscan_cluster[j] && dists[i][j] < min_dist {
                        min_dist = dists[i][j];
                        min_i = i;
                        min_j = j;
                        found_same_cluster = true;
                    }
                }
            }
        }

        // Second pass: if no same-cluster merge found, allow any merge
        if !found_same_cluster {
            min_dist = f64::MAX;
            for i in 0..n {
                if cluster_id[i] < 0 {
                    continue;
                }
                for j in (i + 1)..n {
                    if cluster_id[j] < 0 {
                        continue;
                    }
                    if dists[i][j] < min_dist {
                        min_dist = dists[i][j];
                        min_i = i;
                        min_j = j;
                    }
                }
            }
        }

        let new_cluster_id = (n + merge_idx) as isize;
        let left_id = cluster_id[min_i] as usize;
        let right_id = cluster_id[min_j] as usize;
        let left_size = cluster_sizes[min_i];
        let right_size = cluster_sizes[min_j];
        let new_size = left_size + right_size;

        // Record the merge
        nodes.push(DendrogramNode {
            left: left_id,
            right: right_id,
            height: min_dist / 2.0, // UPGMA uses half the distance as height
            size: new_size,
        });
        max_height = max_height.max(min_dist / 2.0);

        // Merge children lists for leaf ordering
        // Get children of left and right clusters by their cluster IDs
        let left_children = children.get(&left_id).cloned().unwrap_or_default();
        let right_children = children.get(&right_id).cloned().unwrap_or_default();
        let mut new_children = left_children;
        new_children.extend(right_children);
        children.insert(new_cluster_id as usize, new_children);

        // Update distances using UPGMA formula
        for k in 0..n {
            if k == min_i || k == min_j || cluster_id[k] < 0 {
                continue;
            }
            let new_dist = (dists[min_i][k] * left_size as f64
                + dists[min_j][k] * right_size as f64)
                / new_size as f64;
            dists[min_i][k] = new_dist;
            dists[k][min_i] = new_dist;
        }

        // Mark min_j as merged into min_i
        cluster_id[min_j] = -1;
        cluster_id[min_i] = new_cluster_id;
        cluster_sizes[min_i] = new_size;
    }

    // Get leaf order from the root (the last cluster ID created)
    let root_cluster_id = 2 * n - 2;
    let leaf_order = children.get(&root_cluster_id).cloned().unwrap_or_default();

    Dendrogram {
        nodes,
        leaf_order,
        max_height,
    }
}

/// Cut the dendrogram tree at a given height threshold and return cluster assignments.
/// Returns a vector where cluster_ids[i] is the cluster ID for leaf i.
fn cut_dendrogram_at_height(dendrogram: &Dendrogram, threshold: f64) -> Vec<usize> {
    let n_leaves = dendrogram.leaf_order.len();
    if n_leaves == 0 {
        return Vec::new();
    }
    if dendrogram.nodes.is_empty() {
        // Single leaf - one cluster
        return vec![0];
    }

    // Use Union-Find to track which leaves belong to which cluster
    let mut uf = UnionFind::new(n_leaves);

    // Process merges in order (they're already sorted by height due to UPGMA)
    for node in &dendrogram.nodes {
        if node.height <= threshold {
            // This merge happens below the threshold, so merge the clusters
            // Need to find representative leaves from left and right subtrees
            let left_leaf = find_leftmost_leaf(dendrogram, node.left, n_leaves);
            let right_leaf = find_leftmost_leaf(dendrogram, node.right, n_leaves);
            uf.union(left_leaf, right_leaf);
        }
        // If height > threshold, don't merge - these become separate clusters
    }

    // Assign consecutive cluster IDs
    let mut root_to_cluster: FxHashMap<usize, usize> = FxHashMap::default();
    let mut cluster_ids = Vec::with_capacity(n_leaves);
    let mut next_cluster = 0;

    for i in 0..n_leaves {
        let root = uf.find(i);
        let cluster = *root_to_cluster.entry(root).or_insert_with(|| {
            let c = next_cluster;
            next_cluster += 1;
            c
        });
        cluster_ids.push(cluster);
    }

    cluster_ids
}

/// Find the leftmost (smallest index) leaf in a subtree
fn find_leftmost_leaf(dendrogram: &Dendrogram, node_idx: usize, n_leaves: usize) -> usize {
    if node_idx < n_leaves {
        return node_idx;
    }
    let internal_idx = node_idx - n_leaves;
    if internal_idx >= dendrogram.nodes.len() {
        return 0;
    }
    // Recursively find leftmost leaf in left subtree
    find_leftmost_leaf(dendrogram, dendrogram.nodes[internal_idx].left, n_leaves)
}

/// Find optimal threshold for UPGMA tree cutting using the "elbow" method.
/// Looks for the largest gap in merge heights.
fn find_optimal_upgma_threshold(dendrogram: &Dendrogram, max_clusters: Option<usize>) -> f64 {
    if dendrogram.nodes.is_empty() {
        return 0.5;
    }

    let n_leaves = dendrogram.leaf_order.len();
    let max_clusters = max_clusters.unwrap_or_else(|| n_leaves.div_ceil(9)); // ~11% like DBSCAN

    // Collect all merge heights
    let mut heights: Vec<f64> = dendrogram.nodes.iter().map(|n| n.height).collect();
    heights.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Find threshold that gives approximately max_clusters clusters
    // Start from highest threshold (fewer clusters) and work down
    for i in (0..heights.len()).rev() {
        let threshold = heights[i];
        let clusters = cut_dendrogram_at_height(dendrogram, threshold);
        let num_clusters = clusters.iter().max().map(|&m| m + 1).unwrap_or(1);

        if num_clusters >= max_clusters {
            // Found a good threshold
            debug!(
                "UPGMA auto-threshold: {:.4} gives {} clusters (target: {})",
                threshold, num_clusters, max_clusters
            );
            return threshold;
        }
    }

    // If we couldn't find a good threshold, use the smallest height
    let threshold = heights.first().copied().unwrap_or(0.0);
    debug!("UPGMA using minimum threshold: {:.4}", threshold);
    threshold
}

/// Union-Find data structure for DBSCAN clustering
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]); // Path compression
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let px = self.find(x);
        let py = self.find(y);
        if px != py {
            // Union by rank
            if self.rank[px] < self.rank[py] {
                self.parent[px] = py;
            } else if self.rank[px] > self.rank[py] {
                self.parent[py] = px;
            } else {
                self.parent[py] = px;
                self.rank[px] += 1;
            }
        }
    }

    fn count_clusters(&mut self) -> usize {
        let n = self.parent.len();
        let mut roots: FxHashSet<usize> = FxHashSet::default();
        for i in 0..n {
            roots.insert(self.find(i));
        }
        roots.len()
    }
}

/// Run DBSCAN with minPts=1 on distance matrix, return number of clusters
/// With minPts=1, DBSCAN is equivalent to finding connected components
/// where edges exist for distance <= eps
fn dbscan_count_clusters(dist_matrix: &[Vec<f64>], eps: f64) -> usize {
    let n = dist_matrix.len();
    if n == 0 {
        return 0;
    }

    let mut uf = UnionFind::new(n);

    // Connect points within eps distance
    for i in 0..n {
        for j in (i + 1)..n {
            if dist_matrix[i][j] <= eps {
                uf.union(i, j);
            }
        }
    }

    uf.count_clusters()
}

/// Run DBSCAN with minPts=1, return cluster assignments
fn dbscan_cluster(dist_matrix: &[Vec<f64>], eps: f64) -> Vec<usize> {
    let n = dist_matrix.len();
    if n == 0 {
        return Vec::new();
    }

    let mut uf = UnionFind::new(n);

    // Connect points within eps distance
    for i in 0..n {
        for j in (i + 1)..n {
            if dist_matrix[i][j] <= eps {
                uf.union(i, j);
            }
        }
    }

    // Assign cluster IDs (0-indexed, consecutive)
    let mut root_to_cluster: FxHashMap<usize, usize> = FxHashMap::default();
    let mut cluster_ids = Vec::with_capacity(n);
    let mut next_cluster = 0;

    for i in 0..n {
        let root = uf.find(i);
        let cluster = *root_to_cluster.entry(root).or_insert_with(|| {
            let c = next_cluster;
            next_cluster += 1;
            c
        });
        cluster_ids.push(cluster);
    }

    cluster_ids
}

/// Find optimal eps using cosigt's stabilization detection
/// Tests eps from 0.001 to 0.300, finds where cluster count stabilizes
fn find_optimal_eps(
    dist_matrix: &[Vec<f64>],
    n_paths: usize,
    max_clusters_override: Option<usize>,
) -> f64 {
    if dist_matrix.is_empty() {
        return 0.30;
    }

    // Determine max_clusters: use override if provided, otherwise ~11% of paths
    let max_clusters = max_clusters_override.unwrap_or_else(|| n_paths.div_ceil(9));
    debug!(
        "DBSCAN max_clusters: {} ({})",
        max_clusters,
        if max_clusters_override.is_some() {
            "user override"
        } else {
            "automatic"
        }
    );

    // cosigt: pclust <- length(table(dbscan(distanceMatrix, eps = 0, minPts = 1)$cluster))
    let mut prev_clusters = dbscan_count_clusters(dist_matrix, 0.0);
    debug!("DBSCAN eps scan: eps=0.00 -> {} clusters", prev_clusters);

    // eps from 0.005 to 0.300 in steps of 0.005
    for eps_int in 1..=60 {
        let eps = eps_int as f64 * 0.005;
        let curr_clusters = dbscan_count_clusters(dist_matrix, eps);

        // cosigt: if (abs(pclust - cclust) <= 1)
        let change = (prev_clusters as i64 - curr_clusters as i64).abs();
        debug!(
            "DBSCAN eps scan: eps={:.3} -> {} clusters (change={} from prev={}, max_allowed={})",
            eps, curr_clusters, change, prev_clusters, max_clusters
        );

        // Modified stabilization: also accept when we first reach <= max_clusters
        // This captures cases where cluster count jumps directly to the target
        let first_hit_max = prev_clusters > max_clusters && curr_clusters <= max_clusters;

        if (change <= 1 || first_hit_max) && curr_clusters <= max_clusters {
            if first_hit_max && change > 1 {
                debug!(
                    "DBSCAN: first hit max_clusters at eps {:.3} with {} clusters (jumped from {})",
                    eps, curr_clusters, prev_clusters
                );
            } else {
                debug!(
                    "DBSCAN: stabilized at eps {:.3} with {} clusters (max allowed: {})",
                    eps, curr_clusters, max_clusters
                );
            }
            return eps;
        }
        prev_clusters = curr_clusters;
    }

    // cosigt: return(ifelse(eps < 0.3, optimal_eps, 0.3))
    debug!("DBSCAN: no stabilization found, using fallback eps 0.30");
    0.30
}

/// Compute base-pair weighted Jaccard similarity (matching odgi similarity)
/// For each node: add min(bp_a_on_node, bp_b_on_node) to intersection
/// jaccard = intersection / (bp_a + bp_b - intersection)
fn weighted_jaccard_similarity(
    counts_a: &FxHashMap<u64, u64>, // node_id -> total bp on that node for path a
    counts_b: &FxHashMap<u64, u64>, // node_id -> total bp on that node for path b
    bp_a: u64,                      // total bp in path a
    bp_b: u64,                      // total bp in path b
) -> f64 {
    if bp_a == 0 && bp_b == 0 {
        return 1.0;
    }

    // Compute intersection: sum of min(bp_a_on_node, bp_b_on_node) for all nodes
    let mut intersection: u64 = 0;
    for (&node, &bp_a_on_node) in counts_a {
        if let Some(&bp_b_on_node) = counts_b.get(&node) {
            intersection += bp_a_on_node.min(bp_b_on_node);
        }
    }

    let union = bp_a + bp_b - intersection;
    if union == 0 {
        1.0
    } else {
        intersection as f64 / union as f64
    }
}

/// Compute EDR (estimated difference rate) from Jaccard similarity
/// EDR = (1 - jaccard) / (1 + jaccard)
/// This matches odgi similarity's estimated.difference.rate
fn jaccard_to_edr(jaccard: f64) -> f64 {
    (1.0 - jaccard) / (1.0 + jaccard)
}

/// Cluster paths by EDR (estimated difference rate)
/// Uses base-pair weighted Jaccard similarity like odgi similarity
/// If use_upgma is true, uses pure UPGMA hierarchical clustering with tree cutting
/// Otherwise uses DBSCAN (matching cosigt exactly)
fn cluster_paths_by_similarity(
    paths: &[&GfaPath],
    segment_lengths: &[u64], // segment_id -> length (0-indexed by segment_id - 1)
    threshold: Option<f64>,
    use_all_nodes: bool,
    max_clusters: Option<usize>,
    compute_dendrogram: bool,
    use_upgma: bool,
    upgma_threshold: Option<f64>,
) -> ClusteringResult {
    if paths.is_empty() {
        return ClusteringResult {
            ordering: Vec::new(),
            cluster_ids: Vec::new(),
            num_clusters: 0,
            representatives: Vec::new(),
            cluster_sizes: Vec::new(),
            dendrogram: None,
        };
    }

    let n = paths.len();

    // Build bp-weighted node counts for each path (node_id -> total bp on that node)
    // This matches odgi similarity: for each step, add segment length to that node's count
    let path_bp_counts: Vec<FxHashMap<u64, u64>> = paths
        .par_iter()
        .map(|path| {
            let mut counts: FxHashMap<u64, u64> = FxHashMap::default();
            for step in &path.steps {
                let seg_len = segment_lengths
                    .get(step.segment_id as usize)
                    .copied()
                    .unwrap_or(0);
                *counts.entry(step.segment_id).or_insert(0) += seg_len;
            }
            counts
        })
        .collect();

    // Collect all unique nodes
    let mut all_nodes: FxHashSet<u64> = FxHashSet::default();
    for counts in &path_bp_counts {
        for &node in counts.keys() {
            all_nodes.insert(node);
        }
    }
    let total_unique_nodes = all_nodes.len();

    // Determine which nodes to use for clustering (based on bp variation, not just count)
    let nodes_to_use: FxHashSet<u64> = if use_all_nodes {
        debug!(
            "Clustering mode: --cluster-all-nodes (using all {} nodes)",
            total_unique_nodes
        );
        all_nodes
    } else {
        // Find variable nodes: nodes where bp count varies across paths
        let variable_nodes: FxHashSet<u64> = all_nodes
            .into_iter()
            .filter(|&node| {
                let first_bp = path_bp_counts[0].get(&node).copied().unwrap_or(0);
                path_bp_counts
                    .iter()
                    .skip(1)
                    .any(|counts| counts.get(&node).copied().unwrap_or(0) != first_bp)
            })
            .collect();

        let invariant_nodes = total_unique_nodes - variable_nodes.len();
        debug!("Clustering mode: variable nodes only (using {} of {} nodes, {} invariant nodes excluded)",
               variable_nodes.len(), total_unique_nodes, invariant_nodes);
        variable_nodes
    };

    // Build filtered bp counts (only include nodes_to_use)
    let filtered_bp_counts: Vec<FxHashMap<u64, u64>> = path_bp_counts
        .iter()
        .map(|counts| {
            counts
                .iter()
                .filter(|(node, _)| nodes_to_use.contains(node))
                .map(|(&node, &bp)| (node, bp))
                .collect()
        })
        .collect();

    // Compute total bp for each path
    // When using all nodes, use full path lengths (matching odgi)
    // When using variable nodes only, use filtered lengths (consistent intersection/denominator)
    let total_bp: Vec<u64> = if use_all_nodes {
        path_bp_counts
            .iter()
            .map(|counts| counts.values().sum())
            .collect()
    } else {
        filtered_bp_counts
            .iter()
            .map(|counts| counts.values().sum())
            .collect()
    };

    // Build full pairwise EDR matrix (matching cosigt: uses normalized EDR)
    debug!("Computing {}x{} pairwise EDR matrix", n, n);

    // Compute upper triangle in parallel: EDR for each pair
    let filtered_bp_counts_ref = &filtered_bp_counts;
    let total_bp_ref = &total_bp;
    let pairs: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n)
                .map(move |j| {
                    let jaccard = weighted_jaccard_similarity(
                        &filtered_bp_counts_ref[i],
                        &filtered_bp_counts_ref[j],
                        total_bp_ref[i],
                        total_bp_ref[j],
                    );
                    let edr = jaccard_to_edr(jaccard);
                    (i, j, edr)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Find max EDR for normalization (like cosigt: maxD <- max(regularMatrix))
    let max_edr = pairs.iter().map(|(_, _, edr)| *edr).fold(0.0f64, f64::max);
    debug!("Max EDR: {:.6}", max_edr);

    // Debug: print first few EDR values for comparison with odgi
    for (i, j, edr) in pairs.iter().take(5) {
        let jaccard = weighted_jaccard_similarity(
            &filtered_bp_counts[*i],
            &filtered_bp_counts[*j],
            total_bp[*i],
            total_bp[*j],
        );
        debug!(
            "EDR: {} vs {} = {:.6} (jaccard={:.6}, bp_a={}, bp_b={})",
            paths[*i].name, paths[*j].name, edr, jaccard, total_bp[*i], total_bp[*j]
        );
    }

    // Build normalized distance matrix (like cosigt: normRegularMatrix <- regularMatrix / maxD)
    let mut dist_matrix: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for (i, j, edr) in pairs {
        let norm_edr = if max_edr > 0.0 { edr / max_edr } else { 0.0 };
        dist_matrix[i][j] = norm_edr;
        dist_matrix[j][i] = norm_edr;
    }

    // Log distance distribution
    let mut all_dists: Vec<f64> = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            all_dists.push(dist_matrix[i][j]);
        }
    }
    all_dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if !all_dists.is_empty() {
        debug!(
            "Distance range: {:.3} - {:.3}",
            all_dists[0],
            all_dists[all_dists.len() - 1]
        );
        if all_dists.len() >= 4 {
            let q1 = all_dists[all_dists.len() / 4];
            let median = all_dists[all_dists.len() / 2];
            let q3 = all_dists[3 * all_dists.len() / 4];
            debug!(
                "Distance quartiles: Q1={:.3}, median={:.3}, Q3={:.3}",
                q1, median, q3
            );
        }
    }

    // Get cluster assignments using either UPGMA or DBSCAN
    let (cluster_assignments, dendrogram_for_upgma): (Vec<usize>, Option<Dendrogram>) = if use_upgma
    {
        // Pure UPGMA mode: build dendrogram first, then cut at threshold
        debug!("Using UPGMA hierarchical clustering");
        let dg = build_dendrogram(&dist_matrix, None); // No DBSCAN constraint for pure UPGMA

        // Determine cut threshold
        let cut_threshold = match upgma_threshold {
            Some(t) => {
                debug!("Using user-specified UPGMA threshold: {:.4}", t);
                t * dg.max_height // Scale to actual height range
            }
            None => find_optimal_upgma_threshold(&dg, max_clusters),
        };

        let clusters = cut_dendrogram_at_height(&dg, cut_threshold);
        let num_clusters = clusters.iter().max().map(|&m| m + 1).unwrap_or(1);
        debug!(
            "UPGMA cut at height {:.4} gives {} clusters",
            cut_threshold, num_clusters
        );

        (clusters, Some(dg))
    } else {
        // DBSCAN mode (original behavior)
        // Find optimal eps (or convert user threshold to eps)
        let eps = match threshold {
            Some(t) => {
                let e = 1.0 - t; // Convert similarity threshold to distance eps
                debug!("Using user-specified threshold {:.2} (eps = {:.2})", t, e);
                e
            }
            None => find_optimal_eps(&dist_matrix, n, max_clusters),
        };
        debug!("DBSCAN eps: {:.2}", eps);

        // Run DBSCAN to get cluster assignments
        let clusters = dbscan_cluster(&dist_matrix, eps);
        let num_clusters = clusters.iter().max().map(|&m| m + 1).unwrap_or(1);
        debug!("DBSCAN detected {} clusters", num_clusters);

        (clusters, None)
    };

    let num_clusters = cluster_assignments
        .iter()
        .max()
        .map(|&m| m + 1)
        .unwrap_or(1);

    // Group paths by cluster
    let mut cluster_members: Vec<Vec<usize>> = vec![Vec::new(); num_clusters];
    for (i, &cluster) in cluster_assignments.iter().enumerate() {
        cluster_members[cluster].push(i);
    }

    // Sort clusters by size (largest first) for consistent ordering
    cluster_members.sort_by_key(|v| std::cmp::Reverse(v.len()));

    // Compute medoid for each cluster (path with minimum average distance to others)
    let mut representatives: Vec<usize> = Vec::with_capacity(num_clusters);
    let mut cluster_sizes: Vec<usize> = Vec::with_capacity(num_clusters);

    for members in &cluster_members {
        cluster_sizes.push(members.len());

        if members.len() == 1 {
            // Singleton: the single member is the representative
            representatives.push(members[0]);
        } else {
            // Find medoid: member with minimum average distance to all others
            let mut best_medoid = members[0];
            let mut best_avg_dist = f64::MAX;

            for &candidate in members {
                let sum_dist: f64 = members
                    .iter()
                    .filter(|&&m| m != candidate)
                    .map(|&m| dist_matrix[candidate][m])
                    .sum();
                let avg_dist = sum_dist / (members.len() - 1) as f64;

                if avg_dist < best_avg_dist {
                    best_avg_dist = avg_dist;
                    best_medoid = candidate;
                }
            }
            representatives.push(best_medoid);
        }
    }

    // Build final ordering: within each cluster, order by greedy nearest-neighbor
    let mut ordering = Vec::with_capacity(n);
    let mut final_cluster_ids = Vec::with_capacity(n);

    for (cluster_id, members) in cluster_members.iter().enumerate() {
        if members.is_empty() {
            continue;
        }

        if members.len() == 1 {
            ordering.push(members[0]);
            final_cluster_ids.push(cluster_id);
        } else {
            // Greedy nearest-neighbor within cluster
            let mut placed = vec![false; members.len()];

            // Start with the member that has the most base pairs
            let start = members
                .iter()
                .enumerate()
                .max_by_key(|&(_, &idx)| total_bp[idx])
                .map(|(local_idx, _)| local_idx)
                .unwrap();

            placed[start] = true;
            ordering.push(members[start]);
            final_cluster_ids.push(cluster_id);
            let mut current = start;

            while ordering.len() < ordering.capacity() && placed.iter().filter(|&&p| !p).count() > 0
            {
                // Find nearest unplaced member within this cluster
                let current_global = members[current];
                let mut best_local = None;
                let mut best_dist = f64::MAX;

                for (local_idx, &global_idx) in members.iter().enumerate() {
                    if !placed[local_idx] {
                        let dist = dist_matrix[current_global][global_idx];
                        if dist < best_dist {
                            best_dist = dist;
                            best_local = Some(local_idx);
                        }
                    }
                }

                if let Some(local_idx) = best_local {
                    placed[local_idx] = true;
                    ordering.push(members[local_idx]);
                    final_cluster_ids.push(cluster_id);
                    current = local_idx;
                } else {
                    break;
                }
            }
        }
    }

    debug!(
        "Final ordering: {} paths in {} clusters",
        ordering.len(),
        num_clusters
    );

    // Build or reuse dendrogram
    let dendrogram = if use_upgma {
        // For UPGMA mode, we already have the dendrogram
        dendrogram_for_upgma
    } else if compute_dendrogram {
        // For DBSCAN mode, build dendrogram constrained by clusters
        Some(build_dendrogram(&dist_matrix, Some(&cluster_assignments)))
    } else {
        None
    };

    // If dendrogram is available, use its leaf order for visualization
    let (final_ordering, final_cluster_ids) = if let Some(ref dg) = dendrogram {
        // Map dendrogram leaf order to cluster IDs
        let mut dg_ordering = Vec::with_capacity(n);
        let mut dg_cluster_ids = Vec::with_capacity(n);

        for &orig_idx in &dg.leaf_order {
            dg_ordering.push(orig_idx);
            // Find cluster ID for this path from cluster assignments
            dg_cluster_ids.push(cluster_assignments[orig_idx]);
        }
        (dg_ordering, dg_cluster_ids)
    } else {
        (ordering, final_cluster_ids)
    };

    ClusteringResult {
        ordering: final_ordering,
        cluster_ids: final_cluster_ids,
        num_clusters,
        representatives,
        cluster_sizes,
        dendrogram,
    }
}

#[derive(Default, Clone)]
struct BinInfo {
    mean_depth: f64,
    mean_inv: f64,
    mean_pos: f64,      // mean position within path (for darkness gradient)
    mean_uncalled: f64, // proportion of uncalled bases (N's) in bin
    highlighted: bool,  // whether this bin contains highlighted nodes
}

/// Draw a line on the buffer (Bresenham's algorithm)
fn draw_line(
    buffer: &mut [u8],
    width: u32,
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,
    r: u8,
    g: u8,
    b: u8,
) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    let mut x = x0;
    let mut y = y0;

    loop {
        if x >= 0 && y >= 0 {
            let idx = ((y as u32 * width + x as u32) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = r;
                buffer[idx + 1] = g;
                buffer[idx + 2] = b;
                buffer[idx + 3] = 255;
            }
        }
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            err += dx;
            y += sy;
        }
    }
}

/// Render dendrogram recursively, returning the Y coordinate of the node's connection point
fn render_dendrogram_node(
    buffer: &mut [u8],
    width: u32,
    dendrogram: &Dendrogram,
    node_idx: usize,
    n_leaves: usize,
    x_offset: u32,
    dendro_width: u32,
    pix_per_path: u32,
    leaf_y_positions: &[u32],
) -> (u32, f64) {
    // Returns (y_position, height)
    if node_idx < n_leaves {
        // Leaf node - return its Y position (center of the path row)
        let y = leaf_y_positions[node_idx] + pix_per_path / 2;
        return (y, 0.0);
    }

    let internal_idx = node_idx - n_leaves;
    if internal_idx >= dendrogram.nodes.len() {
        // Safety check: invalid node index
        return (0, 0.0);
    }
    let node = &dendrogram.nodes[internal_idx];

    // Recursively get positions of children
    let (left_y, _left_h) = render_dendrogram_node(
        buffer,
        width,
        dendrogram,
        node.left,
        n_leaves,
        x_offset,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
    );
    let (right_y, _right_h) = render_dendrogram_node(
        buffer,
        width,
        dendrogram,
        node.right,
        n_leaves,
        x_offset,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
    );

    // Calculate X position based on merge height
    // X goes from right (leaves at dendro_width) to left (root at 0)
    let x = if dendrogram.max_height > 0.0 {
        x_offset + ((1.0 - node.height / dendrogram.max_height) * (dendro_width - 5) as f64) as u32
    } else {
        x_offset + dendro_width / 2
    };

    // Draw horizontal lines from children to this node's X
    let left_x = if node.left < n_leaves {
        x_offset + dendro_width - 2 // Leaves are at the right edge
    } else {
        let left_internal = node.left - n_leaves;
        if left_internal < dendrogram.nodes.len() {
            let left_node = &dendrogram.nodes[left_internal];
            x_offset
                + ((1.0 - left_node.height / dendrogram.max_height) * (dendro_width - 5) as f64)
                    as u32
        } else {
            x_offset + dendro_width / 2
        }
    };
    let right_x = if node.right < n_leaves {
        x_offset + dendro_width - 2
    } else {
        let right_internal = node.right - n_leaves;
        if right_internal < dendrogram.nodes.len() {
            let right_node = &dendrogram.nodes[right_internal];
            x_offset
                + ((1.0 - right_node.height / dendrogram.max_height) * (dendro_width - 5) as f64)
                    as u32
        } else {
            x_offset + dendro_width / 2
        }
    };

    // Draw lines (dark grey color)
    let line_color = (80u8, 80u8, 80u8);

    // Horizontal line from left child to this X
    draw_line(
        buffer,
        width,
        left_x as i32,
        left_y as i32,
        x as i32,
        left_y as i32,
        line_color.0,
        line_color.1,
        line_color.2,
    );
    // Horizontal line from right child to this X
    draw_line(
        buffer,
        width,
        right_x as i32,
        right_y as i32,
        x as i32,
        right_y as i32,
        line_color.0,
        line_color.1,
        line_color.2,
    );
    // Vertical line connecting the two horizontal lines
    draw_line(
        buffer,
        width,
        x as i32,
        left_y as i32,
        x as i32,
        right_y as i32,
        line_color.0,
        line_color.1,
        line_color.2,
    );

    // Return the midpoint Y
    let mid_y = (left_y + right_y) / 2;
    (mid_y, node.height)
}

/// Render the full dendrogram for PNG output
/// leaf_y_positions: pre-computed Y positions for each leaf (indexed by original path index)
fn render_dendrogram_png(
    buffer: &mut [u8],
    width: u32,
    dendrogram: &Dendrogram,
    dendro_width: u32,
    pix_per_path: u32,
    leaf_y_positions: &[u32],
) {
    if dendrogram.nodes.is_empty() || dendrogram.leaf_order.len() <= 1 {
        return;
    }

    // Number of leaves from the dendrogram
    let n_leaves = dendrogram.leaf_order.len();

    // Root is the last internal node
    let root_idx = n_leaves + dendrogram.nodes.len() - 1;

    render_dendrogram_node(
        buffer,
        width,
        dendrogram,
        root_idx,
        n_leaves,
        0,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
    );
}

/// Render dendrogram node recursively for SVG, returning Y position and collecting path elements
fn render_dendrogram_node_svg(
    dendrogram: &Dendrogram,
    node_idx: usize,
    n_leaves: usize,
    x_offset: f64,
    dendro_width: f64,
    pix_per_path: f64,
    leaf_y_positions: &[f64],
    paths: &mut Vec<String>,
) -> (f64, f64) {
    // Returns (y_position, height)
    if node_idx < n_leaves {
        // Leaf node - return its Y position (center of the path row)
        let y = leaf_y_positions[node_idx] + pix_per_path / 2.0;
        return (y, 0.0);
    }

    let internal_idx = node_idx - n_leaves;
    if internal_idx >= dendrogram.nodes.len() {
        return (0.0, 0.0);
    }
    let node = &dendrogram.nodes[internal_idx];

    // Recursively get positions of children
    let (left_y, _) = render_dendrogram_node_svg(
        dendrogram,
        node.left,
        n_leaves,
        x_offset,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
        paths,
    );
    let (right_y, _) = render_dendrogram_node_svg(
        dendrogram,
        node.right,
        n_leaves,
        x_offset,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
        paths,
    );

    // Calculate X position based on merge height
    let x = if dendrogram.max_height > 0.0 {
        x_offset + (1.0 - node.height / dendrogram.max_height) * (dendro_width - 5.0)
    } else {
        x_offset + dendro_width / 2.0
    };

    // Calculate child X positions
    let left_x = if node.left < n_leaves {
        x_offset + dendro_width - 2.0
    } else {
        let left_internal = node.left - n_leaves;
        if left_internal < dendrogram.nodes.len() {
            let left_node = &dendrogram.nodes[left_internal];
            x_offset + (1.0 - left_node.height / dendrogram.max_height) * (dendro_width - 5.0)
        } else {
            x_offset + dendro_width / 2.0
        }
    };
    let right_x = if node.right < n_leaves {
        x_offset + dendro_width - 2.0
    } else {
        let right_internal = node.right - n_leaves;
        if right_internal < dendrogram.nodes.len() {
            let right_node = &dendrogram.nodes[right_internal];
            x_offset + (1.0 - right_node.height / dendrogram.max_height) * (dendro_width - 5.0)
        } else {
            x_offset + dendro_width / 2.0
        }
    };

    // Add SVG path elements for the lines
    let line_color = "#505050";
    // Horizontal line from left child to this X
    paths.push(format!(
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1"/>"#,
        left_x, left_y, x, left_y, line_color
    ));
    // Horizontal line from right child to this X
    paths.push(format!(
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1"/>"#,
        right_x, right_y, x, right_y, line_color
    ));
    // Vertical line connecting the two horizontal lines
    paths.push(format!(
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1"/>"#,
        x, left_y, x, right_y, line_color
    ));

    let mid_y = (left_y + right_y) / 2.0;
    (mid_y, node.height)
}

/// Render the full dendrogram for SVG output, returns SVG elements as string
/// leaf_y_positions: pre-computed Y positions for each leaf (indexed by original path index)
fn render_dendrogram_svg(
    dendrogram: &Dendrogram,
    dendro_width: f64,
    pix_per_path: f64,
    leaf_y_positions: &[f64],
) -> String {
    if dendrogram.nodes.is_empty() || dendrogram.leaf_order.len() <= 1 {
        return String::new();
    }

    let n_leaves = dendrogram.leaf_order.len();
    let root_idx = n_leaves + dendrogram.nodes.len() - 1;
    let mut paths = Vec::new();

    render_dendrogram_node_svg(
        dendrogram,
        root_idx,
        n_leaves,
        0.0,
        dendro_width,
        pix_per_path,
        leaf_y_positions,
        &mut paths,
    );

    paths.join("\n")
}

/// Render annotation legend at the top of the image (PNG)
fn render_annotation_legend_png(
    buffer: &mut [u8],
    width: u32,
    _left_margin: u32,
    categories: &[String],
    category_colors: &FxHashMap<String, (u8, u8, u8)>,
    legend_height: u32,
    char_size: u32,
) {
    let swatch_size = 12u32;
    let swatch_padding = 8u32;
    let text_padding = 4u32;
    let item_spacing = 12u32;

    // Calculate available width (reserve space for "+N" indicator)
    let available_width = width.saturating_sub(swatch_padding * 2 + 50);

    // Calculate width needed for each category
    let category_widths: Vec<u32> = categories
        .iter()
        .map(|cat| swatch_size + text_padding + (cat.len() as u32 * char_size) + item_spacing)
        .collect();

    // Determine how many categories fit
    let mut total_items_width = 0u32;
    let mut visible_count = 0usize;
    for w in &category_widths {
        if total_items_width + w <= available_width {
            total_items_width += w;
            visible_count += 1;
        } else {
            break;
        }
    }

    // If none fit, show at least one truncated
    if visible_count == 0 && !categories.is_empty() {
        visible_count = 1;
        total_items_width = category_widths.first().copied().unwrap_or(0);
    }

    // Calculate starting x position to center the legend
    let hidden_count = categories.len().saturating_sub(visible_count);
    let indicator_width = if hidden_count > 0 {
        let indicator = format!("+{}", hidden_count);
        indicator.len() as u32 * char_size + item_spacing
    } else {
        0
    };
    let total_legend_width = total_items_width + indicator_width;
    let x_start = (width.saturating_sub(total_legend_width)) / 2;

    let mut x_pos = x_start;
    let y_center = legend_height / 2;
    let swatch_y = y_center.saturating_sub(swatch_size / 2);

    for category in categories.iter().take(visible_count) {
        if let Some(&(r, g, b)) = category_colors.get(category) {
            // Draw color swatch
            for sx in 0..swatch_size {
                for sy in 0..swatch_size {
                    let px = x_pos + sx;
                    let py = swatch_y + sy;
                    if px < width {
                        let idx = ((py * width + px) * 4) as usize;
                        if idx + 3 < buffer.len() {
                            buffer[idx] = r;
                            buffer[idx + 1] = g;
                            buffer[idx + 2] = b;
                            buffer[idx + 3] = 255;
                        }
                    }
                }
            }

            // Draw category label
            let text_x = x_pos + swatch_size + text_padding;
            let text_y = y_center.saturating_sub(char_size / 2);
            for (i, c) in category.chars().enumerate() {
                let char_x = text_x + (i as u32) * char_size;
                if char_x + char_size > width {
                    break;
                }
                let c_byte = c as usize;
                let char_data = if c_byte < 128 {
                    &FONT_5X8[c_byte]
                } else {
                    &FONT_5X8[b'?' as usize]
                };
                write_char(buffer, width, char_x, text_y, char_data, char_size, 0, 0, 0);
            }

            // Move to next item
            x_pos += swatch_size + text_padding + (category.len() as u32 * char_size) + item_spacing;
        }
    }

    // Draw "+N" indicator if there are hidden categories
    let hidden_count = categories.len().saturating_sub(visible_count);
    if hidden_count > 0 {
        let indicator = format!("+{}", hidden_count);
        let text_y = y_center.saturating_sub(char_size / 2);
        for (i, c) in indicator.chars().enumerate() {
            let char_x = x_pos + (i as u32) * char_size;
            if char_x + char_size <= width {
                let c_byte = c as usize;
                let char_data = if c_byte < 128 {
                    &FONT_5X8[c_byte]
                } else {
                    &FONT_5X8[b'?' as usize]
                };
                write_char(buffer, width, char_x, text_y, char_data, char_size, 128, 128, 128);
            }
        }
    }
}

/// Escape XML special characters for SVG text
fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Render annotation legend at the top of the image (SVG)
fn render_annotation_legend_svg(
    categories: &[String],
    category_colors: &FxHashMap<String, (u8, u8, u8)>,
    image_width: f64,
    legend_height: f64,
    font_size: f64,
) -> String {
    let mut svg = String::new();

    let swatch_size = 12.0;
    let text_padding = 4.0;
    let item_spacing = 16.0;

    // Calculate total legend width for centering
    let total_legend_width: f64 = categories
        .iter()
        .filter_map(|cat| category_colors.get(cat).map(|_| cat))
        .map(|cat| {
            let text_width = cat.len() as f64 * font_size * 0.6;
            swatch_size + text_padding + text_width + item_spacing
        })
        .sum();

    // Center the legend
    let x_start = (image_width - total_legend_width).max(0.0) / 2.0;

    let mut x_pos = x_start;
    let y_center = legend_height / 2.0;
    let swatch_y = y_center - swatch_size / 2.0;

    for category in categories {
        if let Some(&(r, g, b)) = category_colors.get(category) {
            // Draw color swatch
            svg.push_str(&format!(
                r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                x_pos, swatch_y, swatch_size, swatch_size, r, g, b
            ));
            svg.push('\n');

            // Draw category label
            let text_x = x_pos + swatch_size + text_padding;
            let text_y = y_center + font_size / 3.0;
            svg.push_str(&format!(
                r#"<text x="{}" y="{}" font-family="'DejaVu Sans Mono', 'Courier New', monospace" font-size="{}" fill="black">{}</text>"#,
                text_x, text_y, font_size, escape_xml(category)
            ));
            svg.push('\n');

            // Estimate text width (approximate: 0.6 * font_size per character)
            let text_width = category.len() as f64 * font_size * 0.6;
            x_pos += swatch_size + text_padding + text_width + item_spacing;
        }
    }

    svg
}

fn write_char(
    buffer: &mut [u8],
    width: u32,
    base_x: u32,
    base_y: u32,
    char_data: &[u8; 8],
    char_size: u32,
    r: u8,
    g: u8,
    b: u8,
) {
    let ratio = char_size / 8;
    for j in 0..8u32 {
        let row = char_data[j as usize];
        let y = base_y + j * ratio;
        for z in (0..8i32).rev() {
            if (row >> z) & 1 == 1 {
                let x = base_x + (7 - z as u32) * ratio;
                for rx in 0..ratio {
                    for ry in 0..ratio {
                        let px = x + rx;
                        let py = y + ry;
                        let idx = ((py * width + px) * 4) as usize;
                        if idx + 3 < buffer.len() {
                            buffer[idx] = r;
                            buffer[idx + 1] = g;
                            buffer[idx + 2] = b;
                            buffer[idx + 3] = 255;
                        }
                    }
                }
            }
        }
    }
}

fn add_path_step(
    buffer: &mut [u8],
    width: u32,
    x: u32,
    y_start: u32,
    pix_per_path: u32,
    r: u8,
    g: u8,
    b: u8,
    no_path_borders: bool,
    black_border: bool,
) {
    let t = y_start;
    if no_path_borders || pix_per_path < 3 {
        let s = t + pix_per_path;
        for y in t..s {
            let idx = ((y * width + x) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = r;
                buffer[idx + 1] = g;
                buffer[idx + 2] = b;
                buffer[idx + 3] = 255;
            }
        }
    } else {
        let s = t + pix_per_path - 1;
        for y in t..s {
            let idx = ((y * width + x) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = r;
                buffer[idx + 1] = g;
                buffer[idx + 2] = b;
                buffer[idx + 3] = 255;
            }
        }
        if black_border {
            let idx = ((s * width + x) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = 0;
                buffer[idx + 1] = 0;
                buffer[idx + 2] = 0;
                buffer[idx + 3] = 255;
            }
        }
    }
}

/// Add a point to the edge visualization area
fn add_edge_point(buffer: &mut [u8], width: u32, x: u32, y: u32, path_space: u32, rgb: u8) {
    let actual_y = y + path_space;
    let idx = ((actual_y * width + x) * 4) as usize;
    if idx + 3 < buffer.len() {
        buffer[idx] = rgb;
        buffer[idx + 1] = rgb;
        buffer[idx + 2] = rgb;
        buffer[idx + 3] = 255;
    }
}

/// Get color for depth using colorbrewer palette (with optional grey for low coverage)
fn get_depth_color(
    mean_depth: f64,
    no_grey_depth: bool,
    palette: Option<&[(u8, u8, u8)]>,
) -> (u8, u8, u8) {
    // Use custom palette if provided, otherwise use Spectral
    if let Some(pal) = palette {
        // Custom palette: distribute colors evenly across depth range
        let n = pal.len();
        if n == 0 {
            return (128, 128, 128);
        }
        // Map mean_depth to palette index (1-based depths, so depth 1 = index 0)
        let idx = if no_grey_depth {
            // Use full palette range
            ((mean_depth - 1.0).max(0.0) / (n as f64)).floor() as usize
        } else {
            // First color for < 0.5x, second for < 1.5x, rest distributed
            if mean_depth < 0.5 {
                return (196, 196, 196); // Grey for very low
            } else if mean_depth < 1.5 {
                return (128, 128, 128); // Grey for low
            }
            ((mean_depth - 1.5) / (n as f64)).floor() as usize
        };
        pal[idx.min(n - 1)]
    } else if no_grey_depth {
        // Use full Spectral 11 range for all depths (skip grey colors at indices 0-1)
        let cuts = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5];
        for (i, &cut) in cuts.iter().enumerate() {
            if mean_depth <= cut {
                return COLORBREWER_SPECTRAL_13[i + 2];
            }
        }
        COLORBREWER_SPECTRAL_13[12]
    } else {
        // Default: use grey for low coverage (indices 0-1), Spectral for rest
        let cuts = [
            0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5,
        ];
        for (i, &cut) in cuts.iter().enumerate() {
            if mean_depth <= cut {
                return COLORBREWER_SPECTRAL_13[i];
            }
        }
        COLORBREWER_SPECTRAL_13[12]
    }
}

fn render(args: &Args, graph: &Graph) -> Vec<u8> {
    // Check for conflicting options
    if args.cluster_paths && args.prefix_merges.is_some() {
        eprintln!("[gfalook] error: -k/--cluster-paths cannot be used with -M/--prefix-merges.");
        std::process::exit(1);
    }
    // Note: compressed_mode conflicts with cluster_paths and prefix_merges are handled by clap

    let mut display_paths: Vec<&GfaPath> = graph.paths.iter().collect();

    if let Some(ref prefix) = args.ignore_prefix {
        display_paths.retain(|p| !p.name.starts_with(prefix));
    }

    if let Some(ref ptd_file) = args.paths_to_display {
        if let Ok(ptd) = load_paths_to_display(ptd_file) {
            let ptd_set: std::collections::HashSet<_> = ptd.iter().collect();
            display_paths.retain(|p| ptd_set.contains(&p.name));
            let path_map: FxHashMap<&String, &GfaPath> =
                display_paths.iter().map(|p| (&p.name, *p)).collect();
            display_paths = ptd
                .iter()
                .filter_map(|name| path_map.get(name).copied())
                .collect();
        }
    }

    let pix_per_path = args.path_height;
    let bottom_padding = 5u32;

    let len_to_visualize = graph.total_length;
    let viz_width = args.width.min(len_to_visualize as u32);

    let bin_width = args
        .bin_width
        .unwrap_or_else(|| len_to_visualize as f64 / viz_width as f64);
    let _scale_x = 1.0; // In binned mode
    let _scale_y = viz_width as f64 / len_to_visualize as f64;

    // Cluster paths by similarity if requested (PNG rendering)
    let cluster_result = if args.cluster_paths {
        debug!(
            "Clustering {} paths by EDR (estimated difference rate)",
            display_paths.len()
        );
        // Build segment lengths vector for EDR computation
        let segment_lengths: Vec<u64> = graph.segments.iter().map(|s| s.sequence_len).collect();
        let original_paths = display_paths.clone(); // Save for medoids TSV
        let result = cluster_paths_by_similarity(
            &display_paths,
            &segment_lengths,
            args.cluster_threshold,
            args.cluster_all_nodes,
            args.max_clusters,
            args.dendrogram || args.use_upgma,
            args.use_upgma,
            args.upgma_threshold,
        );
        display_paths = result.ordering.iter().map(|&i| display_paths[i]).collect();
        // Write cluster assignments to TSV
        write_cluster_tsv(&args.out, &display_paths, &result);
        // Write medoids TSV
        write_medoids_tsv(&args.out, &original_paths, &result);

        // Filter to representatives only if requested (PNG)
        let result = if args.cluster_representatives {
            let rep_set: FxHashSet<usize> = result.representatives.iter().copied().collect();
            let mut filtered_paths = Vec::new();
            let mut filtered_cluster_ids = Vec::new();
            for (pos, &orig_idx) in result.ordering.iter().enumerate() {
                if rep_set.contains(&orig_idx) {
                    filtered_paths.push(display_paths[pos]);
                    filtered_cluster_ids.push(result.cluster_ids[pos]);
                }
            }
            display_paths = filtered_paths;
            ClusteringResult {
                ordering: (0..display_paths.len()).collect(),
                cluster_ids: filtered_cluster_ids,
                num_clusters: result.num_clusters,
                representatives: result.representatives,
                cluster_sizes: result.cluster_sizes,
                dendrogram: result.dendrogram,
            }
        } else {
            result
        };
        Some(result)
    } else {
        None
    };

    // Recalculate path_count after potential filtering by cluster_representatives (PNG)
    let path_count = display_paths.len() as u32;

    // Calculate total gap space needed for cluster separators
    let total_gap = if let Some(ref cr) = cluster_result {
        (cr.num_clusters.saturating_sub(1) as u32) * args.cluster_gap
    } else {
        0
    };

    // Load prefix grouping if specified (PNG) - must be after clustering check
    let path_grouping: Option<PathGrouping> = args.prefix_merges.as_ref().and_then(|p| {
        let paths_vec: Vec<GfaPath> = display_paths.iter().map(|&p| p.clone()).collect();
        match load_prefix_merges(p, &paths_vec) {
            Ok(grouping) => {
                info!(
                    "Read {} valid prefixes for {} groups",
                    grouping.prefixes.len(),
                    grouping.num_groups
                );
                Some(grouping)
            }
            Err(e) => {
                eprintln!("[gfalook] warning: failed to load prefix merges: {}", e);
                None
            }
        }
    });

    // Load annotations if specified
    let annotations: Option<AnnotationData> = args.annotation_file.as_ref().and_then(|p| {
        match load_annotations(p, args.annotation_column) {
            Ok(ann) => {
                info!(
                    "Loaded {} prefixes across {} categories",
                    ann.prefixes.len(),
                    ann.categories.len()
                );
                Some(ann)
            }
            Err(e) => {
                eprintln!("[gfalook] warning: failed to load annotations: {}", e);
                None
            }
        }
    });

    // Effective row count: use num_groups if grouping is enabled, 1 if compressed mode
    let effective_row_count = if args.compressed_mode {
        1
    } else if let Some(ref pg) = path_grouping {
        pg.num_groups as u32
    } else {
        path_count
    };

    debug!("Binned mode");
    debug!("bin width: {:.2e}", bin_width);
    debug!("image width: {}", viz_width);

    // Use prefix names for max_name_len when grouping, "COMPRESSED_MODE" for compressed mode
    let max_name_len = if args.compressed_mode {
        "COMPRESSED_MODE".len()
    } else if let Some(ref pg) = path_grouping {
        pg.prefixes.iter().map(|p| p.len()).max().unwrap_or(10)
    } else if args.cluster_representatives {
        // Account for " (n=X)" suffix: max cluster size determines suffix length
        let max_size = cluster_result
            .as_ref()
            .map(|cr| cr.cluster_sizes.iter().max().copied().unwrap_or(1))
            .unwrap_or(1);
        let suffix_len = format!(" (n={})", max_size).len();
        display_paths
            .iter()
            .map(|p| p.name.len() + suffix_len)
            .max()
            .unwrap_or(10)
    } else {
        display_paths
            .iter()
            .map(|p| p.name.len())
            .max()
            .unwrap_or(10)
    };
    let max_num_of_chars = args.max_num_of_characters.unwrap_or(max_name_len.min(128));
    let char_size = ((pix_per_path / 8) * 8).clamp(8, 64);

    // Cluster bar width (only if clustering is enabled)
    let cluster_bar_width: u32 = if cluster_result.is_some() { 10 } else { 0 };

    // Annotation bar width (only if annotations are loaded)
    let annotation_bar_width: u32 = if annotations.is_some() {
        args.annotation_bar_width
    } else {
        0
    };

    // Gap between cluster bar and annotation bar when both are present
    let bar_gap: u32 = if cluster_result.is_some() && annotations.is_some() { 4 } else { 0 };

    // Legend height (only if annotations are loaded)
    let legend_height: u32 = if annotations.is_some() {
        args.legend_height
    } else {
        0
    };

    // Dendrogram width (only if dendrogram is enabled and we have a dendrogram)
    let dendrogram_width: u32 = if args.dendrogram
        && cluster_result
            .as_ref()
            .map(|cr| cr.dendrogram.is_some())
            .unwrap_or(false)
    {
        args.dendrogram_width
    } else {
        0
    };

    // Disable path names when pack_paths is enabled (they wouldn't make sense)
    let text_only_width = if args.hide_path_names || args.pack_paths {
        0u32
    } else if pix_per_path >= 8 {
        (max_num_of_chars as u32 * char_size) + char_size / 2
    } else {
        0
    };

    // Total left panel width includes dendrogram + cluster bar + gap + annotation bar + path names
    let path_names_width = dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + text_only_width;

    let path_space = effective_row_count * pix_per_path + total_gap;

    // Load colorbrewer palette if specified
    let depth_palette: Option<&[(u8, u8, u8)]> = args.colorbrewer_palette.as_ref().and_then(|arg| {
        if let Some((scheme, _n)) = parse_colorbrewer_arg(arg) {
            get_colorbrewer_palette(&scheme).or_else(|| {
                eprintln!("[gfalook] warning: unknown colorbrewer palette '{}', using default Spectral", scheme);
                None
            })
        } else {
            eprintln!("[gfalook] warning: invalid colorbrewer palette format '{}', expected SCHEME:N", arg);
            None
        }
    });

    // Height for edge visualization area - matches odgi's calculation
    // height = min(len_to_visualize, args.height + bottom_padding)
    // scale_y = height / len_to_visualize
    let height_param = (args.height + bottom_padding) as u64;
    let edge_height = (len_to_visualize.min(height_param)) as u32;
    let scale_y_edges = edge_height as f64 / len_to_visualize as f64;

    let total_width = viz_width + path_names_width;
    // Calculate max axis height for buffer allocation (16 pixels when enabled)
    let max_axis_height: u32 = if args.x_axis.is_some() { 16 } else { 0 };
    // Initial height - will be cropped later based on actual edge rendering (includes legend at top)
    let max_possible_height = legend_height + path_space + max_axis_height + edge_height;

    let mut buffer = vec![255u8; (total_width * max_possible_height * 4) as usize];
    let mut path_names_buffer = if path_names_width > 0 {
        vec![255u8; (path_names_width * max_possible_height * 4) as usize]
    } else {
        Vec::new()
    };

    // Pre-compute leaf Y positions for dendrogram (accounting for cluster gaps and legend)
    let dendrogram_leaf_y_positions: Vec<u32> = if dendrogram_width > 0 {
        if let Some(ref cr) = cluster_result {
            if let Some(ref dg) = cr.dendrogram {
                let n_leaves = dg.leaf_order.len();
                let mut positions = vec![0u32; n_leaves];
                let mut cumulative_gap: u32 = 0;
                let mut prev_cluster_id: Option<usize> = None;

                for (display_pos, &orig_idx) in dg.leaf_order.iter().enumerate() {
                    if orig_idx < n_leaves && display_pos < cr.cluster_ids.len() {
                        // cluster_ids is indexed by display position, not original index
                        let cluster_id = cr.cluster_ids[display_pos];
                        if prev_cluster_id.is_some_and(|prev| prev != cluster_id) {
                            cumulative_gap += args.cluster_gap;
                        }
                        prev_cluster_id = Some(cluster_id);
                        positions[orig_idx] =
                            legend_height + display_pos as u32 * pix_per_path + cumulative_gap;
                    }
                }
                positions
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // Filter annotation categories to only those used by paths in the graph (for legend)
    let filtered_categories: Vec<String> = if let Some(ref ann) = annotations {
        let used_categories: std::collections::HashSet<&String> = display_paths
            .iter()
            .filter_map(|p| ann.get_annotation(&p.name))
            .collect();
        ann.categories
            .iter()
            .filter(|c| used_categories.contains(c))
            .cloned()
            .collect()
    } else {
        Vec::new()
    };

    // Render dendrogram if enabled (PNG)
    if dendrogram_width > 0 && !dendrogram_leaf_y_positions.is_empty() {
        if let Some(ref cr) = cluster_result {
            if let Some(ref dg) = cr.dendrogram {
                render_dendrogram_png(
                    &mut path_names_buffer,
                    path_names_width,
                    dg,
                    dendrogram_width,
                    pix_per_path,
                    &dendrogram_leaf_y_positions,
                );
            }
        }
    }

    // Track maximum y coordinate used (for cropping)
    let mut max_y: u32 = legend_height + path_space + max_axis_height;

    let custom_colors: Option<FxHashMap<String, (u8, u8, u8)>> = args
        .path_colors
        .as_ref()
        .and_then(|p| load_path_colors(p).ok());

    // Load highlight node IDs if specified
    let highlight_nodes: Option<FxHashSet<u64>> = args
        .highlight_node_ids
        .as_ref()
        .and_then(|p| load_highlight_node_ids(p).ok());

    // Track which groups have already been rendered (for path names)
    let mut rendered_groups: FxHashSet<i64> = FxHashSet::default();

    // Calculate max path length for longest-path option
    let max_path_length: u64 = if args.longest_path || args.change_darkness {
        display_paths
            .iter()
            .map(|path| {
                path.steps
                    .iter()
                    .map(|step| {
                        let seg_id = step.segment_id as usize;
                        if seg_id < graph.segments.len() {
                            graph.segments[seg_id].sequence_len
                        } else {
                            0
                        }
                    })
                    .sum::<u64>()
            })
            .max()
            .unwrap_or(1)
    } else {
        1
    };

    // Compressed mode: aggregate bins across all paths and render single row (PNG)
    if args.compressed_mode {
        // Use RdBu palette by default for compressed mode, or user-specified palette
        let compressed_palette = depth_palette.unwrap_or(COLORBREWER_RDBU_11.as_slice());

        // Aggregate bins across all paths
        let mut aggregated_bins: FxHashMap<usize, (f64, u32)> = FxHashMap::default(); // (sum_depth, count)

        for path in display_paths.iter() {
            for step in &path.steps {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    let offset = graph.segment_offsets[seg_id];
                    let seg_len = graph.segments[seg_id].sequence_len;

                    for k in 0..seg_len {
                        let pos = offset + k;
                        let curr_bin = (pos as f64 / bin_width) as usize;
                        let entry = aggregated_bins.entry(curr_bin).or_insert((0.0, 0));
                        entry.0 += 1.0; // Add depth contribution
                        entry.1 += 1; // Count paths covering this position
                    }
                }
            }
        }

        // Normalize aggregated bins to get mean depth across all paths
        let num_paths = display_paths.len() as f64;
        let mut compressed_bins: FxHashMap<usize, f64> = FxHashMap::default();
        for (bin_idx, (sum_depth, _count)) in aggregated_bins.iter() {
            // Normalize: divide by bin_width to get depth, then by num_paths for mean
            let mean_depth = sum_depth / bin_width / num_paths;
            compressed_bins.insert(*bin_idx, mean_depth);
        }

        // Render path name "COMPRESSED_MODE"
        let y_start = legend_height;
        if text_only_width > 0 && pix_per_path >= 8 {
            let display_name = "COMPRESSED_MODE";
            let num_of_chars = display_name.len().min(max_num_of_chars);
            let left_padding = max_num_of_chars - num_of_chars;

            let base_y = y_start + pix_per_path / 2 - char_size / 2;
            for (i, c) in display_name.chars().take(num_of_chars).enumerate() {
                let base_x = (left_padding + i) as u32 * char_size
                    + 3
                    + dendrogram_width
                    + cluster_bar_width
                    + annotation_bar_width;
                let c_byte = c as usize;
                let char_data = if c_byte < 128 {
                    &FONT_5X8[c_byte]
                } else {
                    &FONT_5X8[b'?' as usize]
                };
                write_char(
                    &mut path_names_buffer,
                    path_names_width,
                    base_x,
                    base_y,
                    char_data,
                    char_size,
                    0,
                    0,
                    0,
                );
            }
        }

        // Render aggregated bins (PNG compressed mode)
        for (bin_idx, mean_depth) in &compressed_bins {
            let x = (*bin_idx as u32).min(viz_width - 1);
            let (r, g, b) =
                get_depth_color(*mean_depth, args.no_grey_depth, Some(compressed_palette));
            add_path_step(
                &mut buffer,
                total_width,
                x + path_names_width,
                y_start,
                pix_per_path,
                r,
                g,
                b,
                args.no_path_borders,
                args.black_path_borders,
            );
        }
    }

    // Pack-paths mode: use 2D collision detection to pack paths compactly (PNG)
    if args.pack_paths && !args.compressed_mode {
        // Pre-compute bins for all paths to determine their X ranges
        struct PathBinData {
            min_bin: usize,
            max_bin: usize,
            bins: FxHashMap<usize, BinInfo>,
            color: (u8, u8, u8),
        }

        let mut path_data: Vec<PathBinData> = Vec::with_capacity(display_paths.len());

        for path in display_paths.iter() {
            let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();
            let mut min_bin = usize::MAX;
            let mut max_bin = 0usize;

            let mut path_pos: u64 = 0;
            for step in &path.steps {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    let offset = graph.segment_offsets[seg_id];
                    let seg_len = graph.segments[seg_id].sequence_len;
                    let n_count = graph.segments[seg_id].n_count;
                    let n_proportion = if seg_len > 0 {
                        n_count as f64 / seg_len as f64
                    } else {
                        0.0
                    };
                    let is_highlighted = highlight_nodes
                        .as_ref()
                        .is_some_and(|hn| hn.contains(&step.segment_id));

                    for k in 0..seg_len {
                        let pos = offset + k;
                        let curr_bin = (pos as f64 / bin_width) as usize;
                        min_bin = min_bin.min(curr_bin);
                        max_bin = max_bin.max(curr_bin);

                        let entry = bins.entry(curr_bin).or_default();
                        entry.mean_depth += 1.0;
                        if step.is_reverse {
                            entry.mean_inv += 1.0;
                        }
                        entry.mean_pos += path_pos as f64;
                        entry.mean_uncalled += n_proportion;
                        if is_highlighted {
                            entry.highlighted = true;
                        }
                        path_pos += 1;
                    }
                }
            }

            // Normalize bins
            for (_, v) in bins.iter_mut() {
                if v.mean_depth > 0.0 {
                    v.mean_pos /= v.mean_depth;
                    v.mean_uncalled /= v.mean_depth;
                }
                v.mean_inv /= if v.mean_depth > 0.0 {
                    v.mean_depth
                } else {
                    1.0
                };
                v.mean_depth /= bin_width;
            }

            let color = if let Some(ref colors) = custom_colors {
                colors.get(&path.name).copied().unwrap_or((200, 200, 200)) // Light grey for non-specified paths
            } else {
                compute_path_color(&path.name, args.color_by_prefix)
            };

            path_data.push(PathBinData {
                min_bin,
                max_bin,
                bins,
                color,
            });
        }

        // Layout buffer: for each X position, track the lowest available row
        // Using a simple greedy approach: for each path, find first row where it fits
        let mut occupancy: Vec<Vec<bool>> = vec![vec![false; viz_width as usize]; 1]; // Start with 1 row

        // Assign Y row to each path
        let mut path_rows: Vec<usize> = Vec::with_capacity(path_data.len());

        for pd in &path_data {
            if pd.min_bin == usize::MAX {
                // Empty path, assign to row 0
                path_rows.push(0);
                continue;
            }

            // Find first row where this path fits
            let mut found_row = None;
            for row in 0..occupancy.len() {
                let mut fits = true;
                for bin_idx in pd.min_bin..=pd.max_bin.min(viz_width as usize - 1) {
                    if occupancy[row][bin_idx] {
                        fits = false;
                        break;
                    }
                }
                if fits {
                    found_row = Some(row);
                    break;
                }
            }

            let row = found_row.unwrap_or_else(|| {
                // Need a new row
                occupancy.push(vec![false; viz_width as usize]);
                occupancy.len() - 1
            });

            // Mark bins as occupied
            for bin_idx in pd.min_bin..=pd.max_bin.min(viz_width as usize - 1) {
                occupancy[row][bin_idx] = true;
            }

            path_rows.push(row);
        }

        let packed_rows = occupancy.len() as u32;
        info!(
            "Pack-paths: {} paths packed into {} rows (vs {} original)",
            display_paths.len(),
            packed_rows,
            path_count
        );

        // Resize buffer if packed height is different
        let packed_path_space = packed_rows * pix_per_path;
        let packed_total_height = legend_height + packed_path_space + max_axis_height + edge_height;
        if packed_total_height != max_possible_height {
            buffer = vec![255u8; (total_width * packed_total_height * 4) as usize];
            max_y = legend_height + packed_path_space + max_axis_height;
        }

        // Render each path at its packed Y position
        for (path_idx, (pd, path)) in path_data.iter().zip(display_paths.iter()).enumerate() {
            let y_start = legend_height + path_rows[path_idx] as u32 * pix_per_path;
            let (path_r, path_g, path_b) = pd.color;
            let path_length: u64 = path
                .steps
                .iter()
                .map(|step| {
                    let seg_id = step.segment_id as usize;
                    if seg_id < graph.segments.len() {
                        graph.segments[seg_id].sequence_len
                    } else {
                        0
                    }
                })
                .sum();
            let darkness_length = if args.longest_path {
                max_path_length
            } else {
                path_length
            };

            for (bin_idx, bin_info) in &pd.bins {
                let x = (*bin_idx as u32).min(viz_width - 1);

                // Determine color (same logic as normal rendering)
                let (r, g, b) = if highlight_nodes.is_some() {
                    if bin_info.highlighted {
                        (255, 0, 0)
                    } else {
                        (180, 180, 180)
                    }
                } else if args.color_by_mean_depth {
                    get_depth_color(bin_info.mean_depth, args.no_grey_depth, depth_palette)
                } else if args.color_by_mean_inversion_rate {
                    let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                    (inv_r, 0, 0)
                } else if args.color_by_uncalled_bases {
                    let green = (bin_info.mean_uncalled * 255.0).min(255.0) as u8;
                    (0, green, 0)
                } else if args.show_strand {
                    let apply_strand = args
                        .alignment_prefix
                        .as_ref()
                        .is_none_or(|prefix| path.name.starts_with(prefix));
                    if apply_strand {
                        if bin_info.mean_inv > 0.5 {
                            (200, 50, 50)
                        } else {
                            (50, 50, 200)
                        }
                    } else {
                        (path_r, path_g, path_b)
                    }
                } else {
                    (path_r, path_g, path_b)
                };

                // Apply darkness gradient if enabled
                let (r, g, b) = if args.change_darkness && highlight_nodes.is_none() {
                    let apply_darkness = args
                        .alignment_prefix
                        .as_ref()
                        .is_none_or(|prefix| path.name.starts_with(prefix));
                    if apply_darkness && darkness_length > 0 {
                        let pos_factor = bin_info.mean_pos / darkness_length as f64;
                        let darkness = if bin_info.mean_inv > 0.5 {
                            1.0 - pos_factor
                        } else {
                            pos_factor
                        };
                        if args.white_to_black {
                            let gray = (255.0 * (1.0 - darkness)).round() as u8;
                            (gray, gray, gray)
                        } else {
                            let factor = 1.0 - (darkness * 0.8);
                            (
                                (r as f64 * factor).round() as u8,
                                (g as f64 * factor).round() as u8,
                                (b as f64 * factor).round() as u8,
                            )
                        }
                    } else {
                        (r, g, b)
                    }
                } else {
                    (r, g, b)
                };

                add_path_step(
                    &mut buffer,
                    total_width,
                    x + path_names_width,
                    y_start,
                    pix_per_path,
                    r,
                    g,
                    b,
                    args.no_path_borders,
                    args.black_path_borders,
                );
            }

            // Draw link lines between discontinuous path pieces
            if let Some(link_width) = args.link_path_pieces {
                let mut sorted_bins: Vec<usize> = pd.bins.keys().copied().collect();
                sorted_bins.sort();

                if sorted_bins.len() > 1 {
                    let link_height = ((pix_per_path as f64 * link_width).round() as u32).max(1);
                    let link_y = y_start + pix_per_path / 2 - link_height / 2;

                    for i in 1..sorted_bins.len() {
                        let prev_bin = sorted_bins[i - 1];
                        let curr_bin = sorted_bins[i];

                        if curr_bin > prev_bin + 1 {
                            let x_start =
                                (prev_bin as u32 + 1).min(viz_width - 1) + path_names_width;
                            let x_end = (curr_bin as u32).min(viz_width - 1) + path_names_width;

                            for lx in x_start..x_end {
                                for dy in 0..link_height {
                                    let y = link_y + dy;
                                    let idx = ((y * total_width + lx) * 4) as usize;
                                    if idx + 3 < buffer.len() {
                                        buffer[idx] = path_r;
                                        buffer[idx + 1] = path_g;
                                        buffer[idx + 2] = path_b;
                                        buffer[idx + 3] = 255;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Render each path (PNG) - skip if compressed mode or pack_paths mode
    let mut prev_cluster_id: Option<usize> = None;
    let mut cumulative_gap: u32 = 0;
    let cluster_gap = args.cluster_gap;

    for (path_idx, path) in display_paths.iter().enumerate() {
        // Skip normal rendering in compressed mode or pack_paths mode
        if args.compressed_mode || args.pack_paths {
            break;
        }
        // Check if grouping is enabled and get group index
        let (row_idx, base_name, is_first_in_group) = if let Some(ref pg) = path_grouping {
            let group_idx = pg.path_to_group[path_idx];
            if group_idx < 0 {
                // Skip paths that don't match any prefix
                continue;
            }
            let first = !rendered_groups.contains(&group_idx);
            if first {
                rendered_groups.insert(group_idx);
            }
            (
                group_idx as u32,
                pg.prefixes[group_idx as usize].clone(),
                first,
            )
        } else {
            (path_idx as u32, path.name.clone(), true)
        };

        // Add abundance suffix for cluster representatives
        let display_name = if args.cluster_representatives {
            if let Some(ref cr) = cluster_result {
                let cluster_id = cr.cluster_ids[path_idx];
                let size = cr.cluster_sizes[cluster_id];
                format!("{} (n={})", base_name, size)
            } else {
                base_name
            }
        } else {
            base_name
        };

        // Add gap before new cluster (except first)
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            if prev_cluster_id.is_some_and(|prev| prev != cluster_id) {
                cumulative_gap += cluster_gap;
            }
            prev_cluster_id = Some(cluster_id);
        }

        let y_start = legend_height + row_idx * pix_per_path + cumulative_gap;

        // Render cluster indicator bar on the left (only for first path in group)
        if is_first_in_group {
            if let Some(ref cr) = cluster_result {
                let cluster_id = cr.cluster_ids[path_idx];
                let (cr_r, cr_g, cr_b) = get_cluster_color(cluster_id);
                for x in dendrogram_width..(dendrogram_width + cluster_bar_width) {
                    add_path_step(
                        &mut path_names_buffer,
                        path_names_width,
                        x,
                        y_start,
                        pix_per_path,
                        cr_r,
                        cr_g,
                        cr_b,
                        true,
                        false,
                    ); // no border for cluster bar
                }
            }

            // Render annotation indicator bar (after cluster bar + gap)
            if let Some(ref ann) = annotations {
                if let Some(category) = ann.get_annotation(&path.name) {
                    if let Some(&(ar, ag, ab)) = ann.category_colors.get(category) {
                        let ann_bar_x_start = dendrogram_width + cluster_bar_width + bar_gap;
                        for x in ann_bar_x_start..(ann_bar_x_start + annotation_bar_width) {
                            add_path_step(
                                &mut path_names_buffer,
                                path_names_width,
                                x,
                                y_start,
                                pix_per_path,
                                ar,
                                ag,
                                ab,
                                true,
                                false,
                            );
                        }
                    }
                }
            }
        }

        let (path_r, path_g, path_b) = if let Some(ref colors) = custom_colors {
            colors.get(&path.name).copied().unwrap_or((200, 200, 200)) // Light grey for non-specified paths
        } else {
            compute_path_color(&path.name, args.color_by_prefix)
        };

        // Render path name (only once per group) - PNG normal paths
        if is_first_in_group && text_only_width > 0 && pix_per_path >= 8 {
            let num_of_chars = display_name.len().min(max_num_of_chars);
            let path_name_too_long = display_name.len() > num_of_chars;
            let left_padding = max_num_of_chars - num_of_chars;

            if args.color_path_names_background {
                for x in (left_padding as u32 * char_size + dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width)
                    ..path_names_width
                {
                    add_path_step(
                        &mut path_names_buffer,
                        path_names_width,
                        x,
                        y_start,
                        pix_per_path,
                        path_r,
                        path_g,
                        path_b,
                        args.no_path_borders,
                        args.black_path_borders,
                    );
                }
            }

            let base_y = y_start + pix_per_path / 2 - char_size / 2;
            for (i, c) in display_name.chars().take(num_of_chars).enumerate() {
                // +3 offset to match odgi's text positioning, shifted by dendrogram + cluster_bar + annotation_bar
                let base_x = (left_padding + i) as u32 * char_size
                    + 3
                    + dendrogram_width
                    + cluster_bar_width
                    + annotation_bar_width;
                let char_data = if i == num_of_chars - 1 && path_name_too_long {
                    &TRAILING_DOTS
                } else {
                    let c_byte = c as usize;
                    if c_byte < 128 {
                        &FONT_5X8[c_byte]
                    } else {
                        &FONT_5X8[b'?' as usize]
                    }
                };
                write_char(
                    &mut path_names_buffer,
                    path_names_width,
                    base_x,
                    base_y,
                    char_data,
                    char_size,
                    0,
                    0,
                    0,
                );
            }
        }

        // Compute bins for this path (PNG rendering)
        let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();

        // Calculate current path length for darkness gradient
        let path_length: u64 = path
            .steps
            .iter()
            .map(|step| {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    graph.segments[seg_id].sequence_len
                } else {
                    0
                }
            })
            .sum();
        let darkness_length = if args.longest_path {
            max_path_length
        } else {
            path_length
        };

        let mut path_pos: u64 = 0; // Track position within path
        for step in &path.steps {
            let seg_id = step.segment_id as usize;
            if seg_id < graph.segments.len() {
                let offset = graph.segment_offsets[seg_id];
                let seg_len = graph.segments[seg_id].sequence_len;
                let n_count = graph.segments[seg_id].n_count;
                // Proportion of N's in this segment (for uncalled base coloring)
                let n_proportion = if seg_len > 0 {
                    n_count as f64 / seg_len as f64
                } else {
                    0.0
                };

                // Check if this segment is highlighted
                let is_highlighted = highlight_nodes
                    .as_ref()
                    .is_some_and(|hn| hn.contains(&step.segment_id));

                for k in 0..seg_len {
                    let pos = offset + k;
                    let curr_bin = (pos as f64 / bin_width) as usize;
                    let entry = bins.entry(curr_bin).or_default();
                    entry.mean_depth += 1.0;
                    if step.is_reverse {
                        entry.mean_inv += 1.0;
                    }
                    entry.mean_pos += path_pos as f64;
                    entry.mean_uncalled += n_proportion;
                    if is_highlighted {
                        entry.highlighted = true;
                    }
                    path_pos += 1;
                }
            }
        }

        // Normalize bin values (PNG)
        for (_, v) in bins.iter_mut() {
            if v.mean_depth > 0.0 {
                v.mean_pos /= v.mean_depth;
                v.mean_uncalled /= v.mean_depth; // Normalize uncalled proportion
            }
            v.mean_inv /= if v.mean_depth > 0.0 {
                v.mean_depth
            } else {
                1.0
            };
            v.mean_depth /= bin_width;
        }

        // Render bins (PNG)
        for (bin_idx, bin_info) in &bins {
            let x = (*bin_idx as u32).min(viz_width - 1);

            // Determine color for this bin
            let (r, g, b) = if highlight_nodes.is_some() {
                // Highlighting mode: red for highlighted bins, grey for others
                if bin_info.highlighted {
                    (255, 0, 0)
                } else {
                    (180, 180, 180)
                }
            } else if args.color_by_mean_depth {
                // Use colorbrewer palette based on depth
                get_depth_color(bin_info.mean_depth, args.no_grey_depth, depth_palette)
            } else if args.color_by_mean_inversion_rate {
                // Black to red gradient based on inversion rate
                let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                (inv_r, 0, 0)
            } else if args.color_by_uncalled_bases {
                // Black to green gradient based on proportion of uncalled bases (N's)
                let green = (bin_info.mean_uncalled * 255.0).min(255.0) as u8;
                (0, green, 0)
            } else if args.show_strand {
                // Check if alignment_prefix applies (if set, only apply to matching paths)
                let apply_strand = args
                    .alignment_prefix
                    .as_ref()
                    .is_none_or(|prefix| path.name.starts_with(prefix));

                if apply_strand {
                    if bin_info.mean_inv > 0.5 {
                        (200, 50, 50) // Red for reverse
                    } else {
                        (50, 50, 200) // Blue for forward
                    }
                } else {
                    (path_r, path_g, path_b)
                }
            } else {
                (path_r, path_g, path_b)
            };

            // Apply darkness gradient if enabled
            let (r, g, b) = if args.change_darkness && highlight_nodes.is_none() {
                // Check if alignment_prefix applies
                let apply_darkness = args
                    .alignment_prefix
                    .as_ref()
                    .is_none_or(|prefix| path.name.starts_with(prefix));

                if apply_darkness && darkness_length > 0 {
                    // Calculate darkness factor based on position
                    let pos_factor = bin_info.mean_pos / darkness_length as f64;
                    // In binned mode: inversion rate determines gradient direction
                    let darkness = if bin_info.mean_inv > 0.5 {
                        1.0 - pos_factor // gradient from right for inverted
                    } else {
                        pos_factor // gradient from left for forward
                    };

                    if args.white_to_black {
                        // White to black gradient
                        let gray = (255.0 * (1.0 - darkness)).round() as u8;
                        (gray, gray, gray)
                    } else {
                        // Darken the path color
                        let factor = 1.0 - (darkness * 0.8); // darken up to 80%
                        (
                            (r as f64 * factor).round() as u8,
                            (g as f64 * factor).round() as u8,
                            (b as f64 * factor).round() as u8,
                        )
                    }
                } else {
                    (r, g, b)
                }
            } else {
                (r, g, b)
            };

            add_path_step(
                &mut buffer,
                total_width,
                x + path_names_width,
                y_start,
                pix_per_path,
                r,
                g,
                b,
                args.no_path_borders,
                args.black_path_borders,
            );
        }

        // Draw link lines between discontinuous path pieces
        if let Some(link_width) = args.link_path_pieces {
            let mut sorted_bins: Vec<usize> = bins.keys().copied().collect();
            sorted_bins.sort();

            if sorted_bins.len() > 1 {
                let link_height = ((pix_per_path as f64 * link_width).round() as u32).max(1);
                let link_y = y_start + pix_per_path / 2 - link_height / 2;

                for i in 1..sorted_bins.len() {
                    let prev_bin = sorted_bins[i - 1];
                    let curr_bin = sorted_bins[i];

                    // If there's a gap between bins, draw a connecting line
                    if curr_bin > prev_bin + 1 {
                        let x_start = (prev_bin as u32 + 1).min(viz_width - 1) + path_names_width;
                        let x_end = (curr_bin as u32).min(viz_width - 1) + path_names_width;

                        // Draw thin horizontal line
                        for x in x_start..x_end {
                            for dy in 0..link_height {
                                let y = link_y + dy;
                                let idx = ((y * total_width + x) * 4) as usize;
                                if idx + 3 < buffer.len() {
                                    buffer[idx] = path_r;
                                    buffer[idx + 1] = path_g;
                                    buffer[idx + 2] = path_b;
                                    buffer[idx + 3] = 255;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Calculate x-axis dimensions if enabled
    let axis_char_size = 8u32; // Use native 5x8 font
    let axis_tick_height = 4u32;
    let axis_padding = 2u32;
    let axis_label_height = axis_char_size;
    let axis_total_height = if args.x_axis.is_some() {
        axis_tick_height + axis_label_height + axis_padding * 2
    } else {
        0
    };

    // Render x-axis if requested (between paths and edges)
    if let Some(ref coord_system) = args.x_axis {
        let axis_y = legend_height + path_space + axis_padding;

        // Draw axis label on the left (in path_names_buffer if available)
        // Strip the :start-end range from the label when showing absolute coordinates
        let label_text = if coord_system.to_lowercase() == "pangenomic" {
            "pangenomic".to_string()
        } else if args.x_axis_absolute {
            strip_subpath_range(coord_system).to_string()
        } else {
            coord_system.clone()
        };

        // Draw label text in path_names_buffer (aligned like path names, bold effect)
        if path_names_width > 0 && text_only_width > 0 {
            let max_label_chars = ((text_only_width) / char_size) as usize;
            let display_label: String = if label_text.len() > max_label_chars && max_label_chars > 3
            {
                format!(
                    "{}...",
                    &label_text
                        .chars()
                        .take(max_label_chars.saturating_sub(3))
                        .collect::<String>()
                )
            } else {
                label_text.chars().take(max_label_chars).collect()
            };

            // Center vertically in axis area, similar to path name centering
            let label_y = axis_y + axis_total_height / 2 - char_size / 2;
            let left_padding = max_label_chars.saturating_sub(display_label.len());

            for (i, c) in display_label.chars().enumerate() {
                // +3 offset to match path name positioning, shifted by dendrogram + cluster_bar + annotation_bar
                let char_x = (left_padding + i) as u32 * char_size
                    + 3
                    + dendrogram_width
                    + cluster_bar_width
                    + annotation_bar_width;
                let c_byte = c as usize;
                let char_data = if c_byte < 128 {
                    &FONT_5X8[c_byte]
                } else {
                    &FONT_5X8[b'?' as usize]
                };
                // Draw twice with 1-pixel offset for bold effect
                write_char(
                    &mut path_names_buffer,
                    path_names_width,
                    char_x,
                    label_y,
                    char_data,
                    char_size,
                    0,
                    0,
                    0,
                );
                if char_x + 1 < path_names_width {
                    write_char(
                        &mut path_names_buffer,
                        path_names_width,
                        char_x + 1,
                        label_y,
                        char_data,
                        char_size,
                        0,
                        0,
                        0,
                    );
                }
            }
        }

        // Calculate tick positions and labels
        let num_ticks = args.x_ticks.max(2) as usize;
        let is_pangenomic = coord_system.to_lowercase() == "pangenomic";

        // Warn if --x-axis-absolute is used with pangenomic
        if args.x_axis_absolute && is_pangenomic {
            debug!("--x-axis-absolute has no effect with pangenomic coordinates");
        }

        // For pangenomic coordinates, use total graph length
        // For path-based coordinates, find the path and use its length
        // Also calculate pixel range where the path actually appears
        let (coord_start, coord_end, pixel_start, pixel_end) = if is_pangenomic {
            (0u64, len_to_visualize, 0u32, viz_width)
        } else if let Some(path) = graph.paths.iter().find(|p| p.name == *coord_system) {
            // Calculate path length and pangenomic positions from its steps
            let mut path_len: u64 = 0;
            let mut pangenomic_start: Option<u64> = None;
            let mut pangenomic_end: u64 = 0;

            for step in &path.steps {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    let seg_len = graph.segments[seg_id].sequence_len;
                    let seg_offset = graph.segment_offsets[seg_id];

                    // Track first segment's pangenomic position
                    if pangenomic_start.is_none() {
                        pangenomic_start = Some(seg_offset);
                    }
                    // Track last segment's end position
                    pangenomic_end = seg_offset + seg_len;
                    path_len += seg_len;
                }
            }

            let pangenomic_start = pangenomic_start.unwrap_or(0);

            // Convert pangenomic positions to pixel positions
            let pix_start = ((pangenomic_start as f64 / bin_width) as u32).min(viz_width);
            let pix_end = ((pangenomic_end as f64 / bin_width) as u32).min(viz_width);

            // Add subpath start offset if --x-axis-absolute is enabled
            let offset = if args.x_axis_absolute {
                parse_subpath_start(coord_system)
            } else {
                0
            };
            (offset, offset + path_len, pix_start, pix_end)
        } else {
            debug!(
                "Path '{}' not found, using pangenomic coordinates",
                coord_system
            );
            (0u64, len_to_visualize, 0u32, viz_width)
        };

        // Calculate the pixel width of the path's range
        let path_pixel_width = pixel_end.saturating_sub(pixel_start);

        // Draw horizontal axis line only where the path exists (PNG)
        let axis_line_start = path_names_width + pixel_start;
        let axis_line_end = path_names_width + pixel_end;
        for x in axis_line_start..axis_line_end {
            let idx = ((axis_y * total_width + x) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = 0;
                buffer[idx + 1] = 0;
                buffer[idx + 2] = 0;
                buffer[idx + 3] = 255;
            }
        }

        // Draw ticks and labels only where the path exists
        for i in 0..num_ticks {
            let t = i as f64 / (num_ticks - 1) as f64;
            // Map tick position to the path's pixel range
            let x_pos = path_names_width
                + pixel_start
                + (t * (path_pixel_width as f64 - 1.0).max(0.0)) as u32;
            let coord_value = coord_start as f64 + t * (coord_end - coord_start) as f64;

            // Draw tick mark
            for ty in 0..axis_tick_height {
                let idx = (((axis_y + ty) * total_width + x_pos) * 4) as usize;
                if idx + 3 < buffer.len() {
                    buffer[idx] = 0;
                    buffer[idx + 1] = 0;
                    buffer[idx + 2] = 0;
                    buffer[idx + 3] = 255;
                }
            }

            // Format and draw tick label
            let label = format_coordinate(coord_value as u64);
            let label_y = axis_y + axis_tick_height;

            // Calculate label x position based on tick position
            let label_width = (label.len() as u32) * axis_char_size;
            let label_x = if i == 0 {
                x_pos // Left-aligned for first tick
            } else if i == num_ticks - 1 {
                x_pos.saturating_sub(label_width) // Right-aligned for last tick
            } else {
                x_pos.saturating_sub(label_width / 2) // Center-aligned for middle ticks
            };

            for (j, c) in label.chars().enumerate() {
                let char_x = label_x + (j as u32) * axis_char_size;
                if char_x + axis_char_size <= total_width {
                    let c_byte = c as usize;
                    let char_data = if c_byte < 128 {
                        &FONT_5X8[c_byte]
                    } else {
                        &FONT_5X8[b'?' as usize]
                    };
                    // Draw twice with 1-pixel offset for bold effect
                    write_char(
                        &mut buffer,
                        total_width,
                        char_x,
                        label_y,
                        char_data,
                        axis_char_size,
                        0,
                        0,
                        0,
                    );
                    if char_x + 1 + axis_char_size <= total_width {
                        write_char(
                            &mut buffer,
                            total_width,
                            char_x + 1,
                            label_y,
                            char_data,
                            axis_char_size,
                            0,
                            0,
                            0,
                        );
                    }
                }
            }
        }
    }

    // Adjust path_space to include legend height and axis height for edge rendering
    let path_space_with_axis = legend_height + path_space + axis_total_height;

    // Render edges in the bottom area
    let mut edge_count = 0;
    for edge in &graph.edges {
        let from_id = edge.from_id as usize;
        let to_id = edge.to_id as usize;

        if from_id < graph.segments.len() && to_id < graph.segments.len() {
            // Get positions of from and to segments
            let from_offset = graph.segment_offsets[from_id];
            let from_len = graph.segments[from_id].sequence_len;
            let to_offset = graph.segment_offsets[to_id];

            // Calculate edge endpoints based on orientation
            // For forward orientation, edge exits from end of segment
            // For reverse orientation, edge exits from start of segment
            let a_pos = if edge.from_rev {
                from_offset as f64 / bin_width
            } else {
                (from_offset + from_len) as f64 / bin_width
            };

            let b_pos = if edge.to_rev {
                (to_offset + graph.segments[to_id].sequence_len) as f64 / bin_width
            } else {
                to_offset as f64 / bin_width
            };

            let (a, b) = if a_pos < b_pos {
                (a_pos, b_pos)
            } else {
                (b_pos, a_pos)
            };

            // dist = (b - a) * bin_width (in bp), used for vertical extent
            // odgi calculates this as integer
            let dist = ((b - a) * bin_width) as u64;

            // Use round() for x coordinates to match odgi's std::round()
            let ax = (a.round() as u32).min(viz_width.saturating_sub(1));
            let bx = (b.round() as u32).min(viz_width.saturating_sub(1));

            // Draw vertical line at a - iterate in world coords, scale to pixels
            // odgi: for (; i < dist; i += 1.0 / scale_y) { add_point(a, i, ...) }
            let mut i = 0.0f64;
            while i < dist as f64 {
                let y = (i * scale_y_edges).round() as u32;
                if y < edge_height {
                    add_edge_point(
                        &mut buffer,
                        total_width,
                        ax + path_names_width,
                        y,
                        path_space_with_axis,
                        0,
                    );
                    max_y = max_y.max(path_space_with_axis + y + 1);
                }
                i += 1.0 / scale_y_edges;
            }

            // Draw horizontal line from a to b at height i (where loop ended)
            let h_y = (i * scale_y_edges).round() as u32;
            let h = h_y.min(edge_height.saturating_sub(1));
            let mut x_f = a;
            while x_f <= b {
                let x = (x_f.round() as u32).min(viz_width.saturating_sub(1));
                if x < viz_width {
                    add_edge_point(
                        &mut buffer,
                        total_width,
                        x + path_names_width,
                        h,
                        path_space_with_axis,
                        0,
                    );
                    max_y = max_y.max(path_space_with_axis + h + 1);
                }
                x_f += 1.0; // In binned mode, scale_x is effectively 1
            }

            // Draw vertical line at b
            let mut j = 0.0f64;
            while j < dist as f64 {
                let y = (j * scale_y_edges).round() as u32;
                if y < edge_height {
                    add_edge_point(
                        &mut buffer,
                        total_width,
                        bx + path_names_width,
                        y,
                        path_space_with_axis,
                        0,
                    );
                }
                j += 1.0 / scale_y_edges;
            }

            edge_count += 1;
        }
    }

    debug!("Drew {} edges", edge_count);

    // Apply crop - max_y already includes path_space_with_axis, add padding
    let total_height = (path_space_with_axis + edge_height).min(max_y + bottom_padding);

    // Combine path names and main image
    if path_names_width > 0 {
        for y in 0..total_height {
            for x in 0..path_names_width {
                let src_idx = ((y * path_names_width + x) * 4) as usize;
                let dst_idx = ((y * total_width + x) * 4) as usize;
                if src_idx + 3 < path_names_buffer.len() && dst_idx + 3 < buffer.len() {
                    buffer[dst_idx] = path_names_buffer[src_idx];
                    buffer[dst_idx + 1] = path_names_buffer[src_idx + 1];
                    buffer[dst_idx + 2] = path_names_buffer[src_idx + 2];
                    buffer[dst_idx + 3] = path_names_buffer[src_idx + 3];
                }
            }
        }
    }

    // Render annotation legend at the top using full image width (PNG)
    if let Some(ref ann) = annotations {
        render_annotation_legend_png(
            &mut buffer,
            total_width,
            0, // legend starts at left edge
            &filtered_categories,
            &ann.category_colors,
            legend_height,
            char_size,
        );
    }

    // Return cropped buffer
    let cropped_size = (total_width * total_height * 4) as usize;
    let mut result = Vec::with_capacity(8 + cropped_size);
    result.extend_from_slice(&total_width.to_le_bytes());
    result.extend_from_slice(&total_height.to_le_bytes());
    result.extend_from_slice(&buffer[..cropped_size]);
    result
}

/// Write clustering results to a TSV file
fn write_cluster_tsv(
    output_path: &Path,
    display_paths: &[&GfaPath],
    cluster_result: &ClusteringResult,
) {
    // Derive TSV path from output path: foo.png -> foo.clusters.tsv
    let tsv_path = output_path.with_extension("clusters.tsv");

    let mut content = String::from("path.name\tcluster\n");
    for (path_idx, path) in display_paths.iter().enumerate() {
        let cluster_id = cluster_result.cluster_ids[path_idx];
        content.push_str(&format!("{}\t{}\n", path.name, cluster_id));
    }

    match std::fs::write(&tsv_path, content) {
        Ok(_) => info!("Cluster assignments saved to {:?}", tsv_path),
        Err(e) => eprintln!("Warning: could not write cluster TSV: {}", e),
    }
}

/// Write cluster medoids (representatives) to a TSV file
fn write_medoids_tsv(
    output_path: &Path,
    original_paths: &[&GfaPath],
    cluster_result: &ClusteringResult,
) {
    // Derive TSV path from output path: foo.png -> foo.medoids.tsv
    let tsv_path = output_path.with_extension("medoids.tsv");

    let mut content = String::from("cluster\tmedoid.path\tcluster.size\n");
    for (cluster_id, (&medoid_idx, &size)) in cluster_result
        .representatives
        .iter()
        .zip(cluster_result.cluster_sizes.iter())
        .enumerate()
    {
        let medoid_name = &original_paths[medoid_idx].name;
        content.push_str(&format!("{}\t{}\t{}\n", cluster_id, medoid_name, size));
    }

    match std::fs::write(&tsv_path, content) {
        Ok(_) => info!("Cluster medoids saved to {:?}", tsv_path),
        Err(e) => eprintln!("Warning: could not write medoids TSV: {}", e),
    }
}

/// Format coordinate value with K/M/G suffixes for readability
fn format_coordinate(value: u64) -> String {
    if value >= 1_000_000_000 {
        format!("{:.1}G", value as f64 / 1_000_000_000.0)
    } else if value >= 1_000_000 {
        format!("{:.1}M", value as f64 / 1_000_000.0)
    } else if value >= 1_000 {
        format!("{:.1}K", value as f64 / 1_000.0)
    } else {
        value.to_string()
    }
}

/// Parse the start position from a path name in "name:start-end" format.
/// Returns 0 if the format doesn't match.
fn parse_subpath_start(path_name: &str) -> u64 {
    // Look for the last colon followed by "number-number"
    if let Some(colon_pos) = path_name.rfind(':') {
        let range_part = &path_name[colon_pos + 1..];
        if let Some(dash_pos) = range_part.find('-') {
            if let Ok(start) = range_part[..dash_pos].parse::<u64>() {
                return start;
            }
        }
    }
    0
}

/// Strip the ":start-end" range from a path name if present.
/// Returns the base name without the range.
fn strip_subpath_range(path_name: &str) -> &str {
    // Look for the last colon followed by "number-number"
    if let Some(colon_pos) = path_name.rfind(':') {
        let range_part = &path_name[colon_pos + 1..];
        if let Some(dash_pos) = range_part.find('-') {
            // Verify both parts are numbers
            if range_part[..dash_pos].parse::<u64>().is_ok()
                && range_part[dash_pos + 1..].parse::<u64>().is_ok()
            {
                return &path_name[..colon_pos];
            }
        }
    }
    path_name
}

/// Render graph as SVG with vector fonts
fn render_svg(args: &Args, graph: &Graph) -> String {
    // Check for conflicting options
    if args.cluster_paths && args.prefix_merges.is_some() {
        eprintln!("[gfalook] error: -k/--cluster-paths cannot be used with -M/--prefix-merges.");
        std::process::exit(1);
    }

    let mut display_paths: Vec<&GfaPath> = graph.paths.iter().collect();

    if let Some(ref prefix) = args.ignore_prefix {
        display_paths.retain(|p| !p.name.starts_with(prefix));
    }

    if let Some(ref ptd_file) = args.paths_to_display {
        if let Ok(ptd) = load_paths_to_display(ptd_file) {
            let ptd_set: std::collections::HashSet<_> = ptd.iter().collect();
            display_paths.retain(|p| ptd_set.contains(&p.name));
            let path_map: FxHashMap<&String, &GfaPath> =
                display_paths.iter().map(|p| (&p.name, *p)).collect();
            display_paths = ptd
                .iter()
                .filter_map(|name| path_map.get(name).copied())
                .collect();
        }
    }

    let pix_per_path = args.path_height;

    let len_to_visualize = graph.total_length;

    // Calculate width - if show_all_nodes, ensure smallest segment gets at least node_width pixels
    let viz_width = if args.show_all_nodes {
        // Find the smallest segment length
        let min_seg_len = graph
            .segments
            .iter()
            .map(|s| s.sequence_len)
            .filter(|&len| len > 0)
            .min()
            .unwrap_or(1) as u32;

        // Calculate width needed so smallest segment gets node_width pixels
        // bin_width = total_length / viz_width, we want bin_width <= min_seg_len / node_width
        // So viz_width >= total_length * node_width / min_seg_len
        let min_width = ((len_to_visualize * args.node_width as u64) / min_seg_len as u64) as u32;

        debug!(
            "show_all_nodes: min_seg={}bp, need {}px width for {}px/node",
            min_seg_len, min_width, args.node_width
        );
        min_width.max(args.width)
    } else {
        args.width.min(len_to_visualize as u32)
    };

    let bin_width = args
        .bin_width
        .unwrap_or_else(|| len_to_visualize as f64 / viz_width as f64);

    // Cluster paths by similarity if requested (SVG rendering)
    let cluster_result = if args.cluster_paths {
        debug!(
            "Clustering {} paths by EDR (estimated difference rate)",
            display_paths.len()
        );
        // Build segment lengths vector for EDR computation
        let segment_lengths: Vec<u64> = graph.segments.iter().map(|s| s.sequence_len).collect();
        let original_paths = display_paths.clone(); // Save for medoids TSV
        let result = cluster_paths_by_similarity(
            &display_paths,
            &segment_lengths,
            args.cluster_threshold,
            args.cluster_all_nodes,
            args.max_clusters,
            args.dendrogram || args.use_upgma,
            args.use_upgma,
            args.upgma_threshold,
        );
        display_paths = result.ordering.iter().map(|&i| display_paths[i]).collect();
        // Write cluster assignments to TSV
        write_cluster_tsv(&args.out, &display_paths, &result);
        // Write medoids TSV
        write_medoids_tsv(&args.out, &original_paths, &result);

        // Filter to representatives only if requested (SVG)
        let result = if args.cluster_representatives {
            let rep_set: FxHashSet<usize> = result.representatives.iter().copied().collect();
            let mut filtered_paths = Vec::new();
            let mut filtered_cluster_ids = Vec::new();
            for (pos, &orig_idx) in result.ordering.iter().enumerate() {
                if rep_set.contains(&orig_idx) {
                    filtered_paths.push(display_paths[pos]);
                    filtered_cluster_ids.push(result.cluster_ids[pos]);
                }
            }
            display_paths = filtered_paths;
            ClusteringResult {
                ordering: (0..display_paths.len()).collect(),
                cluster_ids: filtered_cluster_ids,
                num_clusters: result.num_clusters,
                representatives: result.representatives,
                cluster_sizes: result.cluster_sizes,
                dendrogram: result.dendrogram,
            }
        } else {
            result
        };
        Some(result)
    } else {
        None
    };

    // Recalculate path_count after potential filtering by cluster_representatives (SVG)
    let path_count = display_paths.len() as u32;

    // Load prefix grouping if specified (SVG) - must be after clustering check
    let path_grouping: Option<PathGrouping> = args.prefix_merges.as_ref().and_then(|p| {
        let paths_vec: Vec<GfaPath> = display_paths.iter().map(|&p| p.clone()).collect();
        match load_prefix_merges(p, &paths_vec) {
            Ok(grouping) => {
                info!(
                    "Read {} valid prefixes for {} groups",
                    grouping.prefixes.len(),
                    grouping.num_groups
                );
                Some(grouping)
            }
            Err(e) => {
                eprintln!("[gfalook] warning: failed to load prefix merges: {}", e);
                None
            }
        }
    });

    // Load annotations if specified (SVG)
    let annotations: Option<AnnotationData> = args.annotation_file.as_ref().and_then(|p| {
        match load_annotations(p, args.annotation_column) {
            Ok(ann) => {
                info!(
                    "Loaded {} prefixes across {} categories (SVG)",
                    ann.prefixes.len(),
                    ann.categories.len()
                );
                Some(ann)
            }
            Err(e) => {
                eprintln!("[gfalook] warning: failed to load annotations: {}", e);
                None
            }
        }
    });

    // Effective row count: use num_groups if grouping is enabled, 1 if compressed mode
    let effective_row_count = if args.compressed_mode {
        1
    } else if let Some(ref pg) = path_grouping {
        pg.num_groups as u32
    } else {
        path_count
    };

    // Calculate text width based on longest path/prefix name, "COMPRESSED_MODE" for compressed mode
    let max_name_len = if args.compressed_mode {
        "COMPRESSED_MODE".len()
    } else if let Some(ref pg) = path_grouping {
        pg.prefixes.iter().map(|p| p.len()).max().unwrap_or(10)
    } else if args.cluster_representatives {
        // Account for " (n=X)" suffix: max cluster size determines suffix length
        let max_size = cluster_result
            .as_ref()
            .map(|cr| cr.cluster_sizes.iter().max().copied().unwrap_or(1))
            .unwrap_or(1);
        let suffix_len = format!(" (n={})", max_size).len();
        display_paths
            .iter()
            .map(|p| p.name.len() + suffix_len)
            .max()
            .unwrap_or(10)
    } else {
        display_paths
            .iter()
            .map(|p| p.name.len())
            .max()
            .unwrap_or(10)
    };
    let font_size = (pix_per_path as f64 * 0.8).max(8.0);
    let char_width = font_size * 0.6; // Approximate monospace character width
                                      // Disable path names when pack_paths is enabled (they wouldn't make sense)
    let text_width = if args.hide_path_names || args.pack_paths {
        0.0
    } else {
        (max_name_len as f64 * char_width) + 10.0
    };

    // Cluster bar width (only if clustering is enabled)
    let cluster_bar_width = if cluster_result.is_some() { 10.0 } else { 0.0 };

    // Annotation bar width (only if annotations are loaded)
    let annotation_bar_width: f64 = if annotations.is_some() {
        args.annotation_bar_width as f64
    } else {
        0.0
    };

    // Gap between cluster bar and annotation bar when both are present
    let bar_gap: f64 = if cluster_result.is_some() && annotations.is_some() { 4.0 } else { 0.0 };

    // Legend height (only if annotations are loaded)
    let legend_height: f64 = if annotations.is_some() {
        args.legend_height as f64
    } else {
        0.0
    };

    // Dendrogram width (only if dendrogram is enabled and we have a dendrogram)
    let dendrogram_width: f64 = if args.dendrogram
        && cluster_result
            .as_ref()
            .map(|cr| cr.dendrogram.is_some())
            .unwrap_or(false)
    {
        args.dendrogram_width as f64
    } else {
        0.0
    };

    let path_space = effective_row_count * pix_per_path;
    let bottom_padding = 5u32;
    let height_param = (args.height + bottom_padding) as u64;
    let edge_height = (len_to_visualize.min(height_param)) as u32;
    let scale_y_edges = edge_height as f64 / len_to_visualize as f64;

    let total_width = viz_width as f64 + text_width + cluster_bar_width + bar_gap + annotation_bar_width + dendrogram_width;
    let total_height = legend_height as u32 + path_space + edge_height;

    // Load colorbrewer palette if specified (SVG)
    let depth_palette: Option<&[(u8, u8, u8)]> = args.colorbrewer_palette.as_ref().and_then(|arg| {
        if let Some((scheme, _n)) = parse_colorbrewer_arg(arg) {
            get_colorbrewer_palette(&scheme).or_else(|| {
                eprintln!("[gfalook] warning: unknown colorbrewer palette '{}', using default Spectral", scheme);
                None
            })
        } else {
            eprintln!("[gfalook] warning: invalid colorbrewer palette format '{}', expected SCHEME:N", arg);
            None
        }
    });

    let custom_colors: Option<FxHashMap<String, (u8, u8, u8)>> = args
        .path_colors
        .as_ref()
        .and_then(|p| load_path_colors(p).ok());

    // Load highlight node IDs if specified
    let highlight_nodes: Option<FxHashSet<u64>> = args
        .highlight_node_ids
        .as_ref()
        .and_then(|p| load_highlight_node_ids(p).ok());

    // Track which groups have already been rendered (for path names)
    let mut rendered_groups: FxHashSet<i64> = FxHashSet::default();

    // Calculate max path length for longest-path option
    let max_path_length: u64 = if args.longest_path || args.change_darkness {
        display_paths
            .iter()
            .map(|path| {
                path.steps
                    .iter()
                    .map(|step| {
                        let seg_id = step.segment_id as usize;
                        if seg_id < graph.segments.len() {
                            graph.segments[seg_id].sequence_len
                        } else {
                            0
                        }
                    })
                    .sum::<u64>()
            })
            .max()
            .unwrap_or(1)
    } else {
        1
    };

    // Pre-compute leaf Y positions for dendrogram (accounting for cluster gaps) - SVG
    let dendrogram_leaf_y_positions_svg: Vec<f64> = if dendrogram_width > 0.0 {
        if let Some(ref cr) = cluster_result {
            if let Some(ref dg) = cr.dendrogram {
                let n_leaves = dg.leaf_order.len();
                let mut positions = vec![0.0f64; n_leaves];
                let mut cumulative_gap: f64 = 0.0;
                let mut prev_cluster_id: Option<usize> = None;

                for (display_pos, &orig_idx) in dg.leaf_order.iter().enumerate() {
                    if orig_idx < n_leaves && display_pos < cr.cluster_ids.len() {
                        // cluster_ids is indexed by display position, not original index
                        let cluster_id = cr.cluster_ids[display_pos];
                        if prev_cluster_id.is_some_and(|prev| prev != cluster_id) {
                            cumulative_gap += args.cluster_gap as f64;
                        }
                        prev_cluster_id = Some(cluster_id);
                        positions[orig_idx] =
                            legend_height + display_pos as f64 * pix_per_path as f64 + cumulative_gap;
                    }
                }
                positions
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    let mut svg = String::new();

    // SVG header
    svg.push_str(&format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">
<style>
  .path-name {{ font-family: 'DejaVu Sans Mono', 'Courier New', monospace; font-size: {}px; }}
</style>
<rect width="100%" height="100%" fill="white"/>
"#,
        total_width, total_height, total_width, total_height, font_size
    ));

    // Render annotation legend at the top if annotations are loaded (SVG)
    if let Some(ref ann) = annotations {
        // Filter categories to only those used by paths in the graph
        let used_categories: std::collections::HashSet<&String> = display_paths
            .iter()
            .filter_map(|p| ann.get_annotation(&p.name))
            .collect();
        let filtered_categories: Vec<String> = ann
            .categories
            .iter()
            .filter(|c| used_categories.contains(c))
            .cloned()
            .collect();

        let legend_svg = render_annotation_legend_svg(
            &filtered_categories,
            &ann.category_colors,
            total_width,
            legend_height,
            font_size,
        );
        svg.push_str(&legend_svg);
    }

    // Render dendrogram if enabled (SVG)
    if dendrogram_width > 0.0 && !dendrogram_leaf_y_positions_svg.is_empty() {
        if let Some(ref cr) = cluster_result {
            if let Some(ref dg) = cr.dendrogram {
                let dendro_svg = render_dendrogram_svg(
                    dg,
                    dendrogram_width,
                    pix_per_path as f64,
                    &dendrogram_leaf_y_positions_svg,
                );
                svg.push_str(&dendro_svg);
                svg.push('\n');
            }
        }
    }

    // Track max_y for edge rendering
    let mut max_y: f64 = legend_height + path_space as f64;

    // Compressed mode: aggregate bins across all paths and render single row (SVG)
    if args.compressed_mode {
        // Use RdBu palette by default for compressed mode, or user-specified palette
        let compressed_palette = depth_palette.unwrap_or(COLORBREWER_RDBU_11.as_slice());

        // Aggregate bins across all paths
        let mut aggregated_bins: FxHashMap<usize, f64> = FxHashMap::default();

        for path in display_paths.iter() {
            for step in &path.steps {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    let offset = graph.segment_offsets[seg_id];
                    let seg_len = graph.segments[seg_id].sequence_len;

                    for k in 0..seg_len {
                        let pos = offset + k;
                        let curr_bin = (pos as f64 / bin_width) as usize;
                        *aggregated_bins.entry(curr_bin).or_insert(0.0) += 1.0;
                    }
                }
            }
        }

        // Normalize to mean depth
        let num_paths = display_paths.len() as f64;
        let compressed_bins: FxHashMap<usize, f64> = aggregated_bins
            .iter()
            .map(|(k, v)| (*k, v / bin_width / num_paths))
            .collect();

        // Render path name "COMPRESSED_MODE"
        let y_start = legend_height;
        if !args.hide_path_names {
            let text_y = y_start + (pix_per_path as f64 / 2.0) + (font_size / 3.0);
            svg.push_str(&format!(
                r#"<text x="{}" y="{}" class="path-name" fill="black">{}</text>"#,
                dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + 5.0,
                text_y,
                "COMPRESSED_MODE"
            ));
            svg.push('\n');
        }

        // Group consecutive bins with same color for rect merging
        let mut sorted_bins: Vec<(usize, f64)> = compressed_bins.into_iter().collect();
        sorted_bins.sort_by_key(|(k, _)| *k);

        let mut prev_x: Option<usize> = None;
        let mut run_start: usize = 0;
        let mut run_color: (u8, u8, u8) = (0, 0, 0);

        for (bin_idx, mean_depth) in &sorted_bins {
            let (r, g, b) =
                get_depth_color(*mean_depth, args.no_grey_depth, Some(compressed_palette));

            if let Some(px) = prev_x {
                if *bin_idx == px + 1 && (r, g, b) == run_color {
                    // Continue the run
                } else {
                    // Output the previous run
                    let x = dendrogram_width + text_width + cluster_bar_width + bar_gap + annotation_bar_width + run_start as f64;
                    let width = (px - run_start + 1) as f64;
                    svg.push_str(&format!(
                        r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                        x, y_start, width, pix_per_path, run_color.0, run_color.1, run_color.2
                    ));
                    svg.push('\n');
                    // Start new run
                    run_start = *bin_idx;
                    run_color = (r, g, b);
                }
            } else {
                // First bin
                run_start = *bin_idx;
                run_color = (r, g, b);
            }
            prev_x = Some(*bin_idx);
        }
        // Output last run
        if let Some(px) = prev_x {
            let x = dendrogram_width + text_width + cluster_bar_width + bar_gap + annotation_bar_width + run_start as f64;
            let width = (px - run_start + 1) as f64;
            svg.push_str(&format!(
                r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                x, y_start, width, pix_per_path, run_color.0, run_color.1, run_color.2
            ));
            svg.push('\n');
        }
    }

    // Pack-paths mode: use 2D collision detection to pack paths compactly (SVG)
    if args.pack_paths && !args.compressed_mode {
        // Pre-compute bins for all paths to determine their X ranges
        struct PathBinDataSvg {
            min_bin: usize,
            max_bin: usize,
            bins: FxHashMap<usize, BinInfo>,
            color: (u8, u8, u8),
        }

        let mut path_data: Vec<PathBinDataSvg> = Vec::with_capacity(display_paths.len());

        for path in display_paths.iter() {
            let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();
            let mut min_bin = usize::MAX;
            let mut max_bin = 0usize;

            let mut path_pos: u64 = 0;
            for step in &path.steps {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    let offset = graph.segment_offsets[seg_id];
                    let seg_len = graph.segments[seg_id].sequence_len;
                    let n_count = graph.segments[seg_id].n_count;
                    let n_proportion = if seg_len > 0 {
                        n_count as f64 / seg_len as f64
                    } else {
                        0.0
                    };
                    let is_highlighted = highlight_nodes
                        .as_ref()
                        .is_some_and(|hn| hn.contains(&step.segment_id));

                    for k in 0..seg_len {
                        let pos = offset + k;
                        let curr_bin = (pos as f64 / bin_width) as usize;
                        min_bin = min_bin.min(curr_bin);
                        max_bin = max_bin.max(curr_bin);

                        let entry = bins.entry(curr_bin).or_default();
                        entry.mean_depth += 1.0;
                        if step.is_reverse {
                            entry.mean_inv += 1.0;
                        }
                        entry.mean_pos += path_pos as f64;
                        entry.mean_uncalled += n_proportion;
                        if is_highlighted {
                            entry.highlighted = true;
                        }
                        path_pos += 1;
                    }
                }
            }

            // Normalize bins
            for (_, v) in bins.iter_mut() {
                if v.mean_depth > 0.0 {
                    v.mean_pos /= v.mean_depth;
                    v.mean_uncalled /= v.mean_depth;
                }
                v.mean_inv /= if v.mean_depth > 0.0 {
                    v.mean_depth
                } else {
                    1.0
                };
                v.mean_depth /= bin_width;
            }

            let color = if let Some(ref colors) = custom_colors {
                colors.get(&path.name).copied().unwrap_or((200, 200, 200)) // Light grey for non-specified paths
            } else {
                compute_path_color(&path.name, args.color_by_prefix)
            };

            path_data.push(PathBinDataSvg {
                min_bin,
                max_bin,
                bins,
                color,
            });
        }

        // Layout buffer: for each X position, track the lowest available row
        let mut occupancy: Vec<Vec<bool>> = vec![vec![false; viz_width as usize]; 1];
        let mut path_rows: Vec<usize> = Vec::with_capacity(path_data.len());

        for pd in &path_data {
            if pd.min_bin == usize::MAX {
                path_rows.push(0);
                continue;
            }

            let mut found_row = None;
            for row in 0..occupancy.len() {
                let mut fits = true;
                for bin_idx in pd.min_bin..=pd.max_bin.min(viz_width as usize - 1) {
                    if occupancy[row][bin_idx] {
                        fits = false;
                        break;
                    }
                }
                if fits {
                    found_row = Some(row);
                    break;
                }
            }

            let row = found_row.unwrap_or_else(|| {
                occupancy.push(vec![false; viz_width as usize]);
                occupancy.len() - 1
            });

            for bin_idx in pd.min_bin..=pd.max_bin.min(viz_width as usize - 1) {
                occupancy[row][bin_idx] = true;
            }

            path_rows.push(row);
        }

        let packed_rows = occupancy.len() as u32;
        let packed_path_space = packed_rows * pix_per_path;
        max_y = legend_height + packed_path_space as f64;

        // Render each path at its packed Y position
        for (path_idx, (pd, path)) in path_data.iter().zip(display_paths.iter()).enumerate() {
            let y_start = legend_height + path_rows[path_idx] as f64 * pix_per_path as f64;
            let (path_r, path_g, path_b) = pd.color;
            let path_length: u64 = path
                .steps
                .iter()
                .map(|step| {
                    let seg_id = step.segment_id as usize;
                    if seg_id < graph.segments.len() {
                        graph.segments[seg_id].sequence_len
                    } else {
                        0
                    }
                })
                .sum();
            let darkness_length = if args.longest_path {
                max_path_length
            } else {
                path_length
            };

            // Group bins by color for rect merging
            let mut sorted_bins: Vec<(usize, &BinInfo)> =
                pd.bins.iter().map(|(k, v)| (*k, v)).collect();
            sorted_bins.sort_by_key(|(k, _)| *k);

            let mut prev_x: Option<usize> = None;
            let mut run_start: usize = 0;
            let mut run_color: (u8, u8, u8) = (0, 0, 0);

            for (bin_idx, bin_info) in &sorted_bins {
                // Calculate color
                let (r, g, b) = if highlight_nodes.is_some() {
                    if bin_info.highlighted {
                        (255, 0, 0)
                    } else {
                        (180, 180, 180)
                    }
                } else if args.color_by_mean_depth {
                    get_depth_color(bin_info.mean_depth, args.no_grey_depth, depth_palette)
                } else if args.color_by_mean_inversion_rate {
                    let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                    (inv_r, 0, 0)
                } else if args.color_by_uncalled_bases {
                    let green = (bin_info.mean_uncalled * 255.0).min(255.0) as u8;
                    (0, green, 0)
                } else if args.show_strand {
                    let apply_strand = args
                        .alignment_prefix
                        .as_ref()
                        .is_none_or(|prefix| path.name.starts_with(prefix));
                    if apply_strand {
                        if bin_info.mean_inv > 0.5 {
                            (200, 50, 50)
                        } else {
                            (50, 50, 200)
                        }
                    } else {
                        (path_r, path_g, path_b)
                    }
                } else {
                    (path_r, path_g, path_b)
                };

                let (r, g, b) = if args.change_darkness && highlight_nodes.is_none() {
                    let apply_darkness = args
                        .alignment_prefix
                        .as_ref()
                        .is_none_or(|prefix| path.name.starts_with(prefix));
                    if apply_darkness && darkness_length > 0 {
                        let pos_factor = bin_info.mean_pos / darkness_length as f64;
                        let darkness = if bin_info.mean_inv > 0.5 {
                            1.0 - pos_factor
                        } else {
                            pos_factor
                        };
                        if args.white_to_black {
                            let gray = (255.0 * (1.0 - darkness)).round() as u8;
                            (gray, gray, gray)
                        } else {
                            let factor = 1.0 - (darkness * 0.8);
                            (
                                (r as f64 * factor).round() as u8,
                                (g as f64 * factor).round() as u8,
                                (b as f64 * factor).round() as u8,
                            )
                        }
                    } else {
                        (r, g, b)
                    }
                } else {
                    (r, g, b)
                };

                if let Some(px) = prev_x {
                    if *bin_idx == px + 1 && (r, g, b) == run_color {
                        // Continue run
                    } else {
                        // Output run
                        let x =
                            dendrogram_width + text_width + cluster_bar_width + bar_gap + annotation_bar_width + run_start as f64;
                        let width = (px - run_start + 1) as f64;
                        svg.push_str(&format!(
                            r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                            x, y_start, width, pix_per_path, run_color.0, run_color.1, run_color.2
                        ));
                        svg.push('\n');
                        run_start = *bin_idx;
                        run_color = (r, g, b);
                    }
                } else {
                    run_start = *bin_idx;
                    run_color = (r, g, b);
                }
                prev_x = Some(*bin_idx);
            }
            // Output last run
            if let Some(px) = prev_x {
                let x = dendrogram_width + text_width + cluster_bar_width + bar_gap + annotation_bar_width + run_start as f64;
                let width = (px - run_start + 1) as f64;
                svg.push_str(&format!(
                    r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                    x, y_start, width, pix_per_path, run_color.0, run_color.1, run_color.2
                ));
                svg.push('\n');
            }
        }
    }

    // Render each path (SVG) - skip if compressed mode or pack_paths mode
    let mut prev_cluster_id: Option<usize> = None;
    let mut cumulative_gap: f64 = 0.0;
    let cluster_gap = args.cluster_gap as f64;

    for (path_idx, path) in display_paths.iter().enumerate() {
        // Skip normal rendering in compressed mode or pack_paths mode
        if args.compressed_mode || args.pack_paths {
            break;
        }
        // Check if grouping is enabled and get group index
        let (row_idx, base_name, is_first_in_group) = if let Some(ref pg) = path_grouping {
            let group_idx = pg.path_to_group[path_idx];
            if group_idx < 0 {
                // Skip paths that don't match any prefix
                continue;
            }
            let first = !rendered_groups.contains(&group_idx);
            if first {
                rendered_groups.insert(group_idx);
            }
            (
                group_idx as u32,
                pg.prefixes[group_idx as usize].clone(),
                first,
            )
        } else {
            (path_idx as u32, path.name.clone(), true)
        };

        // Add abundance suffix for cluster representatives
        let display_name = if args.cluster_representatives {
            if let Some(ref cr) = cluster_result {
                let cluster_id = cr.cluster_ids[path_idx];
                let size = cr.cluster_sizes[cluster_id];
                format!("{} (n={})", base_name, size)
            } else {
                base_name
            }
        } else {
            base_name
        };

        // Add gap before new cluster (except first)
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            if prev_cluster_id.is_some_and(|prev| prev != cluster_id) {
                cumulative_gap += cluster_gap;
            }
            prev_cluster_id = Some(cluster_id);
        }

        let y_start = legend_height + (row_idx * pix_per_path) as f64 + cumulative_gap;

        // Render cluster indicator bar on the left (only for first path in group)
        if is_first_in_group {
            if let Some(ref cr) = cluster_result {
                let cluster_id = cr.cluster_ids[path_idx];
                let (cr, cg, cb) = get_cluster_color(cluster_id);
                svg.push_str(&format!(
                    r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                    dendrogram_width, y_start, cluster_bar_width, pix_per_path, cr, cg, cb
                ));
                svg.push('\n');
            }

            // Render annotation indicator bar (after cluster bar)
            if let Some(ref ann) = annotations {
                if let Some(category) = ann.get_annotation(&path.name) {
                    if let Some(&(ar, ag, ab)) = ann.category_colors.get(category) {
                        svg.push_str(&format!(
                            r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                            dendrogram_width + cluster_bar_width,
                            y_start,
                            annotation_bar_width,
                            pix_per_path,
                            ar, ag, ab
                        ));
                        svg.push('\n');
                    }
                }
            }
        }

        let (path_r, path_g, path_b) = if let Some(ref colors) = custom_colors {
            colors.get(&path.name).copied().unwrap_or((200, 200, 200)) // Light grey for non-specified paths
        } else {
            compute_path_color(&path.name, args.color_by_prefix)
        };

        // Render path name (full name, vector font) - only once per group
        if is_first_in_group && !args.hide_path_names {
            let text_y = y_start + (pix_per_path as f64 / 2.0) + (font_size / 3.0);
            let text_color = if args.color_path_names_background {
                // White text on colored background
                svg.push_str(&format!(
                    r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                    dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width,
                    y_start,
                    text_width,
                    pix_per_path,
                    path_r,
                    path_g,
                    path_b
                ));
                svg.push('\n');
                "white"
            } else {
                "black"
            };
            svg.push_str(&format!(
                r#"<text x="{}" y="{}" class="path-name" fill="{}">{}</text>"#,
                dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + 5.0,
                text_y,
                text_color,
                escape_xml(&display_name)
            ));
            svg.push('\n');
        }

        // Compute bins for this path (SVG rendering)
        let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();

        // Calculate current path length for darkness gradient
        let path_length: u64 = path
            .steps
            .iter()
            .map(|step| {
                let seg_id = step.segment_id as usize;
                if seg_id < graph.segments.len() {
                    graph.segments[seg_id].sequence_len
                } else {
                    0
                }
            })
            .sum();
        let darkness_length = if args.longest_path {
            max_path_length
        } else {
            path_length
        };

        let mut path_pos: u64 = 0; // Track position within path
        for step in &path.steps {
            let seg_id = step.segment_id as usize;
            if seg_id < graph.segments.len() {
                let offset = graph.segment_offsets[seg_id];
                let seg_len = graph.segments[seg_id].sequence_len;
                let n_count = graph.segments[seg_id].n_count;
                // Proportion of N's in this segment (for uncalled base coloring)
                let n_proportion = if seg_len > 0 {
                    n_count as f64 / seg_len as f64
                } else {
                    0.0
                };

                // Check if this segment is highlighted
                let is_highlighted = highlight_nodes
                    .as_ref()
                    .is_some_and(|hn| hn.contains(&step.segment_id));

                for k in 0..seg_len {
                    let pos = offset + k;
                    let curr_bin = (pos as f64 / bin_width) as usize;
                    let entry = bins.entry(curr_bin).or_default();
                    entry.mean_depth += 1.0;
                    if step.is_reverse {
                        entry.mean_inv += 1.0;
                    }
                    entry.mean_pos += path_pos as f64;
                    entry.mean_uncalled += n_proportion;
                    if is_highlighted {
                        entry.highlighted = true;
                    }
                    path_pos += 1;
                }
            }
        }

        // Normalize bin values (SVG)
        for (_, v) in bins.iter_mut() {
            if v.mean_depth > 0.0 {
                v.mean_pos /= v.mean_depth;
                v.mean_uncalled /= v.mean_depth; // Normalize uncalled proportion
            }
            v.mean_inv /= if v.mean_depth > 0.0 {
                v.mean_depth
            } else {
                1.0
            };
            v.mean_depth /= bin_width;
        }

        // Render bins as rectangles
        let rect_height = if args.no_path_borders || pix_per_path < 3 {
            pix_per_path as f64
        } else {
            (pix_per_path - 1) as f64
        };

        // Merge consecutive bins with same color into single rectangles
        let mut bin_list: Vec<(&usize, &BinInfo)> = bins.iter().collect();
        bin_list.sort_by_key(|(idx, _)| **idx);

        // Helper to get color for a bin
        let get_bin_color = |bin_info: &BinInfo| -> (u8, u8, u8) {
            let (r, g, b) = if highlight_nodes.is_some() {
                // Highlighting mode: red for highlighted bins, grey for others
                if bin_info.highlighted {
                    (255, 0, 0)
                } else {
                    (180, 180, 180)
                }
            } else if args.color_by_mean_depth {
                get_depth_color(bin_info.mean_depth, args.no_grey_depth, depth_palette)
            } else if args.color_by_mean_inversion_rate {
                let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                (inv_r, 0, 0)
            } else if args.color_by_uncalled_bases {
                // Black to green gradient based on proportion of uncalled bases (N's)
                let green = (bin_info.mean_uncalled * 255.0).min(255.0) as u8;
                (0, green, 0)
            } else if args.show_strand {
                // Check if alignment_prefix applies
                let apply_strand = args
                    .alignment_prefix
                    .as_ref()
                    .is_none_or(|prefix| path.name.starts_with(prefix));

                if apply_strand {
                    if bin_info.mean_inv > 0.5 {
                        (200, 50, 50)
                    } else {
                        (50, 50, 200)
                    }
                } else {
                    (path_r, path_g, path_b)
                }
            } else {
                (path_r, path_g, path_b)
            };

            // Apply darkness gradient if enabled
            if args.change_darkness && highlight_nodes.is_none() {
                let apply_darkness = args
                    .alignment_prefix
                    .as_ref()
                    .is_none_or(|prefix| path.name.starts_with(prefix));

                if apply_darkness && darkness_length > 0 {
                    let pos_factor = bin_info.mean_pos / darkness_length as f64;
                    let darkness = if bin_info.mean_inv > 0.5 {
                        1.0 - pos_factor
                    } else {
                        pos_factor
                    };

                    if args.white_to_black {
                        let gray = (255.0 * (1.0 - darkness)).round() as u8;
                        (gray, gray, gray)
                    } else {
                        let factor = 1.0 - (darkness * 0.8);
                        (
                            (r as f64 * factor).round() as u8,
                            (g as f64 * factor).round() as u8,
                            (b as f64 * factor).round() as u8,
                        )
                    }
                } else {
                    (r, g, b)
                }
            } else {
                (r, g, b)
            }
        };

        if !bin_list.is_empty() {
            let mut run_start = *bin_list[0].0;
            let mut run_color = get_bin_color(bin_list[0].1);
            let mut run_end = run_start;

            for i in 1..bin_list.len() {
                let (&bin_idx, bin_info) = bin_list[i];
                let color = get_bin_color(bin_info);

                // Check if this bin continues the run (consecutive and same color)
                if bin_idx == run_end + 1 && color == run_color {
                    run_end = bin_idx;
                } else {
                    // Output the previous run
                    let x = dendrogram_width
                        + cluster_bar_width
                        + text_width
                        + (run_start as f64).min((viz_width - 1) as f64);
                    let width = (run_end - run_start + 1) as f64;
                    svg.push_str(&format!(
                        r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                        x, y_start, width, rect_height, run_color.0, run_color.1, run_color.2
                    ));
                    svg.push('\n');

                    // Start a new run
                    run_start = bin_idx;
                    run_color = color;
                    run_end = bin_idx;
                }
            }

            // Output the final run
            let x = dendrogram_width
                + cluster_bar_width
                + text_width
                + (run_start as f64).min((viz_width - 1) as f64);
            let width = (run_end - run_start + 1) as f64;
            svg.push_str(&format!(
                r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                x, y_start, width, rect_height, run_color.0, run_color.1, run_color.2
            ));
            svg.push('\n');
        }

        // Add border line if needed
        if !args.no_path_borders && pix_per_path >= 3 {
            let border_y = y_start + rect_height;
            let border_color = if args.black_path_borders {
                "black"
            } else {
                "white"
            };
            svg.push_str(&format!(
                r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1"/>"#,
                dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + text_width,
                border_y,
                total_width,
                border_y,
                border_color
            ));
            svg.push('\n');
        }

        // Draw link lines between discontinuous path pieces
        if let Some(link_width) = args.link_path_pieces {
            let mut sorted_bins: Vec<usize> = bins.keys().copied().collect();
            sorted_bins.sort();

            if sorted_bins.len() > 1 {
                let link_height = (pix_per_path as f64 * link_width).max(1.0);
                let link_y = y_start + (pix_per_path as f64 / 2.0) - (link_height / 2.0);

                for i in 1..sorted_bins.len() {
                    let prev_bin = sorted_bins[i - 1];
                    let curr_bin = sorted_bins[i];

                    // If there's a gap between bins, draw a connecting line
                    if curr_bin > prev_bin + 1 {
                        let x_start = dendrogram_width
                            + cluster_bar_width
                            + annotation_bar_width
                            + text_width
                            + (prev_bin as f64 + 1.0).min((viz_width - 1) as f64);
                        let x_end = dendrogram_width
                            + cluster_bar_width
                            + annotation_bar_width
                            + text_width
                            + (curr_bin as f64).min((viz_width - 1) as f64);
                        let line_width = x_end - x_start;

                        if line_width > 0.0 {
                            svg.push_str(&format!(
                                r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                                x_start, link_y, line_width, link_height, path_r, path_g, path_b
                            ));
                            svg.push('\n');
                        }
                    }
                }
            }
        }
    }

    // Update path space to include cumulative gap
    let path_space_with_gap = path_space as f64 + cumulative_gap;
    max_y = max_y.max(path_space_with_gap);

    // Calculate x-axis dimensions if enabled
    let axis_font_size = 10.0;
    let tick_height = 5.0;
    let axis_padding = 3.0;
    let label_height = axis_font_size + 2.0;
    let axis_total_height = if args.x_axis.is_some() {
        tick_height + label_height + axis_padding * 2.0
    } else {
        0.0
    };

    // Render x-axis if requested (between paths and edges)
    if let Some(ref coord_system) = args.x_axis {
        // Y position for the axis line (at the bottom of paths)
        let axis_y = legend_height + path_space_with_gap + axis_padding;

        // X start of the axis line (end will be calculated based on path's pangenomic extent)
        let axis_x_start = dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + text_width;

        // Draw axis label on the left
        // Strip the :start-end range from the label when showing absolute coordinates
        let label_text = if coord_system.to_lowercase() == "pangenomic" {
            "pangenomic".to_string()
        } else if args.x_axis_absolute {
            strip_subpath_range(coord_system).to_string()
        } else {
            coord_system.clone()
        };

        // Truncate label if too long for the left panel (use same char_width as path names)
        let max_label_chars = (text_width / char_width) as usize;
        let display_label = if label_text.len() > max_label_chars && max_label_chars > 3 {
            format!("{}...", &label_text[..max_label_chars.saturating_sub(3)])
        } else {
            label_text.clone()
        };

        // Position label like path names (same x offset, vertically centered in axis area)
        let label_y = axis_y + (axis_total_height / 2.0) + (font_size / 3.0);
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" class="path-name" font-weight="bold" fill="black">{}</text>"#,
            dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + 5.0,
            label_y,
            escape_xml(&display_label)
        ));
        svg.push('\n');

        // Calculate tick positions and labels
        let num_ticks = args.x_ticks.max(2) as usize;
        let is_pangenomic = coord_system.to_lowercase() == "pangenomic";

        // Warn if --x-axis-absolute is used with pangenomic
        if args.x_axis_absolute && is_pangenomic {
            debug!("--x-axis-absolute has no effect with pangenomic coordinates");
        }

        // For pangenomic coordinates, use total graph length
        // For path-based coordinates, find the path and use its length
        // Also calculate pixel range where the path actually appears
        let (coord_start, coord_end, pixel_start, pixel_end) = if is_pangenomic {
            (0u64, len_to_visualize, 0.0f64, viz_width as f64)
        } else {
            // Find the path with the specified name
            if let Some(path) = graph.paths.iter().find(|p| p.name == *coord_system) {
                // Calculate path length and pangenomic positions from its steps
                let mut path_len: u64 = 0;
                let mut pangenomic_start: Option<u64> = None;
                let mut pangenomic_end: u64 = 0;

                for step in &path.steps {
                    let seg_id = step.segment_id as usize;
                    if seg_id < graph.segments.len() {
                        let seg_len = graph.segments[seg_id].sequence_len;
                        let seg_offset = graph.segment_offsets[seg_id];

                        // Track first segment's pangenomic position
                        if pangenomic_start.is_none() {
                            pangenomic_start = Some(seg_offset);
                        }
                        // Track last segment's end position
                        pangenomic_end = seg_offset + seg_len;
                        path_len += seg_len;
                    }
                }

                let pangenomic_start = pangenomic_start.unwrap_or(0);

                // Convert pangenomic positions to pixel positions
                let pix_start = (pangenomic_start as f64 / bin_width).min(viz_width as f64);
                let pix_end = (pangenomic_end as f64 / bin_width).min(viz_width as f64);

                // Add subpath start offset if --x-axis-absolute is enabled
                let offset = if args.x_axis_absolute {
                    parse_subpath_start(coord_system)
                } else {
                    0
                };
                (offset, offset + path_len, pix_start, pix_end)
            } else {
                // Path not found, fall back to pangenomic
                debug!(
                    "Path '{}' not found, using pangenomic coordinates",
                    coord_system
                );
                (0u64, len_to_visualize, 0.0f64, viz_width as f64)
            }
        };

        // Calculate the pixel width of the path's range
        let path_pixel_width = pixel_end - pixel_start;

        // Draw horizontal axis line only where the path exists
        let axis_line_x_start = axis_x_start as f64 + pixel_start;
        let axis_line_x_end = axis_x_start as f64 + pixel_end;
        svg.push_str(&format!(
            r#"<line x1="{:.1}" y1="{}" x2="{:.1}" y2="{}" stroke="black" stroke-width="1"/>"#,
            axis_line_x_start, axis_y, axis_line_x_end, axis_y
        ));
        svg.push('\n');

        // Draw ticks only where the path exists
        for i in 0..num_ticks {
            let t = i as f64 / (num_ticks - 1) as f64;
            // Map tick position to the path's pixel range
            let x_pos = axis_x_start as f64 + pixel_start + t * (path_pixel_width - 1.0).max(0.0);
            let coord_value = coord_start as f64 + t * (coord_end - coord_start) as f64;

            // Draw tick mark
            svg.push_str(&format!(
                r#"<line x1="{:.1}" y1="{}" x2="{:.1}" y2="{}" stroke="black" stroke-width="1"/>"#,
                x_pos,
                axis_y,
                x_pos,
                axis_y + tick_height
            ));
            svg.push('\n');

            // Format coordinate value (use K/M/G suffixes for large numbers)
            let label = format_coordinate(coord_value as u64);

            // Calculate text anchor based on position
            let anchor = if i == 0 {
                "start"
            } else if i == num_ticks - 1 {
                "end"
            } else {
                "middle"
            };

            svg.push_str(&format!(
                r#"<text x="{:.1}" y="{}" font-family="'DejaVu Sans Mono', 'Courier New', monospace" font-size="{}" font-weight="bold" text-anchor="{}" fill="black">{}</text>"#,
                x_pos,
                axis_y + tick_height + axis_font_size,
                axis_font_size,
                anchor,
                label
            ));
            svg.push('\n');
        }

        // Update max_y to include axis
        max_y = path_space_with_gap + axis_total_height;
    }

    // Render edges as SVG paths (offset by x-axis height if present)
    let edge_base_y = path_space_with_gap + axis_total_height;

    for edge in &graph.edges {
        let from_id = edge.from_id as usize;
        let to_id = edge.to_id as usize;

        if from_id < graph.segments.len() && to_id < graph.segments.len() {
            let from_offset = graph.segment_offsets[from_id];
            let from_len = graph.segments[from_id].sequence_len;
            let to_offset = graph.segment_offsets[to_id];

            let a_pos = if edge.from_rev {
                from_offset as f64 / bin_width
            } else {
                (from_offset + from_len) as f64 / bin_width
            };

            let b_pos = if edge.to_rev {
                (to_offset + graph.segments[to_id].sequence_len) as f64 / bin_width
            } else {
                to_offset as f64 / bin_width
            };

            let (a, b) = if a_pos < b_pos {
                (a_pos, b_pos)
            } else {
                (b_pos, a_pos)
            };
            let dist = (b - a) * bin_width;
            let h = (dist * scale_y_edges).min(edge_height as f64 - 1.0);

            let ax = dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + text_width + a.round();
            let bx = dendrogram_width + cluster_bar_width + bar_gap + annotation_bar_width + text_width + b.round();

            // Skip degenerate edges (zero height and minimal width - not visible)
            if h < 1.0 && (bx - ax).abs() < 2.0 {
                continue;
            }

            // Draw U-shaped edge as SVG path
            svg.push_str(&format!(
                r#"<path d="M{:.1},{:.1} L{:.1},{:.1} L{:.1},{:.1} L{:.1},{:.1}" fill="none" stroke="black" stroke-width="1"/>"#,
                ax, edge_base_y,
                ax, edge_base_y + h,
                bx, edge_base_y + h,
                bx, edge_base_y
            ));
            svg.push('\n');

            max_y = max_y.max(edge_base_y + h + 1.0);
        }
    }

    // Close SVG
    svg.push_str("</svg>\n");

    // Update viewBox height to crop to actual content
    let final_height = max_y + bottom_padding as f64;
    svg = svg.replace(
        &format!(
            r#"height="{}" viewBox="0 0 {} {}"#,
            total_height, total_width, total_height
        ),
        &format!(
            r#"height="{}" viewBox="0 0 {} {}"#,
            final_height, total_width, final_height
        ),
    );

    svg
}

fn main() {
    let args = Args::parse();

    // Initialize logger based on verbosity
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();

    info!("Starting visualization...");

    let graph = match parse_gfa(&args.idx) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Error loading GFA file: {}", e);
            std::process::exit(1);
        }
    };

    if graph.paths.is_empty() {
        eprintln!("Warning: No paths found in the GFA file.");
    }

    // Detect output format by file extension
    let is_svg = args
        .out
        .extension()
        .map(|ext| ext.eq_ignore_ascii_case("svg"))
        .unwrap_or(false);

    if is_svg {
        info!("Rendering SVG...");
    } else {
        info!("Rendering image...");
    }

    if is_svg {
        // SVG output
        let svg_content = render_svg(&args, &graph);

        info!("Saving to {:?}...", args.out);

        let mut file = match File::create(&args.out) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error creating file: {}", e);
                std::process::exit(1);
            }
        };

        if let Err(e) = file.write_all(svg_content.as_bytes()) {
            eprintln!("Error writing SVG: {}", e);
            std::process::exit(1);
        }
    } else {
        // PNG output
        let buffer = render(&args, &graph);

        let width = u32::from_le_bytes([buffer[0], buffer[1], buffer[2], buffer[3]]);
        let height = u32::from_le_bytes([buffer[4], buffer[5], buffer[6], buffer[7]]);
        let pixels = &buffer[8..];

        let mut rgb_pixels = Vec::with_capacity((width * height * 3) as usize);
        for chunk in pixels.chunks(4) {
            if chunk.len() >= 3 {
                rgb_pixels.push(chunk[0]);
                rgb_pixels.push(chunk[1]);
                rgb_pixels.push(chunk[2]);
            }
        }

        info!("Saving to {:?}...", args.out);

        let img = image::RgbImage::from_raw(width, height, rgb_pixels)
            .expect("Failed to create image from buffer");

        if let Err(e) = img.save(&args.out) {
            eprintln!("Error saving image: {}", e);
            std::process::exit(1);
        }
    }

    info!("Done.");
}
