use clap::Parser;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "gfalook")]
#[command(about = "Visualize a variation graph in 1D.", long_about = None)]
struct Args {
    // MANDATORY OPTIONS
    /// Load the variation graph in GFA format from this FILE.
    #[arg(short = 'i', long = "idx", value_name = "FILE")]
    idx: PathBuf,

    /// Write the visualization to this FILE (PNG or SVG based on extension).
    #[arg(short = 'o', long = "out", value_name = "FILE")]
    out: PathBuf,

    // Visualization Options
    /// Set the width in pixels of the output image.
    #[arg(short = 'x', long = "width", value_name = "N", default_value_t = 1500)]
    width: u32,

    /// Set the height in pixels of the output image.
    #[arg(short = 'y', long = "height", value_name = "N", default_value_t = 500)]
    height: u32,

    /// The height in pixels for a path.
    #[arg(short = 'a', long = "path-height", value_name = "N", default_value_t = 10)]
    path_height: u32,

    /// The padding in pixels on the x-axis for a path.
    #[arg(short = 'X', long = "path-x-padding", value_name = "N", default_value_t = 0)]
    path_x_padding: u32,

    /// Don't show path borders.
    #[arg(short = 'n', long = "no-path-borders")]
    no_path_borders: bool,

    /// Draw path borders in black (default is white).
    #[arg(short = 'b', long = "black-path-borders")]
    black_path_borders: bool,

    /// Pack all paths rather than displaying a single path per row.
    #[arg(short = 'R', long = "pack-paths")]
    pack_paths: bool,

    /// Show thin links of this relative width to connect path pieces.
    #[arg(short = 'L', long = "link-path-pieces", value_name = "FLOAT")]
    link_path_pieces: Option<f64>,

    /// Apply alignment related visual motifs to paths which have this name prefix.
    #[arg(short = 'A', long = "alignment-prefix", value_name = "STRING")]
    alignment_prefix: Option<String>,

    /// Use red and blue coloring to display forward and reverse alignments.
    #[arg(short = 'S', long = "show-strand")]
    show_strand: bool,

    /// Change the color respect to the node strandness (black for forward, red for reverse).
    #[arg(short = 'z', long = "color-by-mean-inversion-rate")]
    color_by_mean_inversion_rate: bool,

    /// Change the color with respect to the uncalled bases.
    #[arg(short = 'N', long = "color-by-uncalled-bases")]
    color_by_uncalled_bases: bool,

    /// Color paths by their names looking at the prefix before the given character.
    #[arg(short = 's', long = "color-by-prefix", value_name = "CHAR")]
    color_by_prefix: Option<char>,

    /// Read per-path RGB colors from FILE.
    #[arg(short = 'F', long = "path-colors", value_name = "FILE")]
    path_colors: Option<PathBuf>,

    /// Automatically order paths by similarity.
    #[arg(short = 'k', long = "cluster-paths", conflicts_with = "paths_to_display")]
    cluster_paths: bool,

    /// Similarity threshold for cluster detection (automatic if not specified).
    #[arg(long = "cluster-threshold", value_name = "F", requires = "cluster_paths")]
    cluster_threshold: Option<f64>,

    /// Use all nodes for clustering instead of only variable nodes.
    #[arg(long = "cluster-all-nodes", requires = "cluster_paths")]
    cluster_all_nodes: bool,

    /// Gap in pixels between clusters.
    #[arg(long = "cluster-gap", value_name = "N", requires = "cluster_paths", default_value_t = 10)]
    cluster_gap: u32,

    /// Color nodes listed in FILE in red and all other nodes in grey.
    #[arg(short = 'J', long = "highlight-node-ids", value_name = "FILE")]
    highlight_node_ids: Option<PathBuf>,

    /// Merge paths beginning with prefixes listed in FILE.
    #[arg(short = 'M', long = "prefix-merges", value_name = "FILE")]
    prefix_merges: Option<PathBuf>,

    /// Ignore paths starting with the given PREFIX.
    #[arg(short = 'I', long = "ignore-prefix", value_name = "PREFIX")]
    ignore_prefix: Option<String>,

    // Intervals Selection Options
    /// Nucleotide range to visualize: STRING=[PATH:]start-end.
    #[arg(short = 'r', long = "path-range", value_name = "STRING")]
    path_range: Option<String>,

    // Path Selection Options
    /// List of paths to display in the specified order.
    #[arg(short = 'p', long = "paths-to-display", value_name = "FILE")]
    paths_to_display: Option<PathBuf>,

    // Path Names Viz Options
    /// Hide the path names on the left of the generated image.
    #[arg(short = 'H', long = "hide-path-names")]
    hide_path_names: bool,

    /// Color path names background with the same color as paths.
    #[arg(short = 'C', long = "color-path-names-background")]
    color_path_names_background: bool,

    /// Maximum number of characters to display for each path name.
    #[arg(short = 'c', long = "max-num-of-characters", value_name = "N")]
    max_num_of_characters: Option<usize>,

    // Binned Mode Options
    /// The bin width specifies the size of each bin in the binned mode.
    #[arg(short = 'w', long = "bin-width", value_name = "bp")]
    bin_width: Option<f64>,

    /// Change the color with respect to the mean coverage.
    #[arg(short = 'm', long = "color-by-mean-depth")]
    color_by_mean_depth: bool,

    /// Use the colorbrewer palette specified by SCHEME:N.
    #[arg(short = 'B', long = "colorbrewer-palette", value_name = "SCHEME:N")]
    colorbrewer_palette: Option<String>,

    /// Use the colorbrewer palette for <0.5x and ~1x coverage bins.
    #[arg(short = 'G', long = "no-grey-depth")]
    no_grey_depth: bool,

    // Gradient Mode Options
    /// Change the color darkness based on nucleotide position.
    #[arg(short = 'd', long = "change-darkness")]
    change_darkness: bool,

    /// Use the longest path length to change the color darkness.
    #[arg(short = 'l', long = "longest-path")]
    longest_path: bool,

    /// Change the color darkness from white to black.
    #[arg(short = 'u', long = "white-to-black")]
    white_to_black: bool,

    // Compressed Mode Options
    /// Compress the view vertically, summarizing path coverage.
    #[arg(short = 'O', long = "compressed-mode")]
    compressed_mode: bool,

    // Threading
    /// Number of threads to use for parallel operations.
    #[arg(short = 't', long = "threads", value_name = "N")]
    threads: Option<usize>,

    // Logging
    /// Verbosity level (0 = error, 1 = info, 2 = debug).
    #[arg(short = 'v', long = "verbose", value_name = "N", default_value_t = 1)]
    verbose: u8,
}

/// A segment (node) in the graph
#[derive(Debug, Clone)]
struct Segment {
    sequence_len: u64,
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
                let seq_len = parts[2].len() as u64;
                let id = graph.segments.len() as u64;
                graph.segment_name_to_id.insert(name, id);
                graph.segments.push(Segment { sequence_len: seq_len });
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
    let mut edge_set: std::collections::HashSet<(u64, bool, u64, bool)> = std::collections::HashSet::new();

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
                    let (name, is_reverse) = if seg.ends_with('+') {
                        (&seg[..seg.len() - 1], false)
                    } else if seg.ends_with('-') {
                        (&seg[..seg.len() - 1], true)
                    } else {
                        (seg, false)
                    };
                    if let Some(&id) = graph.segment_name_to_id.get(name) {
                        steps.push(PathStep { segment_id: id, is_reverse });
                    }
                }

                graph.paths.push(GfaPath { name: path_name, steps });
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
                                steps.push(PathStep { segment_id: id, is_reverse });
                            }
                        }
                    }
                }

                graph.paths.push(GfaPath { name: path_name, steps });
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
            edge_set.insert(edge_key(from.segment_id, from.is_reverse, to.segment_id, to.is_reverse));
        }
    }

    // Convert edge set to vector
    for (from_id, from_rev, to_id, to_rev) in edge_set {
        graph.edges.push(Edge { from_id, from_rev, to_id, to_rev });
    }

    info!("Found {} paths, {} edges", graph.paths.len(), graph.edges.len());

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

/// Result of path clustering
struct ClusteringResult {
    ordering: Vec<usize>,
    cluster_ids: Vec<usize>,
    num_clusters: usize,
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
/// Tests eps from 0.01 to 0.30, finds where cluster count stabilizes
fn find_optimal_eps(dist_matrix: &[Vec<f64>], n_paths: usize) -> f64 {
    if dist_matrix.is_empty() {
        return 0.30;
    }

    // cosigt: pclust <- length(table(dbscan(distanceMatrix, eps = 0, minPts = 1)$cluster))
    let mut prev_clusters = dbscan_count_clusters(dist_matrix, 0.0);

    // cosigt: for (eps in seq(0.01, 0.30, 0.01))
    for eps_int in 1..=30 {
        let eps = eps_int as f64 / 100.0;
        let curr_clusters = dbscan_count_clusters(dist_matrix, eps);

        // cosigt: if (abs(pclust - cclust) <= 1)
        let change = (prev_clusters as i64 - curr_clusters as i64).abs();
        if change <= 1 {
            // cosigt: if ((cclust <= round(attr(distanceMatrix, "Size") / 10)))
            let max_clusters = n_paths / 10;
            if curr_clusters <= max_clusters {
                debug!("DBSCAN: stabilized at eps {:.2} with {} clusters (max allowed: {})",
                       eps, curr_clusters, max_clusters);
                return eps;
            }
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
    counts_a: &FxHashMap<u64, u64>,  // node_id -> total bp on that node for path a
    counts_b: &FxHashMap<u64, u64>,  // node_id -> total bp on that node for path b
    bp_a: u64,  // total bp in path a
    bp_b: u64,  // total bp in path b
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

/// Cluster paths by EDR (estimated difference rate) using DBSCAN (matching cosigt exactly)
/// Uses base-pair weighted Jaccard similarity like odgi similarity
fn cluster_paths_by_similarity(
    paths: &[&GfaPath],
    segment_lengths: &[u64],  // segment_id -> length (0-indexed by segment_id - 1)
    threshold: Option<f64>,
    use_all_nodes: bool,
) -> ClusteringResult {
    if paths.is_empty() {
        return ClusteringResult { ordering: Vec::new(), cluster_ids: Vec::new(), num_clusters: 0 };
    }

    let n = paths.len();

    // Build bp-weighted node counts for each path (node_id -> total bp on that node)
    // This matches odgi similarity: for each step, add segment length to that node's count
    let path_bp_counts: Vec<FxHashMap<u64, u64>> = paths
        .par_iter()
        .map(|path| {
            let mut counts: FxHashMap<u64, u64> = FxHashMap::default();
            for step in &path.steps {
                let seg_len = segment_lengths.get(step.segment_id as usize - 1).copied().unwrap_or(0);
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
        debug!("Clustering mode: --cluster-all-nodes (using all {} nodes)", total_unique_nodes);
        all_nodes
    } else {
        // Find variable nodes: nodes where bp count varies across paths
        let variable_nodes: FxHashSet<u64> = all_nodes
            .into_iter()
            .filter(|&node| {
                let first_bp = path_bp_counts[0].get(&node).copied().unwrap_or(0);
                path_bp_counts.iter().skip(1).any(|counts| {
                    counts.get(&node).copied().unwrap_or(0) != first_bp
                })
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

    // Compute filtered total bp for each path (only from nodes_to_use)
    let filtered_total_bp: Vec<u64> = filtered_bp_counts
        .iter()
        .map(|counts| counts.values().sum())
        .collect();

    // Build full pairwise EDR matrix (matching cosigt: uses normalized EDR)
    debug!("Computing {}x{} pairwise EDR matrix", n, n);

    // Compute upper triangle in parallel: EDR for each pair
    let filtered_bp_counts_ref = &filtered_bp_counts;
    let filtered_total_bp_ref = &filtered_total_bp;
    let pairs: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n)
                .map(move |j| {
                    let jaccard = weighted_jaccard_similarity(
                        &filtered_bp_counts_ref[i],
                        &filtered_bp_counts_ref[j],
                        filtered_total_bp_ref[i],
                        filtered_total_bp_ref[j],
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
        debug!("Distance range: {:.3} - {:.3}", all_dists[0], all_dists[all_dists.len() - 1]);
        if all_dists.len() >= 4 {
            let q1 = all_dists[all_dists.len() / 4];
            let median = all_dists[all_dists.len() / 2];
            let q3 = all_dists[3 * all_dists.len() / 4];
            debug!("Distance quartiles: Q1={:.3}, median={:.3}, Q3={:.3}", q1, median, q3);
        }
    }

    // Find optimal eps (or convert user threshold to eps)
    let eps = match threshold {
        Some(t) => {
            let e = 1.0 - t; // Convert similarity threshold to distance eps
            debug!("Using user-specified threshold {:.2} (eps = {:.2})", t, e);
            e
        }
        None => find_optimal_eps(&dist_matrix, n),
    };
    debug!("DBSCAN eps: {:.2}", eps);

    // Run DBSCAN to get cluster assignments
    let dbscan_clusters = dbscan_cluster(&dist_matrix, eps);
    let num_clusters = dbscan_clusters.iter().max().map(|&m| m + 1).unwrap_or(1);
    debug!("DBSCAN detected {} clusters", num_clusters);

    // Group paths by cluster
    let mut cluster_members: Vec<Vec<usize>> = vec![Vec::new(); num_clusters];
    for (i, &cluster) in dbscan_clusters.iter().enumerate() {
        cluster_members[cluster].push(i);
    }

    // Sort clusters by size (largest first) for consistent ordering
    cluster_members.sort_by(|a, b| b.len().cmp(&a.len()));

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
                .max_by_key(|&(_, &idx)| filtered_total_bp[idx])
                .map(|(local_idx, _)| local_idx)
                .unwrap();

            placed[start] = true;
            ordering.push(members[start]);
            final_cluster_ids.push(cluster_id);
            let mut current = start;

            while ordering.len() < ordering.capacity() && placed.iter().filter(|&&p| !p).count() > 0 {
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

    debug!("Final ordering: {} paths in {} clusters", ordering.len(), num_clusters);

    ClusteringResult {
        ordering,
        cluster_ids: final_cluster_ids,
        num_clusters,
    }
}

#[derive(Default, Clone)]
struct BinInfo {
    mean_depth: f64,
    mean_inv: f64,
}

fn write_char(
    buffer: &mut [u8],
    width: u32,
    base_x: u32,
    base_y: u32,
    char_data: &[u8; 8],
    char_size: u32,
    r: u8, g: u8, b: u8,
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
    r: u8, g: u8, b: u8,
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

/// Get color for depth using Spectral colorbrewer palette (with grey for low coverage)
fn get_depth_color(mean_depth: f64) -> (u8, u8, u8) {
    // Cuts at 0.5, 1.5, 2.5, ..., 12.5 (13 bins for 13 colors)
    let cuts = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5];

    for (i, &cut) in cuts.iter().enumerate() {
        if mean_depth <= cut {
            return COLORBREWER_SPECTRAL_13[i];
        }
    }
    COLORBREWER_SPECTRAL_13[12]
}

fn render(args: &Args, graph: &Graph) -> Vec<u8> {
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
            display_paths = ptd.iter().filter_map(|name| path_map.get(name).copied()).collect();
        }
    }

    let path_count = display_paths.len() as u32;
    let pix_per_path = args.path_height;
    let bottom_padding = 5u32;

    let len_to_visualize = graph.total_length;
    let viz_width = args.width.min(len_to_visualize as u32);

    let bin_width = args.bin_width.unwrap_or_else(|| len_to_visualize as f64 / viz_width as f64);
    let _scale_x = 1.0; // In binned mode
    let _scale_y = viz_width as f64 / len_to_visualize as f64;

    // Cluster paths by similarity if requested
    let cluster_result = if args.cluster_paths {
        debug!("Clustering {} paths by EDR (estimated difference rate)", display_paths.len());
        // Build segment lengths vector for EDR computation
        let segment_lengths: Vec<u64> = graph.segments.iter().map(|s| s.sequence_len).collect();
        let result = cluster_paths_by_similarity(&display_paths, &segment_lengths, args.cluster_threshold, args.cluster_all_nodes);
        display_paths = result.ordering.iter().map(|&i| display_paths[i]).collect();
        // Write cluster assignments to TSV
        write_cluster_tsv(&args.out, &display_paths, &result);
        Some(result)
    } else {
        None
    };

    // Calculate total gap space needed for cluster separators
    let total_gap = if let Some(ref cr) = cluster_result {
        (cr.num_clusters.saturating_sub(1) as u32) * args.cluster_gap
    } else {
        0
    };

    debug!("Binned mode");
    debug!("bin width: {:.2e}", bin_width);
    debug!("image width: {}", viz_width);

    let max_name_len = display_paths.iter().map(|p| p.name.len()).max().unwrap_or(10);
    let max_num_of_chars = args.max_num_of_characters.unwrap_or(max_name_len.min(128));
    let char_size = ((pix_per_path / 8) * 8).min(64).max(8);

    // Cluster bar width (only if clustering is enabled)
    let cluster_bar_width: u32 = if cluster_result.is_some() { 10 } else { 0 };

    let text_only_width = if args.hide_path_names {
        0u32
    } else if pix_per_path >= 8 {
        (max_num_of_chars as u32 * char_size) + char_size / 2
    } else {
        0
    };

    // Total left panel width includes cluster bar + path names
    let path_names_width = text_only_width + cluster_bar_width;

    let path_space = path_count * pix_per_path + total_gap;

    // Height for edge visualization area - matches odgi's calculation
    // height = min(len_to_visualize, args.height + bottom_padding)
    // scale_y = height / len_to_visualize
    let height_param = (args.height + bottom_padding) as u64;
    let edge_height = (len_to_visualize.min(height_param)) as u32;
    let scale_y_edges = edge_height as f64 / len_to_visualize as f64;

    let total_width = viz_width + path_names_width;
    // Initial height - will be cropped later based on actual edge rendering
    let max_possible_height = path_space + edge_height;

    let mut buffer = vec![255u8; (total_width * max_possible_height * 4) as usize];
    let mut path_names_buffer = if path_names_width > 0 {
        vec![255u8; (path_names_width * max_possible_height * 4) as usize]
    } else {
        Vec::new()
    };

    // Track maximum y coordinate used (for cropping)
    let mut max_y: u32 = path_space;

    let custom_colors: Option<FxHashMap<String, (u8, u8, u8)>> = args
        .path_colors
        .as_ref()
        .and_then(|p| load_path_colors(p).ok());

    // Render each path
    let mut prev_cluster_id: Option<usize> = None;
    let mut cumulative_gap: u32 = 0;
    let cluster_gap = args.cluster_gap;

    for (path_idx, path) in display_paths.iter().enumerate() {
        // Add gap before new cluster (except first)
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            if prev_cluster_id.map_or(false, |prev| prev != cluster_id) {
                cumulative_gap += cluster_gap;
            }
            prev_cluster_id = Some(cluster_id);
        }

        let y_start = path_idx as u32 * pix_per_path + cumulative_gap;

        // Render cluster indicator bar on the left
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            let (cr_r, cr_g, cr_b) = get_cluster_color(cluster_id);
            for x in 0..cluster_bar_width {
                add_path_step(&mut path_names_buffer, path_names_width, x, y_start, pix_per_path,
                    cr_r, cr_g, cr_b, true, false); // no border for cluster bar
            }
        }

        let (path_r, path_g, path_b) = if let Some(ref colors) = custom_colors {
            colors.get(&path.name).copied()
                .unwrap_or_else(|| compute_path_color(&path.name, args.color_by_prefix))
        } else {
            compute_path_color(&path.name, args.color_by_prefix)
        };

        // Render path name
        if text_only_width > 0 && pix_per_path >= 8 {
            let num_of_chars = path.name.len().min(max_num_of_chars);
            let path_name_too_long = path.name.len() > num_of_chars;
            let left_padding = max_num_of_chars - num_of_chars;

            if args.color_path_names_background {
                for x in (left_padding as u32 * char_size + cluster_bar_width)..path_names_width {
                    add_path_step(&mut path_names_buffer, path_names_width, x, y_start, pix_per_path,
                        path_r, path_g, path_b, args.no_path_borders, args.black_path_borders);
                }
            }

            let base_y = y_start + pix_per_path / 2 - char_size / 2;
            for (i, c) in path.name.chars().take(num_of_chars).enumerate() {
                // +3 offset to match odgi's text positioning, shifted by cluster_bar_width
                let base_x = (left_padding + i) as u32 * char_size + 3 + cluster_bar_width;
                let char_data = if i == num_of_chars - 1 && path_name_too_long {
                    &TRAILING_DOTS
                } else {
                    let c_byte = c as usize;
                    if c_byte < 128 { &FONT_5X8[c_byte] } else { &FONT_5X8[b'?' as usize] }
                };
                write_char(&mut path_names_buffer, path_names_width, base_x, base_y, char_data, char_size, 0, 0, 0);
            }
        }

        // Compute bins for this path
        let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();

        for step in &path.steps {
            let seg_id = step.segment_id as usize;
            if seg_id < graph.segments.len() {
                let offset = graph.segment_offsets[seg_id];
                let seg_len = graph.segments[seg_id].sequence_len;

                for k in 0..seg_len {
                    let pos = offset + k;
                    let curr_bin = (pos as f64 / bin_width) as usize;
                    let entry = bins.entry(curr_bin).or_default();
                    entry.mean_depth += 1.0;
                    if step.is_reverse {
                        entry.mean_inv += 1.0;
                    }
                }
            }
        }

        // Normalize bin values
        for (_, v) in bins.iter_mut() {
            v.mean_inv /= if v.mean_depth > 0.0 { v.mean_depth } else { 1.0 };
            v.mean_depth /= bin_width;
        }

        // Render bins
        for (bin_idx, bin_info) in &bins {
            let x = (*bin_idx as u32).min(viz_width - 1);

            // Determine color for this bin
            let (r, g, b) = if args.color_by_mean_depth {
                // Use colorbrewer palette based on depth
                get_depth_color(bin_info.mean_depth)
            } else if args.color_by_mean_inversion_rate {
                // Black to red gradient based on inversion rate
                let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                (inv_r, 0, 0)
            } else if args.show_strand {
                if bin_info.mean_inv > 0.5 {
                    (200, 50, 50) // Red for reverse
                } else {
                    (50, 50, 200) // Blue for forward
                }
            } else {
                (path_r, path_g, path_b)
            };

            add_path_step(&mut buffer, total_width, x + path_names_width, y_start, pix_per_path,
                r, g, b, args.no_path_borders, args.black_path_borders);
        }
    }

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

            let (a, b) = if a_pos < b_pos { (a_pos, b_pos) } else { (b_pos, a_pos) };

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
                    add_edge_point(&mut buffer, total_width, ax + path_names_width, y, path_space, 0);
                    max_y = max_y.max(path_space + y + 1);
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
                    add_edge_point(&mut buffer, total_width, x + path_names_width, h, path_space, 0);
                    max_y = max_y.max(path_space + h + 1);
                }
                x_f += 1.0; // In binned mode, scale_x is effectively 1
            }

            // Draw vertical line at b
            let mut j = 0.0f64;
            while j < dist as f64 {
                let y = (j * scale_y_edges).round() as u32;
                if y < edge_height {
                    add_edge_point(&mut buffer, total_width, bx + path_names_width, y, path_space, 0);
                }
                j += 1.0 / scale_y_edges;
            }

            edge_count += 1;
        }
    }

    debug!("Drew {} edges", edge_count);

    // Apply crop - max_y already includes path_space, add padding
    let total_height = (path_space + edge_height).min(max_y + bottom_padding);

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
    output_path: &PathBuf,
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

/// Escape special XML characters
fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Render graph as SVG with vector fonts
fn render_svg(args: &Args, graph: &Graph) -> String {
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
            display_paths = ptd.iter().filter_map(|name| path_map.get(name).copied()).collect();
        }
    }

    let path_count = display_paths.len() as u32;
    let pix_per_path = args.path_height;

    let len_to_visualize = graph.total_length;
    let viz_width = args.width.min(len_to_visualize as u32);
    let bin_width = args.bin_width.unwrap_or_else(|| len_to_visualize as f64 / viz_width as f64);

    // Cluster paths by similarity if requested
    let cluster_result = if args.cluster_paths {
        debug!("Clustering {} paths by EDR (estimated difference rate)", display_paths.len());
        // Build segment lengths vector for EDR computation
        let segment_lengths: Vec<u64> = graph.segments.iter().map(|s| s.sequence_len).collect();
        let result = cluster_paths_by_similarity(&display_paths, &segment_lengths, args.cluster_threshold, args.cluster_all_nodes);
        display_paths = result.ordering.iter().map(|&i| display_paths[i]).collect();
        // Write cluster assignments to TSV
        write_cluster_tsv(&args.out, &display_paths, &result);
        Some(result)
    } else {
        None
    };

    // Calculate text width based on longest path name
    let max_name_len = display_paths.iter().map(|p| p.name.len()).max().unwrap_or(10);
    let font_size = (pix_per_path as f64 * 0.8).max(8.0);
    let char_width = font_size * 0.6; // Approximate monospace character width
    let text_width = if args.hide_path_names {
        0.0
    } else {
        (max_name_len as f64 * char_width) + 10.0
    };

    // Cluster bar width (only if clustering is enabled)
    let cluster_bar_width = if cluster_result.is_some() { 10.0 } else { 0.0 };

    let path_space = path_count * pix_per_path;
    let bottom_padding = 5u32;
    let height_param = (args.height + bottom_padding) as u64;
    let edge_height = (len_to_visualize.min(height_param)) as u32;
    let scale_y_edges = edge_height as f64 / len_to_visualize as f64;

    let total_width = viz_width as f64 + text_width + cluster_bar_width;
    let total_height = path_space + edge_height;

    let custom_colors: Option<FxHashMap<String, (u8, u8, u8)>> = args
        .path_colors
        .as_ref()
        .and_then(|p| load_path_colors(p).ok());

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

    // Track max_y for edge rendering
    let mut max_y: f64 = path_space as f64;

    // Render each path
    let mut prev_cluster_id: Option<usize> = None;
    let mut cumulative_gap: f64 = 0.0;
    let cluster_gap = args.cluster_gap as f64;

    for (path_idx, path) in display_paths.iter().enumerate() {
        // Add gap before new cluster (except first)
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            if prev_cluster_id.map_or(false, |prev| prev != cluster_id) {
                cumulative_gap += cluster_gap;
            }
            prev_cluster_id = Some(cluster_id);
        }

        let path_y = path_idx as u32;
        let y_start = (path_y * pix_per_path) as f64 + cumulative_gap;

        // Render cluster indicator bar on the left
        if let Some(ref cr) = cluster_result {
            let cluster_id = cr.cluster_ids[path_idx];
            let (cr, cg, cb) = get_cluster_color(cluster_id);
            svg.push_str(&format!(
                r#"<rect x="0" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                y_start, cluster_bar_width, pix_per_path, cr, cg, cb
            ));
            svg.push('\n');
        }

        let (path_r, path_g, path_b) = if let Some(ref colors) = custom_colors {
            colors.get(&path.name).copied()
                .unwrap_or_else(|| compute_path_color(&path.name, args.color_by_prefix))
        } else {
            compute_path_color(&path.name, args.color_by_prefix)
        };

        // Render path name (full name, vector font)
        if !args.hide_path_names {
            let text_y = y_start + (pix_per_path as f64 / 2.0) + (font_size / 3.0);
            let text_color = if args.color_path_names_background {
                // White text on colored background
                svg.push_str(&format!(
                    r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgb({},{},{})"/>"#,
                    cluster_bar_width, y_start, text_width, pix_per_path, path_r, path_g, path_b
                ));
                svg.push('\n');
                "white"
            } else {
                "black"
            };
            svg.push_str(&format!(
                r#"<text x="{}" y="{}" class="path-name" fill="{}">{}</text>"#,
                cluster_bar_width + 5.0, text_y, text_color, escape_xml(&path.name)
            ));
            svg.push('\n');
        }

        // Compute bins for this path
        let mut bins: FxHashMap<usize, BinInfo> = FxHashMap::default();

        for step in &path.steps {
            let seg_id = step.segment_id as usize;
            if seg_id < graph.segments.len() {
                let offset = graph.segment_offsets[seg_id];
                let seg_len = graph.segments[seg_id].sequence_len;

                for k in 0..seg_len {
                    let pos = offset + k;
                    let curr_bin = (pos as f64 / bin_width) as usize;
                    let entry = bins.entry(curr_bin).or_default();
                    entry.mean_depth += 1.0;
                    if step.is_reverse {
                        entry.mean_inv += 1.0;
                    }
                }
            }
        }

        // Normalize bin values
        for (_, v) in bins.iter_mut() {
            v.mean_inv /= if v.mean_depth > 0.0 { v.mean_depth } else { 1.0 };
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
            if args.color_by_mean_depth {
                get_depth_color(bin_info.mean_depth)
            } else if args.color_by_mean_inversion_rate {
                let inv_r = (bin_info.mean_inv * 255.0).min(255.0) as u8;
                (inv_r, 0, 0)
            } else if args.show_strand {
                if bin_info.mean_inv > 0.5 {
                    (200, 50, 50)
                } else {
                    (50, 50, 200)
                }
            } else {
                (path_r, path_g, path_b)
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
                    let x = cluster_bar_width + text_width + (run_start as f64).min((viz_width - 1) as f64);
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
            let x = cluster_bar_width + text_width + (run_start as f64).min((viz_width - 1) as f64);
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
            let border_color = if args.black_path_borders { "black" } else { "white" };
            svg.push_str(&format!(
                r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1"/>"#,
                cluster_bar_width + text_width, border_y, total_width, border_y, border_color
            ));
            svg.push('\n');
        }
    }

    // Update path space to include cumulative gap
    let path_space_with_gap = path_space as f64 + cumulative_gap;
    max_y = max_y.max(path_space_with_gap);

    // Render edges as SVG paths
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

            let (a, b) = if a_pos < b_pos { (a_pos, b_pos) } else { (b_pos, a_pos) };
            let dist = ((b - a) * bin_width) as f64;
            let h = (dist * scale_y_edges).min(edge_height as f64 - 1.0);

            let ax = cluster_bar_width + text_width + a.round();
            let bx = cluster_bar_width + text_width + b.round();
            let base_y = path_space_with_gap;

            // Skip degenerate edges (zero height and minimal width - not visible)
            if h < 1.0 && (bx - ax).abs() < 2.0 {
                continue;
            }

            // Draw U-shaped edge as SVG path
            svg.push_str(&format!(
                r#"<path d="M{:.1},{:.1} L{:.1},{:.1} L{:.1},{:.1} L{:.1},{:.1}" fill="none" stroke="black" stroke-width="1"/>"#,
                ax, base_y,
                ax, base_y + h,
                bx, base_y + h,
                bx, base_y
            ));
            svg.push('\n');

            max_y = max_y.max(base_y + h + 1.0);
        }
    }

    // Close SVG
    svg.push_str("</svg>\n");

    // Update viewBox height to crop to actual content
    let final_height = max_y + bottom_padding as f64;
    svg = svg.replace(
        &format!(r#"height="{}" viewBox="0 0 {} {}"#, total_height, total_width, total_height),
        &format!(r#"height="{}" viewBox="0 0 {} {}"#, final_height, total_width, final_height),
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
    let is_svg = args.out.extension()
        .map(|ext| ext.to_ascii_lowercase() == "svg")
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
