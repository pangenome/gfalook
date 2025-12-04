use clap::Parser;
use rustc_hash::FxHashMap;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "gfalook")]
#[command(about = "Visualize a variation graph in 1D.", long_about = None)]
struct Args {
    // MANDATORY OPTIONS
    /// Load the variation graph in GFA format from this FILE.
    #[arg(short = 'i', long = "idx", value_name = "FILE")]
    idx: PathBuf,

    /// Write the visualization in PNG format to this FILE.
    #[arg(short = 'o', long = "out", value_name = "FILE")]
    out: PathBuf,

    // Visualization Options
    /// Set the width in pixels of the output image (default: 1500).
    #[arg(short = 'x', long = "width", value_name = "N", default_value = "1500")]
    width: u32,

    /// Set the height in pixels of the output image (default: 500).
    #[arg(short = 'y', long = "height", value_name = "N", default_value = "500")]
    height: u32,

    /// The height in pixels for a path.
    #[arg(short = 'a', long = "path-height", value_name = "N")]
    path_height: Option<u32>,

    /// The padding in pixels on the x-axis for a path.
    #[arg(short = 'X', long = "path-x-padding", value_name = "N", default_value = "0")]
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
    #[arg(short = 'k', long = "cluster-paths")]
    cluster_paths: bool,

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

    // Processing Information
    /// Write the current progress to stderr.
    #[arg(short = 'P', long = "progress")]
    progress: bool,
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

/// ColorBrewer RdBu 11-class diverging palette (default for -m)
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

/// Parse a GFA file efficiently
fn parse_gfa(path: &PathBuf, progress: bool) -> std::io::Result<Graph> {
    let mut graph = Graph::new();

    if progress {
        eprintln!("[gfalook::parse] Loading GFA file...");
    }

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

    if progress {
        eprintln!(
            "[gfalook::parse] Found {} segments, total length: {} bp",
            graph.segments.len(),
            graph.total_length
        );
    }

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

    if progress {
        eprintln!("[gfalook::parse] Found {} paths, {} edges", graph.paths.len(), graph.edges.len());
    }

    Ok(graph)
}

/// Compute SHA256-based path color (matching odgi algorithm)
fn compute_path_color(path_name: &str, color_by_prefix: Option<char>) -> (u8, u8, u8) {
    let hash_input = if let Some(sep) = color_by_prefix {
        path_name.split(sep).next().unwrap_or(path_name)
    } else {
        path_name
    };

    let mut hasher = Sha256::new();
    hasher.update(hash_input.as_bytes());
    let result = hasher.finalize();

    let r = result[24];
    let g = result[8];
    let b = result[16];

    let sum = r as f32 + g as f32 + b as f32;
    if sum > 0.0 {
        let r_norm = r as f32 / sum;
        let g_norm = g as f32 / sum;
        let b_norm = b as f32 / sum;

        let factor = 1.5f32;
        let r_out = (r_norm * 255.0 * factor).min(255.0) as u8;
        let g_out = (g_norm * 255.0 * factor).min(255.0) as u8;
        let b_out = (b_norm * 255.0 * factor).min(255.0) as u8;

        (r_out, g_out, b_out)
    } else {
        (128, 128, 128)
    }
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

/// Get color for depth using RdBu colorbrewer palette
fn get_depth_color(mean_depth: f64) -> (u8, u8, u8) {
    // Cuts at 0.5, 1.5, 2.5, ..., 10.5
    let cuts = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5];

    for (i, &cut) in cuts.iter().enumerate() {
        if mean_depth <= cut {
            return COLORBREWER_RDBU_11[i];
        }
    }
    COLORBREWER_RDBU_11[10]
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
    let pix_per_path = args.path_height.unwrap_or(10);
    let bottom_padding = 5u32;

    let len_to_visualize = graph.total_length;
    let viz_width = args.width.min(len_to_visualize as u32);

    let bin_width = args.bin_width.unwrap_or_else(|| len_to_visualize as f64 / viz_width as f64);
    let scale_x = 1.0; // In binned mode
    let scale_y = viz_width as f64 / len_to_visualize as f64;

    if args.progress {
        eprintln!("[gfalook::viz] Binned mode");
        eprintln!("[gfalook::viz] bin width: {:.2e}", bin_width);
        eprintln!("[gfalook::viz] image width: {}", viz_width);
    }

    let max_name_len = display_paths.iter().map(|p| p.name.len()).max().unwrap_or(10);
    let max_num_of_chars = args.max_num_of_characters.unwrap_or(max_name_len.min(128));
    let char_size = ((pix_per_path / 8) * 8).min(64).max(8);

    let path_names_width = if args.hide_path_names {
        0u32
    } else if pix_per_path >= 8 {
        (max_num_of_chars as u32 * char_size) + char_size / 2
    } else {
        0
    };

    let path_space = path_count * pix_per_path;

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
    for (path_idx, path) in display_paths.iter().enumerate() {
        let path_y = path_idx as u32;

        let (path_r, path_g, path_b) = if let Some(ref colors) = custom_colors {
            colors.get(&path.name).copied()
                .unwrap_or_else(|| compute_path_color(&path.name, args.color_by_prefix))
        } else {
            compute_path_color(&path.name, args.color_by_prefix)
        };

        // Render path name
        if path_names_width > 0 && pix_per_path >= 8 {
            let num_of_chars = path.name.len().min(max_num_of_chars);
            let path_name_too_long = path.name.len() > num_of_chars;
            let left_padding = max_num_of_chars - num_of_chars;

            if args.color_path_names_background {
                let y_start = path_y * pix_per_path;
                for x in (left_padding as u32 * char_size)..path_names_width {
                    add_path_step(&mut path_names_buffer, path_names_width, x, y_start, pix_per_path,
                        path_r, path_g, path_b, args.no_path_borders, args.black_path_borders);
                }
            }

            let base_y = path_y * pix_per_path + pix_per_path / 2 - char_size / 2;
            for (i, c) in path.name.chars().take(num_of_chars).enumerate() {
                let base_x = (left_padding + i) as u32 * char_size;
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
        let y_start = path_y * pix_per_path;

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

            // dist = (b - a) in bins * bin_width (to bp) * scale_y_edges (to pixels)
            // This matches odgi's calculation
            let dist_bp = (b - a) * bin_width;
            let dist = (dist_bp * scale_y_edges) as u32;

            // Draw vertical line at a
            let ax = (a as u32).min(viz_width.saturating_sub(1));
            for i in 0..dist.min(edge_height) {
                add_edge_point(&mut buffer, total_width, ax + path_names_width, i, path_space, 0);
                max_y = max_y.max(path_space + i + 1);
            }

            // Draw horizontal line from a to b at height dist
            let bx = (b as u32).min(viz_width.saturating_sub(1));
            let h = dist.min(edge_height.saturating_sub(1));
            for x in ax..=bx {
                if x < viz_width {
                    add_edge_point(&mut buffer, total_width, x + path_names_width, h, path_space, 0);
                    max_y = max_y.max(path_space + h + 1);
                }
            }

            // Draw vertical line at b
            for i in 0..dist.min(edge_height) {
                add_edge_point(&mut buffer, total_width, bx + path_names_width, i, path_space, 0);
            }

            edge_count += 1;
        }
    }

    if args.progress {
        eprintln!("[gfalook::viz] Drew {} edges", edge_count);
    }

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

fn main() {
    let args = Args::parse();

    if args.progress {
        eprintln!("[gfalook] Starting visualization...");
    }

    let graph = match parse_gfa(&args.idx, args.progress) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Error loading GFA file: {}", e);
            std::process::exit(1);
        }
    };

    if graph.paths.is_empty() {
        eprintln!("Warning: No paths found in the GFA file.");
    }

    if args.progress {
        eprintln!("[gfalook] Rendering image...");
    }

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

    if args.progress {
        eprintln!("[gfalook] Saving to {:?}...", args.out);
    }

    let img = image::RgbImage::from_raw(width, height, rgb_pixels)
        .expect("Failed to create image from buffer");

    if let Err(e) = img.save(&args.out) {
        eprintln!("Error saving image: {}", e);
        std::process::exit(1);
    }

    if args.progress {
        eprintln!("[gfalook] Done.");
    }
}
