#!/usr/bin/env python3
"""
UMI-based deduplication for CRISPR screening data.

This script performs genome-free deduplication using:
1. Library-guided matching: each read is matched to the closest library sgRNA (exact or fuzzy)
2. UMI pooling: all reads matching the same sgRNA pool their UMIs together
3. UMI clustering: UMI-tools algorithms deduplicate UMIs within each sgRNA

This approach is specifically designed for CRISPR screening where:
- We have a KNOWN library of sgRNA sequences
- Reads should be matched TO THE LIBRARY (not clustered among themselves)
- UMI sequencing errors must be handled AFTER library matching

Best Practices (UMI Deduplication):
- Default method: 'directional' (umi_tools) - Gold standard for UMI dedup
  - Uses count-based network clustering where high-count UMIs absorb error variants
  - Accounts for PCR amplification bias (A can absorb B if count(A) >= 2*count(B) - 1)
  - Reference: Smith et al. 2017, Genome Research (UMI-tools paper)
- Alternative methods:
  - 'adjacency': Less stringent, good for sparse data
  - 'cluster': Connected components only, maximum error correction
  - 'unique': No clustering, counts distinct UMIs (use when UMI errors negligible)
- Edit distance: Default 1 is appropriate for most 8-12bp UMIs

Author: nf-core/crisprseq contributors
Released under the MIT license.
"""

import argparse
import gzip
import sys
from collections import defaultdict, Counter

# Check for optional dependencies
try:
    from umi_tools import UMIClusterer
    HAS_UMI_TOOLS = True
except ImportError:
    HAS_UMI_TOOLS = False

try:
    import Levenshtein
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False


def levenshtein_distance_pure(s1, s2):
    """Pure Python Levenshtein distance implementation.
    
    Used as fallback when python-levenshtein is not installed.
    """
    if len(s1) < len(s2):
        return levenshtein_distance_pure(s2, s1)
    if len(s2) == 0:
        return len(s1)
    
    prev_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    
    return prev_row[-1]


def edit_distance(s1, s2):
    """Calculate edit distance using best available method."""
    if HAS_LEVENSHTEIN:
        return Levenshtein.distance(s1, s2)
    return levenshtein_distance_pure(s1, s2)


def parse_args():
    parser = argparse.ArgumentParser(
        description="UMI-based deduplication for CRISPR screening"
    )
    # Input/output
    parser.add_argument("--r1", required=True, help="R1 FASTQ file (gzipped)")
    parser.add_argument("--r2", help="R2 FASTQ file (gzipped) - optional")
    parser.add_argument("--library", required=True, help="Library TSV file")
    parser.add_argument("--output", required=True, help="Output count table")
    parser.add_argument("--sample-name", default="sample",
                        help="Sample name for count table header")
    
    # sgRNA sequence matching parameters
    parser.add_argument("--sgRNA-start", type=int, default=0,
                        help="Position where sgRNA starts in read (default: 0)")
    parser.add_argument("--sgRNA-length", type=int, default=20,
                        help="sgRNA length (default: 20)")
    parser.add_argument("--sgrna-edit-dist", type=int, default=1,
                        help="Max edit distance for sgRNA-to-library matching (default: 1). "
                             "Set to 0 for exact matching only.")
    parser.add_argument("--sgrna-exact-only", action="store_true",
                        help="Require exact sgRNA matches only (equivalent to --sgrna-edit-dist 0)")
    
    # UMI clustering parameters  
    parser.add_argument("--umi-edit-dist", type=int, default=1,
                        help="Edit distance for UMI clustering (default: 1). "
                             "Set to 0 to count unique UMIs without clustering.")
    parser.add_argument("--umi-method", default="directional",
                        choices=["directional", "cluster", "adjacency", "unique"],
                        help="UMI clustering method (default: directional). "
                             "'unique' counts each distinct UMI sequence.")
    
    # Legacy parameter (maps to umi-edit-dist for backwards compatibility)
    parser.add_argument("--edit-dist", type=int, default=None,
                        help="[DEPRECATED] Use --umi-edit-dist instead. "
                             "Sets edit distance for both sgRNA and UMI if specific params not set.")
    
    return parser.parse_args()


def read_library(library_file):
    """Read sgRNA library TSV file.
    
    Expected format: sgRNA_id\tsequence\tgene (tab-separated, with header)
    """
    library = {}  # sequence -> (sgRNA_id, gene)
    with open(library_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                sgrna_id, sequence, gene = parts[0], parts[1].upper(), parts[2]
                library[sequence] = (sgrna_id, gene)
    return library


def read_fastq_with_umi(fastq_file, sgrna_start, sgrna_length):
    """Read FASTQ file and extract sgRNA sequence + UMI from header.
    
    Expects UMI to be appended to read name by umi_tools extract (format: @readname_UMISEQUENCE)
    
    Returns: list of (sgrna_seq, umi) tuples
    """
    reads = []
    opener = gzip.open if fastq_file.endswith('.gz') else open
    
    with opener(fastq_file, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()
            
            # Extract UMI from header (umi_tools format: @readname_UMISEQUENCE)
            if '_' in header:
                umi = header.split('_')[-1].split()[0]  # Get UMI, handle any trailing info
            else:
                umi = "NOUMI"
            
            # Extract sgRNA sequence
            sgrna_seq = seq[sgrna_start:sgrna_start + sgrna_length].upper()
            
            if len(sgrna_seq) == sgrna_length:
                reads.append((sgrna_seq, umi))
    
    return reads


def cluster_umis_umi_tools(umis_with_counts, method, edit_distance):
    """Cluster UMIs using UMI-tools algorithms.
    
    Args:
        umis_with_counts: dict of UMI -> count
        method: clustering method (directional, cluster, adjacency, unique)
        edit_distance: maximum edit distance for clustering
    
    Returns: number of unique UMI clusters
    """
    if method == "unique" or not HAS_UMI_TOOLS:
        return len(umis_with_counts)
    
    clusterer = UMIClusterer(cluster_method=method)
    
    # Format: {umi_bytes: count}
    umi_counts = {umi.encode(): count for umi, count in umis_with_counts.items()}
    
    clusters = clusterer(umi_counts, threshold=edit_distance)
    return len(clusters)


def cluster_umis_simple(umis_with_counts, edit_dist):
    """Simple UMI clustering fallback.
    
    Greedy clustering based on edit distance.
    Now uses pure Python Levenshtein as fallback when library not available.
    """
    if edit_dist == 0:
        return len(umis_with_counts)
    
    umi_list = list(umis_with_counts.keys())
    clusters = []
    assigned = set()
    
    # Sort by count (highest first) for directional-like behavior
    umi_list.sort(key=lambda x: umis_with_counts[x], reverse=True)
    
    for umi in umi_list:
        if umi in assigned:
            continue
        
        cluster = [umi]
        assigned.add(umi)
        
        for other in umi_list:
            if other not in assigned:
                if edit_distance(umi, other) <= edit_dist:
                    cluster.append(other)
                    assigned.add(other)
        
        clusters.append(cluster)
    
    return len(clusters)


def match_to_library(centroid, library, max_edit_dist):
    """Match a sequence centroid to the closest library sgRNA.
    
    Returns: (sgRNA_id, gene) or None if no match within max_edit_dist
    """
    # Exact match first
    if centroid in library:
        return library[centroid]
    
    # Fuzzy match using edit distance
    if max_edit_dist > 0:
        best_match = None
        best_dist = max_edit_dist + 1
        
        for lib_seq, (sgrna_id, gene) in library.items():
            dist = edit_distance(centroid, lib_seq)
            if dist < best_dist:
                best_dist = dist
                best_match = (sgrna_id, gene)
        
        if best_match and best_dist <= max_edit_dist:
            return best_match
    
    return None


def main():
    args = parse_args()
    
    # Handle parameter precedence and backwards compatibility
    # If legacy --edit-dist is set but new params aren't, use it for both
    if args.edit_dist is not None:
        print("WARNING: --edit-dist is deprecated. Use --sgrna-edit-dist and --umi-edit-dist instead.", 
              file=sys.stderr)
        sgrna_edit_dist = args.sgrna_edit_dist if args.sgrna_edit_dist != 1 else args.edit_dist
        umi_edit_dist = args.umi_edit_dist if args.umi_edit_dist != 1 else args.edit_dist
    else:
        sgrna_edit_dist = args.sgrna_edit_dist
        umi_edit_dist = args.umi_edit_dist
    
    # Handle --sgrna-exact-only flag
    if args.sgrna_exact_only:
        sgrna_edit_dist = 0
    
    print(f"Reading library from {args.library}...", file=sys.stderr)
    library = read_library(args.library)
    print(f"  Loaded {len(library)} sgRNAs", file=sys.stderr)
    
    print(f"Reading FASTQ file {args.r1}...", file=sys.stderr)
    reads = read_fastq_with_umi(args.r1, args.sgRNA_start, args.sgRNA_length)
    print(f"  Read {len(reads)} reads", file=sys.stderr)
    
    # Count unique read sequences for stats
    unique_seqs = set(seq for seq, umi in reads)
    print(f"  Found {len(unique_seqs)} unique sgRNA sequences", file=sys.stderr)
    
    # STEP 1: Match each read to library FIRST (before any grouping)
    # This ensures error-containing reads get matched to their correct library sgRNA
    if sgrna_edit_dist == 0:
        print(f"Matching reads to library (exact matches only)...", file=sys.stderr)
    else:
        print(f"Matching reads to library (max edit distance {sgrna_edit_dist})...", file=sys.stderr)
    
    sgrna_to_umis = defaultdict(list)  # sgRNA_id -> list of UMIs
    sgrna_genes = {}  # sgRNA_id -> gene
    unmatched_reads = 0
    exact_matches = 0
    fuzzy_matches = 0
    
    for sgrna_seq, umi in reads:
        match = match_to_library(sgrna_seq, library, sgrna_edit_dist)
        if match:
            sgrna_id, gene = match
            sgrna_to_umis[sgrna_id].append(umi)
            sgrna_genes[sgrna_id] = gene
            # Track match type for stats
            if sgrna_seq in library:
                exact_matches += 1
            else:
                fuzzy_matches += 1
        else:
            unmatched_reads += 1
    
    print(f"  Exact matches: {exact_matches}", file=sys.stderr)
    if sgrna_edit_dist > 0:
        print(f"  Fuzzy matches (corrected): {fuzzy_matches}", file=sys.stderr)
    print(f"  Unmatched reads: {unmatched_reads}", file=sys.stderr)
    print(f"  Reads assigned to {len(sgrna_to_umis)} sgRNAs", file=sys.stderr)
    
    # STEP 2: Cluster UMIs for each sgRNA (deduplication)
    # Now all reads matching the same sgRNA have their UMIs pooled together
    # UMI clustering handles UMI sequencing errors
    if umi_edit_dist == 0 or args.umi_method == "unique":
        print(f"Counting unique UMIs per sgRNA (no clustering)...", file=sys.stderr)
    else:
        print(f"Counting unique UMIs per sgRNA ({args.umi_method} method, edit dist {umi_edit_dist})...", 
              file=sys.stderr)
    
    sgrna_counts = {}  # sgRNA_id -> unique UMI count
    for sgrna_id, umis in sgrna_to_umis.items():
        umi_counts = Counter(umis)
        
        if umi_edit_dist == 0 or args.umi_method == "unique":
            # No clustering - count each unique UMI
            unique_count = len(umi_counts)
        elif HAS_UMI_TOOLS:
            unique_count = cluster_umis_umi_tools(
                umi_counts, args.umi_method, umi_edit_dist
            )
        else:
            unique_count = cluster_umis_simple(umi_counts, umi_edit_dist)
        
        sgrna_counts[sgrna_id] = unique_count
    
    total_umis = sum(sgrna_counts.values())
    print(f"  Total unique UMI counts: {total_umis}", file=sys.stderr)
    
    # STEP 3: Write output count table (MAGeCK-compatible format)
    print(f"Writing output to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as f:
        # Header: sgRNA Gene Sample1 [Sample2 ...]
        f.write(f"sgRNA\tGene\t{args.sample_name}\n")
        
        # Write counts for all library sgRNAs (including zeros)
        for lib_seq, (sgrna_id, gene) in library.items():
            count = sgrna_counts.get(sgrna_id, 0)
            f.write(f"{sgrna_id}\t{gene}\t{count}\n")
    
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
