/* tslint:disable */
/* eslint-disable */

/**
 * Main data store - holds all samples and pipeline data
 */
export class GenomeData {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Calculate SNP distance matrix between all samples using polymorphic_sites
     * Returns JSON: { "samples": [...], "matrix": [[...], ...], "comparable": [[...], ...] }
     * pipeline_filter: which pipeline's data to use
     * mode: "vcf_bam" (use BAM bases when no VCF) or "vcf_ref" (use reference when no VCF, matches Snippy)
     */
    calculate_distance_matrix(pipeline_filter: string, mode: string): string;
    /**
     * Calculate distance matrix with quality filters
     * mode: "vcf_ref", "vcf_bam", or "bam_only"
     * min_depth: minimum depth to consider a position
     * min_consensus: minimum consensus percentage (0-100)
     * min_qual: minimum VCF QUAL score (only applies to VCF-sourced alleles)
     */
    calculate_distance_matrix_filtered(pipeline_filter: string, mode: string, min_depth: number, min_consensus: number, min_qual: number): string;
    /**
     * Get consensus SNP statistics (positions where ALL VCF pipelines agree)
     * Returns global and per-sample consensus vs GT comparison
     */
    get_consensus_stats(): string;
    /**
     * Get coverage statistics per sample per pipeline
     */
    get_coverage_stats(): string;
    /**
     * Get file paths for reproducibility (sample -> pipeline -> {vcf_path, bam_path})
     */
    get_file_paths(): string;
    /**
     * Legacy method for backwards compatibility
     */
    get_filtered_positions(samples_json: string, filters: string): string;
    /**
     * Get all SNP positions that match the given filters
     * filters: comma-separated list of pipeline IDs, or special filters:
     * - "consensus": positions where all pipelines agree
     * - "discordant": positions where pipelines disagree
     * - "exclusive:<pipeline>": positions only in that pipeline
     * - "gaps:<pipeline>": positions where pipeline has a gap
     * filter_mode: "and" (all filters must match) or "or" (any filter matches)
     * sample_mode: "any" (at least one sample) or "all" (all samples)
     */
    get_filtered_positions_v2(samples_json: string, filters: string, filter_mode: string, sample_mode: string): string;
    /**
     * Get report generation timestamp
     */
    get_generated_at(): string;
    /**
     * Get ground truth pileup statistics as JSON
     * Returns: { total_snps, per_sample, covered_positions, pipeline_comparison }
     */
    get_ground_truth_pileup(): string;
    /**
     * Get ground truth pipeline ID (if any)
     */
    get_ground_truth_pipeline(): string | undefined;
    /**
     * Get KPI summary as JSON
     */
    get_kpis(): string;
    /**
     * Get MNP (Multi-Nucleotide Polymorphism) statistics per pipeline
     * Returns: { pipeline_id: { mnps_found, snps_from_mnps } }
     */
    get_mnp_stats(): string;
    /**
     * Get per-sample SNP intersection with GT
     * Returns: { sample_id: { pipeline_id: { intersection, pipeline_snps, gt_snps, pct_of_pipeline, pct_of_gt } } }
     */
    get_per_sample_intersection_with_gt(): string;
    /**
     * Get per-sample statistics as JSON
     * Returns: { sample_id: { pipeline_id: { snps, snps_in_gt_gaps, agreement_with_gt, ... } } }
     */
    get_per_sample_stats(): string;
    /**
     * Get command for a pipeline (if any)
     */
    get_pipeline_command(pipeline_id: string): string | undefined;
    /**
     * Get all pipeline IDs as JSON array
     */
    get_pipeline_ids(): string;
    /**
     * Get display label for a pipeline
     */
    get_pipeline_label(pipeline_id: string): string;
    /**
     * Get reference length
     */
    get_ref_length(): number;
    /**
     * Get reference name
     */
    get_ref_name(): string;
    /**
     * Get reference nucleotide at position (0-based index)
     */
    get_ref_nuc(pos: number): string;
    /**
     * Get all sample IDs as JSON array
     */
    get_sample_ids(): string;
    /**
     * Get display label for a sample
     */
    get_sample_label(sample_id: string): string;
    /**
     * Get SNP at position (returns "ref,alt,qual,depth" or empty string)
     */
    get_snp(sample: string, pipeline: string, pos: number): string;
    /**
     * Get SNP alt allele at position (returns empty if no SNP)
     */
    get_snp_alt(sample: string, pipeline: string, pos: number): string;
    /**
     * Get SNP intersection statistics between pipelines
     * Returns: { pipeline_a: { pipeline_b: { intersection, pct_of_a, pct_of_b } } }
     */
    get_snp_intersection(): string;
    /**
     * Get SNPs in ground truth gaps statistics as JSON
     */
    get_snps_in_gt_gaps(): string;
    /**
     * Get pipelines that have VCF data (used for consensus/discordant)
     */
    get_vcf_pipelines(): string;
    /**
     * Get warnings as JSON array
     */
    get_warnings(): string;
    /**
     * Check if position is in a gap for a sample/pipeline
     */
    is_gap(sample: string, pipeline: string, pos: number): boolean;
    /**
     * Check if a pipeline is the ground truth
     */
    is_ground_truth(pipeline_id: string): boolean;
    /**
     * Load data from JSON report (the main entry point)
     */
    load_json(json: string): void;
    /**
     * Create a new empty GenomeData
     */
    constructor();
    /**
     * Render filtered view (compact, only SNP positions)
     */
    render_filtered(samples_json: string, positions_json: string, offset: number, limit: number): string;
    /**
     * Render HTML for a region
     */
    render_region(samples_json: string, start: number, end: number): string;
}

/**
 * Initialize panic hook for better error messages
 */
export function init(): void;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_genomedata_free: (a: number, b: number) => void;
    readonly genomedata_new: () => number;
    readonly genomedata_load_json: (a: number, b: number, c: number) => [number, number];
    readonly genomedata_get_ref_length: (a: number) => number;
    readonly genomedata_get_ref_name: (a: number) => [number, number];
    readonly genomedata_get_sample_ids: (a: number) => [number, number];
    readonly genomedata_get_sample_label: (a: number, b: number, c: number) => [number, number];
    readonly genomedata_get_pipeline_ids: (a: number) => [number, number];
    readonly genomedata_get_pipeline_label: (a: number, b: number, c: number) => [number, number];
    readonly genomedata_get_pipeline_command: (a: number, b: number, c: number) => [number, number];
    readonly genomedata_get_ground_truth_pipeline: (a: number) => [number, number];
    readonly genomedata_is_ground_truth: (a: number, b: number, c: number) => number;
    readonly genomedata_get_vcf_pipelines: (a: number) => [number, number];
    readonly genomedata_get_generated_at: (a: number) => [number, number];
    readonly genomedata_get_warnings: (a: number) => [number, number];
    readonly genomedata_get_snps_in_gt_gaps: (a: number) => [number, number];
    readonly genomedata_get_ground_truth_pileup: (a: number) => [number, number];
    readonly genomedata_get_mnp_stats: (a: number) => [number, number];
    readonly genomedata_get_file_paths: (a: number) => [number, number];
    readonly genomedata_get_ref_nuc: (a: number, b: number) => number;
    readonly genomedata_is_gap: (a: number, b: number, c: number, d: number, e: number, f: number) => number;
    readonly genomedata_get_snp: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly genomedata_get_snp_alt: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly genomedata_render_region: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly genomedata_get_kpis: (a: number) => [number, number];
    readonly genomedata_get_per_sample_stats: (a: number) => [number, number];
    readonly genomedata_get_snp_intersection: (a: number) => [number, number];
    readonly genomedata_get_per_sample_intersection_with_gt: (a: number) => [number, number];
    readonly genomedata_get_consensus_stats: (a: number) => [number, number];
    readonly genomedata_get_filtered_positions_v2: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => [number, number];
    readonly genomedata_get_filtered_positions: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly genomedata_render_filtered: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly genomedata_calculate_distance_matrix: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly genomedata_get_coverage_stats: (a: number) => [number, number];
    readonly genomedata_calculate_distance_matrix_filtered: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly init: () => void;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
