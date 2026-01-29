/* @ts-self-types="./coreguard_wasm.d.ts" */

/**
 * Main data store - holds all samples and pipeline data
 */
export class GenomeData {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        GenomeDataFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_genomedata_free(ptr, 0);
    }
    /**
     * Calculate SNP distance matrix between all samples using polymorphic_sites
     * Returns JSON: { "samples": [...], "matrix": [[...], ...], "comparable": [[...], ...] }
     * pipeline_filter: which pipeline's data to use
     * mode: "vcf_bam" (use BAM bases when no VCF) or "vcf_ref" (use reference when no VCF, matches Snippy)
     * @param {string} pipeline_filter
     * @param {string} mode
     * @returns {string}
     */
    calculate_distance_matrix(pipeline_filter, mode) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(pipeline_filter, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(mode, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_calculate_distance_matrix(this.__wbg_ptr, ptr0, len0, ptr1, len1);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Calculate distance matrix with quality filters
     * mode: "vcf_ref", "vcf_bam", or "bam_only"
     * min_depth: minimum depth to consider a position
     * min_consensus: minimum consensus percentage (0-100)
     * min_qual: minimum VCF QUAL score (only applies to VCF-sourced alleles)
     * @param {string} pipeline_filter
     * @param {string} mode
     * @param {number} min_depth
     * @param {number} min_consensus
     * @param {number} min_qual
     * @returns {string}
     */
    calculate_distance_matrix_filtered(pipeline_filter, mode, min_depth, min_consensus, min_qual) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(pipeline_filter, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(mode, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_calculate_distance_matrix_filtered(this.__wbg_ptr, ptr0, len0, ptr1, len1, min_depth, min_consensus, min_qual);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Get consensus SNP statistics (positions where ALL VCF pipelines agree)
     * Returns global and per-sample consensus vs GT comparison
     * @returns {string}
     */
    get_consensus_stats() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_consensus_stats(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get coverage statistics per sample per pipeline
     * @returns {string}
     */
    get_coverage_stats() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_coverage_stats(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get report description (markdown content)
     * Returns: description string or null if not available
     * @returns {string}
     */
    get_description() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_description(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get file paths for reproducibility (sample -> pipeline -> {vcf_path, bam_path})
     * @returns {string}
     */
    get_file_paths() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_file_paths(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Legacy method for backwards compatibility
     * @param {string} samples_json
     * @param {string} filters
     * @returns {string}
     */
    get_filtered_positions(samples_json, filters) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(samples_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(filters, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_filtered_positions(this.__wbg_ptr, ptr0, len0, ptr1, len1);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Get all SNP positions that match the given filters
     * filters: comma-separated list of pipeline IDs, or special filters:
     * - "consensus": positions where all pipelines agree
     * - "discordant": positions where pipelines disagree
     * - "exclusive:<pipeline>": positions only in that pipeline
     * - "gaps:<pipeline>": positions where pipeline has a gap
     * filter_mode: "and" (all filters must match) or "or" (any filter matches)
     * sample_mode: "any" (at least one sample) or "all" (all samples)
     * @param {string} samples_json
     * @param {string} filters
     * @param {string} filter_mode
     * @param {string} sample_mode
     * @returns {string}
     */
    get_filtered_positions_v2(samples_json, filters, filter_mode, sample_mode) {
        let deferred5_0;
        let deferred5_1;
        try {
            const ptr0 = passStringToWasm0(samples_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(filters, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ptr2 = passStringToWasm0(filter_mode, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len2 = WASM_VECTOR_LEN;
            const ptr3 = passStringToWasm0(sample_mode, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len3 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_filtered_positions_v2(this.__wbg_ptr, ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3);
            deferred5_0 = ret[0];
            deferred5_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
        }
    }
    /**
     * Get report generation timestamp
     * @returns {string}
     */
    get_generated_at() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_generated_at(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get ground truth pileup statistics as JSON
     * Returns: { total_snps, per_sample, covered_positions, pipeline_comparison }
     * @returns {string}
     */
    get_ground_truth_pileup() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_ground_truth_pileup(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get ground truth pipeline ID (if any)
     * @returns {string | undefined}
     */
    get_ground_truth_pipeline() {
        const ret = wasm.genomedata_get_ground_truth_pipeline(this.__wbg_ptr);
        let v1;
        if (ret[0] !== 0) {
            v1 = getStringFromWasm0(ret[0], ret[1]).slice();
            wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        }
        return v1;
    }
    /**
     * Get KPI summary as JSON
     * @returns {string}
     */
    get_kpis() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_kpis(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get MNP (Multi-Nucleotide Polymorphism) statistics per pipeline
     * Returns: { pipeline_id: { mnps_found, snps_from_mnps } }
     * @returns {string}
     */
    get_mnp_stats() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_mnp_stats(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get per-sample SNP intersection with GT
     * Returns: { sample_id: { pipeline_id: { intersection, pipeline_snps, gt_snps, pct_of_pipeline, pct_of_gt } } }
     * @returns {string}
     */
    get_per_sample_intersection_with_gt() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_per_sample_intersection_with_gt(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get per-sample statistics as JSON
     * Returns: { sample_id: { pipeline_id: { snps, snps_in_gt_gaps, agreement_with_gt, ... } } }
     * @returns {string}
     */
    get_per_sample_stats() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_per_sample_stats(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get command for a pipeline (if any)
     * @param {string} pipeline_id
     * @returns {string | undefined}
     */
    get_pipeline_command(pipeline_id) {
        const ptr0 = passStringToWasm0(pipeline_id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_get_pipeline_command(this.__wbg_ptr, ptr0, len0);
        let v2;
        if (ret[0] !== 0) {
            v2 = getStringFromWasm0(ret[0], ret[1]).slice();
            wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        }
        return v2;
    }
    /**
     * Get detailed pipeline concordance statistics
     * Returns 4 metrics for each pipeline pair:
     * - concordance_any: positions where at least 1 sample has SNP in both pipelines
     * - concordance_all: positions where ALL samples have SNP in both pipelines
     * - consensus_any: positions where at least 1 sample has same allele in both pipelines
     * - consensus_all: positions where ALL samples have same allele in both pipelines
     * @returns {string}
     */
    get_pipeline_concordance() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_pipeline_concordance(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get all pipeline IDs as JSON array
     * @returns {string}
     */
    get_pipeline_ids() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_pipeline_ids(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get display label for a pipeline
     * @param {string} pipeline_id
     * @returns {string}
     */
    get_pipeline_label(pipeline_id) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ptr0 = passStringToWasm0(pipeline_id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_pipeline_label(this.__wbg_ptr, ptr0, len0);
            deferred2_0 = ret[0];
            deferred2_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * Get reference length
     * @returns {number}
     */
    get_ref_length() {
        const ret = wasm.genomedata_get_ref_length(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Get reference name
     * @returns {string}
     */
    get_ref_name() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_ref_name(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get reference nucleotide at position (0-based index)
     * @param {number} pos
     * @returns {string}
     */
    get_ref_nuc(pos) {
        const ret = wasm.genomedata_get_ref_nuc(this.__wbg_ptr, pos);
        return String.fromCodePoint(ret);
    }
    /**
     * Get all sample IDs as JSON array
     * @returns {string}
     */
    get_sample_ids() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_sample_ids(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get display label for a sample
     * @param {string} sample_id
     * @returns {string}
     */
    get_sample_label(sample_id) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ptr0 = passStringToWasm0(sample_id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_sample_label(this.__wbg_ptr, ptr0, len0);
            deferred2_0 = ret[0];
            deferred2_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * Get SNP at position (returns "ref,alt,qual,depth" or empty string)
     * @param {string} sample
     * @param {string} pipeline
     * @param {number} pos
     * @returns {string}
     */
    get_snp(sample, pipeline, pos) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(sample, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(pipeline, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_snp(this.__wbg_ptr, ptr0, len0, ptr1, len1, pos);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Get SNP alt allele at position (returns empty if no SNP)
     * @param {string} sample
     * @param {string} pipeline
     * @param {number} pos
     * @returns {string}
     */
    get_snp_alt(sample, pipeline, pos) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(sample, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(pipeline, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_get_snp_alt(this.__wbg_ptr, ptr0, len0, ptr1, len1, pos);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Get SNP intersection statistics between pipelines
     * Returns: { pipeline_a: { pipeline_b: { intersection, pct_of_a, pct_of_b } } }
     * @returns {string}
     */
    get_snp_intersection() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_snp_intersection(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get SNPs in gaps for ALL pipeline pairs as JSON
     * Returns: { gap_pipeline: { snp_pipeline: { total_snps, snps_in_gaps, percentage } } }
     * @returns {string}
     */
    get_snps_in_gaps() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_snps_in_gaps(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get SNPs in ground truth gaps statistics as JSON (DEPRECATED - use get_snps_in_gaps)
     * @returns {string}
     */
    get_snps_in_gt_gaps() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_snps_in_gt_gaps(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get pipelines that have VCF data (used for consensus/discordant)
     * @returns {string}
     */
    get_vcf_pipelines() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_vcf_pipelines(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Get warnings as JSON array
     * @returns {string}
     */
    get_warnings() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.genomedata_get_warnings(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Check if a pipeline's SNPs come from BAM pileup (no variant calling)
     * @param {string} pipeline_id
     * @returns {boolean}
     */
    is_from_bam_pileup(pipeline_id) {
        const ptr0 = passStringToWasm0(pipeline_id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_is_from_bam_pileup(this.__wbg_ptr, ptr0, len0);
        return ret !== 0;
    }
    /**
     * Check if position is in a gap for a sample/pipeline
     * @param {string} sample
     * @param {string} pipeline
     * @param {number} pos
     * @returns {boolean}
     */
    is_gap(sample, pipeline, pos) {
        const ptr0 = passStringToWasm0(sample, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(pipeline, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_is_gap(this.__wbg_ptr, ptr0, len0, ptr1, len1, pos);
        return ret !== 0;
    }
    /**
     * Check if a pipeline is the ground truth
     * @param {string} pipeline_id
     * @returns {boolean}
     */
    is_ground_truth(pipeline_id) {
        const ptr0 = passStringToWasm0(pipeline_id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_is_ground_truth(this.__wbg_ptr, ptr0, len0);
        return ret !== 0;
    }
    /**
     * Load data from binary (bincode) report - faster than JSON
     * @param {Uint8Array} data
     */
    load_binary(data) {
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_load_binary(this.__wbg_ptr, ptr0, len0);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Load data from JSON report (the main entry point)
     * @param {string} json
     */
    load_json(json) {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genomedata_load_json(this.__wbg_ptr, ptr0, len0);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Create a new empty GenomeData
     */
    constructor() {
        const ret = wasm.genomedata_new();
        this.__wbg_ptr = ret >>> 0;
        GenomeDataFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Render filtered view (compact, only SNP positions)
     * @param {string} samples_json
     * @param {string} positions_json
     * @param {number} offset
     * @param {number} limit
     * @returns {string}
     */
    render_filtered(samples_json, positions_json, offset, limit) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(samples_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(positions_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_render_filtered(this.__wbg_ptr, ptr0, len0, ptr1, len1, offset, limit);
            deferred3_0 = ret[0];
            deferred3_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * Render HTML for a region
     * @param {string} samples_json
     * @param {number} start
     * @param {number} end
     * @returns {string}
     */
    render_region(samples_json, start, end) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ptr0 = passStringToWasm0(samples_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.genomedata_render_region(this.__wbg_ptr, ptr0, len0, start, end);
            deferred2_0 = ret[0];
            deferred2_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
}
if (Symbol.dispose) GenomeData.prototype[Symbol.dispose] = GenomeData.prototype.free;

/**
 * Initialize panic hook for better error messages
 */
export function init() {
    wasm.init();
}

function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_be289d5034ed271b: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_error_7534b8e9a36f1ab4: function(arg0, arg1) {
            let deferred0_0;
            let deferred0_1;
            try {
                deferred0_0 = arg0;
                deferred0_1 = arg1;
                console.error(getStringFromWasm0(arg0, arg1));
            } finally {
                wasm.__wbindgen_free(deferred0_0, deferred0_1, 1);
            }
        },
        __wbg_log_6b5ca2e6124b2808: function(arg0) {
            console.log(arg0);
        },
        __wbg_new_8a6f238a6ece86ea: function() {
            const ret = new Error();
            return ret;
        },
        __wbg_now_a3af9a2f4bbaa4d1: function() {
            const ret = Date.now();
            return ret;
        },
        __wbg_stack_0ed75d68575b0f3c: function(arg0, arg1) {
            const ret = arg1.stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbindgen_cast_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./coreguard_wasm_bg.js": import0,
    };
}

const GenomeDataFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_genomedata_free(ptr >>> 0, 1));

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return decodeText(ptr, len);
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1, 1) >>> 0;
    getUint8ArrayMemory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasm;
function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    wasmModule = module;
    cachedDataViewMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('coreguard_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
