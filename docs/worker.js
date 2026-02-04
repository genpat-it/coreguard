// CoreGuard Web Worker - handles all heavy WASM operations off the main thread
// This is a MODULE worker - must be loaded with { type: 'module' }

import init, { GenomeData } from './coreguard_wasm.js';

let genomeData = null;
let wasmReady = false;

// Initialize WASM
async function initWasm() {
    try {
        await init();
        wasmReady = true;
        console.log('[Worker] WASM initialized');
        return true;
    } catch (err) {
        console.error('[Worker] WASM init failed:', err);
        return false;
    }
}

// Handle messages from main thread
self.onmessage = async function(e) {
    const { id, action, payload } = e.data;

    try {
        let result;

        switch (action) {
            case 'init':
                const success = await initWasm();
                result = { success };
                break;

            case 'loadBinary':
                genomeData = new GenomeData();
                genomeData.load_binary(new Uint8Array(payload.data));
                const pIds = JSON.parse(genomeData.get_pipeline_ids());
                const sIds = JSON.parse(genomeData.get_sample_ids());
                result = {
                    refName: genomeData.get_ref_name(),
                    refLength: genomeData.get_ref_length(),
                    sampleIds: sIds,
                    pipelineIds: pIds,
                    generatedAt: genomeData.get_generated_at(),
                    description: genomeData.get_description(),
                    warnings: JSON.parse(genomeData.get_warnings() || '[]'),
                    vcfPipelines: JSON.parse(genomeData.get_vcf_pipelines()),
                    groundTruthPipeline: genomeData.get_ground_truth_pipeline(),
                    pipelineDistanceMatrices: JSON.parse(genomeData.get_pipeline_distance_matrices() || '{}'),
                    // Pre-fetch pipeline info to avoid round-trips
                    pipelineInfo: pIds.reduce((acc, p) => {
                        acc[p] = {
                            label: genomeData.get_pipeline_label(p),
                            command: genomeData.get_pipeline_command(p),
                            isGroundTruth: genomeData.is_ground_truth(p),
                            isFromBamPileup: genomeData.is_from_bam_pileup(p)
                        };
                        return acc;
                    }, {}),
                    // Pre-fetch sample labels
                    sampleLabels: sIds.reduce((acc, s) => {
                        acc[s] = genomeData.get_sample_label(s);
                        return acc;
                    }, {})
                };
                break;

            case 'loadJson':
                genomeData = new GenomeData();
                genomeData.load_json(payload.json);
                const pIdsJ = JSON.parse(genomeData.get_pipeline_ids());
                const sIdsJ = JSON.parse(genomeData.get_sample_ids());
                result = {
                    refName: genomeData.get_ref_name(),
                    refLength: genomeData.get_ref_length(),
                    sampleIds: sIdsJ,
                    pipelineIds: pIdsJ,
                    generatedAt: genomeData.get_generated_at(),
                    description: genomeData.get_description(),
                    warnings: JSON.parse(genomeData.get_warnings() || '[]'),
                    vcfPipelines: JSON.parse(genomeData.get_vcf_pipelines()),
                    groundTruthPipeline: genomeData.get_ground_truth_pipeline(),
                    pipelineDistanceMatrices: JSON.parse(genomeData.get_pipeline_distance_matrices() || '{}'),
                    pipelineInfo: pIdsJ.reduce((acc, p) => {
                        acc[p] = {
                            label: genomeData.get_pipeline_label(p),
                            command: genomeData.get_pipeline_command(p),
                            isGroundTruth: genomeData.is_ground_truth(p),
                            isFromBamPileup: genomeData.is_from_bam_pileup(p)
                        };
                        return acc;
                    }, {}),
                    sampleLabels: sIdsJ.reduce((acc, s) => {
                        acc[s] = genomeData.get_sample_label(s);
                        return acc;
                    }, {})
                };
                break;

            case 'getPipelineInfo':
                result = {
                    pipelines: payload.pipelineIds.map(p => ({
                        id: p,
                        label: genomeData.get_pipeline_label(p),
                        command: genomeData.get_pipeline_command(p),
                        isGroundTruth: genomeData.is_ground_truth(p),
                        isFromBamPileup: genomeData.is_from_bam_pileup(p)
                    }))
                };
                break;

            case 'getSampleLabel':
                result = { label: genomeData.get_sample_label(payload.sampleId) };
                break;

            case 'getPipelineLabel':
                result = { label: genomeData.get_pipeline_label(payload.pipelineId) };
                break;

            case 'getPerSampleStats':
                result = { stats: JSON.parse(genomeData.get_per_sample_stats()) };
                break;

            case 'getCoverageStats':
                result = { stats: JSON.parse(genomeData.get_coverage_stats()) };
                break;

            case 'getFilteredPositions':
                const positionsJson = genomeData.get_filtered_positions_v2(
                    JSON.stringify(payload.samples),
                    payload.filterStr,
                    payload.filterMode,
                    payload.sampleMode
                );
                result = { positions: JSON.parse(positionsJson) };
                break;

            case 'renderFiltered':
                const html = genomeData.render_filtered(
                    JSON.stringify(payload.samples),
                    JSON.stringify(payload.positions),
                    payload.offset,
                    payload.limit
                );
                result = { html };
                break;

            case 'calculateDistanceMatrix':
                const matrixJson = genomeData.calculate_distance_matrix_filtered(
                    payload.pipeline,
                    payload.mode,
                    payload.minDepth,
                    payload.minConsensus,
                    payload.minQual
                );
                result = { matrix: JSON.parse(matrixJson) };
                break;

            case 'getReproducibilityCommands':
                result = {
                    commands: payload.pipelineIds.map(p => ({
                        id: p,
                        label: genomeData.get_pipeline_label(p),
                        command: genomeData.get_pipeline_command(p)
                    }))
                };
                break;

            case 'getConsensusStats':
                result = { stats: JSON.parse(genomeData.get_consensus_stats()) };
                break;

            case 'getPerSampleIntersectionWithGt':
                result = { intersection: JSON.parse(genomeData.get_per_sample_intersection_with_gt()) };
                break;

            case 'getFilePaths':
                result = { paths: JSON.parse(genomeData.get_file_paths()) };
                break;

            case 'getKpis':
                result = { kpis: JSON.parse(genomeData.get_kpis()) };
                break;

            case 'getMnpStats':
                result = { stats: JSON.parse(genomeData.get_mnp_stats()) };
                break;

            case 'getPipelineConcordance':
                result = { concordance: JSON.parse(genomeData.get_pipeline_concordance()) };
                break;

            case 'getSnpsInGaps':
                result = { snpsInGaps: JSON.parse(genomeData.get_snps_in_gaps()) };
                break;

            case 'getSnp':
                result = { snp: genomeData.get_snp(payload.sample, payload.pipeline, payload.pos) };
                break;

            case 'getGlobalStats':
                result = { stats: JSON.parse(genomeData.get_global_stats()) };
                break;

            case 'getPairwiseUsableStats':
                result = { stats: JSON.parse(genomeData.get_pairwise_usable_stats()) };
                break;

            case 'getReviewerPairwiseStats':
                result = { stats: JSON.parse(genomeData.get_reviewer_pairwise_stats()) };
                break;

            case 'getGlobalStatsForPipeline':
                result = { stats: JSON.parse(genomeData.get_global_stats_for_pipeline(payload.pipelineId)) };
                break;

            case 'getPairwiseUsableStatsForPipeline':
                result = { stats: JSON.parse(genomeData.get_pairwise_usable_stats_for_pipeline(payload.pipelineId)) };
                break;

            case 'getGtDiscVsPipelines':
                result = { stats: JSON.parse(genomeData.get_gt_disc_vs_pipelines()) };
                break;

            default:
                throw new Error(`Unknown action: ${action}`);
        }

        self.postMessage({ id, success: true, result });

    } catch (err) {
        console.error(`[Worker] Error in ${action}:`, err);
        self.postMessage({ id, success: false, error: err.message });
    }
};

console.log('[Worker] CoreGuard worker loaded');
