//! Report generation and output

use crate::cluster::ClusterAnalysis;
use crate::gaps::ReferenceGaps;
use crate::gff::GeneZoneAnalysis;
use crate::impact::ImpactAnalysis;
use crate::pairwise::{PairwiseResults, PairwiseGapResult};
use crate::reference::ReferenceResults;
use crate::scoring::ScoredSample;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Execution metadata
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ExecutionMetadata {
    /// Command line used
    #[serde(default)]
    pub command_line: String,
    /// Timestamp of execution
    #[serde(default)]
    pub timestamp: String,
    /// Duration in seconds
    #[serde(default)]
    pub duration_secs: f64,
    /// Input mode (FASTA or FASTQ)
    #[serde(default)]
    pub input_mode: String,
    /// Input directory path
    #[serde(default)]
    pub input_path: String,
    /// Reference file path
    #[serde(default)]
    pub reference_path: String,
}

/// Input file information
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct InputFileInfo {
    /// Sample name
    pub name: String,
    /// File path(s)
    pub paths: Vec<String>,
    /// Total file size in bytes
    pub total_size_bytes: u64,
    /// File size formatted (e.g., "85.2 MB")
    pub size_formatted: String,
}

/// Gene zone analysis summary
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GeneZoneSummary {
    /// GFF file path
    pub gff_path: String,
    /// Total genes in annotation
    pub total_genes: usize,
    /// Total CDS in annotation
    pub total_cds: usize,
    /// Per-sample gene zone analysis
    pub samples: Vec<GeneZoneAnalysis>,
    /// Total unique genes affected across all samples
    pub total_affected_genes: usize,
    /// Genes completely removed in at least one sample
    pub genes_removed_any_sample: usize,
}

/// Complete analysis report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    /// Version of the tool
    pub version: String,
    /// Execution metadata
    #[serde(default)]
    pub metadata: ExecutionMetadata,
    /// Input file information
    #[serde(default)]
    pub input_files: Vec<InputFileInfo>,
    /// Analysis summary
    pub summary: Summary,
    /// Per-sample results
    pub samples: Vec<ScoredSample>,
    /// Pairwise gap quality matrix (intersection/union)
    pub quality_matrix: Vec<Vec<f64>>,
    /// Sample names for matrix indexing
    pub sample_names: Vec<String>,
    /// Reference analysis summary
    pub reference_summary: ReferenceSummary,
    /// Pairwise gap analysis results
    pub pairwise_gap_results: Vec<PairwiseGapResult>,
    /// Gap regions from reference alignment
    pub reference_gaps: Vec<ReferenceGaps>,
    /// Core alignment impact analysis (gap contagion prediction)
    pub impact_analysis: ImpactAnalysis,
    /// Sample clustering and outlier detection
    pub cluster_analysis: ClusterAnalysis,
    /// Gene zone analysis (if GFF provided)
    #[serde(default)]
    pub gene_zone_analysis: Option<GeneZoneSummary>,
    /// Empty column analysis (threshold-based core)
    #[serde(default)]
    pub empty_column_analysis: Option<crate::empty_column::EmptyColumnAnalysis>,
    /// SNP pipeline comparison (if --snippy-dir provided)
    #[serde(default)]
    pub snp_comparison: Option<crate::snp_compare::SnpComparison>,
    /// Snippy-core comparison (if --snippycore-dir provided)
    #[serde(default)]
    pub snippycore_comparison: Option<crate::snp_compare::SnippyCoreComparison>,
    /// CFSAN comparison (if --cfsan-dir provided)
    #[serde(default)]
    pub cfsan_comparison: Option<crate::snp_compare::CfsanComparison>,
    /// VCF comparison between pipelines (if --compare-variants provided)
    #[serde(default)]
    pub vcf_comparison: Option<crate::vcf_compare::VcfComparison>,
}

/// Summary statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Summary {
    pub total_samples: usize,
    pub flagged_samples: usize,
    pub samples_to_review: usize,
    pub samples_to_exclude: usize,
    /// Average reference coverage across samples
    pub avg_ref_coverage: f64,
    /// Minimum reference coverage
    pub min_ref_coverage: f64,
    /// Dataset gap quality score (from pairwise analysis)
    pub dataset_gap_quality: f64,
}

/// Reference alignment summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceSummary {
    pub reference_name: String,
    pub reference_length: usize,
    pub avg_identity: f64,
    pub avg_sample_coverage: f64,
    pub avg_reference_coverage: f64,
}

impl Report {
    /// Create a new report from scored samples
    pub fn new(
        samples: Vec<ScoredSample>,
        pairwise: &PairwiseResults,
        reference: &ReferenceResults,
        impact_analysis: ImpactAnalysis,
        cluster_analysis: ClusterAnalysis,
    ) -> Self {
        let sample_names: Vec<String> = samples.iter().map(|s| s.name.clone()).collect();

        // Build gap quality matrix (intersection/union)
        let n = sample_names.len();
        let mut quality_matrix = vec![vec![1.0; n]; n];

        for (i, name_a) in sample_names.iter().enumerate() {
            for (j, name_b) in sample_names.iter().enumerate() {
                if i != j {
                    if let Some(result) = pairwise.get(name_a, name_b) {
                        quality_matrix[i][j] = result.quality_score;
                    }
                }
            }
        }

        // Calculate summary
        let flagged_samples = samples.iter().filter(|s| s.flagged).count();
        let samples_to_review = samples
            .iter()
            .filter(|s| s.recommendation == crate::scoring::Recommendation::Review)
            .count();
        let samples_to_exclude = samples
            .iter()
            .filter(|s| s.recommendation == crate::scoring::Recommendation::Exclude)
            .count();

        let coverages: Vec<f64> = samples.iter().map(|s| s.metrics.ref_reference_coverage).collect();
        let avg_ref_coverage = if coverages.is_empty() {
            0.0
        } else {
            coverages.iter().sum::<f64>() / coverages.len() as f64
        };
        let min_ref_coverage = coverages.iter().cloned().fold(f64::INFINITY, f64::min);

        let summary = Summary {
            total_samples: samples.len(),
            flagged_samples,
            samples_to_review,
            samples_to_exclude,
            avg_ref_coverage,
            min_ref_coverage: if min_ref_coverage.is_infinite() { 0.0 } else { min_ref_coverage },
            dataset_gap_quality: pairwise.dataset_quality(),
        };

        // Reference summary
        let reference_summary = ReferenceSummary {
            reference_name: reference.reference_name.clone(),
            reference_length: reference.reference_length,
            avg_identity: reference.avg_identity(),
            avg_sample_coverage: reference.avg_sample_coverage(),
            avg_reference_coverage: reference.avg_reference_coverage(),
        };

        // Gap regions
        let pairwise_gap_results = pairwise.results.clone();
        let reference_gaps = reference.get_all_gaps();

        Report {
            version: env!("CARGO_PKG_VERSION").to_string(),
            metadata: ExecutionMetadata::default(),
            input_files: Vec::new(),
            summary,
            samples,
            quality_matrix,
            sample_names,
            reference_summary,
            pairwise_gap_results,
            reference_gaps,
            impact_analysis,
            cluster_analysis,
            gene_zone_analysis: None,
            empty_column_analysis: None,
            snp_comparison: None,
            snippycore_comparison: None,
            cfsan_comparison: None,
            vcf_comparison: None,
        }
    }

    /// Set gene zone analysis
    pub fn with_gene_zone_analysis(mut self, analysis: GeneZoneSummary) -> Self {
        self.gene_zone_analysis = Some(analysis);
        self
    }

    /// Set empty column analysis
    pub fn with_empty_column_analysis(mut self, analysis: crate::empty_column::EmptyColumnAnalysis) -> Self {
        self.empty_column_analysis = Some(analysis);
        self
    }

    /// Set SNP pipeline comparison
    pub fn with_snp_comparison(mut self, comparison: crate::snp_compare::SnpComparison) -> Self {
        self.snp_comparison = Some(comparison);
        self
    }

    /// Set snippy-core comparison
    pub fn with_snippycore_comparison(mut self, comparison: crate::snp_compare::SnippyCoreComparison) -> Self {
        self.snippycore_comparison = Some(comparison);
        self
    }

    /// Set CFSAN comparison
    pub fn with_cfsan_comparison(mut self, comparison: crate::snp_compare::CfsanComparison) -> Self {
        self.cfsan_comparison = Some(comparison);
        self
    }

    /// Set VCF comparison
    pub fn with_vcf_comparison(mut self, comparison: crate::vcf_compare::VcfComparison) -> Self {
        self.vcf_comparison = Some(comparison);
        self
    }

    /// Set execution metadata
    pub fn with_metadata(
        mut self,
        command_line: String,
        timestamp: String,
        duration_secs: f64,
        input_mode: String,
        input_path: String,
        reference_path: String,
    ) -> Self {
        self.metadata = ExecutionMetadata {
            command_line,
            timestamp,
            duration_secs,
            input_mode,
            input_path,
            reference_path,
        };
        self
    }

    /// Set input file information
    pub fn with_input_files(mut self, input_files: Vec<InputFileInfo>) -> Self {
        self.input_files = input_files;
        self
    }
}

/// Write report to file (format determined by extension)
pub fn write_report(report: &Report, path: &Path) -> Result<()> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("json");

    match ext {
        "json" => write_json(report, path),
        "tsv" | "txt" => write_tsv(report, path),
        "html" | "htm" => write_html(report, path),
        _ => write_json(report, path),
    }
}

fn write_json(report: &Report, path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, report)
        .with_context(|| format!("Failed to write JSON to {}", path.display()))?;
    Ok(())
}

fn write_html(report: &Report, path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // Generate HTML with embedded data and interactive visualizations
    let json_data = serde_json::to_string(report)?;

    write!(writer, r##"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>coreguard Report</title>
    <style>
        :root {{
            --low: #28a745;
            --medium: #007bff;
            --high: #ffc107;
            --critical: #dc3545;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #212529;
            --border: #dee2e6;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            padding: 20px;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        h1 {{ color: #2c3e50; margin-bottom: 10px; }}
        h2 {{ color: #34495e; margin: 20px 0 15px; border-bottom: 2px solid var(--border); padding-bottom: 5px; }}
        .subtitle {{ color: #6c757d; margin-bottom: 30px; }}
        .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 30px; }}
        .card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .card-title {{ font-size: 0.85rem; color: #6c757d; text-transform: uppercase; letter-spacing: 0.5px; }}
        .card-value {{ font-size: 2rem; font-weight: bold; margin-top: 5px; }}
        .card-value.low {{ color: var(--low); }}
        .card-value.medium {{ color: var(--medium); }}
        .card-value.high {{ color: var(--high); }}
        .card-value.critical {{ color: var(--critical); }}
        table {{ width: 100%; border-collapse: collapse; background: var(--card-bg); border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        th, td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid var(--border); }}
        th {{ background: #f1f3f4; font-weight: 600; }}
        tr:hover {{ background: #f8f9fa; }}
        .badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.8rem;
            font-weight: 600;
        }}
        .badge-low {{ background: #d4edda; color: #155724; }}
        .badge-medium {{ background: #cce5ff; color: #004085; }}
        .badge-high {{ background: #fff3cd; color: #856404; }}
        .badge-critical {{ background: #f8d7da; color: #721c24; }}
        .badge-flagged {{ background: #f8d7da; color: #721c24; }}
        .badge-ok {{ background: #d4edda; color: #155724; }}
        .recommendations {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px 20px; margin: 20px 0; border-radius: 0 8px 8px 0; }}
        .recommendations li {{ margin: 5px 0; }}
        .chart-container {{ background: var(--card-bg); border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }}
        .heatmap {{ display: grid; gap: 2px; }}
        .heatmap-cell {{
            aspect-ratio: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 0.7rem;
            color: white;
            border-radius: 2px;
            cursor: pointer;
            transition: transform 0.1s;
        }}
        .heatmap-cell:hover {{ transform: scale(1.1); z-index: 1; }}
        .heatmap-label {{ font-size: 0.75rem; color: #666; display: flex; align-items: center; justify-content: center; }}
        .bar-chart {{ display: flex; flex-direction: column; gap: 8px; }}
        .bar-row {{ display: flex; align-items: center; gap: 10px; }}
        .bar-label {{ width: 120px; font-size: 0.85rem; text-align: right; }}
        .bar-container {{ flex: 1; height: 24px; background: #e9ecef; border-radius: 4px; overflow: hidden; }}
        .bar {{ height: 100%; border-radius: 4px; transition: width 0.5s ease; }}
        .bar-value {{ width: 80px; font-size: 0.85rem; font-weight: 600; }}
        .tooltip {{
            position: fixed;
            background: #333;
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 0.85rem;
            pointer-events: none;
            z-index: 1000;
            display: none;
        }}
        .section {{ margin-bottom: 40px; }}
        .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        @media (max-width: 900px) {{ .two-col {{ grid-template-columns: 1fr; }} }}
        footer {{ text-align: center; color: #6c757d; margin-top: 40px; padding-top: 20px; border-top: 1px solid var(--border); }}
    </style>
</head>
<body>
    <div class="container">
        <h1>coreguard Report</h1>
        <p class="subtitle">Pre-alignment QC for SNP Pipeline Analysis | v{version}</p>

        <div class="cards" id="summary-cards"></div>

        <div id="recommendations-section"></div>

        <h2>Sample Risk Assessment</h2>
        <div class="chart-container">
            <div class="bar-chart" id="risk-chart"></div>
        </div>

        <h2>Sample Details</h2>
        <div class="chart-container">
            <table id="samples-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Risk Level</th>
                        <th>Status</th>
                        <th>Core Reduction</th>
                        <th>Unique Gap Bases</th>
                        <th>Avg Identity</th>
                        <th>Recommendation</th>
                    </tr>
                </thead>
                <tbody></tbody>
            </table>
        </div>

        <div class="two-col">
            <div class="section">
                <h2>Identity Heatmap</h2>
                <div class="chart-container">
                    <div class="heatmap" id="identity-heatmap"></div>
                </div>
            </div>
            <div class="section">
                <h2>Core Impact Analysis</h2>
                <div class="chart-container" id="impact-info"></div>
            </div>
        </div>

        <footer>
            Generated by coreguard v{version} | <a href="https://github.com/genpat-it/coreguard">GitHub</a>
        </footer>
    </div>

    <div class="tooltip" id="tooltip"></div>

    <script>
    const data = {json_data};

    // Summary cards
    const summaryCards = document.getElementById('summary-cards');
    const summary = data.summary;
    const flaggedClass = summary.flagged_samples > 0 ? 'critical' : 'low';

    summaryCards.innerHTML = `
        <div class="card">
            <div class="card-title">Total Samples</div>
            <div class="card-value">${{summary.total_samples}}</div>
        </div>
        <div class="card">
            <div class="card-title">Flagged Samples</div>
            <div class="card-value ${{flaggedClass}}">${{summary.flagged_samples}}</div>
        </div>
        <div class="card">
            <div class="card-title">Average Identity</div>
            <div class="card-value">${{(summary.avg_identity * 100).toFixed(1)}}%</div>
        </div>
        <div class="card">
            <div class="card-title">Estimated Core</div>
            <div class="card-value">${{(data.impact_analysis.estimated_core_all / 1e6).toFixed(2)}} Mb</div>
        </div>
    `;

    // Recommendations
    const recsSection = document.getElementById('recommendations-section');
    if (data.cluster_analysis.recommendations.length > 0) {{
        recsSection.innerHTML = `
            <div class="recommendations">
                <strong>Recommendations:</strong>
                <ul>${{data.cluster_analysis.recommendations.map(r => `<li>${{r}}</li>`).join('')}}</ul>
            </div>
        `;
    }}

    // Risk bar chart
    const riskChart = document.getElementById('risk-chart');
    const maxReduction = Math.max(...data.impact_analysis.samples.map(s => s.estimated_core_reduction), 0.01);

    riskChart.innerHTML = data.impact_analysis.samples.map(s => {{
        const pct = (s.estimated_core_reduction / maxReduction) * 100;
        const color = s.risk_level === 'Critical' ? 'var(--critical)' :
                      s.risk_level === 'High' ? 'var(--high)' :
                      s.risk_level === 'Medium' ? 'var(--medium)' : 'var(--low)';
        return `
            <div class="bar-row">
                <div class="bar-label">${{s.sample}}</div>
                <div class="bar-container">
                    <div class="bar" style="width: ${{pct}}%; background: ${{color}};"></div>
                </div>
                <div class="bar-value">${{(s.estimated_core_reduction * 100).toFixed(2)}}%</div>
            </div>
        `;
    }}).join('');

    // Samples table
    const tbody = document.querySelector('#samples-table tbody');
    data.samples.forEach(s => {{
        const impact = data.impact_analysis.samples.find(i => i.sample === s.name) || {{}};
        const riskClass = (impact.risk_level || 'Low').toLowerCase();
        const statusClass = s.flagged ? 'flagged' : 'ok';
        const row = document.createElement('tr');
        row.innerHTML = `
            <td><strong>${{s.name}}</strong></td>
            <td><span class="badge badge-${{riskClass}}">${{impact.risk_level || 'Low'}}</span></td>
            <td><span class="badge badge-${{statusClass}}">${{s.flagged ? 'Flagged' : 'OK'}}</span></td>
            <td>${{((impact.estimated_core_reduction || 0) * 100).toFixed(2)}}%</td>
            <td>${{(impact.unique_gap_bases || 0).toLocaleString()}} bp</td>
            <td>${{(s.metrics.avg_identity * 100).toFixed(1)}}%</td>
            <td>${{s.recommendation}}</td>
        `;
        tbody.appendChild(row);
    }});

    // Identity heatmap
    const heatmap = document.getElementById('identity-heatmap');
    const n = data.sample_names.length;
    heatmap.style.gridTemplateColumns = `80px repeat(${{n}}, 1fr)`;

    // Header row
    heatmap.innerHTML = '<div class="heatmap-label"></div>' +
        data.sample_names.map(name => `<div class="heatmap-label">${{name.substring(0, 8)}}</div>`).join('');

    // Data rows
    for (let i = 0; i < n; i++) {{
        heatmap.innerHTML += `<div class="heatmap-label">${{data.sample_names[i].substring(0, 8)}}</div>`;
        for (let j = 0; j < n; j++) {{
            const val = data.identity_matrix[i][j];
            const pct = ((val - 0.9) / 0.1) * 100;  // Scale 0.9-1.0 to 0-100%
            const hue = Math.max(0, Math.min(120, pct * 1.2));  // Red to green
            const bg = i === j ? '#333' : `hsl(${{hue}}, 70%, 50%)`;
            heatmap.innerHTML += `
                <div class="heatmap-cell" style="background: ${{bg}};"
                     data-tooltip="${{data.sample_names[i]}} vs ${{data.sample_names[j]}}: ${{(val * 100).toFixed(2)}}%">
                    ${{i === j ? '-' : (val * 100).toFixed(1)}}
                </div>
            `;
        }}
    }}

    // Impact info
    const impactInfo = document.getElementById('impact-info');
    const impact = data.impact_analysis;
    impactInfo.innerHTML = `
        <p><strong>Reference Length:</strong> ${{(impact.reference_length / 1e6).toFixed(2)}} Mb</p>
        <p><strong>Estimated Core (all samples):</strong> ${{(impact.estimated_core_all / 1e6).toFixed(2)}} Mb (${{((1 - impact.total_core_reduction) * 100).toFixed(1)}}%)</p>
        <p><strong>Total Core Reduction:</strong> ${{(impact.total_core_reduction * 100).toFixed(2)}}%</p>
        <hr style="margin: 15px 0; border: none; border-top: 1px solid var(--border);">
        <p><strong>Risk Distribution:</strong></p>
        <ul>
            <li style="color: var(--low);">Low: ${{impact.samples.filter(s => s.risk_level === 'Low').length}}</li>
            <li style="color: var(--medium);">Medium: ${{impact.samples.filter(s => s.risk_level === 'Medium').length}}</li>
            <li style="color: var(--high);">High: ${{impact.samples.filter(s => s.risk_level === 'High').length}}</li>
            <li style="color: var(--critical);">Critical: ${{impact.samples.filter(s => s.risk_level === 'Critical').length}}</li>
        </ul>
    `;

    // Tooltip
    const tooltip = document.getElementById('tooltip');
    document.querySelectorAll('[data-tooltip]').forEach(el => {{
        el.addEventListener('mouseenter', e => {{
            tooltip.textContent = e.target.dataset.tooltip;
            tooltip.style.display = 'block';
        }});
        el.addEventListener('mousemove', e => {{
            tooltip.style.left = (e.clientX + 10) + 'px';
            tooltip.style.top = (e.clientY + 10) + 'px';
        }});
        el.addEventListener('mouseleave', () => {{
            tooltip.style.display = 'none';
        }});
    }});
    </script>
</body>
</html>
"##, version = report.version, json_data = json_data)?;

    Ok(())
}

/// Write BED files for gap visualization in IGV
pub fn write_bed_files(report: &Report, output_dir: &Path) -> Result<()> {
    std::fs::create_dir_all(output_dir)?;

    // 1. Pairwise gap quality BED (per sample pair)
    for pg in &report.pairwise_gap_results {
        let bed_path = output_dir.join(format!("pairwise_{}_{}.bed", pg.sample_a, pg.sample_b));
        let file = File::create(&bed_path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "track name=\"{}_vs_{}_quality\" description=\"Gap quality: {:.2}\"",
                 pg.sample_a, pg.sample_b, pg.quality_score)?;

        // Write a single entry summarizing the pair
        writeln!(writer, "# Gap union: {} bp, intersection: {} bp, quality: {:.2}",
                 pg.gap_union_bases, pg.gap_intersection_bases, pg.quality_score)?;
        writeln!(writer, "# Unique gaps in {}: {} bp, in {}: {} bp",
                 pg.sample_a, pg.unique_gap_a_bases, pg.sample_b, pg.unique_gap_b_bases)?;
    }

    // 2. Reference gaps BED (per sample)
    for rg in &report.reference_gaps {
        let bed_path = output_dir.join(format!("reference_{}.bed", rg.sample));
        let file = File::create(&bed_path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "track name=\"{}_vs_ref_gaps\" description=\"Reference alignment gaps\" color=255,0,0",
                 rg.sample)?;

        // Sample regions not aligned to reference
        for region in &rg.sample_unaligned {
            writeln!(writer, "{}\t{}\t{}\tsample_unaligned\t0\t+",
                     rg.sample, region.start, region.end)?;
        }

        // Reference regions not covered by sample (use reference name)
        let ref_bed_path = output_dir.join(format!("reference_uncovered_{}.bed", rg.sample));
        let ref_file = File::create(&ref_bed_path)?;
        let mut ref_writer = BufWriter::new(ref_file);

        writeln!(ref_writer, "track name=\"ref_uncovered_by_{}\" description=\"Reference regions not covered\" color=0,0,255",
                 rg.sample)?;

        let ref_name = &report.reference_summary.reference_name;

        for region in &rg.reference_uncovered {
            writeln!(ref_writer, "{}\t{}\t{}\tuncovered_by_{}\t0\t+",
                     ref_name, region.start, region.end, rg.sample)?;
        }
    }

    // 3. Impact analysis - unique gaps BED (most important!)
    let unique_bed_path = output_dir.join("unique_gaps_all.bed");
    let file = File::create(&unique_bed_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "track name=\"unique_gaps\" description=\"Sample-specific gaps (cause core reduction)\" color=255,0,0")?;

    for sample_impact in &report.impact_analysis.samples {
        for region in &sample_impact.unique_gap_coords {
            writeln!(writer, "{}\t{}\t{}\tunique_gap_{}\t{}\t+",
                     sample_impact.sample,
                     region.start,
                     region.end,
                     sample_impact.sample,
                     sample_impact.unique_gap_bases)?;
        }
    }

    // 4. Summary BED with all unique gaps on reference coordinates
    {
        let ref_summary = &report.reference_summary;
        let ref_bed_path = output_dir.join("unique_gaps_reference.bed");
        let file = File::create(&ref_bed_path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "track name=\"unique_gaps_on_ref\" description=\"Unique gaps mapped to reference\" itemRgb=On")?;

        // Color by risk level
        for sample_impact in &report.impact_analysis.samples {
            let color = match sample_impact.risk_level {
                crate::impact::RiskLevel::Critical => "255,0,0",     // Red
                crate::impact::RiskLevel::High => "255,165,0",      // Orange
                crate::impact::RiskLevel::Medium => "255,255,0",    // Yellow
                crate::impact::RiskLevel::Low => "0,255,0",         // Green
            };
            let _ = ref_summary; // suppress unused warning

            for region in &sample_impact.unique_gap_coords {
                writeln!(writer, "{}\t{}\t{}\t{}\t0\t+\t{}\t{}\t{}",
                         ref_summary.reference_name,
                         region.start,
                         region.end,
                         sample_impact.sample,
                         region.start,
                         region.end,
                         color)?;
            }
        }
    }

    log::info!("BED files written to: {}", output_dir.display());
    Ok(())
}

fn write_tsv(report: &Report, path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // Header
    writeln!(
        writer,
        "sample\tscore\tflagged\trecommendation\tflags\tref_coverage\tgap_quality\tunique_gaps\tn50\tnum_contigs"
    )?;

    // Data rows
    for s in &report.samples {
        let flags_str = s
            .flags
            .iter()
            .map(|f| f.to_string())
            .collect::<Vec<_>>()
            .join(",");

        writeln!(
            writer,
            "{}\t{:.4}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
            s.name,
            s.score,
            s.flagged,
            s.recommendation,
            if flags_str.is_empty() { "-" } else { &flags_str },
            s.metrics.ref_reference_coverage,
            s.metrics.avg_gap_quality,
            s.metrics.total_unique_gaps,
            s.metrics.n50,
            s.metrics.num_contigs,
        )?;
    }

    // Gap quality matrix section
    writeln!(writer)?;
    writeln!(writer, "# Gap Quality Matrix (intersection/union)")?;
    write!(writer, "sample")?;
    for name in &report.sample_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;

    for (i, name) in report.sample_names.iter().enumerate() {
        write!(writer, "{}", name)?;
        for j in 0..report.sample_names.len() {
            write!(writer, "\t{:.4}", report.quality_matrix[i][j])?;
        }
        writeln!(writer)?;
    }

    // Reference summary
    let ref_sum = &report.reference_summary;
    writeln!(writer)?;
    writeln!(writer, "# Reference Summary")?;
    writeln!(writer, "reference_name\t{}", ref_sum.reference_name)?;
    writeln!(writer, "reference_length\t{}", ref_sum.reference_length)?;
    writeln!(writer, "avg_identity\t{:.4}", ref_sum.avg_identity)?;
    writeln!(writer, "avg_sample_coverage\t{:.4}", ref_sum.avg_sample_coverage)?;
    writeln!(writer, "avg_reference_coverage\t{:.4}", ref_sum.avg_reference_coverage)?;

    Ok(())
}
