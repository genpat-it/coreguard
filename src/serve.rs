//! Web server for interactive IGV.js genome browser
//!
//! Serves the coreguard report with embedded IGV.js for visualizing
//! gaps across all samples.

use crate::fasta::Sample;
use crate::gaps::Region;
use crate::output::Report;
use crate::pairwise::PairwiseResults;
use crate::vcf_compare::VcfComparison;
use anyhow::Result;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use tiny_http::{Header, Response, Server};

/// Start the web server with IGV.js genome browser
pub fn start_server(
    report: &Report,
    reference: Option<&Sample>,
    bed_dir: Option<&Path>,
    port: u16,
    pairwise: &PairwiseResults,
    vcf_comparison: Option<&VcfComparison>,
    snippy_gaps: &HashMap<String, Vec<Region>>,
    cfsan_gaps: &HashMap<String, Vec<Region>>,
    bam_dir: Option<&Path>,
) -> Result<()> {
    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr).map_err(|e| anyhow::anyhow!("Failed to start server: {}", e))?;

    let url = format!("http://localhost:{}", port);
    log::info!("Server running at {}", url);
    log::info!("Press Ctrl+C to stop");

    // Open browser
    if let Err(e) = webbrowser::open(&url) {
        log::warn!("Could not open browser: {}. Please open {} manually.", e, url);
    }

    // Prepare data
    let report_json = serde_json::to_string(report)?;

    // Get reference name from report
    let reference_name = report.reference_summary.reference_name.clone();
    let _ = pairwise; // For future use - pairwise gap data

    // Generate FASTA/FAI using just the chromosome ID (not full description)
    // This is critical for BAM files to match - BAM headers use short names like "AL591824.1"
    let (reference_fasta, reference_fai, reference_name_for_igv) = if let Some(r) = reference {
        let contig_name = r.contigs.first()
            .map(|c| c.name.clone())
            .unwrap_or_else(|| reference_name.clone());
        // Extract just chromosome ID (before first space)
        // e.g., "AL591824.1 Listeria monocytogenes..." -> "AL591824.1"
        let chrom_id = contig_name
            .split_whitespace()
            .next()
            .unwrap_or(&contig_name)
            .to_string();
        (
            r.to_fasta_string_with_name(&chrom_id),
            r.to_fai_string_with_name(&chrom_id),
            chrom_id,
        )
    } else {
        (String::new(), String::new(), reference_name.clone())
    };

    // Load BED files if available
    let mut bed_files: HashMap<String, String> = HashMap::new();
    if let Some(dir) = bed_dir {
        if dir.exists() {
            for entry in std::fs::read_dir(dir)? {
                let entry = entry?;
                let path = entry.path();
                if path.extension().map(|e| e == "bed").unwrap_or(false) {
                    let name = path.file_name().unwrap().to_string_lossy().to_string();
                    let content = std::fs::read_to_string(&path)?;
                    bed_files.insert(name, content);
                }
            }
        }
    }

    // Generate BED files from report if not provided
    if bed_files.is_empty() {
        bed_files = generate_bed_files_from_report(report, snippy_gaps, cfsan_gaps);
    }

    // Load BAM files if available (store paths, not content - BAM files are large)
    let mut bam_files: HashMap<String, std::path::PathBuf> = HashMap::new();
    if let Some(dir) = bam_dir {
        if dir.exists() {
            // Recursively find BAM files
            fn find_bam_files(dir: &Path, bam_files: &mut HashMap<String, std::path::PathBuf>) {
                if let Ok(entries) = std::fs::read_dir(dir) {
                    for entry in entries.filter_map(|e| e.ok()) {
                        let path = entry.path();
                        if path.is_dir() {
                            find_bam_files(&path, bam_files);
                        } else if path.extension().map(|e| e == "bam").unwrap_or(false) {
                            // Extract sample name from path
                            // Try patterns like: sample/reads.sorted.bam or sample.bam
                            let sample_name = if let Some(parent) = path.parent() {
                                if let Some(parent_name) = parent.file_name() {
                                    parent_name.to_string_lossy().to_string()
                                } else {
                                    path.file_stem().unwrap_or_default().to_string_lossy().to_string()
                                }
                            } else {
                                path.file_stem().unwrap_or_default().to_string_lossy().to_string()
                            };

                            // Prefer sorted/deduped BAM files
                            let filename = path.file_name().unwrap_or_default().to_string_lossy();
                            let priority = if filename.contains("indelrealigned") {
                                4
                            } else if filename.contains("deduped") {
                                3
                            } else if filename.contains("sorted") && !filename.contains("unsorted") {
                                2
                            } else {
                                1
                            };

                            // Only replace if this file has higher priority
                            let key = sample_name.clone();
                            if let Some(existing) = bam_files.get(&key) {
                                let existing_name = existing.file_name().unwrap_or_default().to_string_lossy();
                                let existing_priority = if existing_name.contains("indelrealigned") {
                                    4
                                } else if existing_name.contains("deduped") {
                                    3
                                } else if existing_name.contains("sorted") && !existing_name.contains("unsorted") {
                                    2
                                } else {
                                    1
                                };
                                if priority > existing_priority {
                                    bam_files.insert(key, path);
                                }
                            } else {
                                bam_files.insert(key, path);
                            }
                        }
                    }
                }
            }
            find_bam_files(dir, &mut bam_files);
            if !bam_files.is_empty() {
                log::info!("Found {} BAM files for IGV visualization", bam_files.len());
                for (name, path) in &bam_files {
                    log::debug!("  {} -> {}", name, path.display());
                }
            }
        }
    }

    let bed_files = Arc::new(bed_files);
    let bam_files = Arc::new(bam_files);
    let report_json = Arc::new(report_json);
    let reference_fasta = Arc::new(reference_fasta);
    let reference_fai = Arc::new(reference_fai);
    let reference_name = Arc::new(reference_name);
    let reference_name_for_igv = Arc::new(reference_name_for_igv);
    let reference_length = report.reference_summary.reference_length;
    let vcf_comparison = vcf_comparison.cloned().map(Arc::new);
    let snippy_gaps = Arc::new(snippy_gaps.clone());
    let cfsan_gaps = Arc::new(cfsan_gaps.clone());

    // Handle requests
    for request in server.incoming_requests() {
        let path = request.url().to_string();

        let response = match path.as_str() {
            "/" | "/index.html" => {
                let html = generate_igv_html(&report_json, &reference_name, reference_length, &bed_files);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            "/report.json" => {
                Response::from_string(report_json.as_str())
                    .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
            }
            "/reference.fasta" | "/reference.fa" => {
                Response::from_string(reference_fasta.as_str())
                    .with_header(Header::from_bytes("Content-Type", "text/plain").unwrap())
            }
            "/reference.fasta.fai" | "/reference.fa.fai" => {
                Response::from_string(reference_fai.as_str())
                    .with_header(Header::from_bytes("Content-Type", "text/plain").unwrap())
            }
            path if path.starts_with("/bed/") => {
                let bed_name = &path[5..]; // Remove "/bed/"
                if let Some(content) = bed_files.get(bed_name) {
                    Response::from_string(content.as_str())
                        .with_header(Header::from_bytes("Content-Type", "text/plain").unwrap())
                } else {
                    Response::from_string("Not found").with_status_code(404)
                }
            }
            "/bed_list" => {
                let list: Vec<&String> = bed_files.keys().collect();
                let json = serde_json::to_string(&list).unwrap_or_default();
                Response::from_string(json)
                    .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
            }
            "/bam_list" => {
                let list: Vec<&String> = bam_files.keys().collect();
                let json = serde_json::to_string(&list).unwrap_or_default();
                Response::from_string(json)
                    .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
            }
            path if path.starts_with("/bam/") => {
                // Serve BAM or BAI files: /bam/sample.bam or /bam/sample.bam.bai
                let file_path = &path[5..]; // Remove "/bam/"
                let is_bai = file_path.ends_with(".bai");
                let sample_name = if is_bai {
                    file_path.trim_end_matches(".bam.bai")
                } else {
                    file_path.trim_end_matches(".bam")
                };

                if let Some(bam_path) = bam_files.get(sample_name) {
                    let file_to_serve = if is_bai {
                        // Try .bam.bai first, then .bai
                        // with_extension replaces from last dot, so we need to append
                        let mut bai_path = bam_path.clone();
                        bai_path.set_extension("bam.bai");
                        if bai_path.exists() {
                            Some(bai_path)
                        } else {
                            // Try just adding .bai
                            let alt_bai = std::path::PathBuf::from(format!("{}.bai", bam_path.display()));
                            if alt_bai.exists() {
                                Some(alt_bai)
                            } else {
                                None
                            }
                        }
                    } else {
                        Some(bam_path.clone())
                    };

                    if let Some(file_to_serve) = file_to_serve {
                        match std::fs::read(&file_to_serve) {
                            Ok(content) => {
                                Response::from_data(content)
                                    .with_header(Header::from_bytes("Content-Type", "application/octet-stream").unwrap())
                                    .with_header(Header::from_bytes("Access-Control-Allow-Origin", "*").unwrap())
                            }
                            Err(e) => {
                                log::error!("Failed to read {}: {}", file_to_serve.display(), e);
                                Response::from_string("File read error").with_status_code(500)
                            }
                        }
                    } else {
                        Response::from_string("BAI index not found").with_status_code(404)
                    }
                } else {
                    Response::from_string("BAM not found").with_status_code(404)
                }
            }
            "/gaps.svg" => {
                let svg = generate_gaps_svg(report, reference_length);
                Response::from_string(svg)
                    .with_header(Header::from_bytes("Content-Type", "image/svg+xml").unwrap())
            }
            "/pairwise_gaps.json" => {
                // Return pairwise gap analysis results
                let json = serde_json::to_string_pretty(&report.pairwise_gap_results)
                    .unwrap_or_else(|_| "[]".to_string());
                Response::from_string(json)
                    .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
            }
            path if path.starts_with("/api/vcf_positions") => {
                // Paginated VCF positions endpoint
                // Parse query params: ?page=0&page_size=50&sample=all&status=all
                if let Some(ref vcf_comp) = vcf_comparison {
                    let params = parse_query_params(path);
                    let page: usize = params.get("page").and_then(|s| s.parse().ok()).unwrap_or(0);
                    let page_size: usize = params.get("page_size").and_then(|s| s.parse().ok()).unwrap_or(50);
                    let sample = params.get("sample").map(|s| s.as_str());
                    let status = params.get("status").map(|s| s.as_str());

                    let paginated = vcf_comp.get_paginated_positions(page, page_size, sample, status);
                    let json = serde_json::to_string(&paginated).unwrap_or_else(|_| "{}".to_string());
                    Response::from_string(json)
                        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
                        .with_header(Header::from_bytes("Access-Control-Allow-Origin", "*").unwrap())
                } else {
                    Response::from_string(r#"{"error": "No VCF comparison data available"}"#)
                        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
                        .with_status_code(404)
                }
            }
            path if path.starts_with("/api/region") => {
                // Region endpoint: returns reference sequence + SNPs for a specific region
                // ?start=1000&end=2000
                let params = parse_query_params(path);
                let start: usize = params.get("start").and_then(|s| s.parse().ok()).unwrap_or(1);
                let end: usize = params.get("end").and_then(|s| s.parse().ok()).unwrap_or(start + 1000);
                let end = end.min(start + 10000); // Max 10kb per request

                // Parse reference sequence from FASTA
                let mut ref_seq = String::new();
                let fasta_lines: Vec<&str> = reference_fasta.lines().collect();
                let mut seq_start = false;
                for line in &fasta_lines {
                    if line.starts_with('>') {
                        seq_start = true;
                        continue;
                    }
                    if seq_start {
                        ref_seq.push_str(line.trim());
                    }
                }

                // Extract the region
                let region_seq: String = if start > 0 && start <= ref_seq.len() {
                    let s = start - 1; // Convert to 0-based
                    let e = end.min(ref_seq.len());
                    ref_seq.chars().skip(s).take(e - s).collect()
                } else {
                    String::new()
                };

                // Get SNPs in this region
                let mut snps_in_region: Vec<serde_json::Value> = Vec::new();
                if let Some(ref vcf_comp) = vcf_comparison {
                    for pos_comp in &vcf_comp.discordant_positions {
                        if pos_comp.pos >= start && pos_comp.pos <= end {
                            let status_str = format!("{:?}", pos_comp.status);
                            let ref_allele = pos_comp.pipeline_a.as_ref().map(|p| p.ref_allele.clone())
                                .or_else(|| pos_comp.pipeline_b.as_ref().map(|p| p.ref_allele.clone()));
                            let alt_allele = pos_comp.pipeline_a.as_ref().map(|p| p.alt_allele.clone())
                                .or_else(|| pos_comp.pipeline_b.as_ref().map(|p| p.alt_allele.clone()));

                            // Build quality info for Snippy (pipeline_a)
                            let snippy_qual = pos_comp.pipeline_a.as_ref().map(|p| {
                                serde_json::json!({
                                    "qual": p.qual,
                                    "depth": p.depth,
                                    "ref_obs": p.ref_obs,
                                    "alt_obs": p.alt_obs,
                                    "alt_qual": p.alt_qual,
                                    "filter": p.filter,
                                    "genotype": p.genotype,
                                    "allele_freq": p.allele_freq
                                })
                            });

                            // Build quality info for CFSAN (pipeline_b)
                            let cfsan_qual = pos_comp.pipeline_b.as_ref().map(|p| {
                                serde_json::json!({
                                    "qual": p.qual,
                                    "depth": p.depth,
                                    "ref_obs": p.ref_obs,
                                    "alt_obs": p.alt_obs,
                                    "genotype_qual": p.genotype_qual,
                                    "avg_base_qual": p.avg_base_qual,
                                    "pvalue": p.pvalue,
                                    "filter": p.filter,
                                    "genotype": p.genotype,
                                    "allele_freq": p.allele_freq
                                })
                            });

                            snps_in_region.push(serde_json::json!({
                                "pos": pos_comp.pos,
                                "sample": pos_comp.sample,
                                "status": status_str,
                                "ref": ref_allele,
                                "alt": alt_allele,
                                "snippy": snippy_qual,
                                "cfsan": cfsan_qual
                            }));
                        }
                    }
                }

                // Get sample names from report
                let sample_names: Vec<String> = report.samples.iter().map(|s| s.name.clone()).collect();

                // Get ALL gap regions for each sample in this region (from mm2 coverage analysis)
                let mut gaps_in_region: Vec<serde_json::Value> = Vec::new();
                for ref_gap in &report.reference_gaps {
                    for gap in &ref_gap.reference_uncovered {
                        // Check if gap overlaps with requested region
                        if gap.end >= start && gap.start <= end {
                            gaps_in_region.push(serde_json::json!({
                                "sample": ref_gap.sample,
                                "start": gap.start.max(start),
                                "end": gap.end.min(end),
                                "type": "mm2"
                            }));
                        }
                    }
                }

                // Add Snippy gaps (N positions from aligned FASTAs)
                for (sample_name, sample_gaps) in snippy_gaps.iter() {
                    for gap in sample_gaps {
                        if gap.end >= start && gap.start <= end {
                            gaps_in_region.push(serde_json::json!({
                                "sample": sample_name,
                                "start": gap.start.max(start),
                                "end": gap.end.min(end),
                                "type": "snippy"
                            }));
                        }
                    }
                }

                // Add CFSAN gaps (zero coverage from pileup)
                for (sample_name, sample_gaps) in cfsan_gaps.iter() {
                    for gap in sample_gaps {
                        if gap.end >= start && gap.start <= end {
                            gaps_in_region.push(serde_json::json!({
                                "sample": sample_name,
                                "start": gap.start.max(start),
                                "end": gap.end.min(end),
                                "type": "cfsan"
                            }));
                        }
                    }
                }

                let response_json = serde_json::json!({
                    "start": start,
                    "end": end,
                    "length": region_seq.len(),
                    "reference": region_seq,
                    "snps": snps_in_region,
                    "samples": sample_names,
                    "gaps": gaps_in_region
                });

                Response::from_string(serde_json::to_string(&response_json).unwrap_or_default())
                    .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
                    .with_header(Header::from_bytes("Access-Control-Allow-Origin", "*").unwrap())
            }
            "/3d" => {
                // 3D visualization page using Three.js
                let html = generate_3d_html(&report_json, &reference_name, reference_length);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            "/igv" => {
                // Full-power standalone IGV.js genome browser
                // Use actual contig name from FASTA (must match BAM chromosome names)
                let html = generate_full_igv_html(&report_json, &reference_name_for_igv, reference_length, &bed_files, &bam_files);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            "/compare" => {
                // 3-way comparison visualization: mm2 gaps vs Snippy vs CFSAN
                let html = generate_compare_html(&report_json, &reference_name_for_igv, reference_length);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            p if p == "/nucleotide" || p.starts_with("/nucleotide?") => {
                // Simple nucleotide-level alignment viewer
                let params = parse_query_params(p);
                let start: usize = params.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
                let width: usize = params.get("width").and_then(|s| s.parse().ok()).unwrap_or(100);
                let html = generate_nucleotide_html(&reference_fasta, report, start, width);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            "/wasm/coreguard_wasm.js" => {
                // Serve WASM glue code
                let js = include_str!("../wasm/pkg/coreguard_wasm.js");
                Response::from_string(js)
                    .with_header(Header::from_bytes("Content-Type", "application/javascript").unwrap())
            }
            "/wasm/coreguard_wasm_bg.wasm" => {
                // Serve WASM binary
                let wasm = include_bytes!("../wasm/pkg/coreguard_wasm_bg.wasm");
                Response::from_data(wasm.to_vec())
                    .with_header(Header::from_bytes("Content-Type", "application/wasm").unwrap())
            }
            p if p == "/nucleotide-wasm" || p.starts_with("/nucleotide-wasm?") => {
                // WASM-powered nucleotide viewer
                let html = generate_nucleotide_wasm_html(&reference_fasta, report);
                Response::from_string(html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            _ => {
                Response::from_string("Not found").with_status_code(404)
            }
        };

        if let Err(e) = request.respond(response) {
            log::error!("Failed to send response: {}", e);
        }
    }

    Ok(())
}

/// Parse query parameters from URL path
fn parse_query_params(path: &str) -> HashMap<String, String> {
    let mut params = HashMap::new();
    if let Some(query) = path.split('?').nth(1) {
        for pair in query.split('&') {
            if let Some((key, value)) = pair.split_once('=') {
                // Simple URL decoding for common cases
                let decoded = value
                    .replace("%20", " ")
                    .replace("%2F", "/")
                    .replace("%3A", ":")
                    .replace("+", " ");
                params.insert(key.to_string(), decoded);
            }
        }
    }
    params
}

/// Generate SVG visualization of gaps distribution across the genome
fn generate_gaps_svg(report: &Report, reference_length: usize) -> String {
    let width = 1400usize;  // Wider to accommodate longer labels
    let margin = 200usize;  // More margin for sample names + gap counts
    let track_height = 25usize;
    let track_gap = 5usize;

    // Sort samples by unique_gap_bases (best first)
    let mut sorted_samples: Vec<_> = report.impact_analysis.samples.iter().collect();
    sorted_samples.sort_by_key(|s| s.unique_gap_bases);

    let num_tracks = sorted_samples.len() + 1; // +1 for aggregated track
    let density_legend_space = 30usize; // Space for density legend below ALL GAPS
    let height = margin * 2 + (track_height + track_gap) * num_tracks + density_legend_space + 80; // +80 for legend and title

    let scale = (width - margin * 2) as f64 / reference_length as f64;

    let ref_name = &report.reference_summary.reference_name;

    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" style="font-family: Arial, sans-serif;">
  <defs>
    <style>
      .title {{ font-size: 16px; font-weight: bold; }}
      .label {{ font-size: 11px; }}
      .axis-label {{ font-size: 10px; fill: #666; }}
      .track-label {{ font-size: 11px; font-weight: 500; }}
    </style>
  </defs>
  <rect width="100%" height="100%" fill="#fafafa"/>
  <text x="{}" y="25" class="title" text-anchor="middle">Unique Gaps Distribution - Core Genome Impact</text>
  <text x="{}" y="42" class="label" text-anchor="middle" fill="#666">Reference: {} ({:.2} Mb) | Core: {:.2} Mb ({:.1}%)</text>
"##,
        width, height,
        width / 2,
        width / 2,
        ref_name,
        reference_length as f64 / 1e6,
        report.impact_analysis.estimated_core_all as f64 / 1e6,
        (1.0 - report.impact_analysis.total_core_reduction) * 100.0
    );

    // Genome axis
    let axis_y = margin + 15;
    svg.push_str(&format!(
        "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"2\"/>\n",
        margin, axis_y, width - margin, axis_y
    ));

    // Axis ticks and labels (every 500kb)
    let tick_interval = 500_000usize;
    let mut pos = 0usize;
    while pos <= reference_length {
        let x = margin as f64 + pos as f64 * scale;
        svg.push_str(&format!(
            "  <line x1=\"{:.1}\" y1=\"{}\" x2=\"{:.1}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1\"/>\n",
            x, axis_y - 5, x, axis_y + 5
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{}\" class=\"axis-label\" text-anchor=\"middle\">{:.1}M</text>\n",
            x, axis_y + 18, pos as f64 / 1e6
        ));
        pos += tick_interval;
    }

    let tracks_start_y = axis_y + 35;

    // Calculate total unique gaps for ALL GAPS track
    let total_gap_bases: usize = sorted_samples.iter().map(|s| s.unique_gap_bases).sum();
    let total_gap_regions: usize = sorted_samples.iter().map(|s| s.unique_gap_coords.len()).sum();
    let total_gap_pct = (total_gap_bases as f64 / reference_length as f64) * 100.0;

    // Aggregated track (ALL GAPS) with total count
    let all_gaps_label = format!("ALL GAPS ({} regions, {:.1}kb, {:.2}%)",
        total_gap_regions, total_gap_bases as f64 / 1000.0, total_gap_pct);
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"track-label\" text-anchor=\"end\">{}</text>\n",
        margin - 5, tracks_start_y + track_height / 2 + 4, all_gaps_label
    ));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#e0e0e0\" rx=\"2\"/>\n",
        margin, tracks_start_y, width - margin * 2, track_height
    ));

    // DENSITY PLOT for ALL GAPS track
    // Divide genome into bins and calculate gap density per bin
    let track_width_px = width - margin * 2;
    let num_bins = track_width_px as usize; // 1 bin per pixel for smooth density
    let bin_size = (reference_length as f64 / num_bins as f64).ceil() as usize;

    // Count gap bases in each bin
    let mut bin_counts: Vec<usize> = vec![0; num_bins];
    for sample in &sorted_samples {
        for region in &sample.unique_gap_coords {
            let start_bin = region.start / bin_size;
            let end_bin = (region.end.saturating_sub(1)) / bin_size;
            for bin in start_bin..=end_bin.min(num_bins - 1) {
                // Calculate overlap of gap with this bin
                let bin_start = bin * bin_size;
                let bin_end = (bin + 1) * bin_size;
                let overlap_start = region.start.max(bin_start);
                let overlap_end = region.end.min(bin_end);
                if overlap_end > overlap_start {
                    bin_counts[bin] += overlap_end - overlap_start;
                }
            }
        }
    }

    // Find max density for normalization
    let max_density = bin_counts.iter().max().copied().unwrap_or(1) as f64;

    // Draw density bars with color gradient (white -> yellow -> orange -> red)
    for (bin, &count) in bin_counts.iter().enumerate() {
        if count == 0 {
            continue; // Skip empty bins
        }
        let density = count as f64 / bin_size as f64; // 0.0 to 1.0
        let normalized = (count as f64 / max_density).min(1.0);

        // Color gradient: low density = yellow, high density = red
        let (r, g, b) = if normalized < 0.5 {
            // Yellow to orange (255,255,0) -> (255,165,0)
            let t = normalized * 2.0;
            (255, (255.0 - 90.0 * t) as u8, 0)
        } else {
            // Orange to red (255,165,0) -> (220,53,69)
            let t = (normalized - 0.5) * 2.0;
            ((255.0 - 35.0 * t) as u8, (165.0 - 112.0 * t) as u8, (69.0 * t) as u8)
        };

        let x = margin as i32 + bin as i32;
        svg.push_str(&format!(
            "  <rect x=\"{}\" y=\"{}\" width=\"1\" height=\"{}\" fill=\"rgb({},{},{})\"/>\n",
            x, tracks_start_y, track_height, r, g, b
        ));
    }

    // Add density scale legend below ALL GAPS track
    let legend_x = width as i32 - margin as i32 - 150;
    let legend_y = tracks_start_y + track_height + 3;
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"axis-label\">Density:</text>\n",
        legend_x - 45, legend_y + 8
    ));
    // Draw mini gradient
    for i in 0..50 {
        let normalized = i as f64 / 49.0;
        let (r, g, b) = if normalized < 0.5 {
            let t = normalized * 2.0;
            (255, (255.0 - 90.0 * t) as u8, 0)
        } else {
            let t = (normalized - 0.5) * 2.0;
            ((255.0 - 35.0 * t) as u8, (165.0 - 112.0 * t) as u8, (69.0 * t) as u8)
        };
        svg.push_str(&format!(
            "  <rect x=\"{}\" y=\"{}\" width=\"2\" height=\"10\" fill=\"rgb({},{},{})\"/>\n",
            legend_x + i * 2, legend_y, r, g, b
        ));
    }
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"axis-label\">Low</text>\n",
        legend_x, legend_y + 20
    ));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"axis-label\">High</text>\n",
        legend_x + 85, legend_y + 20
    ));

    // Per-sample tracks (add extra spacing after ALL GAPS track for density legend)
    let density_legend_height = 30; // Extra space for density legend
    for (idx, sample) in sorted_samples.iter().enumerate() {
        let y = tracks_start_y + track_height + density_legend_height + (track_height + track_gap) * idx;

        let color = match sample.risk_level {
            crate::impact::RiskLevel::Critical => "#d32f2f",
            crate::impact::RiskLevel::High => "#f57c00",
            crate::impact::RiskLevel::Medium => "#1976d2",
            crate::impact::RiskLevel::Low => "#388e3c",
        };

        // Track label with gap count and background
        let gap_count = sample.unique_gap_coords.len();
        let gap_bp = sample.unique_gap_bases;
        let label_text = if gap_count > 0 {
            format!("{} ({} gaps, {}kb)", sample.sample, gap_count, gap_bp / 1000)
        } else {
            format!("{} (no gaps)", sample.sample)
        };
        svg.push_str(&format!(
            "  <text x=\"{}\" y=\"{}\" class=\"track-label\" text-anchor=\"end\" fill=\"{}\">{}</text>\n",
            margin - 5, y + track_height / 2 + 4, color, label_text
        ));
        svg.push_str(&format!(
            "  <rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#f5f5f5\" rx=\"2\"/>\n",
            margin, y, width - margin * 2, track_height
        ));

        // Draw gaps for this sample (0.5px minimum for sub-pixel gaps)
        for region in &sample.unique_gap_coords {
            let x = margin as f64 + region.start as f64 * scale;
            let w = ((region.end - region.start) as f64 * scale).max(0.5);
            svg.push_str(&format!(
                "  <rect x=\"{:.1}\" y=\"{}\" width=\"{:.1}\" height=\"{}\" fill=\"{}\" rx=\"1\"/>\n",
                x, y, w, track_height, color
            ));
        }
    }

    // Legend
    let legend_y = height as i32 - 35;
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"label\" font-weight=\"bold\">Risk Level:</text>\n",
        margin, legend_y
    ));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"15\" height=\"15\" fill=\"#388e3c\" rx=\"2\"/>\n\
         <text x=\"{}\" y=\"{}\" class=\"label\">Low</text>\n",
        margin + 80, legend_y - 12, margin + 100, legend_y
    ));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"15\" height=\"15\" fill=\"#1976d2\" rx=\"2\"/>\n\
         <text x=\"{}\" y=\"{}\" class=\"label\">Medium</text>\n",
        margin + 150, legend_y - 12, margin + 170, legend_y
    ));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"15\" height=\"15\" fill=\"#f57c00\" rx=\"2\"/>\n\
         <text x=\"{}\" y=\"{}\" class=\"label\">High</text>\n",
        margin + 240, legend_y - 12, margin + 260, legend_y
    ));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"15\" height=\"15\" fill=\"#d32f2f\" rx=\"2\"/>\n\
         <text x=\"{}\" y=\"{}\" class=\"label\">Critical</text>\n",
        margin + 310, legend_y - 12, margin + 330, legend_y
    ));

    svg.push_str("</svg>");
    svg
}

/// Generate BED files from report data
fn generate_bed_files_from_report(
    report: &Report,
    snippy_gaps: &HashMap<String, Vec<Region>>,
    cfsan_gaps: &HashMap<String, Vec<Region>>,
) -> HashMap<String, String> {
    let mut beds = HashMap::new();

    // Use short chromosome name (without description) for BED files
    let ref_name_full = report.reference_summary.reference_name.clone();
    let ref_name = ref_name_full.split_whitespace().next().unwrap_or(&ref_name_full);

    // === ALL UNIQUE GAPS (aggregated from all samples) ===
    let mut all_gaps: Vec<(usize, usize, String)> = Vec::new();
    for sample_impact in &report.impact_analysis.samples {
        for region in &sample_impact.unique_gap_coords {
            all_gaps.push((region.start, region.end, sample_impact.sample.clone()));
        }
    }
    all_gaps.sort_by_key(|(start, _, _)| *start);

    let total_gap_bp: usize = all_gaps.iter().map(|(s, e, _)| e - s).sum();
    let mut all_gaps_bed = format!(
        "track name=\"ALL_UNIQUE_GAPS\" description=\"All unique gaps ({} regions, {} bp total)\" color=255,0,0\n",
        all_gaps.len(), total_gap_bp
    );
    for (start, end, sample) in &all_gaps {
        all_gaps_bed.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            ref_name, start, end, sample
        ));
    }
    beds.insert("ALL_UNIQUE_GAPS.bed".to_string(), all_gaps_bed);

    // === MM2 GAPS PER SAMPLE ===
    for ref_gap in &report.reference_gaps {
        let sample_name = &ref_gap.sample;
        let gaps = &ref_gap.reference_uncovered;
        if !gaps.is_empty() {
            let total_bp: usize = gaps.iter().map(|g| g.end.saturating_sub(g.start)).sum();
            let mut bed = format!(
                "track name=\"{}_mm2_gaps\" description=\"{} - mm2 gaps ({} regions, {} bp)\" color=233,69,96\n",
                sample_name, sample_name, gaps.len(), total_bp
            );
            for gap in gaps {
                bed.push_str(&format!(
                    "{}\t{}\t{}\tmm2_gap\n",
                    ref_name, gap.start, gap.end
                ));
            }
            beds.insert(format!("{}_mm2_gaps.bed", sample_name), bed);
        }
    }

    // === SNIPPY GAPS PER SAMPLE ===
    for (sample_name, gaps) in snippy_gaps {
        if !gaps.is_empty() {
            let total_bp: usize = gaps.iter().map(|g| g.end.saturating_sub(g.start)).sum();
            let mut bed = format!(
                "track name=\"{}_snippy_gaps\" description=\"{} - Snippy gaps ({} regions, {} bp)\" color=59,130,246\n",
                sample_name, sample_name, gaps.len(), total_bp
            );
            for gap in gaps {
                bed.push_str(&format!(
                    "{}\t{}\t{}\tsnippy_gap\n",
                    ref_name, gap.start, gap.end
                ));
            }
            beds.insert(format!("{}_snippy_gaps.bed", sample_name), bed);
        }
    }

    // === CFSAN GAPS PER SAMPLE ===
    for (sample_name, gaps) in cfsan_gaps {
        if !gaps.is_empty() {
            let total_bp: usize = gaps.iter().map(|g| g.end.saturating_sub(g.start)).sum();
            let mut bed = format!(
                "track name=\"{}_cfsan_gaps\" description=\"{} - CFSAN gaps ({} regions, {} bp)\" color=16,185,129\n",
                sample_name, sample_name, gaps.len(), total_bp
            );
            for gap in gaps {
                bed.push_str(&format!(
                    "{}\t{}\t{}\tcfsan_gap\n",
                    ref_name, gap.start, gap.end
                ));
            }
            beds.insert(format!("{}_cfsan_gaps.bed", sample_name), bed);
        }
    }

    // === SNP PIPELINE GAPS FROM REPORT (when loading existing report) ===
    // If snippy_gaps was empty but we have snp_comparison data, use snp_only_regions
    if snippy_gaps.is_empty() {
        if let Some(ref snp_cmp) = report.snp_comparison {
            for sample_cmp in &snp_cmp.samples {
                if !sample_cmp.snp_only_regions.is_empty() {
                    let total_bp: usize = sample_cmp.snp_only_regions.iter()
                        .map(|g| g.end.saturating_sub(g.start)).sum();
                    let mut bed = format!(
                        "track name=\"{}_snp_only\" description=\"{} - {} SNP-only gaps ({} regions, {} bp)\" color=234,88,12\n",
                        sample_cmp.sample, sample_cmp.sample, snp_cmp.pipeline_name,
                        sample_cmp.snp_only_regions.len(), total_bp
                    );
                    for gap in &sample_cmp.snp_only_regions {
                        bed.push_str(&format!(
                            "{}\t{}\t{}\tsnp_pipeline_gap\n",
                            ref_name, gap.start, gap.end
                        ));
                    }
                    beds.insert(format!("{}_snp_only.bed", sample_cmp.sample), bed);
                }
            }
        }
    }

    // === CFSAN GAPS FROM REPORT (when loading existing report) ===
    // If cfsan_gaps was empty but we have cfsan_comparison data
    if cfsan_gaps.is_empty() {
        if let Some(ref cfsan_cmp) = report.cfsan_comparison {
            for sample_cmp in &cfsan_cmp.sample_comparisons {
                // CFSAN comparison doesn't have gap regions, but we can indicate SNPs in gaps
                // Log info about gaps in CFSAN data
                if sample_cmp.snps_in_gaps > 0 {
                    log::debug!("{}: {} SNPs are in coreguard gaps",
                              sample_cmp.sample, sample_cmp.snps_in_gaps);
                }
            }
        }
    }

    // === VCF COMPARISON POSITIONS ===
    if let Some(ref vcf) = report.vcf_comparison {
        // Concordant SNPs (sample from sample_comparisons)
        let mut concordant_positions: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut snippy_only_positions: Vec<(usize, String, String)> = Vec::new(); // pos, sample, info
        let mut cfsan_only_positions: Vec<(usize, String, String)> = Vec::new();

        // Get positions from top_missed lists
        for pos in &vcf.top_missed_by_b {
            snippy_only_positions.push((pos.pos, pos.sample.clone(),
                format!("{}>{}",
                    pos.pipeline_a.as_ref().map(|p| p.ref_allele.as_str()).unwrap_or("?"),
                    pos.pipeline_a.as_ref().map(|p| p.alt_allele.as_str()).unwrap_or("?")
                )
            ));
        }

        for pos in &vcf.top_missed_by_a {
            cfsan_only_positions.push((pos.pos, pos.sample.clone(),
                format!("{}>{}",
                    pos.pipeline_b.as_ref().map(|p| p.ref_allele.as_str()).unwrap_or("?"),
                    pos.pipeline_b.as_ref().map(|p| p.alt_allele.as_str()).unwrap_or("?")
                )
            ));
        }

        // Snippy-only SNPs BED
        if !snippy_only_positions.is_empty() {
            snippy_only_positions.sort_by_key(|(pos, _, _)| *pos);
            let mut bed = format!(
                "track name=\"SNP_Snippy_only\" description=\"Snippy-only SNPs ({} positions)\" color=59,130,246\n",
                snippy_only_positions.len()
            );
            for (pos, sample, info) in &snippy_only_positions {
                bed.push_str(&format!(
                    "{}\t{}\t{}\t{}:{}\n",
                    ref_name, pos.saturating_sub(1), pos, sample, info
                ));
            }
            beds.insert("SNP_Snippy_only.bed".to_string(), bed);
        }

        // CFSAN-only SNPs BED
        if !cfsan_only_positions.is_empty() {
            cfsan_only_positions.sort_by_key(|(pos, _, _)| *pos);
            let mut bed = format!(
                "track name=\"SNP_CFSAN_only\" description=\"CFSAN-only SNPs ({} positions)\" color=16,185,129\n",
                cfsan_only_positions.len()
            );
            for (pos, sample, info) in &cfsan_only_positions {
                bed.push_str(&format!(
                    "{}\t{}\t{}\t{}:{}\n",
                    ref_name, pos.saturating_sub(1), pos, sample, info
                ));
            }
            beds.insert("SNP_CFSAN_only.bed".to_string(), bed);
        }

        // Summary info
        log::info!("Generated BED files: {} Snippy-only SNPs, {} CFSAN-only SNPs",
            snippy_only_positions.len(), cfsan_only_positions.len());
    }

    beds
}

/// Generate the main HTML page with IGV.js embedded
fn generate_igv_html(
    report_json: &str,
    reference_name: &str,
    reference_length: usize,
    bed_files: &HashMap<String, String>,
) -> String {
    // Generate track configurations for each BED file
    let mut tracks_js = String::from("[");
    let mut first = true;

    for (name, content) in bed_files {
        if !first {
            tracks_js.push_str(",");
        }
        first = false;

        // Determine color based on name and extract risk from track description
        let (color, display_name): (&str, String) = if name.contains("CORE") {
            ("rgb(255, 215, 0)", "CORE GENOME".to_string())  // Gold/yellow
        } else if content.contains("CRITICAL") {
            ("rgb(255, 0, 0)", format!("{} [CRITICAL]", name.replace("_gaps.bed", "")))
        } else if content.contains("HIGH") {
            ("rgb(255, 165, 0)", format!("{} [HIGH]", name.replace("_gaps.bed", "")))
        } else if content.contains("MEDIUM") {
            ("rgb(255, 255, 0)", format!("{} [MEDIUM]", name.replace("_gaps.bed", "")))
        } else {
            ("rgb(100, 200, 100)", name.replace("_gaps.bed", ""))
        };

        tracks_js.push_str(&format!(
            r#"{{
                name: "{}",
                url: "/bed/{}",
                format: "bed",
                color: "{}",
                displayMode: "EXPANDED",
                height: 40,
                visibilityWindow: -1,
                showLabels: true,
                labelField: "name"
            }}"#,
            display_name, name, color
        ));
    }
    tracks_js.push_str("]");

    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>coreguard - IGV Genome Browser</title>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
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
        }}
        .header {{
            background: linear-gradient(135deg, #2c3e50, #3498db);
            color: white;
            padding: 20px;
            text-align: center;
        }}
        .header h1 {{ margin: 0; font-size: 1.8rem; }}
        .header .subtitle {{ opacity: 0.9; margin-top: 5px; }}
        .tabs {{
            display: flex;
            background: #2c3e50;
            padding: 0 20px;
        }}
        .tab {{
            padding: 12px 24px;
            color: rgba(255,255,255,0.7);
            cursor: pointer;
            border-bottom: 3px solid transparent;
            transition: all 0.2s;
        }}
        .tab:hover {{ color: white; }}
        .tab.active {{
            color: white;
            border-bottom-color: #3498db;
            background: rgba(255,255,255,0.1);
        }}
        .content {{ padding: 20px; }}
        .panel {{ display: none; }}
        .panel.active {{ display: block; }}
        #igv-container {{
            border: 1px solid var(--border);
            border-radius: 8px;
            overflow: hidden;
            margin-bottom: 20px;
        }}
        .cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        .card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .card-title {{ font-size: 0.85rem; color: #6c757d; text-transform: uppercase; }}
        .card-value {{ font-size: 2rem; font-weight: bold; margin-top: 5px; }}
        .card-value.critical {{ color: var(--critical); }}
        .card-value.low {{ color: var(--low); }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: var(--card-bg);
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
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
        .badge-high {{ background: #fff3cd; color: #856404; }}
        .badge-critical {{ background: #f8d7da; color: #721c24; }}
        .locus-input {{
            display: flex;
            gap: 10px;
            margin-bottom: 20px;
            align-items: center;
        }}
        .locus-input input {{
            padding: 8px 12px;
            border: 1px solid var(--border);
            border-radius: 4px;
            font-size: 14px;
            width: 300px;
        }}
        .locus-input button {{
            padding: 8px 16px;
            background: #3498db;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        .locus-input button:hover {{ background: #2980b9; }}
        .quick-nav {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            margin-bottom: 20px;
        }}
        .quick-nav button {{
            padding: 6px 12px;
            background: #e9ecef;
            border: 1px solid var(--border);
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.85rem;
        }}
        .quick-nav button:hover {{ background: #dee2e6; }}
        footer {{
            text-align: center;
            padding: 20px;
            color: #6c757d;
            border-top: 1px solid var(--border);
            margin-top: 40px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>coreguard</h1>
        <div class="subtitle">Pre-alignment QC for SNP Pipeline Analysis | Reference: <strong>{reference_name}</strong></div>
    </div>

    <div class="tabs">
        <div class="tab active" data-panel="home">Home</div>
        <div class="tab" data-panel="gaps-map">Gaps Map</div>
        <div class="tab" data-panel="pairwise">Pairwise Gaps</div>
        <div class="tab" data-panel="gene-zone">Gene Zone</div>
        <div class="tab" data-panel="empty-column">Coverage Thresholds</div>
        <div class="tab" data-panel="snp-compare">SNP Pipeline</div>
        <div class="tab" data-panel="vcf-compare">VCF Compare</div>
        <div class="tab" data-panel="snippy-bugfix">Snippy Bugfix</div>
        <div class="tab" data-panel="snp-tracks">SNP Tracks</div>
        <div class="tab" data-panel="igv">IGV Browser</div>
        <div class="tab" data-panel="summary">Summary</div>
        <div class="tab" data-panel="samples">Sample Details</div>
        <div class="tab" data-panel="methods">Methods</div>
        <div class="tab" data-panel="files">Files</div>
        <div class="tab" data-panel="help">Help</div>
    </div>

    <div class="content">
        <!-- HOME TAB -->
        <div class="panel active" id="home">
            <!-- Workflow Diagram Section -->
            <div style="background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px; padding: 25px; margin-bottom: 25px; font-family: 'Times New Roman', serif;">
                <h3 style="margin-bottom: 20px; text-align: center; font-weight: normal; font-size: 1.3rem; border-bottom: 1px solid #ccc; padding-bottom: 10px;">
                    <span style="font-variant: small-caps;">Motivation</span>
                </h3>
                <div style="display: grid; grid-template-columns: 1fr 60px 1fr; gap: 15px; align-items: start;">
                    <!-- Traditional Approach -->
                    <div style="border: 1px solid #999; padding: 15px; background: white;">
                        <h4 style="margin-bottom: 12px; text-align: center; font-weight: normal; font-size: 1rem; border-bottom: 1px solid #ddd; padding-bottom: 8px;">
                            (A) Traditional Workflow
                        </h4>
                        <div style="font-size: 0.9rem; line-height: 1.8;">
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                <i>n</i> samples &rarr; SNP pipeline
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                Suboptimal results detected
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                Remove suspect sample, re-run
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                Iterate until acceptable
                            </div>
                        </div>
                        <p style="margin-top: 12px; font-size: 0.85rem; text-align: center; color: #666; font-style: italic;">
                            Complexity: O(k) pipeline executions
                        </p>
                    </div>

                    <!-- Arrow -->
                    <div style="display: flex; align-items: center; justify-content: center; height: 100%; font-size: 1.1rem; color: #666;">
                        vs.
                    </div>

                    <!-- coreguard Approach -->
                    <div style="border: 1px solid #999; padding: 15px; background: white;">
                        <h4 style="margin-bottom: 12px; text-align: center; font-weight: normal; font-size: 1rem; border-bottom: 1px solid #ddd; padding-bottom: 8px;">
                            (B) Proposed Workflow
                        </h4>
                        <div style="font-size: 0.9rem; line-height: 1.8;">
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                <i>n</i> samples &rarr; <b>coreguard</b>
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                Identify problematic samples
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                Exclude flagged samples
                            </div>
                            <div style="text-align: center; color: #666;">&darr;</div>
                            <div style="border: 1px solid #ccc; padding: 8px; text-align: center; margin: 8px 0; background: #fafafa;">
                                <i>n'</i> samples &rarr; SNP pipeline
                            </div>
                        </div>
                        <p style="margin-top: 12px; font-size: 0.85rem; text-align: center; color: #666; font-style: italic;">
                            Complexity: O(1) pipeline execution
                        </p>
                    </div>
                </div>
                <p style="text-align: center; margin-top: 20px; font-size: 0.9rem; color: #444; max-width: 700px; margin-left: auto; margin-right: auto; line-height: 1.6;">
                    <b>Rationale:</b> A single divergent sample can reduce the core genome by 5&ndash;10%,
                    introducing alignment gaps that propagate to all samples and compromise SNP resolution.
                    Pre-analysis quality control enables identification of such samples <i>a priori</i>.
                </p>
            </div>

            <!-- Quick Summary Cards -->
            <div class="cards" style="margin-bottom: 25px;">
                <div class="card">
                    <div class="card-title">Total Samples</div>
                    <div class="card-value" id="home-total-samples">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Reference</div>
                    <div class="card-value" style="font-size: 1rem;" id="home-reference">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Estimated Core</div>
                    <div class="card-value" id="home-core">-</div>
                </div>
                <div class="card">
                    <div class="card-title">Flagged Samples</div>
                    <div class="card-value" id="home-flagged">-</div>
                </div>
            </div>

            <!-- Navigation hints -->
            <div style="background: #e3f2fd; border-radius: 8px; padding: 20px;">
                <h4 style="margin-bottom: 15px; color: #1565c0;">Explore the Analysis</h4>
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px;">
                    <div style="background: white; padding: 15px; border-radius: 4px; cursor: pointer;" onclick="document.querySelector('[data-panel=gaps-map]').click()">
                        <strong>Gaps Map</strong>
                        <p style="font-size: 0.85rem; color: #666; margin-top: 5px;">Visual distribution of unique gaps across the genome</p>
                    </div>
                    <div style="background: white; padding: 15px; border-radius: 4px; cursor: pointer;" onclick="document.querySelector('[data-panel=pairwise]').click()">
                        <strong>Pairwise Gaps</strong>
                        <p style="font-size: 0.85rem; color: #666; margin-top: 5px;">Gap overlap analysis between sample pairs</p>
                    </div>
                    <div style="background: white; padding: 15px; border-radius: 4px; cursor: pointer;" onclick="document.querySelector('[data-panel=summary]').click()">
                        <strong>Summary</strong>
                        <p style="font-size: 0.85rem; color: #666; margin-top: 5px;">Interactive sample selector with core simulation</p>
                    </div>
                    <div style="background: white; padding: 15px; border-radius: 4px; cursor: pointer;" onclick="document.querySelector('[data-panel=methods]').click()">
                        <strong>Methods</strong>
                        <p style="font-size: 0.85rem; color: #666; margin-top: 5px;">Detailed explanation of all calculations</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- GAPS MAP TAB -->
        <div class="panel" id="gaps-map">
            <h3 style="margin-bottom: 15px;">Unique Gaps Distribution Map</h3>
            <p style="color: #6c757d; margin-bottom: 15px;">
                Visualization of unique gaps across the genome. Samples sorted by impact (best to worst).
            </p>
            <div style="background: #e3f2fd; padding: 12px 15px; border-radius: 6px; margin-bottom: 15px; font-size: 0.9rem;">
                <strong>ℹ️ Note:</strong> Gaps shown here derive directly from read alignment to the reference.
                With <code>min-depth=1</code>, a gap indicates a position where <strong>no reads map</strong>.
                This reflects actual genomic divergence between samples and reference, not coverage issues.
            </div>
            <div style="background: white; border-radius: 8px; padding: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); overflow-x: auto;">
                <img src="/gaps.svg" alt="Gaps Distribution" style="width: 100%; min-width: 800px;"/>
            </div>

            <div style="margin-top: 20px; background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                <h4 style="margin-bottom: 15px;">Risk Level Legend</h4>
                <p style="color: #6c757d; margin-bottom: 10px;">Risk level is calculated based on <strong>Core Reduction</strong> = unique_gap_bases / reference_length</p>
                <table style="width: 100%; border-collapse: collapse;">
                    <tr style="background: #f8f9fa;">
                        <th style="padding: 10px; text-align: left; border-bottom: 2px solid #dee2e6;">Risk Level</th>
                        <th style="padding: 10px; text-align: left; border-bottom: 2px solid #dee2e6;">Core Reduction</th>
                        <th style="padding: 10px; text-align: left; border-bottom: 2px solid #dee2e6;">Meaning</th>
                    </tr>
                    <tr>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;"><span style="display: inline-block; width: 20px; height: 20px; background: #388e3c; border-radius: 4px; vertical-align: middle; margin-right: 8px;"></span><strong>Low</strong></td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">&lt; 0.5%</td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">Minimal impact on core genome</td>
                    </tr>
                    <tr>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;"><span style="display: inline-block; width: 20px; height: 20px; background: #1976d2; border-radius: 4px; vertical-align: middle; margin-right: 8px;"></span><strong>Medium</strong></td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">0.5% - 2%</td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">Moderate impact, review recommended</td>
                    </tr>
                    <tr>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;"><span style="display: inline-block; width: 20px; height: 20px; background: #f57c00; border-radius: 4px; vertical-align: middle; margin-right: 8px;"></span><strong>High</strong></td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">2% - 5%</td>
                        <td style="padding: 10px; border-bottom: 1px solid #dee2e6;">Significant core reduction, consider excluding</td>
                    </tr>
                    <tr>
                        <td style="padding: 10px;"><span style="display: inline-block; width: 20px; height: 20px; background: #d32f2f; border-radius: 4px; vertical-align: middle; margin-right: 8px;"></span><strong>Critical</strong></td>
                        <td style="padding: 10px;">&gt; 5%</td>
                        <td style="padding: 10px;">Severe impact, strongly recommend excluding</td>
                    </tr>
                </table>
            </div>
        </div>

        <div class="panel" id="igv">
            <h3 style="margin-bottom: 15px;">Interactive Genome Browser (IGV.js)</h3>
            <div class="locus-input">
                <input type="text" id="locus-input" placeholder="Enter locus (e.g., {ref_name}:1-100000)">
                <button onclick="goToLocus()">Go</button>
                <span style="color: #6c757d; margin-left: 10px;">or use quick navigation:</span>
            </div>
            <div class="quick-nav" id="quick-nav"></div>
            <div id="igv-container"></div>
        </div>

        <div class="panel" id="summary">
            <div class="cards" id="summary-cards"></div>
            <h3 style="margin: 20px 0 10px;">Sample Selection - Core Impact Simulator</h3>
            <p style="color: #6c757d; margin-bottom: 15px;">Select samples to include and see how core genome changes. Sorted by impact (best to worst).</p>
            <div id="sample-selector" style="background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);"></div>
            <div id="core-result" style="margin-top: 20px; padding: 20px; background: #e8f5e9; border-radius: 8px; font-size: 1.2rem;"></div>
        </div>

        <div class="panel" id="samples">
            <p style="color: #6c757d; margin-bottom: 15px;">Sorted by impact (best to worst)</p>
            <table id="samples-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Covered (bp)</th>
                        <th>Coverage %</th>
                        <th>Risk Level</th>
                        <th>Unique Gap Bases</th>
                        <th>Core Reduction</th>
                    </tr>
                </thead>
                <tbody></tbody>
            </table>
        </div>

        <div class="panel" id="pairwise">
            <h3 style="margin-bottom: 15px;">Pairwise Gap Quality Analysis</h3>
            <p style="color: #6c757d; margin-bottom: 15px;">
                For each sample pair, shows gap <strong>union</strong> (A ∪ B = gaps in either sample) and <strong>intersection</strong> (A ∩ B = gaps in both).
                <br><strong>Quality Score</strong> = intersection / union. High score (near 1.0) = similar gap patterns, low score (near 0.0) = different gaps.
            </p>
            <div style="margin-bottom: 20px; display: flex; flex-wrap: wrap; gap: 20px; align-items: center;">
                <div>
                    <label for="pairwise-filter" style="margin-right: 10px;">Core reduction:</label>
                    <select id="pairwise-filter" onchange="filterPairwise()" style="padding: 8px; border-radius: 4px; border: 1px solid #dee2e6;">
                        <option value="all">All pairs</option>
                        <option value="with-gaps">&gt;0% (any reduction)</option>
                        <option value="high-gaps">&gt;1%</option>
                        <option value="critical-gaps">&gt;5%</option>
                        <option value="very-high-gaps">&gt;10%</option>
                        <option value="extreme-gaps">&gt;15%</option>
                    </select>
                </div>
                <div>
                    <label for="quality-filter" style="margin-right: 10px;">Quality score:</label>
                    <select id="quality-filter" onchange="filterPairwise()" style="padding: 8px; border-radius: 4px; border: 1px solid #dee2e6;">
                        <option value="all">All</option>
                        <option value="high-quality">&gt;80% (similar gaps)</option>
                        <option value="medium-quality">50-80%</option>
                        <option value="low-quality">&lt;50% (different gaps)</option>
                        <option value="very-low-quality">&lt;20% (very different)</option>
                    </select>
                </div>
                <span id="pair-count" style="color: #6c757d;"></span>
            </div>
            <div id="pairwise-container" style="display: grid; gap: 20px;"></div>
        </div>

        <!-- GENE ZONE TAB -->
        <div class="panel" id="gene-zone">
            <h3 style="margin-bottom: 15px;">Gene Zone Analysis</h3>
            <div id="gene-zone-content">
                <p style="color: #666;">Loading gene zone analysis...</p>
            </div>
        </div>

        <!-- EMPTY COLUMN TAB -->
        <div class="panel" id="empty-column">
            <h3 style="margin-bottom: 15px;">Coverage Thresholds Analysis</h3>
            <p style="color: #666; margin-bottom: 20px;">
                How does the core genome size change at different sample inclusion thresholds?
                A position is included in the "core" if at least X% of samples cover it.
            </p>
            <div id="empty-column-content">
                <p style="color: #666;">Loading threshold analysis...</p>
            </div>
        </div>

        <!-- SNP PIPELINE COMPARISON TAB -->
        <div class="panel" id="snp-compare">
            <h3 style="margin-bottom: 15px;">SNP Pipeline Comparison</h3>
            <p style="color: #666; margin-bottom: 20px;">
                Compares coreguard gap predictions with positions actually filtered by SNP analysis pipelines (Snippy, CFSAN).
            </p>
            <div id="snp-compare-content">
                <p style="color: #666;">Loading SNP pipeline comparison...</p>
            </div>
        </div>

        <!-- VCF COMPARISON TAB -->
        <div class="panel" id="vcf-compare">
            <h3 style="margin-bottom: 15px;">VCF Pipeline Comparison</h3>
            <p style="color: #666; margin-bottom: 20px;">
                Directly compares variant calls between Snippy and CFSAN to identify discordant positions and their likely causes.
            </p>
            <div id="vcf-compare-content">
                <p style="color: #666;">Loading VCF comparison...</p>
            </div>
        </div>

        <!-- SNIPPY BUGFIX TAB -->
        <div class="panel" id="snippy-bugfix">
            <h3 style="margin-bottom: 15px;">🔬 Snippy Artifact Correction</h3>
            <p style="color: #666; margin-bottom: 20px;">
                Validates Snippy-reported polymorphisms against BAM pileup data to detect artifacts caused by complex variant decomposition.
            </p>
            <div id="snippy-bugfix-content">
                <p style="color: #666;">Loading BAM validation results...</p>
            </div>
        </div>

        <!-- SNP TRACKS TAB -->
        <div class="panel" id="snp-tracks">
            <h3 style="margin-bottom: 15px;">📊 SNP Pipeline Comparison Tracks</h3>
            <p style="color: #666; margin-bottom: 20px;">
                Visual comparison of SNP calls across samples. Each row shows a sample, each column a genomic position.
                Colors indicate which pipeline retained or removed the SNP.
            </p>

            <!-- Controls -->
            <div style="background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;">
                <div style="display: flex; gap: 15px; flex-wrap: wrap; align-items: center;">
                    <div>
                        <label style="font-size: 0.85rem; color: #666;">Region:</label>
                        <input type="text" id="snp-region-input" placeholder="e.g., 100000-200000"
                               style="padding: 8px; border: 1px solid #ddd; border-radius: 4px; width: 150px;">
                    </div>
                    <div>
                        <label style="font-size: 0.85rem; color: #666;">Window size:</label>
                        <select id="snp-window-size" style="padding: 8px; border: 1px solid #ddd; border-radius: 4px;">
                            <option value="1000">1 kb</option>
                            <option value="5000">5 kb</option>
                            <option value="10000" selected>10 kb</option>
                            <option value="50000">50 kb</option>
                            <option value="100000">100 kb</option>
                        </select>
                    </div>
                    <button id="snp-tracks-update" style="padding: 8px 20px; background: #2196f3; color: white; border: none; border-radius: 4px; cursor: pointer;">
                        Update View
                    </button>
                    <button id="snp-tracks-prev" style="padding: 8px 15px; background: #607d8b; color: white; border: none; border-radius: 4px; cursor: pointer;">
                        ← Prev
                    </button>
                    <button id="snp-tracks-next" style="padding: 8px 15px; background: #607d8b; color: white; border: none; border-radius: 4px; cursor: pointer;">
                        Next →
                    </button>
                </div>
            </div>

            <!-- Legend -->
            <div style="display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 20px; font-size: 0.85rem;">
                <div style="display: flex; align-items: center; gap: 5px;">
                    <span style="display: inline-block; width: 16px; height: 16px; background: #4caf50; border-radius: 2px;"></span>
                    <span>Concordant (both pipelines)</span>
                </div>
                <div style="display: flex; align-items: center; gap: 5px;">
                    <span style="display: inline-block; width: 16px; height: 16px; background: #2196f3; border-radius: 2px;"></span>
                    <span>CFSAN-only</span>
                </div>
                <div style="display: flex; align-items: center; gap: 5px;">
                    <span style="display: inline-block; width: 16px; height: 16px; background: #ff9800; border-radius: 2px;"></span>
                    <span>Snippy-only (real)</span>
                </div>
                <div style="display: flex; align-items: center; gap: 5px;">
                    <span style="display: inline-block; width: 16px; height: 16px; background: #f44336; border-radius: 2px;"></span>
                    <span>Snippy artifact</span>
                </div>
                <div style="display: flex; align-items: center; gap: 5px;">
                    <span style="display: inline-block; width: 16px; height: 16px; background: #e0e0e0; border-radius: 2px;"></span>
                    <span>No SNP</span>
                </div>
            </div>

            <!-- Track visualization container -->
            <div id="snp-tracks-container" style="background: white; border: 1px solid #ddd; border-radius: 8px; padding: 20px; overflow-x: auto;">
                <p style="color: #666; text-align: center;">Loading SNP track data...</p>
            </div>

            <!-- Position info -->
            <div id="snp-tracks-info" style="margin-top: 15px; padding: 15px; background: #f5f5f5; border-radius: 8px; font-size: 0.9rem;">
                <span id="snp-tracks-region-info">Select a region to view SNP tracks</span>
            </div>
        </div>

        <!-- METHODS TAB -->
        <div class="panel" id="methods">
            <div style="max-width: 1000px;">
                <h2 style="color: #2c3e50; margin-bottom: 25px; border-bottom: 2px solid #3498db; padding-bottom: 10px;">Methods</h2>

                <!-- Step 1: Reference Alignment -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 1: Reference Alignment</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        Each sample assembly is aligned to the reference genome using <strong>minimap2</strong> with the <code>asm5</code> preset
                        (optimized for &gt;95% sequence identity assemblies).
                    </p>
                    <pre style="background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 4px; overflow-x: auto;">
minimap2 -x asm5 -c reference.fasta sample.fasta &gt; alignment.paf</pre>
                    <p style="line-height: 1.8; margin-top: 15px;">
                        From the PAF output, we extract for each sample:
                    </p>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li><strong>Covered regions</strong>: Reference positions covered by aligned contigs</li>
                        <li><strong>Gap regions</strong>: Reference positions NOT covered (missing in sample)</li>
                        <li><strong>Identity</strong>: Sequence identity in aligned regions</li>
                    </ul>
                    <div style="background: #fff3cd; padding: 15px; border-radius: 4px; margin-top: 15px;">
                        <strong>Note:</strong> Insertions in reads (CIGAR <code>I</code>) do NOT affect reference coverage.
                        Deletions in reads (CIGAR <code>D</code>) mean the read is missing bases, but the reference is still &quot;covered&quot;
                        by that read at those positions. Only positions with NO aligned reads are counted as gaps.
                    </div>
                </div>

                <!-- Step 2: Gap Detection -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 2: Gap Detection</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        For each sample, we compute a <strong>coverage bitmap</strong>: a boolean array of size <code>reference_length</code>
                        where each position is <code>true</code> (gap) or <code>false</code> (covered).
                    </p>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
Reference:      |========================| (3 Mb)
                0        1M       2M     3M

Sample A:       |========|    |==========|
                covered   GAP   covered

Coverage bitmap:
  bitmap[0..1M]   = false  (covered)
  bitmap[1M..1.5M] = true   (GAP)
  bitmap[1.5M..3M] = false  (covered)</pre>
                    <p style="line-height: 1.8; margin-top: 15px;">
                        <strong>Unique gaps</strong> for sample A = positions that are gaps in A but covered in ALL other samples.
                    </p>
                    <div style="background: #e3f2fd; padding: 15px; border-radius: 4px; margin-top: 15px;">
                        <strong>ℹ️ Gap detection from alignment:</strong> With <code>min-depth=1</code>, a gap is defined as a position where
                        <strong>no read maps</strong>. This reflects actual genomic divergence between samples and reference (e.g., large deletions,
                        mobile genetic elements), not sequencing depth or coverage issues. Gaps derive directly from minimap2 alignment without additional filtering.
                    </div>
                </div>

                <!-- Step 3: Pairwise Analysis -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 3: Pairwise Gap Analysis</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        For each pair of samples (A, B), we perform set operations on their gap bitmaps:
                    </p>
                    <table style="width: 100%; border-collapse: collapse; margin: 15px 0;">
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Metric</th>
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Formula</th>
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Meaning</th>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>Gap Union</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6; font-family: monospace;">gaps_A &cup; gaps_B</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Positions that are gaps in A <em>or</em> B (or both)</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>Gap Intersection</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6; font-family: monospace;">gaps_A &cap; gaps_B</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Positions that are gaps in <em>both</em> A and B</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>Quality Score</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6; font-family: monospace;">intersection / union</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">1.0 = identical gaps, 0.0 = no overlap</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>Pairwise Core</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6; font-family: monospace;">ref_length - union</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Core genome size if using only A and B</td>
                        </tr>
                    </table>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
Example:
  Gap A:        [====]              (1M-1.5M, 500kb)
  Gap B:             [====]         (1.3M-1.8M, 500kb)

  Union:        [========]          (1M-1.8M, 800kb)
  Intersection:   [==]              (1.3M-1.5M, 200kb)

  Quality = 200kb / 800kb = 0.25 (25% overlap)
  Pairwise core = 3Mb - 800kb = 2.2Mb</pre>

                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Pairwise Filters</h4>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        The pairwise view provides two filters to identify problematic sample pairs:
                    </p>
                    <table style="width: 100%; border-collapse: collapse; margin: 15px 0;">
                        <tr style="background: #e3f2fd;">
                            <th colspan="3" style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Core Reduction Filter</th>
                        </tr>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Filter</th>
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Threshold</th>
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Shows pairs where</th>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&gt;0%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">(ref - core) / ref &gt; 0</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Any core reduction (most pairs)</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&gt;5%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">(ref - core) / ref &gt; 0.05</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Significant core loss</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&gt;10%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">(ref - core) / ref &gt; 0.10</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Major core loss - investigate pair</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&gt;15%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">(ref - core) / ref &gt; 0.15</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Critical - likely assembly issues</td>
                        </tr>
                    </table>
                    <table style="width: 100%; border-collapse: collapse; margin: 15px 0;">
                        <tr style="background: #fff3e0;">
                            <th colspan="3" style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Quality Score Filter</th>
                        </tr>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Filter</th>
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Range</th>
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Interpretation</th>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&gt;80%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">intersection/union &gt; 0.8</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Similar gaps - samples are compatible</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">50-80%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">0.5 &le; intersection/union &le; 0.8</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Moderate overlap - some shared gaps</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&lt;50%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">intersection/union &lt; 0.5</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Different gaps - one sample may be problematic</td>
                        </tr>
                        <tr>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">&lt;20%</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">intersection/union &lt; 0.2</td>
                            <td style="padding: 10px; border: 1px solid #dee2e6;">Very different - investigate immediately</td>
                        </tr>
                    </table>
                    <p style="line-height: 1.8; background: #fff8e1; padding: 10px; border-radius: 4px; border-left: 4px solid #ffc107;">
                        <strong>Tip:</strong> Low quality score + high core reduction = one sample has many unique gaps.
                        Check the "Unique to A/B" values to identify which sample is the outlier.
                    </p>
                </div>

                <!-- Step 4: Core Estimation -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 4: Core Genome Estimation</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        The <strong>estimated core genome</strong> is the set of reference positions covered by ALL samples:
                    </p>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
core_positions = reference_positions - (gaps_1 &cup; gaps_2 &cup; ... &cup; gaps_n)

estimated_core_size = |core_positions|

core_reduction = 1 - (estimated_core_size / reference_length)</pre>
                    <p style="line-height: 1.8; margin-top: 15px;">
                        <strong>Core Reduction per sample</strong>: How much core is lost specifically due to this sample's unique gaps.
                    </p>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
unique_gaps_i = gaps_i - (gaps_1 &cup; ... &cup; gaps_{{i-1}} &cup; gaps_{{i+1}} &cup; ... &cup; gaps_n)

core_reduction_i = |unique_gaps_i| / reference_length</pre>
                </div>

                <!-- Step 5: Risk Levels -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 5: Risk Level Classification</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        Each sample is assigned a risk level based on its core reduction:
                    </p>
                    <table style="width: 100%; border-collapse: collapse;">
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Risk Level</th>
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Core Reduction</th>
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Recommendation</th>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><span style="color: #388e3c; font-weight: bold;">Low</span></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">&lt; 0.5%</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Include in analysis</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><span style="color: #1976d2; font-weight: bold;">Medium</span></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">0.5% - 2%</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Review before including</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><span style="color: #f57c00; font-weight: bold;">High</span></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">2% - 5%</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Consider excluding</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><span style="color: #d32f2f; font-weight: bold;">Critical</span></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">&gt; 5%</td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Strongly recommend excluding</td>
                        </tr>
                    </table>
                </div>

                <!-- Step 6: Coverage Thresholds -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #e74c3c; margin-bottom: 15px;">Step 6: Coverage Threshold Analysis</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        Analyzes how the <strong>core genome size</strong> changes at different sample inclusion thresholds.
                        A position is considered "core" if at least X% of samples have coverage at that position.
                    </p>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
For each reference position P:
  coverage[P] = number of samples with read depth &ge; min_depth at P

At threshold T%:
  min_samples = ceil(T% &times; total_samples)
  core_positions = positions where coverage[P] &ge; min_samples
  core_size = |core_positions|</pre>
                    <p style="line-height: 1.8; margin-top: 15px;">
                        <strong>Example:</strong> With 4 samples and threshold 75%:
                    </p>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
min_samples = ceil(0.75 &times; 4) = 3

Position 1000: 4/4 samples cover &rarr; IN CORE
Position 2000: 3/4 samples cover &rarr; IN CORE
Position 3000: 2/4 samples cover &rarr; NOT in core
Position 4000: 1/4 samples cover &rarr; NOT in core</pre>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Coverage Histogram</h4>
                    <p style="line-height: 1.8;">
                        Shows the <strong>distribution of per-position coverage</strong>:
                    </p>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li>How many positions are covered by ALL samples (100% threshold)</li>
                        <li>How many positions are covered by most samples (80-99%)</li>
                        <li>How many positions are covered by few samples (&lt;50%)</li>
                    </ul>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Fragile Regions</h4>
                    <p style="line-height: 1.8;">
                        Regions covered by <strong>&lt;50% of samples</strong> are flagged as "fragile".
                        These regions:
                    </p>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li>May indicate assembly/sequencing problems in specific samples</li>
                        <li>May represent variable genomic regions (prophages, genomic islands)</li>
                        <li>Will be excluded from core at strict thresholds</li>
                    </ul>
                    <p style="line-height: 1.8; background: #e3f2fd; padding: 10px; border-radius: 4px; border-left: 4px solid #2196f3; margin-top: 15px;">
                        <strong>Use case:</strong> If your SNP pipeline requires 95% sample coverage per position,
                        check the "95%" row to see how much of the reference will be in your core.
                    </p>
                </div>

                <!-- Step 7: SNP Pipeline Comparison -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #9c27b0; margin-bottom: 15px;">Step 7: SNP Pipeline Comparison</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        <strong>Optional validation</strong>: Compare coreguard's gap predictions with actual filtering
                        performed by SNP pipelines like <strong>Snippy</strong> or <strong>CFSAN</strong>.
                    </p>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">What is compared</h4>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li><strong>Coreguard gaps</strong>: Positions where coreguard predicts coverage problems (from pairwise analysis)</li>
                        <li><strong>SNP pipeline filtered</strong>: Positions masked by the SNP pipeline (e.g., 'N' in Snippy's aligned.fa)</li>
                    </ul>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Comparison categories</h4>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
For each sample and position P:

<span style="color: #27ae60;">CONCORDANT</span>:     Both coreguard AND SNP pipeline filter P
<span style="color: #e67e22;">CG-ONLY</span>:        Coreguard flags P, but SNP pipeline KEEPS it
<span style="color: #3498db;">SNP-ONLY</span>:       SNP pipeline filters P, but coreguard KEEPS it

Concordance % = |Concordant| / |Concordant &cup; CG-only &cup; SNP-only| &times; 100</pre>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Interpretation</h4>
                    <table style="width: 100%; border-collapse: collapse; margin-top: 10px;">
                        <tr style="background: #f8f9fa;">
                            <th style="padding: 12px; border: 1px solid #dee2e6; text-align: left;">Category</th>
                            <th style="padding: 12px; border: 1px solid #dee2e6; text-align: left;">Meaning</th>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>High CG-only</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Coreguard is <em>more conservative</em> &mdash; flags regions that SNP pipelines still use</td>
                        </tr>
                        <tr style="background: #f8f9fa;">
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>High SNP-only</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">SNP pipeline applies <em>additional quality filters</em> (e.g., depth, mapping quality)</td>
                        </tr>
                        <tr>
                            <td style="padding: 12px; border: 1px solid #dee2e6;"><strong>High concordance</strong></td>
                            <td style="padding: 12px; border: 1px solid #dee2e6;">Good agreement between coreguard predictions and SNP pipeline filtering</td>
                        </tr>
                    </table>
                    <p style="line-height: 1.8; background: #fff3e0; padding: 10px; border-radius: 4px; border-left: 4px solid #ff9800; margin-top: 15px;">
                        <strong>Note:</strong> Coreguard uses <em>coverage-based</em> gap detection, while SNP pipelines
                        apply <em>quality-based</em> masking (depth, mapping quality, base quality). Some discordance is expected.
                    </p>
                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Input Files Used</h4>
                    <div style="background: #f8f9fa; padding: 15px; border-radius: 4px; margin-bottom: 15px;">
                        <p style="margin-bottom: 10px;"><strong>Snippy-core</strong> (via <code>--snippycore-dir</code>):</p>
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace; width: 40%;">snippycore.tab</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Per-position SNP calls for each sample (CHR, POS, REF, sample alleles)</td>
                            </tr>
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">snippycore.txt</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Core genome statistics and summary</td>
                            </tr>
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">hamming_distances*.tsv</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Distance matrix between samples (if available)</td>
                            </tr>
                        </table>

                        <p style="margin-top: 15px; margin-bottom: 10px;"><strong>CFSAN SNP Pipeline</strong> (via <code>--cfsan-dir</code>):</p>
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace; width: 40%;">snplist.txt</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">List of SNP positions called by VarScan2</td>
                            </tr>
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">metrics.tsv</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Per-sample QC metrics (depth, coverage, etc.)</td>
                            </tr>
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">snp_distance_matrix.tsv</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Pairwise SNP distance matrix</td>
                            </tr>
                        </table>
                    </div>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
# Example directory structure
snippy_core_results/
├── snippycore.tab         # Core SNPs
├── snippycore.txt         # Summary
└── hamming_distances.tsv  # Distances

cfsan_results/
├── snplist.txt            # SNP positions
├── metrics.tsv            # QC metrics
└── snp_distance_matrix.tsv</pre>
                </div>

                <!-- Step 8: VCF Compare -->
                <div style="background: white; border-radius: 8px; padding: 25px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #9c27b0; margin-bottom: 15px;">Step 8: VCF Comparison (Direct Variant Comparison)</h3>
                    <p style="line-height: 1.8; margin-bottom: 15px;">
                        <strong>Direct VCF comparison</strong> between pipelines, comparing the actual variant calls
                        without intermediate core genome filtering. This shows exactly which SNPs each pipeline called.
                    </p>
                    <div style="background: #fff3e0; padding: 15px; border-radius: 4px; margin-bottom: 15px; border-left: 4px solid #ff9800;">
                        <strong>Key difference from SNP Pipeline tab:</strong>
                        <ul style="margin-top: 10px; margin-bottom: 0; margin-left: 20px;">
                            <li><strong>SNP Pipeline tab</strong>: Compares coreguard's gap predictions with pipeline's filtered positions</li>
                            <li><strong>VCF Compare tab</strong>: Directly compares variant calls between two pipelines (e.g., Snippy vs CFSAN)</li>
                        </ul>
                    </div>

                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">What is compared</h4>
                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
For each position P in any sample:

<span style="color: #27ae60;">CONCORDANT</span>:       Both pipelines call the SAME variant (ref→alt match)
<span style="color: #e67e22;">ONLY PIPELINE A</span>: Only Snippy called this variant
<span style="color: #ad1457;">ONLY PIPELINE B</span>: Only CFSAN called this variant
<span style="color: #7b1fa2;">DISCORDANT</span>:      Both call a variant but with DIFFERENT alleles</pre>

                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Input Files Used</h4>
                    <div style="background: #f8f9fa; padding: 15px; border-radius: 4px; margin-bottom: 15px;">
                        <p style="margin-bottom: 10px;"><strong>Snippy VCFs</strong> (via <code>--vcf-snippy</code>):</p>
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace; width: 50%;">*_snippy/*.vcf</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">VCF files from Snippy (FreeBayes variant calls)</td>
                            </tr>
                        </table>

                        <p style="margin-top: 15px; margin-bottom: 10px;"><strong>CFSAN VCFs</strong> (via <code>--vcf-cfsan</code>):</p>
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                            <tr>
                                <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace; width: 50%;">samples/*/var.flt.vcf</td>
                                <td style="padding: 8px; border: 1px solid #dee2e6;">Filtered VCF files from VarScan2</td>
                            </tr>
                        </table>
                    </div>

                    <pre style="background: #f8f9fa; padding: 15px; border-radius: 4px; font-family: monospace; line-height: 1.6;">
# Example directory structure
snippy_results/
├── TE15676_snippy/
│   └── TE15676_snippy.vcf    # FreeBayes variants
├── TE15677_snippy/
│   └── TE15677_snippy.vcf
└── ...

cfsan_results/
└── samples/
    ├── TE15676/
    │   └── var.flt.vcf       # VarScan2 variants
    ├── TE15677/
    │   └── var.flt.vcf
    └── ...</pre>

                    <h4 style="color: #555; margin-top: 20px; margin-bottom: 10px;">Correlation with Gaps</h4>
                    <p style="line-height: 1.8;">
                        For each pipeline-unique variant, we check if it falls within a coreguard gap region.
                        This helps explain <em>why</em> the other pipeline might have missed it:
                    </p>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li><strong>In gap region</strong>: Low coverage may explain why the other pipeline didn't call the variant</li>
                        <li><strong>Not in gap</strong>: Difference is due to variant calling algorithm (FreeBayes vs VarScan2 sensitivity)</li>
                    </ul>
                </div>

                <!-- Complexity -->
                <div style="background: #e8f5e9; border-radius: 8px; padding: 20px;">
                    <h4 style="color: #2e7d32; margin-bottom: 10px;">Computational Complexity</h4>
                    <ul style="line-height: 2; margin-left: 20px;">
                        <li><strong>Reference alignment</strong>: O(n) minimap2 calls (n = number of samples)</li>
                        <li><strong>Pairwise analysis</strong>: O(n&sup2;) set operations (fast, in-memory bitmaps)</li>
                        <li><strong>Total</strong>: Dominated by O(n) minimap2 calls</li>
                    </ul>
                    <p style="margin-top: 10px; font-style: italic; color: #555;">
                        Compare to traditional SNP pipelines: O(k &times; pipeline_cost) where k = number of iteration cycles.
                    </p>
                </div>
            </div>
        </div>

        <div class="panel" id="files">
            <!-- Run Info Section -->
            <div id="run-info" style="background: #e8f4fd; border-radius: 8px; padding: 20px; margin-bottom: 25px; border-left: 4px solid #2196f3;">
                <h3 style="margin-bottom: 15px; color: #1565c0;">Run Information</h3>
                <div id="run-info-content"></div>
            </div>

            <h3 style="margin-bottom: 15px;">Download Files</h3>
            <p style="color: #6c757d; margin-bottom: 20px;">
                Download generated files for use in external tools (IGV desktop, custom analysis scripts, etc.)
            </p>
            <div style="display: grid; gap: 15px; max-width: 800px;">
                <div style="background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="color: #2c3e50; margin-bottom: 15px;">Report Files</h4>
                    <div style="display: flex; flex-wrap: wrap; gap: 10px;">
                        <a href="/report.json" download="coreguard_report.json"
                           style="display: inline-flex; align-items: center; padding: 10px 20px; background: #3498db; color: white; text-decoration: none; border-radius: 4px;">
                            <span style="margin-right: 8px;">&#128196;</span> report.json
                        </a>
                        <a href="/gaps.svg" download="gaps_distribution.svg"
                           style="display: inline-flex; align-items: center; padding: 10px 20px; background: #27ae60; color: white; text-decoration: none; border-radius: 4px;">
                            <span style="margin-right: 8px;">&#128200;</span> gaps_distribution.svg
                        </a>
                    </div>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="color: #2c3e50; margin-bottom: 15px;">Pairwise Analysis</h4>
                    <p style="color: #666; font-size: 0.9rem; margin-bottom: 15px;">
                        Gap union/intersection metrics for each sample pair.
                    </p>
                    <div style="display: flex; flex-wrap: wrap; gap: 10px;">
                        <a href="/pairwise_gaps.json" download="pairwise_gaps.json"
                           style="display: inline-flex; align-items: center; padding: 10px 20px; background: #e67e22; color: white; text-decoration: none; border-radius: 4px;">
                            <span style="margin-right: 8px;">&#128203;</span> pairwise_gaps.json
                        </a>
                    </div>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="color: #2c3e50; margin-bottom: 15px;">BED Files (for IGV)</h4>
                    <p style="color: #666; font-size: 0.9rem; margin-bottom: 15px;">
                        BED files can be loaded into IGV desktop for detailed visualization.
                    </p>
                    <div id="bed-file-list" style="display: flex; flex-wrap: wrap; gap: 10px;"></div>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="color: #2c3e50; margin-bottom: 15px;">Reference Files</h4>
                    <div style="display: flex; flex-wrap: wrap; gap: 10px;">
                        <a href="/reference.fasta" download="reference.fasta"
                           style="display: inline-flex; align-items: center; padding: 10px 20px; background: #9b59b6; color: white; text-decoration: none; border-radius: 4px;">
                            <span style="margin-right: 8px;">&#129516;</span> reference.fasta
                        </a>
                        <a href="/reference.fasta.fai" download="reference.fasta.fai"
                           style="display: inline-flex; align-items: center; padding: 10px 20px; background: #8e44ad; color: white; text-decoration: none; border-radius: 4px;">
                            <span style="margin-right: 8px;">&#128209;</span> reference.fasta.fai
                        </a>
                    </div>
                </div>

                <div style="background: #fff3cd; border-radius: 8px; padding: 15px; border-left: 4px solid #ffc107;">
                    <strong>Tip:</strong> To use these files in IGV desktop:
                    <ol style="margin: 10px 0 0 20px;">
                        <li>Download the reference FASTA and FAI files</li>
                        <li>Load the reference genome in IGV (Genomes &gt; Load Genome from File)</li>
                        <li>Download and load the BED file(s) (File &gt; Load from File)</li>
                    </ol>
                </div>
            </div>
        </div>

        <div class="panel" id="help">
            <div style="max-width: 900px; margin: 0 auto;">
                <h2 style="color: #2c3e50; margin-bottom: 20px;">coreguard - Help &amp; Documentation</h2>

                <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #3498db; margin-bottom: 15px;">What is coreguard?</h3>
                    <p style="line-height: 1.8;">
                        <strong>coreguard</strong> is a pre-alignment QC tool designed to identify problematic samples
                        <em>before</em> running SNP pipelines. It helps prevent the inclusion of samples that would
                        introduce gaps, artifacts, or inconsistencies in multi-sample SNP analysis.
                    </p>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #3498db; margin-bottom: 15px;">Key Concepts</h3>
                    <dl style="line-height: 1.8;">
                        <dt style="font-weight: bold; margin-top: 10px;">Gap Contagion</dt>
                        <dd style="margin-left: 20px; color: #555;">When a sample has unique gaps (regions not aligned to reference),
                        those positions are excluded from the "core genome" for ALL samples, reducing SNP resolution.</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">Core Genome</dt>
                        <dd style="margin-left: 20px; color: #555;">The set of genomic positions shared by all samples.
                        SNP analysis is performed only on core positions.</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">Risk Level</dt>
                        <dd style="margin-left: 20px; color: #555;">Based on how much core genome reduction a sample causes:
                            <ul style="margin-top: 10px;">
                                <li><span style="color: #388e3c; font-weight: bold;">Low</span>: &lt;0.5% reduction - minimal impact</li>
                                <li><span style="color: #1976d2; font-weight: bold;">Medium</span>: 0.5-2% reduction - review recommended</li>
                                <li><span style="color: #f57c00; font-weight: bold;">High</span>: 2-5% reduction - consider excluding</li>
                                <li><span style="color: #d32f2f; font-weight: bold;">Critical</span>: &gt;5% reduction - strongly recommend excluding</li>
                            </ul>
                        </dd>
                    </dl>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #3498db; margin-bottom: 15px;">Understanding the Tabs</h3>
                    <dl style="line-height: 1.8;">
                        <dt style="font-weight: bold; margin-top: 10px;">Gaps Map</dt>
                        <dd style="margin-left: 20px; color: #555;">SVG visualization showing unique gaps for each sample
                        mapped to reference coordinates. Samples sorted by impact (best first).</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">Pairwise Gaps</dt>
                        <dd style="margin-left: 20px; color: #555;">Shows gap overlap metrics for each sample pair.
                        Quality Score = intersection/union. High = similar gaps, Low = different gaps (more core reduction when combined).</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">IGV Browser</dt>
                        <dd style="margin-left: 20px; color: #555;">Interactive genome browser (IGV.js) for exploring gap regions
                        in detail. Navigate to specific coordinates or use quick navigation buttons.</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">Summary</dt>
                        <dd style="margin-left: 20px; color: #555;">Interactive sample selector. Check/uncheck samples to see
                        how the estimated core genome changes dynamically.</dd>

                        <dt style="font-weight: bold; margin-top: 15px;">Sample Details</dt>
                        <dd style="margin-left: 20px; color: #555;">Table with per-sample metrics including risk level,
                        unique gap bases, and core reduction percentage.</dd>
                    </dl>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #3498db; margin-bottom: 15px;">Command Line Usage</h3>
                    <pre style="background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 4px; overflow-x: auto; font-size: 0.9rem;">
# Basic analysis
coreguard --input samples/ --reference ref.fasta --output report.json

# With web server (opens browser automatically)
coreguard --input samples/ --reference ref.fasta --output report.json --serve

# Serve existing report without re-running analysis
coreguard --serve-only report.json --reference ref.fasta

# Fast mode (reference-only, skips pairwise - much faster for large datasets)
coreguard --input samples/ --reference ref.fasta --fast --output report.json

# Fast mode but still compute identity matrix
coreguard --input samples/ --reference ref.fasta --fast --with-pairwise --output report.json
                    </pre>
                </div>

                <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #3498db; margin-bottom: 15px;">Command Line Options Reference</h3>
                    <pre style="background: #1a252f; color: #ecf0f1; padding: 15px; border-radius: 4px; overflow-x: auto; font-size: 0.85rem; line-height: 1.6;">
<span style="color: #3498db; font-weight: bold;">INPUT/OUTPUT:</span>
  -i, --input &lt;DIR&gt;           Input directory containing FASTA files
  -s, --samples &lt;FILES&gt;...    Individual FASTA files to analyze
  -r, --reference &lt;FILE&gt;      Reference genome (recommended for accurate gap prediction)
  -o, --output &lt;FILE&gt;         Output file (JSON or TSV based on extension) [default: report.json]
      --bed-output &lt;DIR&gt;      Output BED files for IGV visualization

<span style="color: #3498db; font-weight: bold;">ANALYSIS OPTIONS:</span>
  -t, --threads &lt;N&gt;           Number of threads for parallel processing [default: all CPUs]
      --min-identity &lt;FLOAT&gt;  Minimum identity threshold (flag samples below) [default: 0.9]
      --minimap2-path &lt;PATH&gt;  Path to minimap2 binary [default: minimap2 in PATH]
      --fast                  Fast mode: skip pairwise, use reference-only (N vs N² alignments)
      --with-pairwise         Force pairwise analysis even with --fast (for identity matrix)
      --no-dedup              Disable deduplication (by default, identical sequences analyzed once)

<span style="color: #3498db; font-weight: bold;">WEB SERVER:</span>
      --serve                 Start web server with IGV.js genome browser (opens browser)
      --port &lt;N&gt;              Port for the web server [default: 8765]
      --serve-only &lt;FILE&gt;     Serve existing report without re-running analysis

<span style="color: #3498db; font-weight: bold;">GENERAL:</span>
  -v, --verbose               Verbose output
  -h, --help                  Print help
  -V, --version               Print version
                    </pre>
                </div>

                <div style="background: #e3f2fd; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h3 style="color: #1565c0; margin-bottom: 15px;">Citation</h3>
                    <p style="line-height: 1.8;">
                        If you use coreguard in your research, please cite:<br>
                        <em>TBW (To Be Written)</em>
                    </p>
                    <p style="margin-top: 10px;">
                        <a href="https://github.com/genpat-it/coreguard" target="_blank" style="color: #1565c0;">
                            GitHub Repository
                        </a>
                    </p>
                </div>
            </div>
        </div>
    </div>

    <footer>
        Generated by coreguard | <a href="https://github.com/genpat-it/coreguard">GitHub</a>
    </footer>

    <!-- Back to top button -->
    <button id="back-to-top" onclick="window.scrollTo({{top: 0, behavior: 'smooth'}})"
            style="position: fixed; bottom: 30px; right: 30px; width: 50px; height: 50px;
                   border-radius: 50%; border: none; background: #3498db; color: white;
                   font-size: 24px; cursor: pointer; box-shadow: 0 4px 12px rgba(0,0,0,0.3);
                   display: none; z-index: 1000; transition: all 0.3s ease;"
            title="Torna in alto">
        &#8593;
    </button>

    <script>
    const data = {report_json};
    let browser = null;

    // Suppress IGV alerts and console warnings
    const originalAlert = window.alert;
    window.alert = function(msg) {{
        console.log("IGV alert suppressed:", msg);
    }};
    const originalWarn = console.warn;
    console.warn = function(...args) {{
        const msg = args.join(' ');
        if (msg.includes('range header') || msg.includes('Range header')) {{
            return; // Silently ignore range header warnings
        }}
        originalWarn.apply(console, args);
    }};

    // Initialize IGV (lazy - only when tab is clicked)
    let igvInitialized = false;
    async function initIGV() {{
        if (igvInitialized) return;
        igvInitialized = true;

        const options = {{
            genome: {{
                id: "custom",
                name: "{reference_name}",
                fastaURL: "/reference.fasta",
                indexURL: "/reference.fasta.fai",
                tracks: {tracks_js}
            }},
            locus: "{reference_name}:1-{half_len}",
            showNavigation: true,
            showRuler: true,
            showCenterGuide: true
        }};

        try {{
            browser = await igv.createBrowser(document.getElementById('igv-container'), options);
            console.log("IGV browser initialized");
        }} catch (e) {{
            console.error("IGV init error:", e);
            document.getElementById('igv-container').innerHTML =
                '<div style="padding: 40px; text-align: center; color: #dc3545;">' +
                '<h3>IGV.js could not load</h3>' +
                '<p>Reference genome may be too large or missing. ' +
                'Try loading the BED files manually in IGV desktop.</p></div>';
        }}
    }}

    function goToLocus() {{
        const locus = document.getElementById('locus-input').value;
        if (browser && locus) {{
            browser.search(locus);
        }}
    }}

    // Quick navigation for high-impact regions
    function setupQuickNav() {{
        const nav = document.getElementById('quick-nav');
        const hotspots = [];

        // Find first significant gap region for each sample (one button per sample)
        data.impact_analysis.samples.forEach(s => {{
            if (s.unique_gap_coords && s.unique_gap_coords.length > 0) {{
                // Find the largest gap for this sample
                const largestGap = s.unique_gap_coords.reduce((max, coord) =>
                    (coord.end - coord.start) > (max.end - max.start) ? coord : max
                );
                hotspots.push({{
                    label: `${{s.sample}} (${{s.unique_gap_coords.length}} gaps)`,
                    start: largestGap.start,
                    end: largestGap.end,
                    gapSize: largestGap.end - largestGap.start
                }});
            }}
        }});

        // Sort by gap size (largest first) and limit
        hotspots.sort((a, b) => b.gapSize - a.gapSize);
        const limited = hotspots.slice(0, 10);

        // Add "All" button
        nav.innerHTML = '<button onclick="if(browser) browser.search(&#39;all&#39;)">View All</button>';

        limited.forEach(h => {{
            const btn = document.createElement('button');
            // Extract sample name (remove " (N gaps)" suffix)
            const sampleName = h.label.split(' (')[0];
            btn.textContent = h.label;
            btn.onclick = () => {{
                if (browser) {{
                    // Navigate to the gap region
                    const region = `{reference_name}:${{Math.max(1, h.start - 1000)}}-${{h.end + 1000}}`;
                    browser.search(region);

                    // Filter tracks: show only this sample's tracks + reference
                    browser.trackViews.forEach(tv => {{
                        const trackName = tv.track.name || '';
                        // Show reference tracks (sequence, ruler, etc.) and the selected sample
                        const isReference = trackName.toLowerCase().includes('sequence') ||
                                          trackName.toLowerCase().includes('ruler') ||
                                          trackName.toLowerCase().includes('reference') ||
                                          trackName.toLowerCase().includes('gene') ||
                                          tv.track.type === 'sequence' ||
                                          tv.track.type === 'ruler';
                        const isSelectedSample = trackName.includes(sampleName);
                        const shouldShow = isReference || isSelectedSample;

                        // Use display style for reliable visibility control
                        if (tv.viewportDiv) {{
                            tv.viewportDiv.style.display = shouldShow ? '' : 'none';
                        }}
                        if (tv.trackDiv) {{
                            tv.trackDiv.style.display = shouldShow ? '' : 'none';
                        }}
                    }});
                }}
            }};
            nav.appendChild(btn);
        }});

        // Add "Show All Tracks" button after sample buttons
        const showAllBtn = document.createElement('button');
        showAllBtn.textContent = '↻ All Tracks';
        showAllBtn.style.background = '#d4edda';
        showAllBtn.onclick = () => {{
            if (browser) {{
                // Show all tracks by removing display:none from all track-related elements
                document.querySelectorAll('.igv-track-container, .igv-viewport, [class*="igv"]').forEach(el => {{
                    if (el.style.display === 'none') {{
                        el.style.display = '';
                    }}
                }});
                // Also iterate trackViews
                browser.trackViews.forEach(tv => {{
                    if (tv.viewportDiv) tv.viewportDiv.style.display = '';
                    if (tv.trackDiv) tv.trackDiv.style.display = '';
                    if (tv.outerDiv) tv.outerDiv.style.display = '';
                    // Try to find parent containers
                    let parent = tv.viewportDiv ? tv.viewportDiv.parentElement : null;
                    while (parent && parent.id !== 'igv-container') {{
                        if (parent.style.display === 'none') parent.style.display = '';
                        parent = parent.parentElement;
                    }}
                }});
                console.log('All tracks restored:', browser.trackViews.length);
            }}
        }};
        nav.appendChild(showAllBtn);
    }}

    // Tab switching
    document.querySelectorAll('.tab').forEach(tab => {{
        tab.addEventListener('click', () => {{
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById(tab.dataset.panel).classList.add('active');

            // Initialize IGV only when tab is clicked (lazy loading)
            if (tab.dataset.panel === 'igv') {{
                initIGV();
            }}
        }});
    }});

    // Sort samples by unique_gap_bases (best = lowest gaps first)
    const sortedSamples = [...data.impact_analysis.samples].sort((a, b) => a.unique_gap_bases - b.unique_gap_bases);
    const refLen = data.impact_analysis.reference_length;

    // Track selected samples
    let selectedSamples = new Set(sortedSamples.map(s => s.sample));

    // Calculate core for selected samples
    function calculateCore() {{
        // Collect all unique gap positions from selected samples
        const allGapPositions = new Set();
        sortedSamples.forEach(s => {{
            if (selectedSamples.has(s.sample)) {{
                s.unique_gap_coords.forEach(coord => {{
                    for (let i = coord.start; i < coord.end; i++) {{
                        allGapPositions.add(i);
                    }}
                }});
            }}
        }});
        const totalGapBases = allGapPositions.size;
        const coreSize = refLen - totalGapBases;
        const corePct = (coreSize / refLen) * 100;
        return {{ coreSize, corePct, totalGapBases }};
    }}

    function updateCoreResult() {{
        const result = calculateCore();
        const resultDiv = document.getElementById('core-result');
        const selectedCount = selectedSamples.size;
        const color = result.corePct > 95 ? '#28a745' : result.corePct > 90 ? '#ffc107' : '#dc3545';

        // Update result box
        resultDiv.innerHTML = `
            <strong>Selected:</strong> ${{selectedCount}} / ${{sortedSamples.length}} samples<br>
            <strong style="color: ${{color}};">Estimated Core: ${{(result.coreSize / 1e6).toFixed(3)}} Mb (${{result.corePct.toFixed(3)}}%)</strong><br>
            <span style="color: #6c757d;">Gap bases: ${{result.totalGapBases.toLocaleString()}} bp</span>
        `;
        resultDiv.style.background = result.corePct > 95 ? '#e8f5e9' : result.corePct > 90 ? '#fff8e1' : '#ffebee';

        // Update KPI cards
        const cards = document.getElementById('summary-cards');
        const coreColorClass = result.corePct > 95 ? 'low' : result.corePct > 90 ? 'medium' : 'critical';
        cards.innerHTML = `
            <div class="card">
                <div class="card-title">Selected Samples</div>
                <div class="card-value">${{selectedCount}} / ${{sortedSamples.length}}</div>
            </div>
            <div class="card">
                <div class="card-title">Reference Length</div>
                <div class="card-value">${{(refLen / 1e6).toFixed(3)}} Mb</div>
            </div>
            <div class="card">
                <div class="card-title">Estimated Core</div>
                <div class="card-value ${{coreColorClass}}">${{(result.coreSize / 1e6).toFixed(3)}} Mb</div>
            </div>
            <div class="card">
                <div class="card-title">Core %</div>
                <div class="card-value ${{coreColorClass}}">${{result.corePct.toFixed(3)}}%</div>
            </div>
        `;
    }}

    // Populate summary with sample selector
    function populateSummary() {{
        const cards = document.getElementById('summary-cards');
        const s = data.summary;
        const flaggedClass = s.flagged_samples > 0 ? 'critical' : 'low';

        cards.innerHTML = `
            <div class="card">
                <div class="card-title">Total Samples</div>
                <div class="card-value">${{s.total_samples}}</div>
            </div>
            <div class="card">
                <div class="card-title">Reference Length</div>
                <div class="card-value">${{(refLen / 1e6).toFixed(3)}} Mb</div>
            </div>
            <div class="card">
                <div class="card-title">Core (all samples)</div>
                <div class="card-value">${{(data.impact_analysis.estimated_core_all / 1e6).toFixed(3)}} Mb</div>
            </div>
            <div class="card">
                <div class="card-title">Core %</div>
                <div class="card-value">${{((1 - data.impact_analysis.total_core_reduction) * 100).toFixed(1)}}%</div>
            </div>
        `;

        // Sample selector with checkboxes
        const selector = document.getElementById('sample-selector');
        selector.innerHTML = `
            <div style="margin-bottom: 10px;">
                <button onclick="selectAll()" style="margin-right: 10px; padding: 5px 15px;">Select All</button>
                <button onclick="selectNone()" style="padding: 5px 15px;">Select None</button>
            </div>
        ` + sortedSamples.map((s, idx) => {{
            const color = s.risk_level === 'Critical' ? '#dc3545' :
                          s.risk_level === 'High' ? '#ffc107' :
                          s.risk_level === 'Medium' ? '#17a2b8' : '#28a745';
            const barWidth = Math.max(5, (s.unique_gap_bases / Math.max(...sortedSamples.map(x => x.unique_gap_bases))) * 100);
            return `
                <div style="display: flex; align-items: center; margin: 8px 0; padding: 8px; background: #f8f9fa; border-radius: 4px;">
                    <input type="checkbox" id="chk_${{s.sample}}" checked onchange="toggleSample('${{s.sample}}')" style="margin-right: 10px; width: 18px; height: 18px;">
                    <div style="width: 100px; font-weight: 500;">${{s.sample}}</div>
                    <div style="flex: 1; height: 20px; background: #e9ecef; border-radius: 4px; overflow: hidden; margin: 0 10px;">
                        <div style="width: ${{barWidth}}%; height: 100%; background: ${{color}};"></div>
                    </div>
                    <div style="width: 100px; font-size: 0.85rem; color: ${{color}};">${{s.unique_gap_bases.toLocaleString()}} bp</div>
                    <span class="badge badge-${{s.risk_level.toLowerCase()}}" style="width: 70px; text-align: center;">${{s.risk_level}}</span>
                </div>
            `;
        }}).join('');

        updateCoreResult();
    }}

    function toggleSample(sample) {{
        if (selectedSamples.has(sample)) {{
            selectedSamples.delete(sample);
        }} else {{
            selectedSamples.add(sample);
        }}
        updateCoreResult();
    }}

    function selectAll() {{
        sortedSamples.forEach(s => {{
            selectedSamples.add(s.sample);
            document.getElementById('chk_' + s.sample).checked = true;
        }});
        updateCoreResult();
    }}

    function selectNone() {{
        selectedSamples.clear();
        sortedSamples.forEach(s => {{
            document.getElementById('chk_' + s.sample).checked = false;
        }});
        updateCoreResult();
    }}

    // Populate samples table (sorted by impact)
    function populateSamples() {{
        const tbody = document.querySelector('#samples-table tbody');

        // Build a map of sample -> uncovered bases from reference_gaps
        const uncoveredMap = {{}};
        (data.reference_gaps || []).forEach(rg => {{
            const totalUncovered = (rg.reference_uncovered || []).reduce((sum, r) => sum + (r.end - r.start), 0);
            uncoveredMap[rg.sample] = totalUncovered;
        }});

        sortedSamples.forEach(s => {{
            const riskClass = s.risk_level.toLowerCase();
            const uncovered = uncoveredMap[s.sample] || 0;
            const covered = refLen - uncovered;
            const coveragePct = (covered / refLen * 100).toFixed(1);

            const row = document.createElement('tr');
            row.innerHTML = `
                <td><strong>${{s.sample}}</strong></td>
                <td>${{covered.toLocaleString()}}</td>
                <td>${{coveragePct}}%</td>
                <td><span class="badge badge-${{riskClass}}">${{s.risk_level}}</span></td>
                <td>${{s.unique_gap_bases.toLocaleString()}} bp</td>
                <td>${{(s.estimated_core_reduction * 100).toFixed(3)}}%</td>
            `;
            tbody.appendChild(row);
        }});
    }}

    // Pairwise gaps visualization (NEW: uses gap union/intersection metrics)
    function populatePairwise() {{
        const container = document.getElementById('pairwise-container');
        const pairs = data.pairwise_gap_results || [];

        console.log('Pairwise pairs:', pairs.length);

        // Sort pairs by quality_score (worst first = lowest score)
        const sortedPairs = [...pairs].sort((a, b) => a.quality_score - b.quality_score);

        container.innerHTML = sortedPairs.map((pair, idx) => {{
            const unionBp = pair.gap_union_bases || 0;
            const intersectionBp = pair.gap_intersection_bases || 0;
            const qualityScore = pair.quality_score || 0;
            const uniqueA = pair.unique_gap_a_bases || 0;
            const uniqueB = pair.unique_gap_b_bases || 0;
            const coreSize = pair.pairwise_core_size || refLen;
            const coreReduction = ((refLen - coreSize) / refLen * 100).toFixed(3);

            // Quality interpretation
            const qualityClass = qualityScore > 0.8 ? 'good' :
                                 qualityScore > 0.5 ? 'moderate' :
                                 qualityScore > 0.2 ? 'poor' : 'bad';

            const qualityColor = qualityScore > 0.8 ? '#28a745' :
                                 qualityScore > 0.5 ? '#17a2b8' :
                                 qualityScore > 0.2 ? '#ffc107' : '#dc3545';

            const qualityBg = qualityScore > 0.8 ? '#d4edda' :
                              qualityScore > 0.5 ? '#d1ecf1' :
                              qualityScore > 0.2 ? '#fff3cd' : '#f8d7da';

            const qualityText = qualityScore > 0.8 ? 'High (similar gaps)' :
                                qualityScore > 0.5 ? 'Moderate' :
                                qualityScore > 0.2 ? 'Low (different gaps)' : 'Very Low (no overlap)';

            // Visualize union vs intersection
            const unionWidth = unionBp > 0 ? 100 : 0;
            const intersectionWidth = unionBp > 0 ? (intersectionBp / unionBp * 100) : 0;

            return `
                <div class="pair-card ${{qualityClass}}" data-quality="${{qualityScore}}" data-core-reduction="${{coreReduction}}"
                     style="background: white; border-radius: 8px; padding: 15px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                        <h4 style="margin: 0;">${{pair.sample_a}} vs ${{pair.sample_b}}</h4>
                        <span style="padding: 4px 12px; border-radius: 20px; font-size: 0.8rem;
                                     background: ${{qualityBg}}; color: ${{qualityColor}}; font-weight: bold;">
                            Quality: ${{(qualityScore * 100).toFixed(1)}}%
                        </span>
                    </div>

                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 10px;">
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 4px;">
                            <div style="font-size: 0.8rem; color: #666;">Gap Union (A ∪ B)</div>
                            <div style="font-size: 1.2rem; font-weight: bold;">${{unionBp.toLocaleString()}} bp</div>
                            <div style="font-size: 0.75rem; color: #999;">Core reduction: ${{coreReduction}}%</div>
                        </div>
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 4px;">
                            <div style="font-size: 0.8rem; color: #666;">Gap Intersection (A ∩ B)</div>
                            <div style="font-size: 1.2rem; font-weight: bold;">${{intersectionBp.toLocaleString()}} bp</div>
                            <div style="font-size: 0.75rem; color: #999;">${{qualityText}}</div>
                        </div>
                    </div>

                    <div style="margin-bottom: 10px;">
                        <div style="font-size: 0.8rem; color: #666; margin-bottom: 5px;">Gap Overlap Visualization</div>
                        <div style="position: relative; height: 30px; background: #e9ecef; border-radius: 4px; overflow: hidden;">
                            <div style="position: absolute; width: ${{unionWidth}}%; height: 100%; background: #dc3545; opacity: 0.5;" title="Union"></div>
                            <div style="position: absolute; width: ${{intersectionWidth}}%; height: 100%; background: #28a745;" title="Intersection"></div>
                            <div style="position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%); font-size: 0.8rem; font-weight: bold; color: white; text-shadow: 0 0 2px black;">
                                ${{intersectionWidth.toFixed(0)}}% overlap
                            </div>
                        </div>
                    </div>

                    <div style="font-size: 0.8rem; color: #666; display: flex; gap: 15px;">
                        <span><span style="color: #d32f2f;">■</span> Unique to ${{pair.sample_a}}: ${{uniqueA.toLocaleString()}} bp</span>
                        <span><span style="color: #1976d2;">■</span> Unique to ${{pair.sample_b}}: ${{uniqueB.toLocaleString()}} bp</span>
                        <span><span style="color: #28a745;">■</span> Pairwise core: ${{(coreSize / 1e6).toFixed(3)}} Mb</span>
                    </div>
                </div>
            `;
        }}).join('');

        document.getElementById('pair-count').textContent = `${{pairs.length}} pairs`;
    }}

    function filterPairwise() {{
        const reductionFilter = document.getElementById('pairwise-filter').value;
        const qualityFilter = document.getElementById('quality-filter').value;
        const cards = document.querySelectorAll('.pair-card');
        let visible = 0;
        let total = cards.length;

        cards.forEach(card => {{
            const quality = parseFloat(card.dataset.quality) || 0;
            const coreReduction = parseFloat(card.dataset.coreReduction) || 0;

            // Core reduction filter
            let passReduction = false;
            if (reductionFilter === 'all') passReduction = true;
            else if (reductionFilter === 'with-gaps') passReduction = coreReduction > 0;
            else if (reductionFilter === 'high-gaps') passReduction = coreReduction > 1;
            else if (reductionFilter === 'critical-gaps') passReduction = coreReduction > 5;
            else if (reductionFilter === 'very-high-gaps') passReduction = coreReduction > 10;
            else if (reductionFilter === 'extreme-gaps') passReduction = coreReduction > 15;

            // Quality score filter (quality is 0-1, convert to percentage for comparison)
            let passQuality = false;
            const qualityPct = quality * 100;
            if (qualityFilter === 'all') passQuality = true;
            else if (qualityFilter === 'high-quality') passQuality = qualityPct > 80;
            else if (qualityFilter === 'medium-quality') passQuality = qualityPct >= 50 && qualityPct <= 80;
            else if (qualityFilter === 'low-quality') passQuality = qualityPct < 50;
            else if (qualityFilter === 'very-low-quality') passQuality = qualityPct < 20;

            const show = passReduction && passQuality;
            card.style.display = show ? 'block' : 'none';
            if (show) visible++;
        }});

        document.getElementById('pair-count').textContent = `${{visible}} / ${{total}} pairs shown`;
    }}

    // Populate Home page with summary stats
    function populateHome() {{
        document.getElementById('home-total-samples').textContent = data.summary.total_samples;
        document.getElementById('home-reference').textContent = data.reference_summary.reference_name;

        const coreMb = (data.impact_analysis.estimated_core_all / 1e6).toFixed(3);
        const corePct = ((1 - data.impact_analysis.total_core_reduction) * 100).toFixed(3);
        document.getElementById('home-core').innerHTML = `${{coreMb}} Mb<br><span style="font-size: 0.8rem; color: #666;">(${{corePct}}%)</span>`;

        const flagged = data.summary.flagged_samples;
        const flaggedEl = document.getElementById('home-flagged');
        flaggedEl.textContent = flagged;
        if (flagged > 0) {{
            flaggedEl.classList.add('critical');
        }} else {{
            flaggedEl.classList.add('low');
        }}
    }}

    // Populate Run Information
    function populateRunInfo() {{
        const container = document.getElementById('run-info-content');
        const meta = data.metadata || {{}};

        const formatDuration = (secs) => {{
            if (!secs) return 'N/A';
            if (secs < 60) return `${{secs.toFixed(1)}}s`;
            const mins = Math.floor(secs / 60);
            const remainSecs = (secs % 60).toFixed(1);
            return `${{mins}}m ${{remainSecs}}s`;
        }};

        let html = `
            <table style="width: 100%; font-size: 0.9rem; border-collapse: collapse;">
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666; white-space: nowrap; vertical-align: top;"><strong>Command:</strong></td>
                    <td style="padding: 8px 0; font-family: monospace; word-break: break-all; background: #f8f9fa; padding: 8px; border-radius: 4px;">
                        ${{meta.command_line || 'N/A'}}
                    </td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Timestamp:</strong></td>
                    <td style="padding: 8px 0;">${{meta.timestamp || 'N/A'}}</td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Duration:</strong></td>
                    <td style="padding: 8px 0;">${{formatDuration(meta.duration_secs)}}</td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Input Mode:</strong></td>
                    <td style="padding: 8px 0;">${{meta.input_mode || 'FASTA'}}</td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Input Path:</strong></td>
                    <td style="padding: 8px 0; font-family: monospace;">${{meta.input_path || 'N/A'}}</td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Reference:</strong></td>
                    <td style="padding: 8px 0; font-family: monospace;">${{meta.reference_path || 'N/A'}}</td>
                </tr>
                <tr>
                    <td style="padding: 8px 15px 8px 0; color: #666;"><strong>Version:</strong></td>
                    <td style="padding: 8px 0;">${{data.version || 'N/A'}}</td>
                </tr>
            </table>
        `;

        container.innerHTML = html;
    }}

    // Populate BED file download list
    async function populateBedFiles() {{
        const container = document.getElementById('bed-file-list');
        try {{
            const response = await fetch('/bed_list');
            const files = await response.json();

            if (files.length === 0) {{
                container.innerHTML = '<p style="color: #666;">No BED files available</p>';
                return;
            }}

            container.innerHTML = files.map(f => `
                <a href="/bed/${{f}}" download="${{f}}"
                   style="display: inline-flex; align-items: center; padding: 10px 20px; background: #e74c3c; color: white; text-decoration: none; border-radius: 4px;">
                    <span style="margin-right: 8px;">&#128203;</span> ${{f}}
                </a>
            `).join('');
        }} catch (e) {{
            container.innerHTML = '<p style="color: #dc3545;">Error loading BED files</p>';
        }}
    }}

    // Populate Gene Zone analysis
    function populateGeneZone() {{
        const container = document.getElementById('gene-zone-content');
        const gz = data.gene_zone_analysis;

        if (!gz) {{
            container.innerHTML = `
                <div style="background: #fff3cd; padding: 20px; border-radius: 8px;">
                    <h4 style="color: #856404;">No GFF Annotation Provided</h4>
                    <p>To see which genes are affected by alignment gaps, run coreguard with the <code>--gff</code> option:</p>
                    <pre style="background: #2c3e50; color: #ecf0f1; padding: 10px; border-radius: 4px; margin-top: 10px;">
coreguard --input samples/ --reference ref.fasta --gff annotation.gff --output report.json</pre>
                </div>
            `;
            return;
        }}

        let html = `
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 25px;">
                <div class="card">
                    <div class="card-title">Total Genes</div>
                    <div class="card-value">${{gz.total_genes.toLocaleString()}}</div>
                </div>
                <div class="card">
                    <div class="card-title">Genes Affected</div>
                    <div class="card-value" style="color: #f57c00;">${{gz.total_affected_genes.toLocaleString()}}</div>
                    <div style="font-size: 0.8rem; color: #666;">${{(gz.total_affected_genes / gz.total_genes * 100).toFixed(1)}}% of total</div>
                </div>
                <div class="card">
                    <div class="card-title">Partially Affected</div>
                    <div class="card-value" style="color: #ff9800;">${{(gz.total_affected_genes - gz.genes_removed_any_sample).toLocaleString()}}</div>
                    <div style="font-size: 0.8rem; color: #666;">Have gaps but not 100%</div>
                </div>
                <div class="card">
                    <div class="card-title">Genes Removed</div>
                    <div class="card-value" style="color: #d32f2f;">${{gz.genes_removed_any_sample.toLocaleString()}}</div>
                    <div style="font-size: 0.8rem; color: #666;">100% in gap in &ge;1 sample</div>
                </div>
                <div class="card">
                    <div class="card-title">GFF File</div>
                    <div style="font-size: 0.8rem; word-break: break-all;">${{gz.gff_path.split('/').pop()}}</div>
                </div>
            </div>

            <h4 style="margin-bottom: 15px;">Per-Sample Gene Impact</h4>
            <div style="overflow-x: auto;">
                <table style="width: 100%; border-collapse: collapse; margin-bottom: 25px;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 12px; text-align: right; border: 1px solid #dee2e6;">Genes Affected</th>
                            <th style="padding: 12px; text-align: right; border: 1px solid #dee2e6;">Partially Affected</th>
                            <th style="padding: 12px; text-align: right; border: 1px solid #dee2e6;">Significantly (&gt;50%)</th>
                            <th style="padding: 12px; text-align: right; border: 1px solid #dee2e6;">Completely Removed</th>
                            <th style="padding: 12px; text-align: right; border: 1px solid #dee2e6;">CDS Affected</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const s of gz.samples) {{
            const pct = (s.affected_genes / s.total_genes * 100).toFixed(1);
            const partiallyAffected = s.affected_genes - s.removed_genes;
            html += `
                <tr>
                    <td style="padding: 12px; border: 1px solid #dee2e6; font-weight: bold;">${{s.sample}}</td>
                    <td style="padding: 12px; border: 1px solid #dee2e6; text-align: right;">${{s.affected_genes}} (${{pct}}%)</td>
                    <td style="padding: 12px; border: 1px solid #dee2e6; text-align: right; color: #ff9800;">${{partiallyAffected}}</td>
                    <td style="padding: 12px; border: 1px solid #dee2e6; text-align: right; color: #f57c00;">${{s.significantly_affected}}</td>
                    <td style="padding: 12px; border: 1px solid #dee2e6; text-align: right; color: #d32f2f;">${{s.removed_genes}}</td>
                    <td style="padding: 12px; border: 1px solid #dee2e6; text-align: right;">${{s.affected_cds}}</td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <h4 style="margin-bottom: 15px;">Affected Genes</h4>
            <div id="gene-table-container"></div>
        `;

        container.innerHTML = html;

        // Collect all affected genes across samples
        const allAffected = [];
        for (const s of gz.samples) {{
            for (const af of s.affected_features) {{
                if (af.feature.feature_type === 'gene') {{
                    allAffected.push({{ ...af, sample: s.sample }});
                }}
            }}
        }}

        // Sort by overlap percentage
        allAffected.sort((a, b) => b.overlap_pct - a.overlap_pct);

        // Store globally for pagination and filtering
        window.geneTableDataAll = allAffected;  // Original unfiltered
        window.geneTableData = allAffected;      // Current (possibly filtered)
        window.geneTablePage = 1;
        window.geneTablePerPage = 25;
        window.geneTableSampleFilter = 'all';

        renderGeneTable();
    }}

    // Render paginated gene table
    function renderGeneTable() {{
        const data = window.geneTableData || [];
        const page = window.geneTablePage || 1;
        const perPage = window.geneTablePerPage || 25;
        const totalPages = Math.ceil(data.length / perPage);
        const start = (page - 1) * perPage;
        const end = start + perPage;
        const pageData = data.slice(start, end);

        // Get unique samples for filter
        const allSamples = [...new Set((window.geneTableDataAll || []).map(af => af.sample))].sort();
        const currentFilter = window.geneTableSampleFilter || 'all';

        const geneSearch = window.geneTableGeneSearch || '';

        let html = `
            <div style="display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 15px; margin-bottom: 15px;">
                <div style="display: flex; flex-wrap: wrap; gap: 15px; align-items: center;">
                    <div>
                        <label>Sample: </label>
                        <select id="gene-sample-filter" onchange="filterGeneBySample(this.value)" style="padding: 5px; border-radius: 4px;">
                            <option value="all" ${{currentFilter === 'all' ? 'selected' : ''}}>All samples</option>
                            ${{allSamples.map(s => `<option value="${{s}}" ${{currentFilter === s ? 'selected' : ''}}>${{s}}</option>`).join('')}}
                        </select>
                    </div>
                    <div>
                        <label>Gene: </label>
                        <input type="text" id="gene-search" placeholder="Search gene/product..."
                               value="${{geneSearch}}"
                               onkeyup="searchGene(this.value)"
                               style="padding: 5px; border-radius: 4px; border: 1px solid #ccc; width: 180px;">
                        ${{geneSearch ? `<button onclick="clearGeneSearch()" style="padding: 5px 8px; margin-left: 5px; cursor: pointer;">✕</button>` : ''}}
                    </div>
                    <div>
                        <label>Per page: </label>
                        <select id="gene-per-page" onchange="changeGenePerPage(this.value)" style="padding: 5px; border-radius: 4px;">
                            <option value="25" ${{perPage === 25 ? 'selected' : ''}}>25</option>
                            <option value="50" ${{perPage === 50 ? 'selected' : ''}}>50</option>
                            <option value="100" ${{perPage === 100 ? 'selected' : ''}}>100</option>
                            <option value="all" ${{perPage >= data.length ? 'selected' : ''}}>All (${{data.length}})</option>
                        </select>
                    </div>
                </div>
                <div style="color: #666;">
                    Showing ${{data.length > 0 ? start + 1 : 0}}-${{Math.min(end, data.length)}} of ${{data.length}} genes
                    ${{(currentFilter !== 'all' || geneSearch) ? `(filtered from ${{(window.geneTableDataAll || []).length}})` : ''}}
                </div>
            </div>
            <div style="overflow-x: auto;">
                <table style="width: 100%; border-collapse: collapse;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">#</th>
                            <th onclick="sortGeneTable('gene')" style="padding: 10px; text-align: left; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                Gene ${{getSortIndicator('gene')}}
                            </th>
                            <th onclick="sortGeneTable('product')" style="padding: 10px; text-align: left; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                Product ${{getSortIndicator('product')}}
                            </th>
                            <th onclick="sortGeneTable('position')" style="padding: 10px; text-align: right; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                Position ${{getSortIndicator('position')}}
                            </th>
                            <th onclick="sortGeneTable('overlap')" style="padding: 10px; text-align: right; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                Overlap % ${{getSortIndicator('overlap')}}
                            </th>
                            <th onclick="sortGeneTable('sample')" style="padding: 10px; text-align: left; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                Sample ${{getSortIndicator('sample')}}
                            </th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        pageData.forEach((af, idx) => {{
            const gene = af.feature.gene_name || af.feature.locus_tag || af.feature.id || 'Unknown';
            const product = af.feature.product || '-';
            const pos = `${{af.feature.start.toLocaleString()}}-${{af.feature.end.toLocaleString()}}`;
            const pctClass = af.overlap_pct >= 100 ? 'color: #d32f2f; font-weight: bold;' :
                            af.overlap_pct >= 50 ? 'color: #f57c00;' : '';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; color: #999;">${{start + idx + 1}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-family: monospace;">${{gene}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-size: 0.85rem; max-width: 300px; overflow: hidden; text-overflow: ellipsis;" title="${{product}}">${{product.substring(0, 50)}}${{product.length > 50 ? '...' : ''}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; font-family: monospace;">${{pos}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; ${{pctClass}}">${{af.overlap_pct.toFixed(1)}}%</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6;">${{af.sample}}</td>
                </tr>
            `;
        }});

        html += `
                    </tbody>
                </table>
            </div>
        `;

        // Pagination controls
        if (totalPages > 1) {{
            html += `
                <div style="display: flex; justify-content: center; align-items: center; gap: 10px; margin-top: 15px;">
                    <button onclick="geneTableGoTo(1)" ${{page === 1 ? 'disabled' : ''}}
                            style="padding: 8px 12px; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer; ${{page === 1 ? 'opacity: 0.5;' : ''}}">
                        &laquo; First
                    </button>
                    <button onclick="geneTableGoTo(${{page - 1}})" ${{page === 1 ? 'disabled' : ''}}
                            style="padding: 8px 12px; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer; ${{page === 1 ? 'opacity: 0.5;' : ''}}">
                        &lsaquo; Prev
                    </button>
                    <span style="padding: 8px 15px; background: #f1f3f4; border-radius: 4px;">
                        Page ${{page}} of ${{totalPages}}
                    </span>
                    <button onclick="geneTableGoTo(${{page + 1}})" ${{page === totalPages ? 'disabled' : ''}}
                            style="padding: 8px 12px; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer; ${{page === totalPages ? 'opacity: 0.5;' : ''}}">
                        Next &rsaquo;
                    </button>
                    <button onclick="geneTableGoTo(${{totalPages}})" ${{page === totalPages ? 'disabled' : ''}}
                            style="padding: 8px 12px; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer; ${{page === totalPages ? 'opacity: 0.5;' : ''}}">
                        Last &raquo;
                    </button>
                </div>
            `;
        }}

        document.getElementById('gene-table-container').innerHTML = html;
    }}

    function geneTableGoTo(page) {{
        const totalPages = Math.ceil((window.geneTableData || []).length / window.geneTablePerPage);
        if (page < 1) page = 1;
        if (page > totalPages) page = totalPages;
        window.geneTablePage = page;
        renderGeneTable();
    }}

    function changeGenePerPage(value) {{
        if (value === 'all') {{
            window.geneTablePerPage = (window.geneTableData || []).length;
        }} else {{
            window.geneTablePerPage = parseInt(value);
        }}
        window.geneTablePage = 1;
        renderGeneTable();
    }}

    function filterGeneBySample(sample) {{
        window.geneTableSampleFilter = sample;
        applyGeneFilters();
    }}

    function searchGene(query) {{
        window.geneTableGeneSearch = query;
        applyGeneFilters();
    }}

    function clearGeneSearch() {{
        window.geneTableGeneSearch = '';
        document.getElementById('gene-search').value = '';
        applyGeneFilters();
    }}

    function applyGeneFilters() {{
        const allData = window.geneTableDataAll || [];
        const sampleFilter = window.geneTableSampleFilter || 'all';
        const geneSearch = (window.geneTableGeneSearch || '').toLowerCase().trim();

        let filtered = allData;

        // Filter by sample
        if (sampleFilter !== 'all') {{
            filtered = filtered.filter(af => af.sample === sampleFilter);
        }}

        // Filter by gene search
        if (geneSearch) {{
            filtered = filtered.filter(af => {{
                const gene = (af.feature.gene_name || af.feature.locus_tag || af.feature.id || '').toLowerCase();
                const product = (af.feature.product || '').toLowerCase();
                return gene.includes(geneSearch) || product.includes(geneSearch);
            }});
        }}

        window.geneTableData = filtered;

        // Re-apply current sort
        sortGeneTableData();
        window.geneTablePage = 1;
        renderGeneTable();
    }}

    function sortGeneTableData() {{
        const col = window.geneTableSortCol;
        const data = window.geneTableData || [];

        data.sort((a, b) => {{
            let valA, valB;

            switch(col) {{
                case 'gene':
                    valA = (a.feature.gene_name || a.feature.locus_tag || a.feature.id || '').toLowerCase();
                    valB = (b.feature.gene_name || b.feature.locus_tag || b.feature.id || '').toLowerCase();
                    break;
                case 'product':
                    valA = (a.feature.product || '').toLowerCase();
                    valB = (b.feature.product || '').toLowerCase();
                    break;
                case 'position':
                    valA = a.feature.start;
                    valB = b.feature.start;
                    break;
                case 'overlap':
                    valA = a.overlap_pct;
                    valB = b.overlap_pct;
                    break;
                case 'sample':
                    valA = a.sample.toLowerCase();
                    valB = b.sample.toLowerCase();
                    break;
                default:
                    return 0;
            }}

            let cmp = 0;
            if (typeof valA === 'string') {{
                cmp = valA.localeCompare(valB);
            }} else {{
                cmp = valA - valB;
            }}

            return window.geneTableSortDir === 'asc' ? cmp : -cmp;
        }});
    }}

    // Sort state
    window.geneTableSortCol = 'overlap';
    window.geneTableSortDir = 'desc';

    function getSortIndicator(col) {{
        if (window.geneTableSortCol !== col) return '';
        return window.geneTableSortDir === 'asc' ? '&#9650;' : '&#9660;';
    }}

    function sortGeneTable(col) {{
        // Toggle direction if same column, otherwise default to desc for numbers, asc for strings
        if (window.geneTableSortCol === col) {{
            window.geneTableSortDir = window.geneTableSortDir === 'asc' ? 'desc' : 'asc';
        }} else {{
            window.geneTableSortCol = col;
            window.geneTableSortDir = (col === 'overlap' || col === 'position') ? 'desc' : 'asc';
        }}

        sortGeneTableData();
        window.geneTablePage = 1;
        renderGeneTable();
    }}

    // Populate Empty Column Analysis
    function populateEmptyColumn() {{
        const container = document.getElementById('empty-column-content');
        const analysis = data.empty_column_analysis;

        if (!analysis) {{
            container.innerHTML = '<p style="color: #999;">No coverage threshold analysis available.</p>';
            return;
        }}

        const totalSamples = analysis.total_samples;
        const refLen = analysis.reference_length;

        // Threshold comparison table
        let html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 25px; margin-bottom: 25px;">
                <!-- Threshold Table -->
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="margin-bottom: 15px;">Core Size at Different Thresholds</h4>
                    <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                        <thead>
                            <tr style="background: #f1f3f4;">
                                <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Threshold</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Min Samples</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Core Size (bp)</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Core %</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Empty Regions</th>
                            </tr>
                        </thead>
                        <tbody>
        `;

        for (const t of analysis.threshold_results) {{
            const coreMb = (t.core_positions / 1e6).toFixed(3);
            const coreColor = t.core_pct > 90 ? '#28a745' : t.core_pct > 80 ? '#17a2b8' : t.core_pct > 70 ? '#ffc107' : '#dc3545';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{t.threshold_pct}}%</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{t.min_samples}} / ${{totalSamples}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{t.core_positions.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{coreColor}}; font-weight: bold;">${{t.core_pct.toFixed(1)}}%</span>
                    </td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{t.empty_regions.length}}</td>
                </tr>
            `;
        }}

        html += `
                        </tbody>
                    </table>
                    <p style="margin-top: 10px; font-size: 0.8rem; color: #666;">
                        Threshold = minimum % of samples that must cover a position for it to be "core"
                    </p>
                </div>

                <!-- Coverage Histogram -->
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="margin-bottom: 15px;">Position Coverage Distribution</h4>
                    <div style="display: flex; flex-direction: column; gap: 5px;">
        `;

        // Coverage histogram as bars
        const maxCount = Math.max(...analysis.coverage_histogram.map(h => h.position_count));
        for (let i = totalSamples; i >= 0; i--) {{
            const h = analysis.coverage_histogram[i];
            if (!h) continue;
            const barWidth = maxCount > 0 ? (h.position_count / maxCount * 100) : 0;
            const barColor = i === totalSamples ? '#28a745' : i >= totalSamples * 0.8 ? '#17a2b8' : i >= totalSamples * 0.5 ? '#ffc107' : '#dc3545';
            html += `
                <div style="display: flex; align-items: center; gap: 10px;">
                    <div style="width: 80px; text-align: right; font-size: 0.8rem;">${{i}} samples:</div>
                    <div style="flex: 1; background: #f1f3f4; height: 20px; border-radius: 2px; overflow: hidden;">
                        <div style="width: ${{barWidth}}%; height: 100%; background: ${{barColor}};"></div>
                    </div>
                    <div style="width: 100px; text-align: right; font-size: 0.8rem;">
                        ${{h.position_count.toLocaleString()}} (${{h.position_pct.toFixed(1)}}%)
                    </div>
                </div>
            `;
        }}

        html += `
                    </div>
                    <p style="margin-top: 10px; font-size: 0.8rem; color: #666;">
                        Shows how many positions are covered by exactly N samples
                    </p>
                </div>
            </div>
        `;

        // Fragile regions section with sorting and pagination
        if (analysis.fragile_regions && analysis.fragile_regions.length > 0) {{
            // Store data globally for sorting/pagination
            window.fragileRegions = analysis.fragile_regions.map(fr => ({{
                start: fr.region.start,
                end: fr.region.end,
                length: fr.region.end - fr.region.start,
                max_coverage: fr.max_coverage,
                max_coverage_pct: fr.max_coverage_pct,
                missing_samples: fr.missing_samples
            }}));
            window.fragilePageSize = 50;
            window.fragilePage = 0;
            window.fragileSortCol = 'start';
            window.fragileSortAsc = true;
            window.totalSamples = totalSamples;

            html += `
                <div style="background: #fff8e1; padding: 20px; border-radius: 8px; border-left: 4px solid #ffc107; margin-bottom: 20px;">
                    <h4 style="margin-bottom: 15px; color: #856404;">
                        Fragile Regions (covered by &lt;50% of samples)
                        <span style="font-weight: normal; font-size: 0.9rem; margin-left: 10px;">
                            Total: ${{analysis.fragile_regions.length}} regions
                        </span>
                    </h4>
                    <div style="margin-bottom: 10px; display: flex; gap: 15px; align-items: center;">
                        <label>
                            Rows per page:
                            <select id="fragile-page-size" onchange="fragileChangePageSize(this.value)" style="padding: 4px 8px;">
                                <option value="25">25</option>
                                <option value="50" selected>50</option>
                                <option value="100">100</option>
                                <option value="all">All</option>
                            </select>
                        </label>
                        <span id="fragile-page-info" style="color: #666;"></span>
                        <div style="margin-left: auto;">
                            <button onclick="fragileGoPage(-1)" style="padding: 4px 10px; cursor: pointer;">← Prev</button>
                            <button onclick="fragileGoPage(1)" style="padding: 4px 10px; cursor: pointer;">Next →</button>
                        </div>
                    </div>
                    <div style="max-height: 600px; overflow-y: auto;">
                        <table id="fragile-table" style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                            <thead style="position: sticky; top: 0; background: #f1f3f4;">
                                <tr>
                                    <th onclick="fragileSort('start')" style="padding: 8px; text-align: left; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                        Region <span id="sort-start">▲</span>
                                    </th>
                                    <th onclick="fragileSort('length')" style="padding: 8px; text-align: right; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                        Length <span id="sort-length"></span>
                                    </th>
                                    <th onclick="fragileSort('max_coverage')" style="padding: 8px; text-align: right; border: 1px solid #dee2e6; cursor: pointer; user-select: none;">
                                        Max Coverage <span id="sort-max_coverage"></span>
                                    </th>
                                    <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">
                                        Missing Samples
                                    </th>
                                </tr>
                            </thead>
                            <tbody id="fragile-tbody"></tbody>
                        </table>
                    </div>
                </div>
            `;
        }}

        container.innerHTML = html;

        // Initialize fragile regions table if present
        if (window.fragileRegions && window.fragileRegions.length > 0) {{
            fragileRenderTable();
        }}
    }}

    // Fragile regions sorting function
    function fragileSort(col) {{
        if (window.fragileSortCol === col) {{
            window.fragileSortAsc = !window.fragileSortAsc;
        }} else {{
            window.fragileSortCol = col;
            window.fragileSortAsc = true;
        }}
        window.fragilePage = 0;
        fragileRenderTable();
    }}

    // Fragile regions page navigation
    function fragileGoPage(delta) {{
        const maxPage = window.fragilePageSize === 'all' ? 0 :
            Math.floor((window.fragileRegions.length - 1) / window.fragilePageSize);
        window.fragilePage = Math.max(0, Math.min(maxPage, window.fragilePage + delta));
        fragileRenderTable();
    }}

    // Fragile regions change page size
    function fragileChangePageSize(val) {{
        window.fragilePageSize = val === 'all' ? 'all' : parseInt(val);
        window.fragilePage = 0;
        fragileRenderTable();
    }}

    // Fragile regions render table
    function fragileRenderTable() {{
        const tbody = document.getElementById('fragile-tbody');
        const pageInfo = document.getElementById('fragile-page-info');
        if (!tbody) return;

        // Sort data
        const sorted = [...window.fragileRegions].sort((a, b) => {{
            const col = window.fragileSortCol;
            const aVal = a[col];
            const bVal = b[col];
            const cmp = typeof aVal === 'number' ? aVal - bVal : String(aVal).localeCompare(String(bVal));
            return window.fragileSortAsc ? cmp : -cmp;
        }});

        // Paginate
        const pageSize = window.fragilePageSize === 'all' ? sorted.length : window.fragilePageSize;
        const start = window.fragilePage * pageSize;
        const end = Math.min(start + pageSize, sorted.length);
        const pageData = sorted.slice(start, end);

        // Update sort indicators
        ['start', 'length', 'max_coverage'].forEach(col => {{
            const el = document.getElementById('sort-' + col);
            if (el) {{
                if (col === window.fragileSortCol) {{
                    el.textContent = window.fragileSortAsc ? '▲' : '▼';
                }} else {{
                    el.textContent = '';
                }}
            }}
        }});

        // Update page info
        if (pageInfo) {{
            pageInfo.textContent = `Showing ${{start + 1}}-${{end}} of ${{sorted.length}}`;
        }}

        // Render rows
        let html = '';
        for (const fr of pageData) {{
            html += `
                <tr>
                    <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">
                        ${{fr.start.toLocaleString()}} - ${{fr.end.toLocaleString()}}
                    </td>
                    <td style="padding: 8px; border: 1px solid #dee2e6; text-align: right;">
                        ${{fr.length.toLocaleString()}} bp
                    </td>
                    <td style="padding: 8px; border: 1px solid #dee2e6; text-align: right;">
                        ${{fr.max_coverage}} / ${{window.totalSamples}} (${{fr.max_coverage_pct.toFixed(0)}}%)
                    </td>
                    <td style="padding: 8px; border: 1px solid #dee2e6;">
                        ${{fr.missing_samples.join(', ')}}
                    </td>
                </tr>
            `;
        }}
        tbody.innerHTML = html;
    }}

    // Populate SNP Pipeline Comparison
    function populateSnpCompare() {{
        const container = document.getElementById('snp-compare-content');
        const comparison = data.snp_comparison;
        const scComparison = data.snippycore_comparison;
        const cfsanComparison = data.cfsan_comparison;

        let html = '';

        // Show snippy-core comparison if available
        if (scComparison) {{
            html += `<div id="snippycore-section"></div>`;
        }}

        // Show CFSAN comparison if available
        if (cfsanComparison) {{
            if (scComparison) {{
                html += `<hr style="margin: 30px 0; border: none; border-top: 2px solid #dee2e6;">`;
            }}
            html += `<div id="cfsan-section"></div>`;
        }}

        // Show legacy Snippy comparison if no other comparison available
        if (!scComparison && !cfsanComparison && comparison) {{
            html += `<div id="snippy-legacy-section"></div>`;
        }}

        // Show placeholder if nothing available
        if (!scComparison && !cfsanComparison && !comparison) {{
            container.innerHTML = `
                <div style="background: #f8f9fa; padding: 30px; border-radius: 8px; text-align: center;">
                    <p style="font-size: 1.2rem; color: #666; margin-bottom: 15px;">No SNP pipeline comparison available</p>
                    <p style="color: #999;">Run with <code>--snippycore-dir &lt;PATH&gt;</code> to compare with snippy-core results</p>
                    <p style="color: #999; font-size: 0.85rem; margin-top: 10px;">or <code>--cfsan-dir &lt;PATH&gt;</code> for CFSAN comparison</p>
                </div>
            `;
            return;
        }}

        container.innerHTML = html;

        // Populate sections
        if (scComparison) {{
            populateSnippyCoreComparison(document.getElementById('snippycore-section'), scComparison);
        }}
        if (cfsanComparison) {{
            populateCfsanComparison(document.getElementById('cfsan-section'), cfsanComparison);
        }}
        if (!scComparison && !cfsanComparison && comparison) {{
            populateLegacySnippyComparison(document.getElementById('snippy-legacy-section'), comparison);
        }}
    }}

    // Populate VCF Comparison tab
    function populateVcfCompare() {{
        const container = document.getElementById('vcf-compare-content');
        const vcfComp = data.vcf_comparison;

        if (!vcfComp) {{
            container.innerHTML = `
                <div style="background: #f8f9fa; padding: 30px; border-radius: 8px; text-align: center;">
                    <p style="font-size: 1.2rem; color: #666; margin-bottom: 15px;">No VCF comparison available</p>
                    <p style="color: #999;">Run with <code>--compare-variants --vcf-snippy &lt;DIR&gt; --vcf-cfsan &lt;DIR&gt;</code></p>
                </div>
            `;
            return;
        }}

        const summary = vcfComp.summary;
        const concordanceColor = summary.concordance_rate > 90 ? '#4caf50' :
                                summary.concordance_rate > 70 ? '#ff9800' : '#f44336';

        let html = `
            <div style="background: #e8eaf6; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 0.9rem;">
                <strong>What is VCF Comparison?</strong>
                <p style="margin: 8px 0 0 0;">This directly compares variant calls from ${{summary.pipeline_a}} and ${{summary.pipeline_b}} VCF files.
                Unlike the SNP Pipeline tab (which compares coreguard gaps with pipeline filtering), this shows which SNPs each pipeline actually called.</p>
            </div>

            <!-- Summary Cards -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 15px; margin-bottom: 25px;">
                <div style="background: #e8f5e9; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #388e3c; text-transform: uppercase;">Concordance</div>
                    <div style="font-size: 2rem; font-weight: bold; color: ${{concordanceColor}};">${{summary.concordance_rate.toFixed(1)}}%</div>
                </div>
                <div style="background: #e3f2fd; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #1976d2; text-transform: uppercase;">Concordant</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #1565c0;">${{summary.concordant.toLocaleString()}}</div>
                </div>
                <div style="background: #fff3e0; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #e65100; text-transform: uppercase;">${{summary.pipeline_a}} Only</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #ef6c00;">${{summary.only_pipeline_a.toLocaleString()}}</div>
                </div>
                <div style="background: #fce4ec; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #c2185b; text-transform: uppercase;">${{summary.pipeline_b}} Only</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #ad1457;">${{summary.only_pipeline_b.toLocaleString()}}</div>
                </div>
            </div>

            <!-- Interpretation -->
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px; border-left: 4px solid ${{concordanceColor}}; margin-bottom: 25px;">
                <p style="margin: 0;"><strong>Interpretation:</strong> ${{summary.interpretation}}</p>
            </div>

            <!-- Gap correlation -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">Gap Region Correlation</h4>
                <p style="font-size: 0.9rem; color: #666; margin-bottom: 15px;">
                    How many pipeline-unique SNPs fall within coreguard's gap regions?
                </p>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                    <tr>
                        <td style="padding: 10px; border: 1px solid #dee2e6;">${{summary.pipeline_a}}-only SNPs in gaps</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; font-weight: bold;">${{summary.a_only_in_gaps.toLocaleString()}}</td>
                    </tr>
                    <tr>
                        <td style="padding: 10px; border: 1px solid #dee2e6;">${{summary.pipeline_b}}-only SNPs in gaps</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; font-weight: bold;">${{summary.b_only_in_gaps.toLocaleString()}}</td>
                    </tr>
                </table>
            </div>
        `;

        // BAM Validation section (if available)
        if (summary.bam_validation) {{
            const bv = summary.bam_validation;
            const artifactColor = bv.artifacts > 0 ? '#f44336' : '#4caf50';
            const correctionPct = bv.original_snippy_only > 0 ?
                ((bv.original_snippy_only - bv.corrected_snippy_only) / bv.original_snippy_only * 100).toFixed(1) : 0;

            html += `
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px; border-left: 4px solid ${{artifactColor}};">
                    <h4 style="margin-bottom: 15px;">🔬 BAM Validation (Snippy Artifact Detection)</h4>
                    <p style="font-size: 0.85rem; color: #666; margin-bottom: 15px;">
                        Validates ${{summary.pipeline_a}}-only positions against BAM pileup data to detect artifacts caused by complex variant decomposition.
                        When Snippy reports a polymorphism but the BAM consensus shows a different base, it's likely an artifact.
                    </p>

                    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 15px; margin-bottom: 20px;">
                        <div style="background: #f5f5f5; padding: 15px; border-radius: 8px; text-align: center;">
                            <div style="font-size: 0.75rem; color: #666; text-transform: uppercase;">Checked</div>
                            <div style="font-size: 1.5rem; font-weight: bold; color: #333;">${{bv.total_checked.toLocaleString()}}</div>
                        </div>
                        <div style="background: #e8f5e9; padding: 15px; border-radius: 8px; text-align: center;">
                            <div style="font-size: 0.75rem; color: #388e3c; text-transform: uppercase;">Confirmed Real</div>
                            <div style="font-size: 1.5rem; font-weight: bold; color: #4caf50;">${{bv.real_variants.toLocaleString()}}</div>
                        </div>
                        <div style="background: #ffebee; padding: 15px; border-radius: 8px; text-align: center;">
                            <div style="font-size: 0.75rem; color: #c62828; text-transform: uppercase;">Artifacts</div>
                            <div style="font-size: 1.5rem; font-weight: bold; color: ${{artifactColor}};">${{bv.artifacts.toLocaleString()}}</div>
                        </div>
                        <div style="background: #fff8e1; padding: 15px; border-radius: 8px; text-align: center;">
                            <div style="font-size: 0.75rem; color: #f57c00; text-transform: uppercase;">No BAM Data</div>
                            <div style="font-size: 1.5rem; font-weight: bold; color: #ff9800;">${{bv.no_data.toLocaleString()}}</div>
                        </div>
                    </div>

                    <div style="background: #f5f5f5; padding: 15px; border-radius: 8px;">
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                            <tr>
                                <td style="padding: 8px 0;">Artifact Rate</td>
                                <td style="padding: 8px 0; text-align: right; font-weight: bold; color: ${{artifactColor}};">${{bv.artifact_rate.toFixed(1)}}%</td>
                            </tr>
                            <tr style="border-top: 1px solid #ddd;">
                                <td style="padding: 8px 0;">Original ${{summary.pipeline_a}}-only count</td>
                                <td style="padding: 8px 0; text-align: right; font-weight: bold; color: #ef6c00;">${{bv.original_snippy_only.toLocaleString()}}</td>
                            </tr>
                            <tr>
                                <td style="padding: 8px 0;">Corrected ${{summary.pipeline_a}}-only count</td>
                                <td style="padding: 8px 0; text-align: right; font-weight: bold; color: #4caf50;">${{bv.corrected_snippy_only.toLocaleString()}}</td>
                            </tr>
                            <tr style="border-top: 1px solid #ddd; background: #e8f5e9;">
                                <td style="padding: 8px 0;"><strong>Reduction</strong></td>
                                <td style="padding: 8px 0; text-align: right; font-weight: bold; color: #388e3c;">-${{correctionPct}}%</td>
                            </tr>
                        </table>
                    </div>
                </div>
            `;
        }}

        // Per-sample breakdown
        if (vcfComp.sample_comparisons && vcfComp.sample_comparisons.length > 0) {{
            html += `
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                    <h4 style="margin-bottom: 15px;">Per-Sample Comparison</h4>
                    <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                        <thead>
                            <tr style="background: #f1f3f4;">
                                <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{summary.pipeline_a}} SNPs</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{summary.pipeline_b}} SNPs</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Concordant</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{summary.pipeline_a}}-only</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{summary.pipeline_b}}-only</th>
                                <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Concordance</th>
                            </tr>
                        </thead>
                        <tbody>
            `;

            for (const s of vcfComp.sample_comparisons) {{
                const concColor = s.concordance_rate > 90 ? '#4caf50' : s.concordance_rate > 70 ? '#ff9800' : '#f44336';
                html += `
                    <tr>
                        <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{s.sample}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.snps_a.toLocaleString()}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.snps_b.toLocaleString()}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #4caf50;">${{s.concordant.toLocaleString()}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #ef6c00;">${{s.only_a.toLocaleString()}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #ad1457;">${{s.only_b.toLocaleString()}}</td>
                        <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                            <span style="color: ${{concColor}}; font-weight: bold;">${{s.concordance_rate.toFixed(1)}}%</span>
                        </td>
                    </tr>
                `;
            }}
            html += `</tbody></table></div>`;
        }}

        // Top missed SNPs by pipeline B (high quality in A)
        if (vcfComp.top_missed_by_b && vcfComp.top_missed_by_b.length > 0) {{
            html += `
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                    <h4 style="margin-bottom: 15px; color: #e65100;">Top SNPs Missed by ${{summary.pipeline_b}}</h4>
                    <p style="font-size: 0.85rem; color: #666; margin-bottom: 15px;">
                        These are high-quality variants called by ${{summary.pipeline_a}} (QUAL>100, DP>10) that ${{summary.pipeline_b}} did not call.
                    </p>
                    <table style="width: 100%; border-collapse: collapse; font-size: 0.8rem;">
                        <thead>
                            <tr style="background: #fff3e0;">
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Position</th>
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">Ref→Alt</th>
                                <th style="padding: 8px; text-align: right; border: 1px solid #dee2e6;">QUAL</th>
                                <th style="padding: 8px; text-align: right; border: 1px solid #dee2e6;">DP</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">In Gap?</th>
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Likely Cause</th>
                            </tr>
                        </thead>
                        <tbody>
            `;

            for (const p of vcfComp.top_missed_by_b.slice(0, 20)) {{
                const varA = p.pipeline_a || {{}};
                html += `
                    <tr>
                        <td style="padding: 8px; border: 1px solid #dee2e6; font-family: monospace;">${{p.chrom}}:${{p.pos}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6;">${{p.sample}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6; text-align: center; font-family: monospace;">${{varA.ref_allele || '?'}}→${{varA.alt_allele || '?'}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6; text-align: right;">${{varA.qual ? varA.qual.toFixed(0) : '-'}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6; text-align: right;">${{varA.depth || '-'}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6; text-align: center;">${{p.in_gap_region ? '⚠️ Yes' : '✓ No'}}</td>
                        <td style="padding: 8px; border: 1px solid #dee2e6; font-size: 0.75rem;">${{p.likely_cause || '-'}}</td>
                    </tr>
                `;
            }}
            html += `</tbody></table></div>`;
        }}

        // All Discordant Positions (lazy loading)
        html += `
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">All Discordant Positions</h4>
                <p style="font-size: 0.85rem; color: #666; margin-bottom: 15px;">
                    Browse all variant positions with lazy loading. Use the filters to narrow down the results.
                </p>

                <!-- Filters -->
                <div style="display: flex; gap: 15px; margin-bottom: 15px; flex-wrap: wrap; align-items: center;">
                    <div>
                        <label style="font-size: 0.85rem; color: #666;">Sample:</label>
                        <select id="vcf-sample-filter" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; margin-left: 5px;">
                            <option value="all">All Samples</option>
                        </select>
                    </div>
                    <div>
                        <label style="font-size: 0.85rem; color: #666;">Status:</label>
                        <select id="vcf-status-filter" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; margin-left: 5px;">
                            <option value="all">All Statuses</option>
                            <option value="OnlyPipelineA">${{summary.pipeline_a}} Only</option>
                            <option value="OnlyPipelineB">${{summary.pipeline_b}} Only</option>
                            <option value="DiscordantAllele">Discordant Allele</option>
                        </select>
                    </div>
                    <button id="vcf-apply-filter" style="padding: 5px 15px; background: #1976d2; color: white; border: none; border-radius: 4px; cursor: pointer;">
                        Apply Filter
                    </button>
                    <span id="vcf-total-count" style="font-size: 0.85rem; color: #666;"></span>
                </div>

                <!-- Table -->
                <div style="overflow-x: auto;">
                    <table id="vcf-positions-table" style="width: 100%; border-collapse: collapse; font-size: 0.8rem;">
                        <thead>
                            <tr style="background: #f5f5f5;">
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Position</th>
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">Status</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">${{summary.pipeline_a}}</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">${{summary.pipeline_b}}</th>
                                <th style="padding: 8px; text-align: center; border: 1px solid #dee2e6;">In Gap?</th>
                                <th style="padding: 8px; text-align: left; border: 1px solid #dee2e6;">Likely Cause</th>
                            </tr>
                        </thead>
                        <tbody id="vcf-positions-tbody">
                            <tr><td colspan="7" style="padding: 20px; text-align: center; color: #666;">Loading...</td></tr>
                        </tbody>
                    </table>
                </div>

                <!-- Pagination -->
                <div id="vcf-pagination" style="display: flex; gap: 10px; justify-content: center; align-items: center; margin-top: 15px;">
                </div>
            </div>
        `;

        container.innerHTML = html;

        // Initialize lazy loading
        initVcfLazyLoading(summary);
    }}

    // Initialize VCF positions lazy loading
    function initVcfLazyLoading(summary) {{
        let currentPage = 0;
        const pageSize = 50;

        async function loadPositions(page, sample, status) {{
            const tbody = document.getElementById('vcf-positions-tbody');
            const pagination = document.getElementById('vcf-pagination');
            const totalCount = document.getElementById('vcf-total-count');

            tbody.innerHTML = '<tr><td colspan="7" style="padding: 20px; text-align: center; color: #666;">Loading...</td></tr>';

            try {{
                const params = new URLSearchParams({{
                    page: page,
                    page_size: pageSize,
                    sample: sample || 'all',
                    status: status || 'all'
                }});

                const response = await fetch(`/api/vcf_positions?${{params}}`);
                const data = await response.json();

                if (data.error) {{
                    tbody.innerHTML = `<tr><td colspan="7" style="padding: 20px; text-align: center; color: #f44336;">API not available (run with --compare-variants)</td></tr>`;
                    return;
                }}

                // Update sample filter options
                const sampleFilter = document.getElementById('vcf-sample-filter');
                if (data.available_samples && sampleFilter.options.length <= 1) {{
                    for (const s of data.available_samples) {{
                        const opt = document.createElement('option');
                        opt.value = s;
                        opt.textContent = s;
                        sampleFilter.appendChild(opt);
                    }}
                }}

                // Update total count
                totalCount.textContent = `${{data.total.toLocaleString()}} positions`;

                // Render positions
                if (data.positions.length === 0) {{
                    tbody.innerHTML = '<tr><td colspan="7" style="padding: 20px; text-align: center; color: #666;">No positions match the filter</td></tr>';
                }} else {{
                    let rows = '';
                    for (const p of data.positions) {{
                        const varA = p.pipeline_a || {{}};
                        const varB = p.pipeline_b || {{}};
                        const statusColor = p.status === 'OnlyPipelineA' ? '#ef6c00' :
                                          p.status === 'OnlyPipelineB' ? '#ad1457' :
                                          p.status === 'DiscordantAllele' ? '#7b1fa2' : '#4caf50';
                        rows += `
                            <tr>
                                <td style="padding: 6px; border: 1px solid #dee2e6; font-family: monospace;">${{p.chrom}}:${{p.pos}}</td>
                                <td style="padding: 6px; border: 1px solid #dee2e6;">${{p.sample}}</td>
                                <td style="padding: 6px; border: 1px solid #dee2e6; text-align: center;">
                                    <span style="background: ${{statusColor}}; color: white; padding: 2px 6px; border-radius: 3px; font-size: 0.75rem;">
                                        ${{p.status.replace('Pipeline', '')}}
                                    </span>
                                </td>
                                <td style="padding: 6px; border: 1px solid #dee2e6; text-align: center; font-family: monospace;">
                                    ${{varA.ref_allele ? varA.ref_allele + '→' + varA.alt_allele : '-'}}
                                </td>
                                <td style="padding: 6px; border: 1px solid #dee2e6; text-align: center; font-family: monospace;">
                                    ${{varB.ref_allele ? varB.ref_allele + '→' + varB.alt_allele : '-'}}
                                </td>
                                <td style="padding: 6px; border: 1px solid #dee2e6; text-align: center;">
                                    ${{p.in_gap_region ? '⚠️' : '✓'}}
                                </td>
                                <td style="padding: 6px; border: 1px solid #dee2e6; font-size: 0.75rem;">${{p.likely_cause || '-'}}</td>
                            </tr>
                        `;
                    }}
                    tbody.innerHTML = rows;
                }}

                // Render pagination
                let paginationHtml = '';
                if (data.total_pages > 1) {{
                    paginationHtml += `<button onclick="vcfLoadPage(${{Math.max(0, page - 1)}})" ${{page === 0 ? 'disabled' : ''}}
                        style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">◀ Prev</button>`;
                    paginationHtml += `<span style="font-size: 0.9rem;">Page ${{page + 1}} of ${{data.total_pages}}</span>`;
                    paginationHtml += `<button onclick="vcfLoadPage(${{Math.min(data.total_pages - 1, page + 1)}})" ${{page >= data.total_pages - 1 ? 'disabled' : ''}}
                        style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">Next ▶</button>`;
                }}
                pagination.innerHTML = paginationHtml;
                currentPage = page;

            }} catch (e) {{
                tbody.innerHTML = `<tr><td colspan="7" style="padding: 20px; text-align: center; color: #f44336;">Error loading data: ${{e.message}}</td></tr>`;
            }}
        }}

        // Global function for pagination buttons
        window.vcfLoadPage = (page) => {{
            const sample = document.getElementById('vcf-sample-filter').value;
            const status = document.getElementById('vcf-status-filter').value;
            loadPositions(page, sample, status);
        }};

        // Apply filter button
        document.getElementById('vcf-apply-filter').addEventListener('click', () => {{
            window.vcfLoadPage(0);
        }});

        // Initial load
        loadPositions(0, 'all', 'all');
    }}

    // Populate Snippy Bugfix tab
    function populateSnippyBugfix() {{
        const container = document.getElementById('snippy-bugfix-content');
        const vcfComp = data.vcf_comparison;

        if (!vcfComp || !vcfComp.summary || !vcfComp.summary.bam_validation) {{
            container.innerHTML = `
                <div style="background: #f8f9fa; padding: 30px; border-radius: 8px; text-align: center;">
                    <p style="font-size: 1.2rem; color: #666; margin-bottom: 15px;">No BAM validation available</p>
                    <p style="color: #999;">BAM validation requires:</p>
                    <ul style="text-align: left; max-width: 500px; margin: 10px auto; color: #666;">
                        <li>Run with <code>--compare-variants --vcf-snippy &lt;DIR&gt; --vcf-cfsan &lt;DIR&gt;</code></li>
                        <li>BAM files must exist alongside Snippy VCF files (in *_snippy/ directories)</li>
                        <li><code>samtools</code> must be available in PATH</li>
                    </ul>
                </div>
            `;
            return;
        }}

        const summary = vcfComp.summary;
        const bv = summary.bam_validation;
        const artifactColor = bv.artifacts > 0 ? '#f44336' : '#4caf50';
        const correctionPct = bv.original_snippy_only > 0 ?
            ((bv.original_snippy_only - bv.corrected_snippy_only) / bv.original_snippy_only * 100).toFixed(1) : 0;

        // Calculate corrected metrics
        const originalSnippyOnly = bv.original_snippy_only;
        const correctedSnippyOnly = bv.corrected_snippy_only;
        const artifactsRemoved = originalSnippyOnly - correctedSnippyOnly;

        // Recalculate concordance
        const originalTotal = summary.concordant + summary.only_pipeline_a + summary.only_pipeline_b;
        const correctedTotal = summary.concordant + correctedSnippyOnly + summary.only_pipeline_b;
        const originalConcordance = originalTotal > 0 ? (summary.concordant / originalTotal * 100) : 0;
        const correctedConcordance = correctedTotal > 0 ? (summary.concordant / correctedTotal * 100) : 0;

        let html = `
            <!-- Explanation -->
            <div style="background: #e8eaf6; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 0.9rem;">
                <strong>What is Snippy Artifact Detection?</strong>
                <p style="margin: 8px 0 0 0;">
                    Snippy uses FreeBayes which can generate complex variants (MNPs) that are then decomposed into individual SNPs.
                    This decomposition can create <strong>artifact polymorphisms</strong> - positions reported as variants that don't actually
                    differ from the reference when examining the raw BAM reads.
                </p>
                <p style="margin: 8px 0 0 0;">
                    This analysis validates each Snippy-only position against the BAM pileup. If the consensus base from the BAM
                    differs from what Snippy reported, it's flagged as an artifact.
                </p>
            </div>

            <!-- BEFORE/AFTER Comparison Hero -->
            <div style="display: grid; grid-template-columns: 1fr auto 1fr; gap: 20px; margin-bottom: 30px; align-items: stretch;">
                <!-- BEFORE -->
                <div style="background: linear-gradient(135deg, #ffebee 0%, #ffcdd2 100%); padding: 25px; border-radius: 12px; border: 2px solid #ef5350;">
                    <div style="text-align: center; margin-bottom: 20px;">
                        <span style="background: #f44336; color: white; padding: 5px 15px; border-radius: 20px; font-size: 0.85rem; font-weight: bold;">ORIGINAL (Uncorrected)</span>
                    </div>
                    <div style="text-align: center;">
                        <div style="font-size: 0.9rem; color: #c62828; margin-bottom: 5px;">${{summary.pipeline_a}}-only SNPs</div>
                        <div style="font-size: 3rem; font-weight: bold; color: #d32f2f;">${{originalSnippyOnly.toLocaleString()}}</div>
                    </div>
                    <div style="margin-top: 20px; padding-top: 15px; border-top: 1px solid #ef9a9a;">
                        <table style="width: 100%; font-size: 0.85rem; color: #c62828;">
                            <tr>
                                <td>Concordant</td>
                                <td style="text-align: right; font-weight: bold;">${{summary.concordant.toLocaleString()}}</td>
                            </tr>
                            <tr>
                                <td>${{summary.pipeline_b}}-only</td>
                                <td style="text-align: right; font-weight: bold;">${{summary.only_pipeline_b.toLocaleString()}}</td>
                            </tr>
                            <tr style="border-top: 1px solid #ef9a9a;">
                                <td><strong>Concordance Rate</strong></td>
                                <td style="text-align: right; font-weight: bold;">${{originalConcordance.toFixed(1)}}%</td>
                            </tr>
                        </table>
                    </div>
                </div>

                <!-- ARROW -->
                <div style="display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 20px;">
                    <div style="font-size: 3rem; color: #4caf50;">→</div>
                    <div style="background: #4caf50; color: white; padding: 8px 15px; border-radius: 20px; font-size: 0.85rem; font-weight: bold; margin-top: 10px;">
                        -${{artifactsRemoved.toLocaleString()}} artifacts
                    </div>
                    <div style="color: #666; font-size: 0.8rem; margin-top: 5px;">(${{correctionPct}}% reduction)</div>
                </div>

                <!-- AFTER -->
                <div style="background: linear-gradient(135deg, #e8f5e9 0%, #c8e6c9 100%); padding: 25px; border-radius: 12px; border: 2px solid #66bb6a;">
                    <div style="text-align: center; margin-bottom: 20px;">
                        <span style="background: #4caf50; color: white; padding: 5px 15px; border-radius: 20px; font-size: 0.85rem; font-weight: bold;">CORRECTED (Bugfix Applied)</span>
                    </div>
                    <div style="text-align: center;">
                        <div style="font-size: 0.9rem; color: #2e7d32; margin-bottom: 5px;">${{summary.pipeline_a}}-only SNPs</div>
                        <div style="font-size: 3rem; font-weight: bold; color: #388e3c;">${{correctedSnippyOnly.toLocaleString()}}</div>
                    </div>
                    <div style="margin-top: 20px; padding-top: 15px; border-top: 1px solid #a5d6a7;">
                        <table style="width: 100%; font-size: 0.85rem; color: #2e7d32;">
                            <tr>
                                <td>Concordant</td>
                                <td style="text-align: right; font-weight: bold;">${{summary.concordant.toLocaleString()}}</td>
                            </tr>
                            <tr>
                                <td>${{summary.pipeline_b}}-only</td>
                                <td style="text-align: right; font-weight: bold;">${{summary.only_pipeline_b.toLocaleString()}}</td>
                            </tr>
                            <tr style="border-top: 1px solid #a5d6a7;">
                                <td><strong>Concordance Rate</strong></td>
                                <td style="text-align: right; font-weight: bold;">${{correctedConcordance.toFixed(1)}}%</td>
                            </tr>
                        </table>
                    </div>
                </div>
            </div>

            <!-- Validation Stats -->
            <div style="background: white; padding: 25px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 20px;">BAM Validation Details</h4>
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px;">
                    <div style="background: #f5f5f5; padding: 20px; border-radius: 8px; text-align: center;">
                        <div style="font-size: 0.75rem; color: #666; text-transform: uppercase;">Positions Checked</div>
                        <div style="font-size: 2rem; font-weight: bold; color: #333;">${{bv.total_checked.toLocaleString()}}</div>
                    </div>
                    <div style="background: #e8f5e9; padding: 20px; border-radius: 8px; text-align: center;">
                        <div style="font-size: 0.75rem; color: #388e3c; text-transform: uppercase;">Confirmed Real</div>
                        <div style="font-size: 2rem; font-weight: bold; color: #4caf50;">${{bv.real_variants.toLocaleString()}}</div>
                    </div>
                    <div style="background: #ffebee; padding: 20px; border-radius: 8px; text-align: center;">
                        <div style="font-size: 0.75rem; color: #c62828; text-transform: uppercase;">Artifacts Detected</div>
                        <div style="font-size: 2rem; font-weight: bold; color: ${{artifactColor}};">${{bv.artifacts.toLocaleString()}}</div>
                    </div>
                    <div style="background: #fff8e1; padding: 20px; border-radius: 8px; text-align: center;">
                        <div style="font-size: 0.75rem; color: #f57c00; text-transform: uppercase;">No BAM Data</div>
                        <div style="font-size: 2rem; font-weight: bold; color: #ff9800;">${{bv.no_data.toLocaleString()}}</div>
                    </div>
                </div>

                <div style="margin-top: 20px; background: #f5f5f5; padding: 15px; border-radius: 8px;">
                    <div style="display: flex; justify-content: space-between; align-items: center;">
                        <span style="font-size: 1.1rem;">Artifact Rate:</span>
                        <span style="font-size: 1.5rem; font-weight: bold; color: ${{artifactColor}};">${{bv.artifact_rate.toFixed(1)}}%</span>
                    </div>
                </div>
            </div>

            <!-- Interpretation -->
            <div style="background: ${{bv.artifacts > 0 ? '#fff3e0' : '#e8f5e9'}}; padding: 20px; border-radius: 8px; border-left: 4px solid ${{bv.artifacts > 0 ? '#ff9800' : '#4caf50'}};">
                <h4 style="margin-bottom: 10px;">${{bv.artifacts > 0 ? '⚠️ Artifacts Detected' : '✓ No Artifacts Found'}}</h4>
                <p style="margin: 0; line-height: 1.6;">
                    ${{bv.artifacts > 0 ?
                        `BAM validation detected <strong>${{bv.artifacts.toLocaleString()}} artifacts</strong> (${{bv.artifact_rate.toFixed(1)}}% of Snippy-only positions).
                        These are positions where Snippy reported a variant but the BAM reads show a different consensus base.
                        The corrected Snippy-only count is <strong>${{correctedSnippyOnly.toLocaleString()}}</strong>, which improves the
                        concordance rate from ${{originalConcordance.toFixed(1)}}% to <strong>${{correctedConcordance.toFixed(1)}}%</strong>.` :
                        `All ${{bv.real_variants.toLocaleString()}} Snippy-only positions were confirmed as real variants by BAM validation.
                        No artifacts were detected in this dataset.`
                    }}
                </p>
            </div>
        `;

        // Add corrected distance matrix if available
        if (bv.corrected_distance_matrix && bv.corrected_distance_matrix.pairwise && bv.corrected_distance_matrix.pairwise.length > 0) {{
            const dm = bv.corrected_distance_matrix;
            html += `
                <!-- Corrected Distance Matrix -->
                <div style="background: white; padding: 25px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 25px;">
                    <h4 style="margin-bottom: 15px;">📊 Corrected Pairwise Distance Matrix</h4>
                    <p style="color: #666; font-size: 0.9rem; margin-bottom: 20px;">
                        Distances calculated from <strong>${{bv.real_variants}}</strong> confirmed real variant positions only (artifacts excluded).
                        Values represent SNP differences between samples based on BAM consensus.
                    </p>
            `;

            // Create matrix table
            const samples = dm.samples;
            const matrix = dm.matrix;

            if (samples.length <= 10) {{
                // Full matrix view for small datasets
                html += `<div style="overflow-x: auto;"><table style="border-collapse: collapse; font-size: 0.85rem; width: 100%;">`;

                // Header row
                html += `<tr><th style="padding: 10px; background: #f5f5f5; border: 1px solid #ddd;"></th>`;
                for (const s of samples) {{
                    html += `<th style="padding: 10px; background: #f5f5f5; border: 1px solid #ddd; font-weight: bold;">${{s}}</th>`;
                }}
                html += `</tr>`;

                // Data rows
                for (let i = 0; i < samples.length; i++) {{
                    html += `<tr><td style="padding: 10px; background: #f5f5f5; border: 1px solid #ddd; font-weight: bold;">${{samples[i]}}</td>`;
                    for (let j = 0; j < samples.length; j++) {{
                        const val = matrix[i][j];
                        const bgColor = i === j ? '#e0e0e0' : (val === 0 ? '#c8e6c9' : (val < 5 ? '#e8f5e9' : (val < 20 ? '#fff8e1' : '#ffebee')));
                        html += `<td style="padding: 10px; text-align: center; border: 1px solid #ddd; background: ${{bgColor}}; font-weight: ${{val > 0 ? 'bold' : 'normal'}};">${{val}}</td>`;
                    }}
                    html += `</tr>`;
                }}
                html += `</table></div>`;
            }}

            // Pairwise list (always show)
            html += `
                <h5 style="margin-top: 20px; margin-bottom: 10px;">Pairwise Distances</h5>
                <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 10px; max-height: 300px; overflow-y: auto;">
            `;

            // Sort pairwise by distance
            const sortedPairwise = [...dm.pairwise].sort((a, b) => a.distance - b.distance);

            for (const pw of sortedPairwise) {{
                const distColor = pw.distance === 0 ? '#4caf50' : (pw.distance < 5 ? '#8bc34a' : (pw.distance < 20 ? '#ff9800' : '#f44336'));
                html += `
                    <div style="background: #f8f9fa; padding: 12px; border-radius: 6px; display: flex; justify-content: space-between; align-items: center; border-left: 4px solid ${{distColor}};">
                        <div>
                            <span style="font-weight: 500;">${{pw.sample_a}}</span>
                            <span style="color: #999; margin: 0 8px;">↔</span>
                            <span style="font-weight: 500;">${{pw.sample_b}}</span>
                        </div>
                        <div style="background: ${{distColor}}; color: white; padding: 4px 12px; border-radius: 15px; font-weight: bold; font-size: 0.9rem;">
                            ${{pw.distance}} SNP${{pw.distance !== 1 ? 's' : ''}}
                        </div>
                    </div>
                `;
            }}

            html += `</div></div>`;
        }}

        container.innerHTML = html;
    }}

    // SNP Tracks visualization - 1bp resolution with reference sequence
    function populateSnpTracks() {{
        const container = document.getElementById('snp-tracks-container');
        const infoDiv = document.getElementById('snp-tracks-info');

        if (!data.vcf_comparison) {{
            container.innerHTML = `<p style="color: #666; text-align: center;">No VCF comparison data available</p>`;
            return;
        }}

        // State
        let currentStart = 1;
        let windowSize = 100; // Show 100bp at a time
        let regionData = null;
        const refLength = data.reference_length || 3000000;

        // Color map for nucleotides
        const baseColors = {{
            'A': '#e8f5e9', // light green
            'T': '#e3f2fd', // light blue
            'G': '#fff3e0', // light orange
            'C': '#fce4ec', // light pink
            '-': '#9e9e9e'  // gray for gaps
        }};

        async function loadRegion(start, end) {{
            try {{
                const response = await fetch(`/api/region?start=${{start}}&end=${{end}}`);
                regionData = await response.json();
                renderRegion();
            }} catch (e) {{
                console.error('Error loading region:', e);
                container.innerHTML = `<p style="color: #d32f2f;">Error loading region data</p>`;
            }}
        }}

        function renderRegion() {{
            if (!regionData) return;

            const {{ start, end, reference, snps, samples, gaps }} = regionData;

            // Build SNP lookup: position -> sample -> snp data
            const snpMap = new Map();
            snps.forEach(snp => {{
                if (!snpMap.has(snp.pos)) snpMap.set(snp.pos, new Map());
                snpMap.get(snp.pos).set(snp.sample, snp);
            }});

            // Build gap lookup: sample -> position -> gap type
            const gapMap = new Map();
            samples.forEach(s => gapMap.set(s, new Map()));
            (gaps || []).forEach(gap => {{
                const sampleGaps = gapMap.get(gap.sample);
                if (sampleGaps) {{
                    for (let p = gap.start; p <= gap.end; p++) {{
                        sampleGaps.set(p, gap.type || 'mm2');
                    }}
                }}
            }});

            // Gap colors by type
            const gapColors = {{
                'mm2': '#9e9e9e',      // Gray - minimap2 coverage gaps
                'snippy': '#e65100',   // Dark orange - Snippy gaps (N in aligned FASTA)
                'cfsan': '#0d47a1'     // Dark blue - CFSAN gaps
            }};

            let html = `
                <div style="margin-bottom: 15px; display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;">
                    <div style="font-size: 0.85rem; color: #666;">
                        Region: <strong>${{start.toLocaleString()}}</strong> - <strong>${{end.toLocaleString()}}</strong> (${{reference.length}} bp)
                    </div>
                    <div style="display: flex; gap: 5px;">
                        <button onclick="window.snpTracksNav(-1000)" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">⏪ -1kb</button>
                        <button onclick="window.snpTracksNav(-100)" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">◀ -100</button>
                        <button onclick="window.snpTracksNav(100)" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">+100 ▶</button>
                        <button onclick="window.snpTracksNav(1000)" style="padding: 5px 10px; border: 1px solid #ddd; border-radius: 4px; cursor: pointer;">+1kb ⏩</button>
                    </div>
                </div>
            `;

            // Create table - one row per sample + reference row
            html += `<div style="overflow-x: auto;"><table style="border-collapse: collapse; font-family: monospace; font-size: 12px;">`;

            // Header row with positions
            html += `<thead><tr><th style="padding: 4px 8px; background: #f5f5f5; border: 1px solid #ddd; position: sticky; left: 0; z-index: 10;">Pos</th>`;
            for (let i = 0; i < reference.length; i++) {{
                const pos = start + i;
                const showPos = i % 10 === 0;
                html += `<th style="padding: 2px; background: #f5f5f5; border: 1px solid #eee; min-width: 20px; font-size: 9px; writing-mode: vertical-rl; height: ${{showPos ? '50px' : '20px'}};">${{showPos ? pos : ''}}</th>`;
            }}
            html += `</tr></thead><tbody>`;

            // Reference row
            html += `<tr><td style="padding: 4px 8px; background: #e0e0e0; border: 1px solid #ddd; font-weight: bold; position: sticky; left: 0; z-index: 5;">REF</td>`;
            for (let i = 0; i < reference.length; i++) {{
                const base = reference[i].toUpperCase();
                const bgColor = baseColors[base] || '#fff';
                html += `<td style="padding: 2px; text-align: center; border: 1px solid #eee; background: ${{bgColor}}; font-weight: bold;">${{base}}</td>`;
            }}
            html += `</tr>`;

            // Sample rows
            samples.forEach(sample => {{
                html += `<tr><td style="padding: 4px 8px; background: #f5f5f5; border: 1px solid #ddd; font-weight: bold; position: sticky; left: 0; z-index: 5;">${{sample}}</td>`;
                const sampleGaps = gapMap.get(sample) || new Set();

                for (let i = 0; i < reference.length; i++) {{
                    const pos = start + i;
                    const refBase = reference[i].toUpperCase();
                    const snpData = snpMap.get(pos)?.get(sample);
                    const gapType = sampleGaps.get(pos);

                    let cellContent = refBase;
                    let bgColor = baseColors[refBase] || '#fff';
                    let textColor = '#333';
                    let title = `${{sample}} @ ${{pos}}: ${{refBase}}`;
                    let fontWeight = 'normal';

                    // Check for gap first - gaps take priority
                    if (gapType) {{
                        cellContent = '-';
                        bgColor = gapColors[gapType] || '#9e9e9e';
                        textColor = '#fff';
                        title = `${{sample}} @ ${{pos}}: GAP (${{gapType.toUpperCase()}})`;
                        fontWeight = 'bold';
                    }} else if (snpData) {{
                        cellContent = snpData.alt || '?';
                        fontWeight = 'bold';
                        title = `${{sample}} @ ${{pos}}: ${{snpData.ref}}→${{snpData.alt}} (${{snpData.status}})`;

                        if (snpData.status === 'Concordant') {{
                            bgColor = '#4caf50';
                            textColor = '#fff';
                        }} else if (snpData.status === 'OnlyPipelineB') {{
                            bgColor = '#2196f3';
                            textColor = '#fff';
                        }} else if (snpData.status === 'OnlyPipelineA') {{
                            bgColor = '#ff9800';
                            textColor = '#fff';
                        }} else if (snpData.status === 'DiscordantAllele') {{
                            bgColor = '#9c27b0';
                            textColor = '#fff';
                        }}
                    }}

                    // Add onclick for SNP cells to show quality popup
                    const onclick = snpData ? `onclick="window.showQualityPopup(${{pos}}, '${{sample}}')"` : '';
                    html += `<td title="${{title}}" ${{onclick}} style="padding: 2px; text-align: center; border: 1px solid #eee; background: ${{bgColor}}; color: ${{textColor}}; font-weight: ${{fontWeight}}; cursor: ${{snpData ? 'pointer' : 'default'}};">${{cellContent}}</td>`;
                }}

                html += `</tr>`;
            }});

            html += `</tbody></table></div>`;

            container.innerHTML = html;

            // Update info
            const snpCount = snps.length;
            const statusCounts = {{}};
            snps.forEach(s => {{
                statusCounts[s.status] = (statusCounts[s.status] || 0) + 1;
            }});

            // Count gap positions by type
            const gapCounts = {{ mm2: 0, snippy: 0, cfsan: 0 }};
            gapMap.forEach((positions, sample) => {{
                positions.forEach((gapType, pos) => {{
                    gapCounts[gapType] = (gapCounts[gapType] || 0) + 1;
                }});
            }});
            const totalGaps = gapCounts.mm2 + gapCounts.snippy + gapCounts.cfsan;

            infoDiv.innerHTML = `
                <strong>SNPs:</strong> ${{snpCount}}
                &nbsp;|&nbsp;
                <span style="color: #4caf50;">■</span> Concordant: ${{statusCounts['Concordant'] || 0}}
                &nbsp;|&nbsp;
                <span style="color: #ff9800;">■</span> Snippy-only: ${{statusCounts['OnlyPipelineA'] || 0}}
                &nbsp;|&nbsp;
                <span style="color: #2196f3;">■</span> CFSAN-only: ${{statusCounts['OnlyPipelineB'] || 0}}
                &nbsp;|&nbsp;
                <strong>Gaps:</strong>
                <span style="color: #9e9e9e;">■</span> mm2: ${{gapCounts.mm2}}
                ${{gapCounts.snippy > 0 ? `<span style="color: #e65100;">■</span> Snippy: ${{gapCounts.snippy}}` : ''}}
                ${{gapCounts.cfsan > 0 ? `<span style="color: #0d47a1;">■</span> CFSAN: ${{gapCounts.cfsan}}` : ''}}
            `;

            document.getElementById('snp-region-input').value = `${{start}}-${{end}}`;
        }}

        // Navigation function
        window.snpTracksNav = (delta) => {{
            currentStart = Math.max(1, Math.min(refLength - windowSize, currentStart + delta));
            loadRegion(currentStart, currentStart + windowSize);
        }};

        // Event handlers
        document.getElementById('snp-tracks-update').onclick = () => {{
            const input = document.getElementById('snp-region-input').value;
            const match = input.match(/(\d+)\s*-?\s*(\d*)/);
            if (match) {{
                const start = parseInt(match[1]) || 1;
                const end = parseInt(match[2]) || (start + windowSize);
                windowSize = Math.min(500, end - start);
                currentStart = start;
                loadRegion(start, end);
            }}
        }};

        document.getElementById('snp-window-size').onchange = (e) => {{
            windowSize = parseInt(e.target.value);
            loadRegion(currentStart, currentStart + windowSize);
        }};

        document.getElementById('snp-tracks-prev').onclick = () => window.snpTracksNav(-windowSize);
        document.getElementById('snp-tracks-next').onclick = () => window.snpTracksNav(windowSize);

        // Quality popup function
        window.showQualityPopup = (pos, sample) => {{
            if (!regionData) return;

            const snp = regionData.snps.find(s => s.pos === pos && s.sample === sample);
            if (!snp) return;

            // Build popup content
            let content = `<div style="font-family: monospace; font-size: 13px;">`;
            content += `<h3 style="margin: 0 0 15px 0; color: #333;">${{sample}} @ pos ${{pos.toLocaleString()}}</h3>`;
            content += `<div style="margin-bottom: 10px;"><strong>Variant:</strong> ${{snp.ref}} → ${{snp.alt}}</div>`;
            content += `<div style="margin-bottom: 15px;"><strong>Status:</strong> <span style="background: ${{
                snp.status === 'Concordant' ? '#4caf50' :
                snp.status === 'OnlyPipelineA' ? '#ff9800' :
                snp.status === 'OnlyPipelineB' ? '#2196f3' : '#9c27b0'
            }}; color: white; padding: 2px 8px; border-radius: 4px;">${{snp.status}}</span></div>`;

            // Snippy quality
            content += `<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">`;

            if (snp.snippy) {{
                content += `<div style="background: #fff3e0; padding: 12px; border-radius: 6px; border-left: 4px solid #ff9800;">`;
                content += `<div style="font-weight: bold; margin-bottom: 8px; color: #e65100;">Snippy (FreeBayes)</div>`;
                content += `<table style="width: 100%; font-size: 12px;">`;
                if (snp.snippy.depth != null) content += `<tr><td>Depth (DP):</td><td><strong>${{snp.snippy.depth}}</strong></td></tr>`;
                if (snp.snippy.qual != null) content += `<tr><td>Quality:</td><td><strong>${{snp.snippy.qual.toFixed(1)}}</strong></td></tr>`;
                if (snp.snippy.ref_obs != null) content += `<tr><td>Ref reads (RO):</td><td>${{snp.snippy.ref_obs}}</td></tr>`;
                if (snp.snippy.alt_obs != null) content += `<tr><td>Alt reads (AO):</td><td><strong>${{snp.snippy.alt_obs}}</strong></td></tr>`;
                if (snp.snippy.alt_qual != null) content += `<tr><td>Alt quality (QA):</td><td>${{snp.snippy.alt_qual}}</td></tr>`;
                if (snp.snippy.allele_freq != null) content += `<tr><td>Allele freq:</td><td>${{(snp.snippy.allele_freq * 100).toFixed(1)}}%</td></tr>`;
                if (snp.snippy.genotype) content += `<tr><td>Genotype:</td><td>${{snp.snippy.genotype}}</td></tr>`;
                if (snp.snippy.filter) content += `<tr><td>Filter:</td><td>${{snp.snippy.filter}}</td></tr>`;
                content += `</table></div>`;
            }} else {{
                content += `<div style="background: #f5f5f5; padding: 12px; border-radius: 6px; color: #999;">`;
                content += `<div style="font-weight: bold; margin-bottom: 8px;">Snippy</div>`;
                content += `<em>No variant called</em></div>`;
            }}

            // CFSAN quality
            if (snp.cfsan) {{
                content += `<div style="background: #e3f2fd; padding: 12px; border-radius: 6px; border-left: 4px solid #2196f3;">`;
                content += `<div style="font-weight: bold; margin-bottom: 8px; color: #1565c0;">CFSAN (VarScan2)</div>`;
                content += `<table style="width: 100%; font-size: 12px;">`;
                if (snp.cfsan.depth != null) content += `<tr><td>Depth (DP):</td><td><strong>${{snp.cfsan.depth}}</strong></td></tr>`;
                if (snp.cfsan.genotype_qual != null) content += `<tr><td>Genotype qual (GQ):</td><td><strong>${{snp.cfsan.genotype_qual}}</strong></td></tr>`;
                if (snp.cfsan.ref_obs != null) content += `<tr><td>Ref reads (RD):</td><td>${{snp.cfsan.ref_obs}}</td></tr>`;
                if (snp.cfsan.alt_obs != null) content += `<tr><td>Alt reads (AD):</td><td><strong>${{snp.cfsan.alt_obs}}</strong></td></tr>`;
                if (snp.cfsan.avg_base_qual != null) content += `<tr><td>Avg base qual (ABQ):</td><td>${{snp.cfsan.avg_base_qual}}</td></tr>`;
                if (snp.cfsan.pvalue != null) content += `<tr><td>P-value:</td><td>${{snp.cfsan.pvalue.toExponential(2)}}</td></tr>`;
                if (snp.cfsan.allele_freq != null) content += `<tr><td>Allele freq:</td><td>${{(snp.cfsan.allele_freq * 100).toFixed(1)}}%</td></tr>`;
                if (snp.cfsan.genotype) content += `<tr><td>Genotype:</td><td>${{snp.cfsan.genotype}}</td></tr>`;
                if (snp.cfsan.filter) content += `<tr><td>Filter:</td><td>${{snp.cfsan.filter}}</td></tr>`;
                content += `</table></div>`;
            }} else {{
                content += `<div style="background: #f5f5f5; padding: 12px; border-radius: 6px; color: #999;">`;
                content += `<div style="font-weight: bold; margin-bottom: 8px;">CFSAN</div>`;
                content += `<em>No variant called</em></div>`;
            }}

            content += `</div></div>`;

            // Create and show modal
            let modal = document.getElementById('quality-modal');
            if (!modal) {{
                modal = document.createElement('div');
                modal.id = 'quality-modal';
                modal.style.cssText = 'position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.5); display: flex; justify-content: center; align-items: center; z-index: 10000;';
                modal.onclick = (e) => {{ if (e.target === modal) modal.remove(); }};
                document.body.appendChild(modal);
            }}

            modal.innerHTML = `
                <div style="background: white; padding: 25px; border-radius: 12px; max-width: 600px; max-height: 80vh; overflow: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
                    ${{content}}
                    <div style="text-align: right; margin-top: 15px;">
                        <button onclick="document.getElementById('quality-modal').remove()" style="padding: 8px 20px; background: #666; color: white; border: none; border-radius: 4px; cursor: pointer;">Close</button>
                    </div>
                </div>
            `;
        }};

        // Initial load - find first region with SNPs
        loadRegion(1, windowSize + 1);
    }}

    // Legacy Snippy comparison (old format)
    function populateLegacySnippyComparison(container, comparison) {{
        if (!comparison) return;

        const summary = comparison.summary;
        const totalUnion = summary.total_concordant_bases + summary.total_coreguard_only_bases + summary.total_snp_only_bases;

        let html = `
            <!-- Summary Cards -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 25px;">
                <div style="background: #e8f5e9; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.85rem; color: #388e3c;">Pipeline</div>
                    <div style="font-size: 1.5rem; font-weight: bold; color: #2e7d32;">${{comparison.pipeline_name}}</div>
                </div>
                <div style="background: #e3f2fd; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.85rem; color: #1976d2;">Concordance</div>
                    <div style="font-size: 1.5rem; font-weight: bold; color: #1565c0;">${{summary.avg_concordance_pct.toFixed(1)}}%</div>
                </div>
                <div style="background: #fff8e1; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.85rem; color: #f57c00;">Concordant</div>
                    <div style="font-size: 1.5rem; font-weight: bold; color: #ef6c00;">${{summary.total_concordant_bases.toLocaleString()}} bp</div>
                </div>
            </div>

            <!-- Interpretation -->
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9e9e9e; margin-bottom: 25px;">
                <p style="margin: 0; font-style: italic;">${{summary.interpretation}}</p>
            </div>

            <!-- Comparison Breakdown -->
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 25px; margin-bottom: 25px;">
                <!-- Visual Breakdown -->
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="margin-bottom: 15px;">Overall Breakdown</h4>
                    <div style="display: flex; height: 30px; border-radius: 4px; overflow: hidden; margin-bottom: 15px;">
                        <div style="width: ${{(summary.total_concordant_bases/totalUnion*100).toFixed(1)}}%; background: #4caf50;" title="Concordant"></div>
                        <div style="width: ${{(summary.total_coreguard_only_bases/totalUnion*100).toFixed(1)}}%; background: #2196f3;" title="coreguard-only"></div>
                        <div style="width: ${{(summary.total_snp_only_bases/totalUnion*100).toFixed(1)}}%; background: #ff9800;" title="SNP-only"></div>
                    </div>
                    <div style="display: flex; gap: 20px; font-size: 0.85rem;">
                        <div><span style="display: inline-block; width: 12px; height: 12px; background: #4caf50; border-radius: 2px; margin-right: 5px;"></span>Concordant</div>
                        <div><span style="display: inline-block; width: 12px; height: 12px; background: #2196f3; border-radius: 2px; margin-right: 5px;"></span>coreguard-only</div>
                        <div><span style="display: inline-block; width: 12px; height: 12px; background: #ff9800; border-radius: 2px; margin-right: 5px;"></span>${{comparison.pipeline_name}}-only</div>
                    </div>
                </div>

                <!-- Numbers -->
                <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                    <h4 style="margin-bottom: 15px;">Total Bases</h4>
                    <table style="width: 100%; font-size: 0.9rem;">
                        <tr>
                            <td style="padding: 8px 0;"><span style="color: #4caf50; font-weight: bold;">&#9632;</span> Concordant (both filter)</td>
                            <td style="text-align: right; font-weight: bold;">${{summary.total_concordant_bases.toLocaleString()}} bp</td>
                        </tr>
                        <tr>
                            <td style="padding: 8px 0;"><span style="color: #2196f3; font-weight: bold;">&#9632;</span> coreguard-only (CG flags, SNP keeps)</td>
                            <td style="text-align: right; font-weight: bold;">${{summary.total_coreguard_only_bases.toLocaleString()}} bp</td>
                        </tr>
                        <tr>
                            <td style="padding: 8px 0;"><span style="color: #ff9800; font-weight: bold;">&#9632;</span> ${{comparison.pipeline_name}}-only (SNP filters, CG keeps)</td>
                            <td style="text-align: right; font-weight: bold;">${{summary.total_snp_only_bases.toLocaleString()}} bp</td>
                        </tr>
                    </table>
                </div>
            </div>

            <!-- Per-Sample Table -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                <h4 style="margin-bottom: 15px;">Per-Sample Comparison</h4>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">CG Gaps</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{comparison.pipeline_name}} Filtered</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Concordant</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">CG-only</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">${{comparison.pipeline_name}}-only</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Concordance</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const s of comparison.samples) {{
            const concordanceColor = s.concordance_pct > 50 ? '#4caf50' : s.concordance_pct > 20 ? '#ff9800' : '#f44336';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{s.sample}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.coreguard_gap_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.snp_filtered_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #4caf50;">${{s.concordant_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #2196f3;">${{s.coreguard_only_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #ff9800;">${{s.snp_only_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{concordanceColor}}; font-weight: bold;">${{s.concordance_pct.toFixed(1)}}%</span>
                    </td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <!-- Gap Distribution Overlay -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 25px;">
                <h4 style="margin-bottom: 15px;">Gap Distribution Overlay</h4>
                <p style="color: #666; font-size: 0.85rem; margin-bottom: 15px;">
                    Overlapping view of coreguard gaps (blue) and ${{comparison.pipeline_name}} masked regions (orange) along the reference genome.
                    Green regions indicate concordant filtering.
                </p>
                <div id="gap-overlay-chart"></div>
            </div>

            <!-- Explanation -->
            <div style="margin-top: 20px; padding: 15px; background: #e8f4fd; border-radius: 8px; font-size: 0.9rem;">
                <strong>Interpretation:</strong>
                <ul style="margin: 10px 0 0 20px;">
                    <li><strong>coreguard-only</strong>: Positions that coreguard flags as gaps (coverage &lt; min-depth) but ${{comparison.pipeline_name}} uses in the alignment. These may be low-coverage regions that still have usable data.</li>
                    <li><strong>${{comparison.pipeline_name}}-only</strong>: Positions that ${{comparison.pipeline_name}} masks but coreguard keeps. These are typically filtered for quality reasons (mapping quality, base quality, strand bias) rather than coverage.</li>
                    <li><strong>Concordant</strong>: Both tools agree to filter these positions - highest confidence problematic regions.</li>
                </ul>
            </div>
        `;

        container.innerHTML = html;

        // Generate gap overlay chart
        const chartContainer = document.getElementById('gap-overlay-chart');
        if (chartContainer && comparison.samples.length > 0) {{
            const refLen = comparison.reference_length;
            const width = 900;
            const sampleHeight = 60;
            const margin = {{ left: 100, right: 20, top: 10, bottom: 30 }};
            const chartWidth = width - margin.left - margin.right;
            const height = margin.top + margin.bottom + comparison.samples.length * sampleHeight;

            let svg = `<svg width="${{width}}" height="${{height}}" style="font-family: sans-serif; font-size: 11px;">`;

            // Scale function
            const scale = (pos) => margin.left + (pos / refLen) * chartWidth;

            // Draw samples
            comparison.samples.forEach((s, i) => {{
                const y = margin.top + i * sampleHeight;

                // Sample label
                svg += `<text x="${{margin.left - 10}}" y="${{y + 25}}" text-anchor="end" font-weight="bold">${{s.sample}}</text>`;

                // Background track
                svg += `<rect x="${{margin.left}}" y="${{y + 5}}" width="${{chartWidth}}" height="40" fill="#f5f5f5" stroke="#ddd"/>`;

                // coreguard gaps (blue, top half)
                if (s.coreguard_only_regions) {{
                    s.coreguard_only_regions.forEach(r => {{
                        const x = scale(r.start);
                        const w = Math.max(1, scale(r.end) - scale(r.start));
                        svg += `<rect x="${{x}}" y="${{y + 5}}" width="${{w}}" height="18" fill="#2196f3" opacity="0.7"/>`;
                    }});
                }}

                // Snippy gaps (orange, bottom half)
                if (s.snp_only_regions) {{
                    s.snp_only_regions.forEach(r => {{
                        const x = scale(r.start);
                        const w = Math.max(1, scale(r.end) - scale(r.start));
                        svg += `<rect x="${{x}}" y="${{y + 27}}" width="${{w}}" height="18" fill="#ff9800" opacity="0.7"/>`;
                    }});
                }}

                // Track labels
                svg += `<text x="${{margin.left + 5}}" y="${{y + 18}}" fill="#1565c0" font-size="9">CG</text>`;
                svg += `<text x="${{margin.left + 5}}" y="${{y + 40}}" fill="#e65100" font-size="9">${{comparison.pipeline_name.substring(0,3)}}</text>`;
            }});

            // X axis
            const axisY = height - margin.bottom + 15;
            svg += `<line x1="${{margin.left}}" y1="${{axisY - 10}}" x2="${{margin.left + chartWidth}}" y2="${{axisY - 10}}" stroke="#333"/>`;

            // Tick marks
            for (let i = 0; i <= 10; i++) {{
                const pos = (refLen / 10) * i;
                const x = scale(pos);
                svg += `<line x1="${{x}}" y1="${{axisY - 10}}" x2="${{x}}" y2="${{axisY - 5}}" stroke="#333"/>`;
                const label = pos >= 1000000 ? (pos/1000000).toFixed(1) + 'M' : pos >= 1000 ? (pos/1000).toFixed(0) + 'k' : pos;
                svg += `<text x="${{x}}" y="${{axisY + 5}}" text-anchor="middle" font-size="10">${{label}}</text>`;
            }}

            // Legend
            svg += `<rect x="${{margin.left}}" y="${{height - 12}}" width="12" height="12" fill="#2196f3" opacity="0.7"/>`;
            svg += `<text x="${{margin.left + 16}}" y="${{height - 2}}" font-size="10">coreguard-only</text>`;
            svg += `<rect x="${{margin.left + 110}}" y="${{height - 12}}" width="12" height="12" fill="#ff9800" opacity="0.7"/>`;
            svg += `<text x="${{margin.left + 126}}" y="${{height - 2}}" font-size="10">${{comparison.pipeline_name}}-only</text>`;

            svg += '</svg>';
            chartContainer.innerHTML = svg;
        }}
    }}

    // Populate Snippy-core Comparison (new redesigned view)
    function populateSnippyCoreComparison(container, sc) {{
        const summary = sc.summary;
        // Reference length available for future use: sc.reference_length

        // Determine status colors
        const snpRetentionColor = summary.avg_snp_retention_pct > 99 ? '#4caf50' :
                                  summary.avg_snp_retention_pct > 95 ? '#8bc34a' :
                                  summary.avg_snp_retention_pct > 90 ? '#ff9800' : '#f44336';
        const coreDiffColor = summary.avg_core_diff_pct < 1 ? '#4caf50' :
                              summary.avg_core_diff_pct < 3 ? '#ff9800' : '#f44336';

        let html = `
            <h4 style="margin-bottom: 20px;">Snippy-core Validation</h4>

            <!-- About This Section -->
            <div style="background: #e3f2fd; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 0.9rem;">
                <strong>What is Snippy?</strong>
                <p style="margin: 8px 0 0 0;">Snippy is a reference-based SNP caller that maps reads to a reference genome and identifies variants.
                It produces a multi-sample alignment (snippy-core) showing SNP positions across all samples.</p>
                <p style="margin: 8px 0 0 0;"><strong>This comparison shows:</strong> How coreguard's gap predictions (based on coverage alone)
                align with snippy-core's actual SNP calling. High SNP retention means coreguard's gaps don't affect the variants that matter.</p>
            </div>

            <!-- Files Used -->
            <details style="margin-bottom: 20px; background: #eceff1; border-radius: 8px; padding: 15px;">
                <summary style="cursor: pointer; font-weight: bold; color: #455a64;">Files Used for This Analysis</summary>
                <div style="margin-top: 15px; font-size: 0.85rem;">
                    <table style="width: 100%; border-collapse: collapse;">
                        <tr>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;"><code>snippycore.txt</code></td>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;">Per-sample statistics (aligned bases, lowcov, variant count)</td>
                        </tr>
                        <tr>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;"><code>snippycore.tab</code></td>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;">SNP positions with per-sample alleles (A/T/G/C/N/-)</td>
                        </tr>
                        <tr>
                            <td style="padding: 6px;"><code>hamming_distances_*.tsv</code></td>
                            <td style="padding: 6px;">Pairwise SNP distance matrix between samples</td>
                        </tr>
                    </table>
                </div>
            </details>

            <!-- Command Line Reference -->
            <details style="margin-bottom: 20px; background: #f5f5f5; border-radius: 8px; padding: 15px;">
                <summary style="cursor: pointer; font-weight: bold; color: #333;">Example Command Lines (Reference)</summary>
                <div style="margin-top: 15px;">
                    <div style="background: #fff3e0; padding: 10px; border-radius: 4px; margin-bottom: 15px; font-size: 0.8rem;">
                        <strong>Note:</strong> These are typical command lines for Snippy. The actual commands used to generate your data
                        may vary. Check your pipeline logs or documentation for the exact parameters used.
                    </div>
                    <p style="font-size: 0.85rem; margin-bottom: 10px;"><strong>1. Per-sample Snippy:</strong></p>
                    <pre style="background: #263238; color: #aed581; padding: 12px; border-radius: 4px; overflow-x: auto; font-size: 0.8rem;">snippy --outdir results/\${{SAMPLE}}_snippy \\
  --ref reference.fasta \\
  --R1 reads/\${{SAMPLE}}_R1.fq.gz \\
  --R2 reads/\${{SAMPLE}}_R2.fq.gz \\
  --cpus 4</pre>
                    <p style="font-size: 0.85rem; margin: 15px 0 10px 0;"><strong>2. Snippy-core (multi-sample alignment):</strong></p>
                    <pre style="background: #263238; color: #aed581; padding: 12px; border-radius: 4px; overflow-x: auto; font-size: 0.8rem;">snippy-core --ref reference.fasta --prefix core results/*_snippy</pre>
                    <p style="font-size: 0.75rem; color: #666; margin-top: 10px;">
                        Key parameters that affect results: <code>--mincov</code> (min coverage, default 10),
                        <code>--minfrac</code> (min variant fraction, default 0.9), <code>--mapqual</code> (min mapping quality, default 60)
                    </p>
                </div>
            </details>

            <!-- Summary Cards -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin-bottom: 25px;">
                <div style="background: #e8f5e9; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #388e3c; text-transform: uppercase;">SNP Retention</div>
                    <div style="font-size: 2rem; font-weight: bold; color: ${{snpRetentionColor}};">${{summary.avg_snp_retention_pct.toFixed(1)}}%</div>
                    <div style="font-size: 0.75rem; color: #666;">of SNPs preserved</div>
                </div>
                <div style="background: #e3f2fd; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #1976d2; text-transform: uppercase;">Core Diff</div>
                    <div style="font-size: 2rem; font-weight: bold; color: ${{coreDiffColor}}">${{summary.avg_core_diff_pct.toFixed(2)}}%</div>
                    <div style="font-size: 0.75rem; color: #666;">avg size difference</div>
                </div>
                <div style="background: #fff8e1; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #f57c00; text-transform: uppercase;">Total SNPs</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #ef6c00;">${{summary.total_snps.toLocaleString()}}</div>
                    <div style="font-size: 0.75rem; color: #666;">across all samples</div>
                </div>
                <div style="background: #fce4ec; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #c2185b; text-transform: uppercase;">SNPs in Gaps</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #ad1457;">${{summary.total_snps_in_gaps.toLocaleString()}}</div>
                    <div style="font-size: 0.75rem; color: #666;">would be lost</div>
                </div>
            </div>

            <!-- Interpretation -->
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px; border-left: 4px solid ${{snpRetentionColor}}; margin-bottom: 25px;">
                <p style="margin: 0;"><strong>Interpretation:</strong> ${{summary.interpretation}}</p>
            </div>

            <!-- Core Size Comparison -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">Table 1: Core Size Comparison</h4>
                <div style="background: #fff8e1; padding: 12px; border-radius: 6px; margin-bottom: 15px; font-size: 0.85rem;">
                    <strong>What this table shows:</strong> Compares how much of the reference genome each tool considers "usable."
                    <ul style="margin: 8px 0 0 15px; padding: 0;">
                        <li><strong>coreguard Core</strong> = Reference length - gap bases (positions with zero coverage)</li>
                        <li><strong>Snippy-core Aligned</strong> = Actual aligned bases from snippy-core.txt output</li>
                        <li><strong>Difference</strong> = How much the two estimates differ (small = good agreement)</li>
                    </ul>
                </div>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">coreguard Core</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">coreguard %</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Snippy-core Aligned</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Snippy-core %</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Difference</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const c of sc.core_comparisons) {{
            const diffColor = Math.abs(c.core_diff_pct) < 1 ? '#4caf50' :
                              Math.abs(c.core_diff_pct) < 3 ? '#ff9800' : '#f44336';
            const diffSign = c.core_diff > 0 ? '+' : '';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{c.sample}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{c.cg_core_size.toLocaleString()}} bp</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{c.cg_core_pct.toFixed(1)}}%</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{c.sc_aligned.toLocaleString()}} bp</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{c.sc_aligned_pct.toFixed(1)}}%</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{diffColor}}; font-weight: bold;">${{diffSign}}${{c.core_diff_pct.toFixed(2)}}%</span>
                    </td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <!-- SNP Impact Analysis -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">Table 2: SNP Impact Analysis</h4>
                <div style="background: #e8f5e9; padding: 12px; border-radius: 6px; margin-bottom: 15px; font-size: 0.85rem;">
                    <strong>What this table shows:</strong> How many SNPs would be "lost" if you excluded coreguard's gap regions.
                    <ul style="margin: 8px 0 0 15px; padding: 0;">
                        <li><strong>Total SNPs</strong> = Number of positions in snippy-core.tab where this sample has a valid call (not N or -)</li>
                        <li><strong>SNPs Retained</strong> = SNPs that are NOT in coreguard's gap regions (these are safe)</li>
                        <li><strong>SNPs in Gaps</strong> = SNPs that fall within gap regions (would be lost if you filter by coverage)</li>
                        <li><strong>Retention %</strong> = SNPs Retained / Total SNPs × 100 (higher is better)</li>
                    </ul>
                    <p style="margin: 8px 0 0 0;"><em>Note: "Total SNPs" here refers to positions with valid calls, not total variants found by Snippy.</em></p>
                </div>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Total SNPs</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">SNPs Retained</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">SNPs in Gaps</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Retention %</th>
                            <th style="padding: 10px; text-align: center; border: 1px solid #dee2e6;">Status</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const s of sc.snp_impacts) {{
            const retColor = s.retention_pct > 99 ? '#4caf50' :
                            s.retention_pct > 95 ? '#8bc34a' :
                            s.retention_pct > 90 ? '#ff9800' : '#f44336';
            const status = s.retention_pct > 99 ? '&#10004; Excellent' :
                          s.retention_pct > 95 ? '&#10004; Good' :
                          s.retention_pct > 90 ? '&#9888; Moderate' : '&#10006; Review';
            const statusColor = s.retention_pct > 95 ? '#4caf50' :
                               s.retention_pct > 90 ? '#ff9800' : '#f44336';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{s.sample}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.total_snps.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #4caf50;">${{s.snps_retained.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #f44336;">${{s.snps_in_gaps.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{retColor}}; font-weight: bold;">${{s.retention_pct.toFixed(1)}}%</span>
                    </td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: center; color: ${{statusColor}};">${{status}}</td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <!-- Explanation -->
            <div style="margin-top: 20px; padding: 15px; background: #e8f4fd; border-radius: 8px; font-size: 0.9rem;">
                <strong>What this comparison tells you:</strong>
                <ul style="margin: 10px 0 0 20px;">
                    <li><strong>SNP Retention %</strong>: The percentage of actual variants that are NOT in coreguard's gap regions.
                        High values (>99%) mean coreguard's gaps won't affect your variant calling.</li>
                    <li><strong>Core Size Difference</strong>: How coreguard's estimated core compares to snippy-core's actual aligned bases.
                        Small differences (<1%) indicate good agreement.</li>
                    <li><strong>SNPs in Gaps</strong>: These are real variants that would be lost if you excluded coreguard's gap regions.
                        If this number is high, consider adjusting your min-depth threshold.</li>
                </ul>
                <div style="margin-top: 15px; padding: 10px; background: #fff3e0; border-radius: 4px;">
                    <strong>Note:</strong> coreguard uses coverage-based gap detection, while snippy-core applies additional quality filters.
                    Some difference is expected - the key metric is SNP retention.
                </div>
            </div>
        `;

        // Add distance matrix if available
        if (sc.distance_matrix) {{
            html += renderDistanceMatrix(sc.distance_matrix, 'Snippy-core', '#1565c0');
        }}

        container.innerHTML = html;
    }}

    // Render a distance matrix as a heatmap table
    function renderDistanceMatrix(matrix, title, headerColor) {{
        if (!matrix || !matrix.samples || matrix.samples.length === 0) {{
            return '';
        }}

        const samples = matrix.samples;
        const distances = matrix.distances;
        const n = samples.length;

        // Find max distance for color scaling
        let maxDist = 0;
        for (const row of distances) {{
            for (const d of row) {{
                if (d > maxDist) maxDist = d;
            }}
        }}

        // Generate color based on distance
        const getColor = (d) => {{
            if (d === 0) return '#e8f5e9'; // Light green for zero
            if (maxDist === 0) return '#fff8e1';
            const ratio = d / maxDist;
            if (ratio < 0.25) return '#c8e6c9'; // Green
            if (ratio < 0.5) return '#fff9c4'; // Yellow
            if (ratio < 0.75) return '#ffcc80'; // Orange
            return '#ef9a9a'; // Red
        }};

        let html = `
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 25px;">
                <h4 style="margin-bottom: 15px;">SNP Distance Matrix (${{title}})</h4>
                <div style="background: #fff8e1; padding: 12px; border-radius: 6px; margin-bottom: 15px; font-size: 0.85rem;">
                    <strong>What this shows:</strong> Pairwise SNP distances between samples.
                    Colors indicate distance magnitude: <span style="background:#e8f5e9;padding:2px 6px;border-radius:3px;">0</span>
                    <span style="background:#c8e6c9;padding:2px 6px;border-radius:3px;">low</span>
                    <span style="background:#fff9c4;padding:2px 6px;border-radius:3px;">medium</span>
                    <span style="background:#ffcc80;padding:2px 6px;border-radius:3px;">high</span>
                    <span style="background:#ef9a9a;padding:2px 6px;border-radius:3px;">max</span>
                </div>
                <div style="overflow-x: auto;">
                    <table style="border-collapse: collapse; font-size: 0.85rem; min-width: 100%;">
                        <thead>
                            <tr>
                                <th style="padding: 8px; background: ${{headerColor}}; color: white; border: 1px solid #dee2e6;"></th>
        `;

        // Header row
        for (const s of samples) {{
            html += `<th style="padding: 8px; background: ${{headerColor}}; color: white; border: 1px solid #dee2e6; font-weight: bold;">${{s}}</th>`;
        }}
        html += `</tr></thead><tbody>`;

        // Data rows
        for (let i = 0; i < n; i++) {{
            html += `<tr>
                <td style="padding: 8px; background: ${{headerColor}}; color: white; border: 1px solid #dee2e6; font-weight: bold;">${{samples[i]}}</td>`;
            for (let j = 0; j < n; j++) {{
                const d = distances[i] ? distances[i][j] : 0;
                const bg = getColor(d);
                const fontWeight = i === j ? 'normal' : 'bold';
                html += `<td style="padding: 8px; text-align: center; border: 1px solid #dee2e6; background: ${{bg}}; font-weight: ${{fontWeight}};">${{d}}</td>`;
            }}
            html += `</tr>`;
        }}

        html += `</tbody></table></div>`;

        // Add max distance note
        if (maxDist > 0) {{
            html += `<p style="margin-top: 10px; font-size: 0.8rem; color: #666;">Max distance: ${{maxDist}} SNPs</p>`;
        }} else {{
            html += `<p style="margin-top: 10px; font-size: 0.8rem; color: #f44336;"><strong>Warning:</strong> All distances are zero - no polymorphic SNPs detected!</p>`;
        }}

        html += `</div>`;
        return html;
    }}

    // Populate CFSAN Comparison
    function populateCfsanComparison(container, cfsan) {{
        const summary = cfsan.summary;

        // Determine status colors
        const snpRetentionColor = summary.avg_snp_retention_pct > 99 ? '#4caf50' :
                                  summary.avg_snp_retention_pct > 95 ? '#8bc34a' :
                                  summary.avg_snp_retention_pct > 90 ? '#ff9800' : '#f44336';

        let html = `
            <h4 style="margin-bottom: 20px; color: #7b1fa2;">CFSAN SNP Pipeline Validation</h4>

            <!-- About This Section -->
            <div style="background: #f3e5f5; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 0.9rem;">
                <strong>What is CFSAN SNP Pipeline?</strong>
                <p style="margin: 8px 0 0 0;">CFSAN is the CDC/FDA SNP Pipeline, a reference-based SNP caller commonly used for foodborne pathogen surveillance.
                It generates detailed metrics including coverage depth, mapping rates, and missing positions.</p>
                <p style="margin: 8px 0 0 0;"><strong>This comparison shows:</strong> How coreguard's gap predictions compare with CFSAN's missing positions.
                Unlike Snippy, CFSAN provides explicit coverage metrics (Avg Depth) because it outputs a detailed metrics.tsv file.</p>
            </div>

            <!-- Files Used -->
            <details style="margin-bottom: 20px; background: #ede7f6; border-radius: 8px; padding: 15px;">
                <summary style="cursor: pointer; font-weight: bold; color: #512da8;">Files Used for This Analysis</summary>
                <div style="margin-top: 15px; font-size: 0.85rem;">
                    <table style="width: 100%; border-collapse: collapse;">
                        <tr>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;"><code>metrics.tsv</code></td>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;">Per-sample metrics (% reads mapped, avg depth, phase2 SNPs, missing positions)</td>
                        </tr>
                        <tr>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;"><code>snplist.txt</code></td>
                            <td style="padding: 6px; border-bottom: 1px solid #ccc;">List of SNP positions with per-sample alleles</td>
                        </tr>
                        <tr>
                            <td style="padding: 6px;"><code>snp_distance_matrix.tsv</code></td>
                            <td style="padding: 6px;">Pairwise SNP distance matrix between samples</td>
                        </tr>
                    </table>
                </div>
            </details>

            <!-- Command Line Reference -->
            <details style="margin-bottom: 20px; background: #f5f5f5; border-radius: 8px; padding: 15px;">
                <summary style="cursor: pointer; font-weight: bold; color: #333;">Example Command Line (Reference)</summary>
                <div style="margin-top: 15px;">
                    <div style="background: #fff3e0; padding: 10px; border-radius: 4px; margin-bottom: 15px; font-size: 0.8rem;">
                        <strong>Note:</strong> This is a typical CFSAN command line. The actual commands used to generate your data
                        may vary. Check your pipeline configuration files for the exact parameters used.
                    </div>
                    <pre style="background: #263238; color: #ce93d8; padding: 12px; border-radius: 4px; overflow-x: auto; font-size: 0.8rem;">cfsan_snp_pipeline run \\
  -o results/cfsan \\
  -s samples_dir \\
  -m soft \\
  reference.fasta</pre>
                    <p style="font-size: 0.75rem; color: #666; margin-top: 10px;">
                        Key parameters: <code>-c config.ini</code> for custom config (minBaseQual, minMapQual, minConsFreq).
                        Default: minBaseQual=0, minMapQual=30, minConsFreq=0.9
                    </p>
                </div>
            </details>

            <!-- Summary Cards -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin-bottom: 25px;">
                <div style="background: #f3e5f5; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #7b1fa2; text-transform: uppercase;">SNP Retention</div>
                    <div style="font-size: 2rem; font-weight: bold; color: ${{snpRetentionColor}};">${{summary.avg_snp_retention_pct.toFixed(1)}}%</div>
                    <div style="font-size: 0.75rem; color: #666;">of SNPs preserved</div>
                </div>
                <div style="background: #ede7f6; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #512da8; text-transform: uppercase;">Total SNPs</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #673ab7;">${{summary.total_snps.toLocaleString()}}</div>
                    <div style="font-size: 0.75rem; color: #666;">across all samples</div>
                </div>
                <div style="background: #fce4ec; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #c2185b; text-transform: uppercase;">SNPs in Gaps</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #ad1457;">${{summary.total_snps_in_gaps.toLocaleString()}}</div>
                    <div style="font-size: 0.75rem; color: #666;">would be lost</div>
                </div>
                <div style="background: #e8f5e9; padding: 20px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 0.8rem; color: #388e3c; text-transform: uppercase;">CFSAN Positions</div>
                    <div style="font-size: 2rem; font-weight: bold; color: #2e7d32;">${{cfsan.total_cfsan_snps.toLocaleString()}}</div>
                    <div style="font-size: 0.75rem; color: #666;">total SNP positions</div>
                </div>
            </div>

            <!-- Interpretation -->
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px; border-left: 4px solid ${{snpRetentionColor}}; margin-bottom: 25px;">
                <p style="margin: 0;"><strong>Interpretation:</strong> ${{summary.interpretation}}</p>
            </div>

            <!-- CFSAN Metrics -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">Table 1: CFSAN Sample Metrics</h4>
                <div style="background: #fff8e1; padding: 12px; border-radius: 6px; margin-bottom: 15px; font-size: 0.85rem;">
                    <strong>What this table shows:</strong> Mapping and coverage statistics from CFSAN's <code>metrics.tsv</code> output.
                    <ul style="margin: 8px 0 0 15px; padding: 0;">
                        <li><strong>Reads Mapped</strong> = Percentage of input reads that aligned to the reference</li>
                        <li><strong>Avg Depth</strong> = Mean coverage depth across the genome (CFSAN provides this; Snippy does not)</li>
                        <li><strong>Phase2 SNPs</strong> = Final high-quality SNPs after all filtering (Phase 1 = initial, Phase 2 = filtered)</li>
                        <li><strong>Missing Pos</strong> = Positions excluded from the SNP matrix due to low coverage or quality</li>
                    </ul>
                    <p style="margin: 8px 0 0 0; color: #666;"><em>Note: Snippy doesn't output coverage depth in a summary file, which is why it doesn't appear in the Snippy comparison.</em></p>
                </div>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Reads Mapped</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Avg Depth</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Phase2 SNPs</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Missing Pos</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const m of cfsan.cfsan_metrics) {{
            const depthColor = m.average_depth > 20 ? '#4caf50' : m.average_depth > 10 ? '#ff9800' : '#f44336';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{m.sample}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{m.percent_reads_mapped.toFixed(1)}}%</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{depthColor}};">${{m.average_depth.toFixed(1)}}x</span>
                    </td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{m.phase2_snps.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{m.missing_positions.toLocaleString()}}</td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <!-- SNP Impact Analysis -->
            <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px;">
                <h4 style="margin-bottom: 15px;">Table 2: SNP Impact Analysis</h4>
                <div style="background: #e8f5e9; padding: 12px; border-radius: 6px; margin-bottom: 15px; font-size: 0.85rem;">
                    <strong>What this table shows:</strong> How many CFSAN SNP positions fall within coreguard's gap regions.
                    <ul style="margin: 8px 0 0 15px; padding: 0;">
                        <li><strong>CG Gap Bases</strong> = Total bases that coreguard flagged as gaps (zero coverage)</li>
                        <li><strong>CFSAN SNPs</strong> = SNP positions reported by CFSAN for this sample</li>
                        <li><strong>SNPs in Gaps</strong> = CFSAN SNPs that fall within coreguard's gap regions</li>
                        <li><strong>Retention %</strong> = Percentage of SNPs NOT in gap regions (higher is better)</li>
                    </ul>
                </div>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                    <thead>
                        <tr style="background: #f1f3f4;">
                            <th style="padding: 10px; text-align: left; border: 1px solid #dee2e6;">Sample</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">CG Gap Bases</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">CFSAN SNPs</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">SNPs in Gaps</th>
                            <th style="padding: 10px; text-align: right; border: 1px solid #dee2e6;">Retention %</th>
                            <th style="padding: 10px; text-align: center; border: 1px solid #dee2e6;">Status</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        for (const s of cfsan.sample_comparisons) {{
            const retColor = s.snp_retention_pct > 99 ? '#4caf50' :
                            s.snp_retention_pct > 95 ? '#8bc34a' :
                            s.snp_retention_pct > 90 ? '#ff9800' : '#f44336';
            const status = s.snp_retention_pct > 99 ? '&#10004; Excellent' :
                          s.snp_retention_pct > 95 ? '&#10004; Good' :
                          s.snp_retention_pct > 90 ? '&#9888; Moderate' : '&#10006; Review';
            const statusColor = s.snp_retention_pct > 95 ? '#4caf50' :
                               s.snp_retention_pct > 90 ? '#ff9800' : '#f44336';
            html += `
                <tr>
                    <td style="padding: 10px; border: 1px solid #dee2e6; font-weight: bold;">${{s.sample}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.cg_gap_bases.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">${{s.cfsan_snps.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right; color: #f44336;">${{s.snps_in_gaps.toLocaleString()}}</td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: right;">
                        <span style="color: ${{retColor}}; font-weight: bold;">${{s.snp_retention_pct.toFixed(1)}}%</span>
                    </td>
                    <td style="padding: 10px; border: 1px solid #dee2e6; text-align: center; color: ${{statusColor}};">${{status}}</td>
                </tr>
            `;
        }}

        html += `
                    </tbody>
                </table>
            </div>

            <!-- Explanation -->
            <div style="margin-top: 20px; padding: 15px; background: #f3e5f5; border-radius: 8px; font-size: 0.9rem;">
                <strong>About CFSAN SNP Pipeline:</strong>
                <ul style="margin: 10px 0 0 20px;">
                    <li><strong>CFSAN</strong> is the CDC/FDA SNP Pipeline, commonly used for foodborne pathogen surveillance.</li>
                    <li><strong>Phase2 SNPs</strong>: Final set of high-quality SNPs after all quality filtering steps.</li>
                    <li><strong>Missing Positions</strong>: Positions excluded from the SNP matrix for this sample (usually due to low coverage or quality).</li>
                    <li>High SNP retention (&gt;99%) indicates that coreguard's gap regions don't significantly impact CFSAN's variant calling.</li>
                </ul>
            </div>
        `;

        // Add distance matrix if available
        if (cfsan.distance_matrix) {{
            html += renderDistanceMatrix(cfsan.distance_matrix, 'CFSAN', '#7b1fa2');
        }}

        container.innerHTML = html;
    }}

    // Back to top button visibility
    window.addEventListener('scroll', () => {{
        const btn = document.getElementById('back-to-top');
        if (window.scrollY > 300) {{
            btn.style.display = 'block';
        }} else {{
            btn.style.display = 'none';
        }}
    }});

    // Initialize (IGV loaded lazily when tab clicked)
    document.addEventListener('DOMContentLoaded', () => {{
        populateHome();
        setupQuickNav();
        populateSummary();
        populateSamples();
        populatePairwise();
        populateGeneZone();
        populateEmptyColumn();
        populateSnpCompare();
        populateVcfCompare();
        populateSnippyBugfix();
        populateSnpTracks();
        populateRunInfo();
        populateBedFiles();
    }});
    </script>
</body>
</html>
"##,
        ref_name = reference_name,
        report_json = report_json,
        reference_name = reference_name,
        tracks_js = tracks_js,
        half_len = reference_length / 2
    )
}

impl Sample {
    /// Line width for FASTA output (standard is 60 or 80)
    const FASTA_LINE_WIDTH: usize = 60;

    /// Convert sample to FASTA string with proper line wrapping
    #[allow(dead_code)]
    pub fn to_fasta_string(&self) -> String {
        self.to_fasta_string_with_name(&self.contigs.first().map(|c| c.name.clone()).unwrap_or_default())
    }

    /// Convert sample to FASTA string with custom contig name
    /// Used to match BED file chromosome names
    pub fn to_fasta_string_with_name(&self, name: &str) -> String {
        let mut fasta = String::new();
        // Use provided name for the first/main contig, concatenate all sequences
        fasta.push_str(&format!(">{}\n", name));
        for contig in &self.contigs {
            let seq_str = String::from_utf8_lossy(&contig.sequence);
            // Wrap sequence at FASTA_LINE_WIDTH characters per line
            for chunk in seq_str.as_bytes().chunks(Self::FASTA_LINE_WIDTH) {
                fasta.push_str(&String::from_utf8_lossy(chunk));
                fasta.push('\n');
            }
        }
        fasta
    }

    /// Generate FASTA index (.fai) string
    /// Format: name\tlength\toffset\tline_bases\tline_width
    #[allow(dead_code)]
    pub fn to_fai_string(&self) -> String {
        self.to_fai_string_with_name(&self.contigs.first().map(|c| c.name.clone()).unwrap_or_default())
    }

    /// Generate FASTA index with custom contig name
    pub fn to_fai_string_with_name(&self, name: &str) -> String {
        let line_bases = Self::FASTA_LINE_WIDTH;
        let line_width = line_bases + 1; // bases + newline

        // Calculate total sequence length (all contigs concatenated)
        let total_len: usize = self.contigs.iter().map(|c| c.sequence.len()).sum();

        // Header line length
        let header_line = format!(">{}\n", name);
        let offset = header_line.len();

        // Single FAI entry for the concatenated sequence
        format!(
            "{}\t{}\t{}\t{}\t{}\n",
            name, total_len, offset, line_bases, line_width
        )
    }
}

/// Generate full-power standalone IGV.js genome browser page
fn generate_full_igv_html(
    report_json: &str,
    reference_name: &str,
    reference_length: usize,
    bed_files: &HashMap<String, String>,
    bam_files: &HashMap<String, std::path::PathBuf>,
) -> String {
    // Generate track configurations for each BAM file (alignment tracks)
    let mut tracks_js = String::new();

    // BAM tracks first (alignment/coverage)
    let colors = ["#e94560", "#00d2d3", "#ff9f43", "#10ac84", "#5f27cd", "#ee5a24"];
    for (idx, (sample_name, _path)) in bam_files.iter().enumerate() {
        let color = colors[idx % colors.len()];
        tracks_js.push_str(&format!(
            r#"{{
                name: "{} (Alignment)",
                url: "/bam/{}.bam",
                indexURL: "/bam/{}.bam.bai",
                format: "bam",
                type: "alignment",
                height: 200,
                colorBy: "strand",
                showCoverage: true,
                showAlignments: true,
                viewAsPairs: true,
                color: "{}",
                coverageColor: "{}",
                visibilityWindow: 100000
            }},"#,
            sample_name, sample_name, sample_name, color, color
        ));
    }

    // BED tracks (gap regions)
    for (name, content) in bed_files {
        // Determine color based on name and extract risk from track description
        let (color, display_name): (&str, String) = if name.contains("CORE") {
            ("rgb(255, 215, 0)", "CORE GENOME".to_string())
        } else if content.contains("CRITICAL") {
            ("rgb(255, 0, 0)", format!("{} [CRITICAL]", name.replace("_gaps.bed", "")))
        } else if content.contains("HIGH") {
            ("rgb(255, 165, 0)", format!("{} [HIGH]", name.replace("_gaps.bed", "")))
        } else if content.contains("MEDIUM") {
            ("rgb(255, 255, 0)", format!("{} [MEDIUM]", name.replace("_gaps.bed", "")))
        } else {
            ("rgb(100, 200, 100)", name.replace("_gaps.bed", ""))
        };

        tracks_js.push_str(&format!(
            r#"{{
                name: "{}",
                url: "/bed/{}",
                format: "bed",
                color: "{}",
                displayMode: "EXPANDED",
                height: 50,
                visibilityWindow: -1,
                showLabels: true,
                labelField: "name"
            }},"#,
            display_name, name, color
        ));
    }

    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CoreGuard - IGV.js Full Browser</title>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e;
            color: #eee;
            height: 100vh;
            display: flex;
            flex-direction: column;
        }}
        .header {{
            background: linear-gradient(135deg, #16213e, #0f3460);
            padding: 12px 20px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            border-bottom: 2px solid #e94560;
        }}
        .header h1 {{
            font-size: 1.4rem;
            color: #e94560;
        }}
        .header .info {{
            font-size: 0.85rem;
            color: #aaa;
        }}
        .toolbar {{
            background: #16213e;
            padding: 10px 20px;
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
            border-bottom: 1px solid #333;
        }}
        .toolbar label {{
            color: #aaa;
            font-size: 0.85rem;
        }}
        .toolbar input, .toolbar select, .toolbar button {{
            padding: 8px 12px;
            border: 1px solid #444;
            border-radius: 4px;
            background: #0f3460;
            color: #eee;
            font-size: 0.9rem;
        }}
        .toolbar input:focus, .toolbar select:focus {{
            outline: none;
            border-color: #e94560;
        }}
        .toolbar button {{
            background: #e94560;
            border: none;
            cursor: pointer;
            font-weight: 600;
            transition: background 0.2s;
        }}
        .toolbar button:hover {{
            background: #ff6b6b;
        }}
        .toolbar button.secondary {{
            background: #333;
        }}
        .toolbar button.secondary:hover {{
            background: #444;
        }}
        .quick-nav {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        .quick-nav button {{
            background: #333;
            font-size: 0.8rem;
            padding: 5px 10px;
        }}
        .quick-nav button:hover {{
            background: #e94560;
        }}
        #igv-container {{
            flex: 1;
            overflow: hidden;
            background: white;
        }}
        .track-controls {{
            background: #16213e;
            padding: 10px 20px;
            border-top: 1px solid #333;
        }}
        .track-controls h4 {{
            color: #e94560;
            margin-bottom: 10px;
            font-size: 0.9rem;
        }}
        .track-list {{
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
        }}
        .track-item {{
            display: flex;
            align-items: center;
            gap: 6px;
            padding: 6px 12px;
            background: #0f3460;
            border-radius: 4px;
            font-size: 0.85rem;
        }}
        .track-item input[type="checkbox"] {{
            width: 16px;
            height: 16px;
        }}
        .stats-bar {{
            background: #0f3460;
            padding: 8px 20px;
            display: flex;
            gap: 20px;
            font-size: 0.85rem;
            border-top: 1px solid #333;
        }}
        .stat {{
            display: flex;
            gap: 8px;
        }}
        .stat-label {{
            color: #888;
        }}
        .stat-value {{
            color: #e94560;
            font-weight: 600;
        }}
        .modal {{
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.8);
            z-index: 1000;
            justify-content: center;
            align-items: center;
        }}
        .modal.active {{
            display: flex;
        }}
        .modal-content {{
            background: #16213e;
            padding: 30px;
            border-radius: 8px;
            max-width: 600px;
            width: 90%;
            max-height: 80vh;
            overflow-y: auto;
        }}
        .modal-content h3 {{
            color: #e94560;
            margin-bottom: 20px;
        }}
        .modal-content table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .modal-content th, .modal-content td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #333;
        }}
        .modal-content th {{
            color: #e94560;
        }}
        .close-modal {{
            position: absolute;
            top: 15px;
            right: 20px;
            font-size: 1.5rem;
            cursor: pointer;
            color: #888;
        }}
        .close-modal:hover {{
            color: #e94560;
        }}
        .bookmarks {{
            display: flex;
            gap: 5px;
            flex-wrap: wrap;
        }}
        .bookmark {{
            background: #333;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 0.8rem;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .bookmark:hover {{
            background: #e94560;
        }}
        .bookmark .remove {{
            font-size: 0.7rem;
            opacity: 0.7;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>IGV.js Genome Browser</h1>
        <div class="info">
            Reference: <strong>{reference_name}</strong> |
            Length: <strong>{ref_len_fmt}</strong> bp |
            <a href="/" style="color: #e94560; text-decoration: none;">← Back to Dashboard</a>
        </div>
    </div>

    <div class="toolbar">
        <div>
            <label>Search: </label>
            <input type="text" id="locus-input" placeholder="{reference_name}:1-10000 or gene name" style="width: 280px;">
            <button onclick="goToLocus()">Go</button>
        </div>
        <div>
            <label>Zoom: </label>
            <button class="secondary" onclick="zoomOut()">−</button>
            <button class="secondary" onclick="zoomIn()">+</button>
            <button class="secondary" onclick="zoomToAll()">Fit All</button>
        </div>
        <div class="quick-nav">
            <button onclick="goToRegion(1, 50000)">Start</button>
            <button onclick="goToRegion({ref_mid}, {ref_mid_end})">Middle</button>
            <button onclick="goToRegion({ref_end_start}, {reference_length})">End</button>
            <button onclick="showBookmarks()">Bookmarks</button>
            <button onclick="addBookmark()">+ Bookmark</button>
        </div>
    </div>

    <div id="igv-container"></div>

    <div class="stats-bar">
        <div class="stat">
            <span class="stat-label">Current Region:</span>
            <span class="stat-value" id="current-region">-</span>
        </div>
        <div class="stat">
            <span class="stat-label">Region Size:</span>
            <span class="stat-value" id="region-size">-</span>
        </div>
        <div class="stat">
            <span class="stat-label">Tracks:</span>
            <span class="stat-value" id="track-count">-</span>
        </div>
    </div>

    <div class="track-controls">
        <h4>Track Visibility</h4>
        <div class="track-list" id="track-list"></div>
    </div>

    <!-- Bookmarks Modal -->
    <div class="modal" id="bookmarks-modal">
        <div class="modal-content" style="position: relative;">
            <span class="close-modal" onclick="closeModal()">&times;</span>
            <h3>Saved Bookmarks</h3>
            <div class="bookmarks" id="bookmarks-list"></div>
            <p style="color: #888; margin-top: 15px; font-size: 0.85rem;">
                Click a bookmark to navigate to that region. Bookmarks are stored in your browser.
            </p>
        </div>
    </div>

    <script>
        const reportData = {report_json};
        const refName = '{reference_name}';
        const refLength = {reference_length};

        let browser = null;
        let bookmarks = JSON.parse(localStorage.getItem('coreguard_bookmarks') || '[]');

        // Initialize IGV
        document.addEventListener('DOMContentLoaded', async () => {{
            // Suppress IGV performance alerts by overriding alert for IGV messages
            const originalAlert = window.alert;
            window.alert = function(msg) {{
                if (msg && (msg.includes('performance') || msg.includes('zoom') || msg.includes('visibility'))) {{
                    console.log('IGV alert suppressed:', msg);
                    return;
                }}
                originalAlert.call(window, msg);
            }};

            // Suppress performance alerts
            igv.setGoogleOauthToken && igv.setGoogleOauthToken(null);

            const options = {{
                genome: {{
                    name: refName,
                    fastaURL: '/reference.fasta',
                    indexURL: '/reference.fasta.fai'
                }},
                locus: refName + ':1-50000',
                showNavigation: true,
                showRuler: true,
                showCenterGuide: true,
                showCursorGuide: true,
                showSequenceTrack: true,
                showSampleNames: true,
                showControls: true,
                // Disable performance warnings
                apiKey: null,
                queryParametersSupported: false,
                tracks: [
                    {{
                        name: "Reference Sequence",
                        type: "sequence",
                        order: -Number.MAX_SAFE_INTEGER,
                        height: 50,
                        removable: false
                    }},
                    {tracks_js}
                ]
            }};

            try {{
                browser = await igv.createBrowser(document.getElementById('igv-container'), options);

                // Update stats on location change
                browser.on('locuschange', (referenceFrameList) => {{
                    if (referenceFrameList && referenceFrameList.length > 0) {{
                        const frame = referenceFrameList[0];
                        const start = Math.round(frame.start);
                        const end = Math.round(frame.end);
                        const size = end - start;

                        document.getElementById('current-region').textContent =
                            `${{refName}}:${{start.toLocaleString()}}-${{end.toLocaleString()}}`;
                        document.getElementById('region-size').textContent =
                            size > 1000000 ? `${{(size/1000000).toFixed(2)}} Mb` :
                            size > 1000 ? `${{(size/1000).toFixed(2)}} kb` :
                            `${{size}} bp`;
                    }}
                }});

                // Build track controls
                buildTrackControls();

            }} catch (error) {{
                console.error('IGV initialization error:', error);
                document.getElementById('igv-container').innerHTML =
                    '<div style="padding: 50px; color: #ff6b6b; text-align: center;">' +
                    '<h3>Error initializing IGV.js</h3><p>' + error.message + '</p></div>';
            }}
        }});

        function goToLocus() {{
            const input = document.getElementById('locus-input').value.trim();
            if (input && browser) {{
                browser.search(input);
            }}
        }}

        function goToRegion(start, end) {{
            if (browser) {{
                browser.search(`${{refName}}:${{start}}-${{end}}`);
            }}
        }}

        function zoomIn() {{
            if (browser) browser.zoomIn();
        }}

        function zoomOut() {{
            if (browser) browser.zoomOut();
        }}

        function zoomToAll() {{
            if (browser) {{
                browser.search(`${{refName}}:1-${{refLength}}`);
            }}
        }}

        function buildTrackControls() {{
            if (!browser) return;

            const trackList = document.getElementById('track-list');
            trackList.innerHTML = '';

            const tracks = browser.trackViews.map(tv => tv.track).filter(t => t.name);
            document.getElementById('track-count').textContent = tracks.length;

            tracks.forEach((track, idx) => {{
                const item = document.createElement('div');
                item.className = 'track-item';
                item.innerHTML = `
                    <input type="checkbox" id="track-${{idx}}" checked
                           onchange="toggleTrack(${{idx}}, this.checked)">
                    <label for="track-${{idx}}">${{track.name}}</label>
                `;
                trackList.appendChild(item);
            }});
        }}

        function toggleTrack(idx, visible) {{
            if (!browser) return;
            const track = browser.trackViews[idx]?.track;
            if (track) {{
                track.visible = visible;
                browser.updateViews();
            }}
        }}

        function addBookmark() {{
            if (!browser) return;

            const frames = browser.referenceFrameList;
            if (frames && frames.length > 0) {{
                const frame = frames[0];
                const start = Math.round(frame.start);
                const end = Math.round(frame.end);
                const name = prompt('Bookmark name:', `Region ${{start}}-${{end}}`);

                if (name) {{
                    bookmarks.push({{
                        name: name,
                        locus: `${{refName}}:${{start}}-${{end}}`
                    }});
                    localStorage.setItem('coreguard_bookmarks', JSON.stringify(bookmarks));
                    alert('Bookmark saved!');
                }}
            }}
        }}

        function showBookmarks() {{
            const list = document.getElementById('bookmarks-list');
            list.innerHTML = '';

            if (bookmarks.length === 0) {{
                list.innerHTML = '<p style="color: #888;">No bookmarks saved yet.</p>';
            }} else {{
                bookmarks.forEach((b, idx) => {{
                    const el = document.createElement('div');
                    el.className = 'bookmark';
                    el.innerHTML = `
                        <span onclick="goToBookmark(${{idx}})">${{b.name}}</span>
                        <span class="remove" onclick="removeBookmark(${{idx}})">&times;</span>
                    `;
                    list.appendChild(el);
                }});
            }}

            document.getElementById('bookmarks-modal').classList.add('active');
        }}

        function goToBookmark(idx) {{
            const bookmark = bookmarks[idx];
            if (bookmark && browser) {{
                browser.search(bookmark.locus);
                closeModal();
            }}
        }}

        function removeBookmark(idx) {{
            bookmarks.splice(idx, 1);
            localStorage.setItem('coreguard_bookmarks', JSON.stringify(bookmarks));
            showBookmarks();
        }}

        function closeModal() {{
            document.querySelectorAll('.modal').forEach(m => m.classList.remove('active'));
        }}

        // Keyboard shortcuts
        document.addEventListener('keydown', (e) => {{
            if (e.target.tagName === 'INPUT') return;

            switch (e.key) {{
                case '+':
                case '=':
                    zoomIn();
                    break;
                case '-':
                    zoomOut();
                    break;
                case 'Home':
                    goToRegion(1, 50000);
                    break;
                case 'End':
                    goToRegion({ref_end_start}, {reference_length});
                    break;
                case 'b':
                    showBookmarks();
                    break;
                case 'Escape':
                    closeModal();
                    break;
            }}
        }});

        // Enter key for search
        document.getElementById('locus-input').addEventListener('keypress', (e) => {{
            if (e.key === 'Enter') goToLocus();
        }});
    </script>
</body>
</html>"##,
        reference_name = reference_name,
        ref_len_fmt = reference_length.to_string().chars()
            .rev()
            .collect::<Vec<_>>()
            .chunks(3)
            .map(|c| c.iter().collect::<String>())
            .collect::<Vec<_>>()
            .join(",")
            .chars()
            .rev()
            .collect::<String>(),
        ref_mid = reference_length / 2 - 25000,
        ref_mid_end = reference_length / 2 + 25000,
        ref_end_start = reference_length.saturating_sub(50000),
        report_json = report_json,
        tracks_js = tracks_js,
    )
}

/// Generate 3-way comparison visualization page
fn generate_compare_html(report_json: &str, reference_name: &str, reference_length: usize) -> String {
    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CoreGuard - Pipeline Comparison</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #0a0a1a;
            color: #e0e0e0;
            min-height: 100vh;
        }}
        .header {{
            background: linear-gradient(135deg, #1a1a2e, #16213e);
            padding: 20px 30px;
            border-bottom: 2px solid #e94560;
        }}
        .header h1 {{ color: #e94560; font-size: 1.8rem; }}
        .header .subtitle {{ color: #888; margin-top: 5px; }}
        .container {{ padding: 30px; max-width: 1600px; margin: 0 auto; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 25px; }}
        .card {{
            background: #1a1a2e;
            border-radius: 12px;
            padding: 25px;
            border: 1px solid #333;
        }}
        .card h3 {{
            color: #e94560;
            margin-bottom: 20px;
            font-size: 1.2rem;
            border-bottom: 1px solid #333;
            padding-bottom: 10px;
        }}
        .stat-row {{
            display: flex;
            justify-content: space-between;
            padding: 12px 0;
            border-bottom: 1px solid #222;
        }}
        .stat-row:last-child {{ border-bottom: none; }}
        .stat-label {{ color: #888; }}
        .stat-value {{ font-weight: bold; font-size: 1.1rem; }}
        .stat-value.good {{ color: #4ade80; }}
        .stat-value.bad {{ color: #f87171; }}
        .stat-value.warn {{ color: #fbbf24; }}
        .stat-value.info {{ color: #60a5fa; }}
        .venn-container {{
            display: flex;
            justify-content: center;
            align-items: center;
            height: 300px;
            position: relative;
        }}
        .venn-circle {{
            position: absolute;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: bold;
            font-size: 1.2rem;
        }}
        .genome-track {{
            height: 30px;
            background: #222;
            border-radius: 4px;
            margin: 10px 0;
            position: relative;
            overflow: hidden;
        }}
        .genome-track .label {{
            position: absolute;
            left: 10px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 0.8rem;
            color: #888;
            z-index: 10;
        }}
        .track-region {{
            position: absolute;
            height: 100%;
        }}
        .legend {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            margin-top: 15px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 0.85rem;
        }}
        .legend-color {{
            width: 16px;
            height: 16px;
            border-radius: 3px;
        }}
        .matrix {{
            display: grid;
            gap: 2px;
            font-size: 0.85rem;
        }}
        .matrix-cell {{
            padding: 8px;
            text-align: center;
            background: #222;
            border-radius: 4px;
        }}
        .matrix-header {{
            background: #333;
            font-weight: bold;
        }}
        .interpretation {{
            background: #1e3a5f;
            padding: 20px;
            border-radius: 8px;
            margin-top: 20px;
            border-left: 4px solid #60a5fa;
        }}
        .interpretation h4 {{
            color: #60a5fa;
            margin-bottom: 10px;
        }}
        a {{ color: #e94560; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Pipeline Comparison: mm2 vs Snippy vs CFSAN</h1>
        <div class="subtitle">
            Reference: <strong>{reference_name}</strong> ({ref_len_fmt} bp) |
            <a href="/">Dashboard</a> |
            <a href="/igv">IGV Browser</a>
        </div>
    </div>

    <div class="container">
        <div class="grid">
            <!-- VCF Concordance Summary -->
            <div class="card">
                <h3>SNP Concordance (Snippy vs CFSAN)</h3>
                <div id="vcf-summary">Loading...</div>
            </div>

            <!-- BAM Validation Results -->
            <div class="card">
                <h3>BAM Validation (Artifact Detection)</h3>
                <div id="bam-validation">Loading...</div>
            </div>

            <!-- Venn Diagram -->
            <div class="card">
                <h3>SNP Distribution</h3>
                <div class="venn-container" id="venn-diagram">
                    <div class="venn-circle" style="width: 200px; height: 200px; background: rgba(59, 130, 246, 0.3); border: 2px solid #3b82f6; left: 80px;" id="venn-snippy">
                        <span style="position: absolute; left: 20px;">Snippy</span>
                    </div>
                    <div class="venn-circle" style="width: 200px; height: 200px; background: rgba(16, 185, 129, 0.3); border: 2px solid #10b981; left: 180px;" id="venn-cfsan">
                        <span style="position: absolute; right: 20px;">CFSAN</span>
                    </div>
                    <div style="position: absolute; z-index: 20; text-align: center;" id="venn-overlap"></div>
                </div>
                <div id="venn-legend" class="legend"></div>
            </div>

            <!-- Corrected Distance Matrix -->
            <div class="card">
                <h3>Corrected SNP Distance Matrix</h3>
                <p style="color: #888; margin-bottom: 15px; font-size: 0.9rem;">
                    After removing artifacts, real genetic distances between samples
                </p>
                <div id="distance-matrix">Loading...</div>
            </div>

            <!-- Gap Comparison -->
            <div class="card" style="grid-column: span 2;">
                <h3>Genome-wide Gap & SNP Distribution</h3>
                <p style="color: #888; margin-bottom: 15px; font-size: 0.9rem;">
                    Showing where gaps (mm2) and SNP calls fall across the genome
                </p>
                <div id="genome-tracks">Loading...</div>
                <div class="legend">
                    <div class="legend-item"><div class="legend-color" style="background: #e94560;"></div>mm2 Gaps (no coverage)</div>
                    <div class="legend-item"><div class="legend-color" style="background: #3b82f6;"></div>Snippy-only SNPs</div>
                    <div class="legend-item"><div class="legend-color" style="background: #10b981;"></div>CFSAN-only SNPs</div>
                    <div class="legend-item"><div class="legend-color" style="background: #fbbf24;"></div>Concordant SNPs</div>
                    <div class="legend-item"><div class="legend-color" style="background: #8b5cf6;"></div>Artifacts (detected by BAM)</div>
                </div>
            </div>
        </div>

        <!-- Interpretation -->
        <div class="interpretation" id="interpretation">
            <h4>Interpretation</h4>
            <p id="interpretation-text">Loading analysis...</p>
        </div>
    </div>

    <script>
        const data = {report_json};
        const refLength = {reference_length};
        const refName = '{reference_name}';

        document.addEventListener('DOMContentLoaded', () => {{
            renderVcfSummary();
            renderBamValidation();
            renderVennDiagram();
            renderDistanceMatrix();
            renderGenomeTracks();
            renderInterpretation();
        }});

        function renderVcfSummary() {{
            const vcf = data.vcf_comparison;
            if (!vcf || !vcf.summary) {{
                document.getElementById('vcf-summary').innerHTML = '<p style="color: #888;">No VCF comparison data available</p>';
                return;
            }}
            const s = vcf.summary;
            const html = `
                <div class="stat-row">
                    <span class="stat-label">Total Positions</span>
                    <span class="stat-value info">${{s.total_positions?.toLocaleString() || 'N/A'}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Concordant (both agree)</span>
                    <span class="stat-value good">${{s.concordant?.toLocaleString() || 0}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Snippy-only</span>
                    <span class="stat-value warn">${{s.only_pipeline_a?.toLocaleString() || 0}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">CFSAN-only</span>
                    <span class="stat-value warn">${{s.only_pipeline_b?.toLocaleString() || 0}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Concordance Rate</span>
                    <span class="stat-value ${{s.concordance_rate > 80 ? 'good' : s.concordance_rate > 50 ? 'warn' : 'bad'}}">${{(s.concordance_rate || 0).toFixed(1)}}%</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Snippy-only in gaps</span>
                    <span class="stat-value">${{s.a_only_in_gaps?.toLocaleString() || 0}}</span>
                </div>
            `;
            document.getElementById('vcf-summary').innerHTML = html;
        }}

        function renderBamValidation() {{
            const bam = data.vcf_comparison?.summary?.bam_validation;
            if (!bam) {{
                document.getElementById('bam-validation').innerHTML = '<p style="color: #888;">No BAM validation data</p>';
                return;
            }}
            const html = `
                <div class="stat-row">
                    <span class="stat-label">Positions Checked</span>
                    <span class="stat-value info">${{bam.total_checked?.toLocaleString()}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Real Variants</span>
                    <span class="stat-value good">${{bam.real_variants}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Artifacts Detected</span>
                    <span class="stat-value bad">${{bam.artifacts?.toLocaleString()}}</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Artifact Rate</span>
                    <span class="stat-value bad">${{(bam.artifact_rate || 0).toFixed(2)}}%</span>
                </div>
                <div class="stat-row">
                    <span class="stat-label">Corrected Snippy-only</span>
                    <span class="stat-value">${{bam.original_snippy_only?.toLocaleString()}} → <strong style="color: #4ade80;">${{bam.corrected_snippy_only}}</strong></span>
                </div>
            `;
            document.getElementById('bam-validation').innerHTML = html;
        }}

        function renderVennDiagram() {{
            const vcf = data.vcf_comparison?.summary;
            if (!vcf) return;

            const snippyOnly = vcf.only_pipeline_a || 0;
            const cfsanOnly = vcf.only_pipeline_b || 0;
            const concordant = vcf.concordant || 0;

            document.getElementById('venn-overlap').innerHTML = `
                <div style="font-size: 1.5rem; color: #fbbf24;">${{concordant.toLocaleString()}}</div>
                <div style="font-size: 0.8rem; color: #888;">concordant</div>
            `;

            document.getElementById('venn-legend').innerHTML = `
                <div class="legend-item"><div class="legend-color" style="background: #3b82f6;"></div>Snippy-only: ${{snippyOnly.toLocaleString()}}</div>
                <div class="legend-item"><div class="legend-color" style="background: #10b981;"></div>CFSAN-only: ${{cfsanOnly.toLocaleString()}}</div>
                <div class="legend-item"><div class="legend-color" style="background: #fbbf24;"></div>Both: ${{concordant.toLocaleString()}}</div>
            `;
        }}

        function renderDistanceMatrix() {{
            const bam = data.vcf_comparison?.summary?.bam_validation;
            if (!bam || !bam.corrected_distance_matrix) {{
                document.getElementById('distance-matrix').innerHTML = '<p style="color: #888;">No distance matrix available</p>';
                return;
            }}
            const dm = bam.corrected_distance_matrix;
            const samples = dm.samples || [];
            const matrix = dm.matrix || [];

            let html = `<div class="matrix" style="grid-template-columns: repeat(${{samples.length + 1}}, 1fr);">`;
            html += `<div class="matrix-cell matrix-header"></div>`;
            for (const s of samples) {{
                html += `<div class="matrix-cell matrix-header">${{s}}</div>`;
            }}
            for (let i = 0; i < samples.length; i++) {{
                html += `<div class="matrix-cell matrix-header">${{samples[i]}}</div>`;
                for (let j = 0; j < samples.length; j++) {{
                    const val = matrix[i][j];
                    const color = i === j ? '#333' : val === 0 ? '#065f46' : val <= 2 ? '#166534' : val <= 5 ? '#ca8a04' : '#dc2626';
                    html += `<div class="matrix-cell" style="background: ${{color}};">${{val}}</div>`;
                }}
            }}
            html += `</div>`;
            document.getElementById('distance-matrix').innerHTML = html;
        }}

        function renderGenomeTracks() {{
            // Use impact_analysis for coverage data
            const impactSamples = data.impact_analysis?.samples || [];
            const refSummary = data.reference_summary || {{}};
            const avgCoverage = (refSummary.avg_sample_coverage || 0.85) * 100;

            let html = `
                <div style="margin-bottom: 15px; padding: 10px; background: #222; border-radius: 6px;">
                    <div style="display: flex; justify-content: space-between; margin-bottom: 10px;">
                        <span>Average coverage: <strong style="color: #4ade80;">${{avgCoverage.toFixed(1)}}%</strong></span>
                        <span>Reference: <strong>${{refLength.toLocaleString()}} bp</strong></span>
                    </div>
                </div>
            `;

            // Create sample tracks
            for (const sample of impactSamples) {{
                const name = sample.sample || sample.name || 'Unknown';
                const sharedGaps = sample.shared_gap_bases || 0;
                const uniqueGaps = sample.unique_gap_bases || 0;
                const totalGaps = sharedGaps + uniqueGaps;
                const gapPercent = (totalGaps / refLength) * 100;
                const coveredPercent = 100 - gapPercent;
                const riskLevel = sample.risk_level || 'Unknown';
                const riskColor = riskLevel === 'Low' ? '#4ade80' : riskLevel === 'Medium' ? '#fbbf24' : '#f87171';

                html += `
                    <div class="genome-track" style="height: 40px; margin: 8px 0;">
                        <span class="label" style="width: 120px;">${{name}}</span>
                        <div style="position: absolute; left: 130px; right: 10px; top: 5px; bottom: 5px; background: #333; border-radius: 4px; overflow: hidden;">
                            <div style="position: absolute; left: 0; top: 0; bottom: 0; width: ${{coveredPercent}}%; background: linear-gradient(90deg, #166534, #22c55e);"></div>
                            <div style="position: absolute; right: 0; top: 0; bottom: 0; width: ${{gapPercent}}%; background: #e94560;"></div>
                        </div>
                        <span style="position: absolute; right: 15px; top: 50%; transform: translateY(-50%); font-size: 0.75rem; color: ${{riskColor}}; background: rgba(0,0,0,0.5); padding: 2px 6px; border-radius: 3px;">
                            ${{coveredPercent.toFixed(1)}}% | ${{riskLevel}}
                        </span>
                    </div>
                `;
            }}

            // Add VCF comparison summary bar
            const vcf = data.vcf_comparison?.summary;
            if (vcf) {{
                const concordant = vcf.concordant || 0;
                const snippyOnly = vcf.only_pipeline_a || 0;
                const cfsanOnly = vcf.only_pipeline_b || 0;
                const total = concordant + snippyOnly + cfsanOnly;

                html += `
                    <div style="margin-top: 20px; padding-top: 15px; border-top: 1px solid #333;">
                        <div style="font-size: 0.9rem; color: #888; margin-bottom: 8px;">SNP Distribution (total: ${{total.toLocaleString()}} positions)</div>
                        <div style="height: 30px; background: #222; border-radius: 4px; overflow: hidden; display: flex;">
                            <div style="width: ${{(concordant/total)*100}}%; background: #fbbf24;" title="Concordant: ${{concordant.toLocaleString()}}"></div>
                            <div style="width: ${{(snippyOnly/total)*100}}%; background: #3b82f6;" title="Snippy-only: ${{snippyOnly.toLocaleString()}}"></div>
                            <div style="width: ${{(cfsanOnly/total)*100}}%; background: #10b981;" title="CFSAN-only: ${{cfsanOnly.toLocaleString()}}"></div>
                        </div>
                        <div style="display: flex; justify-content: space-between; margin-top: 5px; font-size: 0.8rem;">
                            <span style="color: #fbbf24;">Concordant: ${{(concordant/total*100).toFixed(1)}}%</span>
                            <span style="color: #3b82f6;">Snippy: ${{(snippyOnly/total*100).toFixed(1)}}%</span>
                            <span style="color: #10b981;">CFSAN: ${{(cfsanOnly/total*100).toFixed(1)}}%</span>
                        </div>
                    </div>
                `;
            }}

            document.getElementById('genome-tracks').innerHTML = html || '<p style="color: #888;">No sample data</p>';
        }}

        function renderInterpretation() {{
            const vcf = data.vcf_comparison?.summary;
            const bam = vcf?.bam_validation;

            let text = '';
            if (vcf) {{
                text += vcf.interpretation || '';
                if (bam && bam.artifact_rate > 90) {{
                    text += ` <strong>Critical finding:</strong> ${{bam.artifact_rate.toFixed(1)}}% of Snippy-unique SNPs are artifacts. `;
                    text += `After BAM validation, only ${{bam.corrected_snippy_only}} real Snippy-unique variants remain. `;
                    text += `The corrected distance matrix shows the true genetic relationships between samples.`;
                }}
            }} else {{
                text = 'No VCF comparison data available. Run with --vcf-snippy and --vcf-cfsan --compare-variants flags.';
            }}

            document.getElementById('interpretation-text').innerHTML = text;
        }}
    </script>
</body>
</html>"##,
        reference_name = reference_name,
        ref_len_fmt = reference_length.to_string().chars()
            .rev()
            .collect::<Vec<_>>()
            .chunks(3)
            .map(|c| c.iter().collect::<String>())
            .collect::<Vec<_>>()
            .join(",")
            .chars()
            .rev()
            .collect::<String>(),
        report_json = report_json,
    )
}

/// Generate 3D visualization page using Three.js
fn generate_3d_html(report_json: &str, reference_name: &str, reference_length: usize) -> String {
    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CoreGuard 3D Genome Visualization</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #0a0a0f;
            color: #e0e0e0;
            overflow: hidden;
        }}
        #container {{
            width: 100vw;
            height: 100vh;
            position: relative;
        }}
        #controls {{
            position: absolute;
            top: 20px;
            left: 20px;
            background: rgba(20, 20, 30, 0.95);
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.5);
            z-index: 100;
            max-width: 320px;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }}
        #controls h2 {{
            font-size: 1.2rem;
            margin-bottom: 15px;
            color: #fff;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        #controls h2::before {{
            content: '';
            display: inline-block;
            width: 8px;
            height: 8px;
            background: #4CAF50;
            border-radius: 50%;
            animation: pulse 2s infinite;
        }}
        @keyframes pulse {{
            0%, 100% {{ opacity: 1; }}
            50% {{ opacity: 0.5; }}
        }}
        .control-group {{
            margin-bottom: 15px;
        }}
        .control-group label {{
            display: block;
            font-size: 0.85rem;
            color: #aaa;
            margin-bottom: 6px;
        }}
        .control-group input[type="range"] {{
            width: 100%;
            height: 6px;
            background: #333;
            border-radius: 3px;
            appearance: none;
            outline: none;
        }}
        .control-group input[type="range"]::-webkit-slider-thumb {{
            appearance: none;
            width: 16px;
            height: 16px;
            background: #4CAF50;
            border-radius: 50%;
            cursor: pointer;
        }}
        .control-group select, .control-group input[type="number"] {{
            width: 100%;
            padding: 8px 10px;
            background: #1a1a24;
            border: 1px solid #333;
            border-radius: 6px;
            color: #fff;
            font-size: 0.9rem;
        }}
        #region-display {{
            font-size: 0.85rem;
            color: #888;
            margin-top: 5px;
            font-family: monospace;
        }}
        .btn {{
            padding: 10px 16px;
            background: linear-gradient(135deg, #4CAF50, #45a049);
            border: none;
            border-radius: 6px;
            color: #fff;
            font-size: 0.9rem;
            cursor: pointer;
            transition: all 0.2s;
            margin-right: 8px;
            margin-bottom: 8px;
        }}
        .btn:hover {{
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(76, 175, 80, 0.4);
        }}
        .btn.secondary {{
            background: linear-gradient(135deg, #555, #444);
        }}
        #legend {{
            position: absolute;
            bottom: 20px;
            left: 360px;
            background: rgba(20, 20, 30, 0.95);
            padding: 15px 20px;
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.5);
            z-index: 100;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }}
        #legend h3 {{
            font-size: 0.9rem;
            margin-bottom: 10px;
            color: #fff;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 6px 0;
            font-size: 0.85rem;
            padding: 4px 8px;
            border-radius: 4px;
            transition: all 0.2s;
        }}
        .legend-item.clickable {{
            cursor: pointer;
        }}
        .legend-item.clickable:hover {{
            background: rgba(255, 255, 255, 0.1);
        }}
        .legend-item.filtered {{
            opacity: 0.4;
            text-decoration: line-through;
        }}
        .legend-color {{
            width: 20px;
            height: 12px;
            margin-right: 10px;
            border-radius: 2px;
            flex-shrink: 0;
        }}
        .filter-status {{
            margin-left: auto;
            font-size: 0.9rem;
            color: #4CAF50;
        }}
        .filter-status.inactive {{
            color: #666;
        }}
        .legend-item.filtered .filter-status {{
            color: #666;
        }}
        #info-panel {{
            position: absolute;
            top: 20px;
            right: 20px;
            background: rgba(20, 20, 30, 0.95);
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.5);
            z-index: 100;
            min-width: 280px;
            display: none;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }}
        #info-panel h3 {{
            font-size: 1rem;
            margin-bottom: 12px;
            color: #fff;
        }}
        #info-panel .info-row {{
            display: flex;
            justify-content: space-between;
            padding: 6px 0;
            border-bottom: 1px solid rgba(255, 255, 255, 0.05);
            font-size: 0.85rem;
        }}
        #info-panel .info-label {{
            color: #888;
        }}
        #info-panel .info-value {{
            color: #fff;
            font-family: monospace;
        }}
        #info-panel .close-btn {{
            position: absolute;
            top: 10px;
            right: 12px;
            background: none;
            border: none;
            color: #666;
            font-size: 1.2rem;
            cursor: pointer;
        }}
        #loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            text-align: center;
            z-index: 200;
        }}
        #loading .spinner {{
            width: 50px;
            height: 50px;
            border: 3px solid #333;
            border-top-color: #4CAF50;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 15px;
        }}
        @keyframes spin {{
            to {{ transform: rotate(360deg); }}
        }}
        #stats {{
            position: absolute;
            bottom: 20px;
            right: 20px;
            background: rgba(20, 20, 30, 0.8);
            padding: 10px 15px;
            border-radius: 8px;
            font-size: 0.8rem;
            font-family: monospace;
            color: #888;
        }}
        .sample-label {{
            position: absolute;
            pointer-events: none;
            color: #fff;
            font-size: 12px;
            text-shadow: 0 1px 3px rgba(0,0,0,0.8);
            white-space: nowrap;
        }}
    </style>
</head>
<body>
    <div id="container">
        <canvas id="canvas3d"></canvas>
    </div>

    <div id="controls">
        <h2>CoreGuard 3D</h2>
        <div class="control-group">
            <label>Reference: {reference_name}</label>
            <div id="region-display">Length: {reference_length} bp</div>
        </div>
        <div class="control-group">
            <label>Region Start (bp)</label>
            <input type="number" id="region-start" value="1" min="1" max="{reference_length}">
        </div>
        <div class="control-group">
            <label>Region Size (bp)</label>
            <select id="region-size" onchange="loadRegion()">
                <option value="50">50 bp (nucleotides)</option>
                <option value="100">100 bp (nucleotides)</option>
                <option value="200">200 bp</option>
                <option value="500">500 bp</option>
                <option value="1000">1,000 bp</option>
                <option value="5000" selected>5,000 bp</option>
                <option value="10000">10,000 bp</option>
                <option value="50000">50,000 bp</option>
                <option value="100000">100,000 bp</option>
            </select>
        </div>
        <div class="control-group">
            <label>SNP Height Scale: <span id="snp-scale-value">2.0</span></label>
            <input type="range" id="snp-scale" min="0.5" max="5" step="0.1" value="2" oninput="updateSlider('snp-scale')">
        </div>
        <div class="control-group">
            <label>Track Spacing: <span id="track-spacing-value">2.0</span></label>
            <input type="range" id="track-spacing" min="1" max="5" step="0.5" value="2" oninput="updateSlider('track-spacing')">
        </div>
        <div style="margin-top: 15px;">
            <button class="btn" onclick="loadRegion()">Load Region</button>
            <button class="btn secondary" onclick="resetCamera()">Reset View</button>
        </div>
        <div style="margin-top: 10px;">
            <button class="btn secondary" onclick="toggleAutoRotate()">Toggle Rotation</button>
            <button class="btn secondary" onclick="jumpToGaps()">Find Gaps</button>
        </div>
        <div style="margin-top: 10px; display: flex; gap: 5px;">
            <button class="btn secondary" onclick="navigate(-1)" style="padding: 5px 10px;">&lt;&lt;</button>
            <button class="btn secondary" onclick="navigate(1)" style="padding: 5px 10px;">&gt;&gt;</button>
            <span style="font-size: 0.8rem; color: #888; align-self: center; margin-left: 5px;">Navigate</span>
        </div>
        <div style="margin-top: 15px;">
            <button class="btn" onclick="loadSnakeGenome()" style="background: linear-gradient(135deg, #9C27B0, #7B1FA2);">
                🍝 Full Genome (Snake)
            </button>
            <div style="font-size: 0.7rem; color: #666; margin-top: 3px;">Same view, whole genome as spaghetti</div>
        </div>
        <div id="genome-overview" style="margin-top: 15px; background: #1a1a24; border-radius: 6px; padding: 8px;">
            <div style="font-size: 0.75rem; color: #888; margin-bottom: 5px;">Genome Overview (click to jump)</div>
            <canvas id="overview-canvas" width="280" height="30" style="cursor: pointer; border-radius: 4px;"></canvas>
        </div>
    </div>

    <div id="legend">
        <h3>Legend <span style="font-size: 0.7rem; color: #888;">(click to filter)</span></h3>
        <div class="legend-item clickable" data-filter="reference" onclick="toggleFilter('reference')">
            <div class="legend-color" style="background: #4CAF50;"></div>Reference
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="snippy-only" onclick="toggleFilter('snippy-only')">
            <div class="legend-color" style="background: #7B1FA2;"></div>SNP Snippy-only
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="cfsan-only" onclick="toggleFilter('cfsan-only')">
            <div class="legend-color" style="background: #0288D1;"></div>SNP CFSAN-only
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="concordant" onclick="toggleFilter('concordant')">
            <div class="legend-color" style="background: #4CAF50;"></div>SNP Concordant
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="discordant" onclick="toggleFilter('discordant')">
            <div class="legend-color" style="background: #ff5722;"></div>SNP Discordant
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="gap-mm2" onclick="toggleFilter('gap-mm2')">
            <div class="legend-color" style="background: #9e9e9e;"></div>Gap (mm2)
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="gap-snippy" onclick="toggleFilter('gap-snippy')">
            <div class="legend-color" style="background: #e65100;"></div>Gap (Snippy)
            <span class="filter-status active">&#10003;</span>
        </div>
        <div class="legend-item clickable" data-filter="gap-cfsan" onclick="toggleFilter('gap-cfsan')">
            <div class="legend-color" style="background: #0d47a1;"></div>Gap (CFSAN)
            <span class="filter-status active">&#10003;</span>
        </div>
    </div>

    <div id="info-panel">
        <button class="close-btn" onclick="closeInfoPanel()">&times;</button>
        <h3 id="info-title">Details</h3>
        <div id="info-content"></div>
    </div>

    <div id="loading">
        <div class="spinner"></div>
        <div>Loading 3D visualization...</div>
    </div>

    <div id="stats">
        <span id="fps">FPS: --</span> |
        <span id="objects">Objects: --</span>
    </div>

    <!-- Genome Map Modal -->
    <div id="genome-map-modal" style="display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.95); z-index: 1000; overflow: auto;">
        <div style="padding: 20px; color: #fff;">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
                <h2 style="margin: 0;">🗺️ Full Genome Map - {reference_length} bp</h2>
                <button onclick="closeGenomeMap()" style="background: #f44336; border: none; color: #fff; padding: 10px 20px; border-radius: 6px; cursor: pointer; font-size: 1rem;">✕ Close</button>
            </div>
            <div style="display: flex; gap: 20px; margin-bottom: 15px; flex-wrap: wrap;">
                <div style="background: rgba(255,255,255,0.1); padding: 10px 15px; border-radius: 6px;">
                    <span style="color: #4CAF50;">■</span> Covered |
                    <span style="color: #9e9e9e;">■</span> Gap (mm2) |
                    <span style="color: #e65100;">■</span> Gap (Snippy) |
                    <span style="color: #0d47a1;">■</span> Gap (CFSAN)
                </div>
                <div style="background: rgba(255,255,255,0.1); padding: 10px 15px; border-radius: 6px;">
                    <span style="color: #7B1FA2;">●</span> SNP Snippy-only |
                    <span style="color: #0288D1;">●</span> SNP CFSAN-only |
                    <span style="color: #4CAF50;">●</span> Concordant |
                    <span style="color: #ff5722;">●</span> Discordant
                </div>
                <div style="background: rgba(255,255,255,0.1); padding: 10px 15px; border-radius: 6px;">
                    Zoom: <select id="map-zoom" onchange="renderGenomeMap()" style="background: #333; color: #fff; border: 1px solid #555; padding: 5px; border-radius: 4px;">
                        <option value="1">1 bp/px (full detail)</option>
                        <option value="10">10 bp/px</option>
                        <option value="100" selected>100 bp/px</option>
                        <option value="1000">1000 bp/px</option>
                    </select>
                </div>
            </div>
            <div id="map-status" style="margin-bottom: 10px; color: #888;">Loading genome data...</div>
            <div id="map-container" style="overflow: auto; max-height: calc(100vh - 200px); background: #0a0a10; border-radius: 8px; padding: 10px;">
                <canvas id="genome-map-canvas" style="display: block;"></canvas>
            </div>
            <div id="map-tooltip" style="position: fixed; display: none; background: rgba(0,0,0,0.9); color: #fff; padding: 10px 15px; border-radius: 6px; font-size: 0.9rem; pointer-events: none; z-index: 1001; border: 1px solid #444;"></div>
        </div>
    </div>

    <!-- Three.js from CDN -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>

    <script>
        // Global state
        let scene, camera, renderer, controls;
        let trackGroup, snpGroup, gapGroup, labelGroup;
        let samples = [];
        let referenceLength = {reference_length};
        let currentStart = 1;
        let currentEnd = 5001;
        let autoRotate = false;
        let frameCount = 0;
        let lastTime = performance.now();
        let lastLoadedData = null;  // Cache last loaded data for re-rendering

        // Filter state - all enabled by default
        const filters = {{
            'reference': true,
            'snippy-only': true,
            'cfsan-only': true,
            'concordant': true,
            'discordant': true,
            'gap-mm2': true,
            'gap-snippy': true,
            'gap-cfsan': true
        }};

        // Colors
        const COLORS = {{
            reference: 0x4CAF50,
            // SNP colors by source
            snpSnippyOnly: 0x7B1FA2,   // Purple - Snippy only
            snpCfsanOnly: 0x0288D1,    // Light blue - CFSAN only
            snpConcordant: 0x4CAF50,   // Green - Both pipelines agree
            snpDiscordant: 0xff5722,   // Orange/Red - Both but different
            // Gap colors
            gapMm2: 0x9e9e9e,
            gapSnippy: 0xe65100,
            gapCfsan: 0x0d47a1,
            // Track colors
            trackBase: 0x37474F,
            trackColors: [0x5C6BC0, 0x26A69A, 0xFFCA28, 0xEF5350, 0xAB47BC, 0x66BB6A, 0xFF7043, 0x29B6F6, 0xD4E157, 0x8D6E63],
            grid: 0x333340
        }};

        // Parse report data
        const reportData = {report_json};

        // Initialize Three.js scene
        function init() {{
            const container = document.getElementById('container');
            const canvas = document.getElementById('canvas3d');

            // Scene
            scene = new THREE.Scene();
            scene.background = new THREE.Color(0x0a0a0f);
            scene.fog = new THREE.Fog(0x0a0a0f, 100, 500);

            // Camera
            camera = new THREE.PerspectiveCamera(
                60,
                window.innerWidth / window.innerHeight,
                0.1,
                1000
            );
            camera.position.set(0, 30, 50);

            // Renderer
            renderer = new THREE.WebGLRenderer({{
                canvas: canvas,
                antialias: true,
                alpha: true
            }});
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
            renderer.shadowMap.enabled = true;
            renderer.shadowMap.type = THREE.PCFSoftShadowMap;

            // Controls
            controls = new THREE.OrbitControls(camera, renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;
            controls.minDistance = 10;
            controls.maxDistance = 200;
            controls.maxPolarAngle = Math.PI / 2;

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.4);
            scene.add(ambientLight);

            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(50, 100, 50);
            directionalLight.castShadow = true;
            directionalLight.shadow.mapSize.width = 2048;
            directionalLight.shadow.mapSize.height = 2048;
            scene.add(directionalLight);

            const pointLight = new THREE.PointLight(0x4CAF50, 0.5, 100);
            pointLight.position.set(0, 20, 0);
            scene.add(pointLight);

            // Groups for organization
            trackGroup = new THREE.Group();
            snpGroup = new THREE.Group();
            gapGroup = new THREE.Group();
            labelGroup = new THREE.Group();
            scene.add(trackGroup);
            scene.add(snpGroup);
            scene.add(gapGroup);
            scene.add(labelGroup);

            // Grid
            const gridHelper = new THREE.GridHelper(200, 50, COLORS.grid, COLORS.grid);
            gridHelper.position.y = -2;
            scene.add(gridHelper);

            // Extract samples from report
            if (reportData.pairwise_results && reportData.pairwise_results.length > 0) {{
                const sampleSet = new Set();
                reportData.pairwise_results.forEach(p => {{
                    sampleSet.add(p.sample_a_name);
                    sampleSet.add(p.sample_b_name);
                }});
                samples = Array.from(sampleSet).slice(0, 10); // Limit to 10 samples
            }} else if (reportData.samples) {{
                samples = reportData.samples.map(s => s.name).slice(0, 10);
            }}

            // Raycaster for mouse interaction
            setupInteraction();

            // Window resize
            window.addEventListener('resize', onWindowResize);

            // Hide loading
            document.getElementById('loading').style.display = 'none';

            // Initial load
            loadRegion();

            // Start animation loop
            animate();
        }}

        function createTracks(data) {{
            // Cache data for re-rendering
            lastLoadedData = data;

            // Clear existing
            while(trackGroup.children.length) trackGroup.remove(trackGroup.children[0]);
            while(snpGroup.children.length) snpGroup.remove(snpGroup.children[0]);
            while(gapGroup.children.length) gapGroup.remove(gapGroup.children[0]);

            const trackSpacing = parseFloat(document.getElementById('track-spacing').value);
            const snpScale = parseFloat(document.getElementById('snp-scale').value);
            const regionSize = data.length || (currentEnd - currentStart);

            // For nucleotide regions (< 500bp), make each nucleotide ~0.8 units wide for readable letters
            // For larger regions, use standard scaling
            const trackLength = regionSize <= 500
                ? regionSize * 0.8  // Wide enough for letters
                : Math.min(regionSize / 100, 100);
            const scaleFactor = trackLength / regionSize;

            // Reference track (bottom) - only if filter enabled
            if (filters['reference']) {{
                const refGeometry = new THREE.BoxGeometry(trackLength, 0.3, 1);
                const refMaterial = new THREE.MeshPhongMaterial({{
                    color: COLORS.reference,
                    shininess: 60
                }});
                const refTrack = new THREE.Mesh(refGeometry, refMaterial);
                refTrack.position.set(0, -1, 0);
                refTrack.userData = {{ type: 'reference', name: 'Reference' }};
                trackGroup.add(refTrack);

                // Show nucleotides if region is small enough (< 500bp)
                const refSeq = data.reference || '';
                if (refSeq.length > 0 && refSeq.length <= 500) {{
                    const ntColors = {{
                        'A': 0x4CAF50, // Green
                        'T': 0xf44336, // Red
                        'G': 0xFFEB3B, // Yellow
                        'C': 0x2196F3, // Blue
                        'N': 0x9e9e9e  // Gray
                    }};
                    const ntWidth = trackLength / refSeq.length;

                    for (let i = 0; i < refSeq.length; i++) {{
                        const nt = refSeq[i].toUpperCase();
                        const color = ntColors[nt] || 0x666666;
                        const xPos = (i * ntWidth) - trackLength/2 + ntWidth/2;

                        // Nucleotide box
                        const ntGeom = new THREE.BoxGeometry(ntWidth * 0.9, 0.25, 0.9);
                        const ntMat = new THREE.MeshPhongMaterial({{
                            color: color,
                            shininess: 80
                        }});
                        const ntMesh = new THREE.Mesh(ntGeom, ntMat);
                        ntMesh.position.set(xPos, -1.4, 0);
                        ntMesh.userData = {{ type: 'nucleotide', base: nt, position: currentStart + i }};
                        trackGroup.add(ntMesh);

                        // Add letter label for all nucleotide regions
                        const canvas = document.createElement('canvas');
                        const ctx = canvas.getContext('2d');
                        canvas.width = 64;
                        canvas.height = 64;
                        ctx.fillStyle = '#ffffff';
                        ctx.font = 'bold 48px monospace';
                        ctx.textAlign = 'center';
                        ctx.fillText(nt, 32, 48);
                        const texture = new THREE.CanvasTexture(canvas);
                        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture, transparent: true }}));
                        sprite.scale.set(ntWidth * 0.9, ntWidth * 0.9, 1);
                        sprite.position.set(xPos, -1.4, 0.6);
                        trackGroup.add(sprite);
                    }}
                    console.log(`Rendered ${{refSeq.length}} nucleotides`);

                    // Add position markers every 50bp
                    for (let i = 0; i < refSeq.length; i += 50) {{
                        const genomePos = currentStart + i;
                        const xPos = (i * ntWidth) - trackLength/2 + ntWidth/2;

                        // Tick mark
                        const tickGeom = new THREE.BoxGeometry(0.1, 0.4, 0.1);
                        const tickMat = new THREE.MeshBasicMaterial({{ color: 0xffffff }});
                        const tickMesh = new THREE.Mesh(tickGeom, tickMat);
                        tickMesh.position.set(xPos, -1.8, 0);
                        trackGroup.add(tickMesh);

                        // Position label
                        const posCanvas = document.createElement('canvas');
                        const posCtx = posCanvas.getContext('2d');
                        posCanvas.width = 128;
                        posCanvas.height = 48;
                        posCtx.fillStyle = '#ffffff';
                        posCtx.font = 'bold 20px Arial';
                        posCtx.textAlign = 'center';
                        posCtx.fillText(genomePos.toLocaleString(), 64, 32);
                        const posTex = new THREE.CanvasTexture(posCanvas);
                        const posSprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: posTex, transparent: true }}));
                        posSprite.scale.set(ntWidth * 6, ntWidth * 2, 1);
                        posSprite.position.set(xPos, -2.3, 0);
                        trackGroup.add(posSprite);
                    }}
                }}

                // Add start position label
                const startCanvas = document.createElement('canvas');
                const startCtx = startCanvas.getContext('2d');
                startCanvas.width = 256;
                startCanvas.height = 64;
                startCtx.fillStyle = '#4CAF50';
                startCtx.font = 'bold 24px Arial';
                startCtx.textAlign = 'right';
                startCtx.fillText(`${{currentStart.toLocaleString()}} bp`, 250, 40);
                const startTex = new THREE.CanvasTexture(startCanvas);
                const startSprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: startTex, transparent: true }}));
                startSprite.scale.set(8, 2, 1);
                startSprite.position.set(-trackLength/2 - 5, -2, 0);
                trackGroup.add(startSprite);

                // Add end position label
                const endCanvas = document.createElement('canvas');
                const endCtx = endCanvas.getContext('2d');
                endCanvas.width = 256;
                endCanvas.height = 64;
                endCtx.fillStyle = '#4CAF50';
                endCtx.font = 'bold 24px Arial';
                endCtx.textAlign = 'left';
                endCtx.fillText(`${{currentEnd.toLocaleString()}} bp`, 5, 40);
                const endTex = new THREE.CanvasTexture(endCanvas);
                const endSprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: endTex, transparent: true }}));
                endSprite.scale.set(8, 2, 1);
                endSprite.position.set(trackLength/2 + 5, -2, 0);
                trackGroup.add(endSprite);
            }}

            // Build lookup maps for SNPs and gaps per sample
            const snpMap = {{}};  // snpMap[sample][position] = alt allele
            const gapSet = {{}};  // gapSet[sample] = Set of positions in gaps

            data.samples?.forEach(sample => {{
                snpMap[sample] = {{}};
                gapSet[sample] = new Set();
            }});

            (data.snps || []).forEach(snp => {{
                if (snp.sample && snp.pos) {{
                    snpMap[snp.sample] = snpMap[snp.sample] || {{}};
                    snpMap[snp.sample][snp.pos] = snp.alt || '?';
                }}
            }});

            (data.gaps || []).forEach(gap => {{
                if (gap.sample && gap.start && gap.end) {{
                    gapSet[gap.sample] = gapSet[gap.sample] || new Set();
                    for (let p = gap.start; p <= gap.end; p++) {{
                        gapSet[gap.sample].add(p);
                    }}
                }}
            }});

            // Sample tracks with distinct colors and labels
            const refSeq = data.reference || '';
            const showNucleotides = refSeq.length > 0 && refSeq.length <= 500;

            data.samples.forEach((sampleName, idx) => {{
                const yPos = idx * trackSpacing + 1;
                const trackColor = COLORS.trackColors[idx % COLORS.trackColors.length];

                // Base track - more visible (thinner if showing nucleotides)
                const trackGeometry = new THREE.BoxGeometry(trackLength, showNucleotides ? 0.15 : 0.3, 1.0);
                const trackMaterial = new THREE.MeshPhongMaterial({{
                    color: trackColor,
                    transparent: true,
                    opacity: showNucleotides ? 0.5 : 0.85,
                    shininess: 30
                }});
                const track = new THREE.Mesh(trackGeometry, trackMaterial);
                track.position.set(0, yPos, 0);
                track.userData = {{ type: 'track', sample: sampleName }};
                trackGroup.add(track);

                // Show nucleotides on sample tracks
                if (showNucleotides) {{
                    const ntColors = {{
                        'A': 0x4CAF50, 'T': 0xf44336, 'G': 0xFFEB3B, 'C': 0x2196F3, 'N': 0x9e9e9e
                    }};
                    const ntWidth = trackLength / refSeq.length;

                    for (let i = 0; i < refSeq.length; i++) {{
                        const genomePos = currentStart + i;
                        const isGap = gapSet[sampleName]?.has(genomePos);
                        const altAllele = snpMap[sampleName]?.[genomePos];

                        // Determine nucleotide to show
                        let nt, isVariant = false;
                        if (isGap) {{
                            nt = 'N';  // Gap
                        }} else if (altAllele) {{
                            nt = altAllele.charAt(0).toUpperCase();  // SNP - show alt
                            isVariant = true;
                        }} else {{
                            nt = refSeq[i].toUpperCase();  // Reference
                        }}

                        const color = ntColors[nt] || 0x666666;
                        const xPos = (i * ntWidth) - trackLength/2 + ntWidth/2;

                        // Nucleotide box - brighter for variants
                        const ntGeom = new THREE.BoxGeometry(ntWidth * 0.85, 0.2, 0.7);
                        const ntMat = new THREE.MeshPhongMaterial({{
                            color: color,
                            emissive: isVariant ? color : 0x000000,
                            emissiveIntensity: isVariant ? 0.5 : 0,
                            transparent: !isVariant,
                            opacity: isVariant ? 1.0 : 0.6,
                            shininess: 80
                        }});
                        const ntMesh = new THREE.Mesh(ntGeom, ntMat);
                        ntMesh.position.set(xPos, yPos - 0.25, 0);
                        trackGroup.add(ntMesh);

                        // Letter sprite
                        const canvas = document.createElement('canvas');
                        const ctx = canvas.getContext('2d');
                        canvas.width = 64;
                        canvas.height = 64;
                        ctx.fillStyle = isVariant ? '#ffffff' : 'rgba(255,255,255,0.7)';
                        ctx.font = isVariant ? 'bold 48px monospace' : '40px monospace';
                        ctx.textAlign = 'center';
                        ctx.fillText(nt, 32, 48);
                        const texture = new THREE.CanvasTexture(canvas);
                        const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture, transparent: true }}));
                        sprite.scale.set(ntWidth * 0.8, ntWidth * 0.8, 1);
                        sprite.position.set(xPos, yPos - 0.25, 0.5);
                        trackGroup.add(sprite);
                    }}
                }}

                // Add sample label as a sprite
                const canvas = document.createElement('canvas');
                const ctx = canvas.getContext('2d');
                canvas.width = 256;
                canvas.height = 64;
                ctx.fillStyle = '#ffffff';
                ctx.font = 'bold 24px Arial';
                ctx.fillText(sampleName, 10, 40);

                const texture = new THREE.CanvasTexture(canvas);
                const spriteMaterial = new THREE.SpriteMaterial({{ map: texture, transparent: true }});
                const sprite = new THREE.Sprite(spriteMaterial);
                sprite.scale.set(8, 2, 1);
                sprite.position.set(-trackLength/2 - 5, yPos, 0);
                trackGroup.add(sprite);
            }});

            // Log what we're rendering
            console.log('Rendering', data.samples.length, 'samples:', data.samples);
            console.log('Gaps:', data.gaps?.length || 0);
            console.log('SNPs:', data.snps?.length || 0);

            // Add SNPs as vertical pillars
            if (data.snps && data.snps.length > 0) {{
                // Group SNPs by position for efficiency
                const snpsByPos = {{}};
                data.snps.forEach(snp => {{
                    if (!snpsByPos[snp.pos]) snpsByPos[snp.pos] = [];
                    snpsByPos[snp.pos].push(snp);
                }});

                Object.entries(snpsByPos).forEach(([pos, snpsAtPos]) => {{
                    const xPos = ((parseInt(pos) - currentStart) * scaleFactor) - trackLength/2;

                    snpsAtPos.forEach(snp => {{
                        const sampleIdx = data.samples.indexOf(snp.sample);
                        if (sampleIdx === -1) return;

                        // Check filter based on status
                        const status = snp.status || 'unknown';
                        let filterKey;
                        if (status === 'OnlyPipelineA' || status === 'SnippyOnly') {{
                            filterKey = 'snippy-only';
                        }} else if (status === 'OnlyPipelineB' || status === 'CfsanOnly') {{
                            filterKey = 'cfsan-only';
                        }} else if (status === 'Concordant' || status === 'Both') {{
                            filterKey = 'concordant';
                        }} else {{
                            filterKey = 'discordant';
                        }}
                        if (!filters[filterKey]) return;  // Skip if filtered out

                        const yBase = sampleIdx * trackSpacing + 1;

                        // Determine color based on status
                        let snpColor, emissiveColor;
                        if (status === 'OnlyPipelineA' || status === 'SnippyOnly') {{
                            snpColor = COLORS.snpSnippyOnly;
                            emissiveColor = 0x4A148C;
                        }} else if (status === 'OnlyPipelineB' || status === 'CfsanOnly') {{
                            snpColor = COLORS.snpCfsanOnly;
                            emissiveColor = 0x01579B;
                        }} else if (status === 'Concordant' || status === 'Both') {{
                            snpColor = COLORS.snpConcordant;
                            emissiveColor = 0x1B5E20;
                        }} else {{
                            // Discordant or unknown
                            snpColor = COLORS.snpDiscordant;
                            emissiveColor = 0xBF360C;
                        }}

                        // SNP pillar height based on quality (normalized)
                        const quality = snp.snippy?.qual || snp.cfsan?.qual || 100;
                        const height = (quality / 1000) * snpScale + 0.5;

                        const snpGeometry = new THREE.CylinderGeometry(0.08, 0.08, height, 8);
                        const snpMaterial = new THREE.MeshPhongMaterial({{
                            color: snpColor,
                            shininess: 80,
                            emissive: emissiveColor,
                            emissiveIntensity: 0.2
                        }});
                        const snpMesh = new THREE.Mesh(snpGeometry, snpMaterial);
                        snpMesh.position.set(xPos, yBase + height/2, 0);
                        snpMesh.userData = {{
                            type: 'snp',
                            position: parseInt(pos),
                            sample: snp.sample,
                            ref: snp.ref,
                            alt: snp.alt,
                            status: snp.status,
                            snippy: snp.snippy,
                            cfsan: snp.cfsan,
                            quality: quality
                        }};
                        snpGroup.add(snpMesh);
                    }});
                }});
            }}

            // Add gaps as colored blocks
            if (data.gaps && data.gaps.length > 0) {{
                data.gaps.forEach(gap => {{
                    const sampleIdx = data.samples.indexOf(gap.sample);
                    if (sampleIdx === -1) return;

                    // Check filter based on gap type
                    const gapFilterKey = 'gap-' + (gap.type || 'mm2');
                    if (!filters[gapFilterKey]) return;  // Skip if filtered out

                    const gapStart = ((gap.start - currentStart) * scaleFactor) - trackLength/2;
                    const gapWidth = (gap.end - gap.start) * scaleFactor;
                    const yPos = sampleIdx * trackSpacing + 1;

                    let color;
                    switch(gap.type) {{
                        case 'snippy': color = COLORS.gapSnippy; break;
                        case 'cfsan': color = COLORS.gapCfsan; break;
                        default: color = COLORS.gapMm2;
                    }}

                    const gapGeometry = new THREE.BoxGeometry(Math.max(gapWidth, 0.2), 0.4, 0.6);
                    const gapMaterial = new THREE.MeshPhongMaterial({{
                        color: color,
                        transparent: true,
                        opacity: 0.7
                    }});
                    const gapMesh = new THREE.Mesh(gapGeometry, gapMaterial);
                    gapMesh.position.set(gapStart + gapWidth/2, yPos, 0);
                    gapMesh.userData = {{
                        type: 'gap',
                        gapType: gap.type,
                        sample: gap.sample,
                        start: gap.start,
                        end: gap.end,
                        size: gap.end - gap.start
                    }};
                    gapGroup.add(gapMesh);
                }});
            }}

            // Update stats
            document.getElementById('objects').textContent =
                'Objects: ' + (trackGroup.children.length + snpGroup.children.length + gapGroup.children.length);
        }}

        async function loadRegion() {{
            const start = parseInt(document.getElementById('region-start').value) || 1;
            const size = parseInt(document.getElementById('region-size').value) || 5000;
            currentStart = start;
            currentEnd = Math.min(start + size, referenceLength);

            document.getElementById('region-display').textContent =
                `Region: ${{currentStart.toLocaleString()}} - ${{currentEnd.toLocaleString()}} bp`;

            try {{
                const response = await fetch(`/api/region?start=${{currentStart}}&end=${{currentEnd}}`);
                const data = await response.json();
                createTracks(data);

                // Update genome overview
                drawGenomeOverview();

                // Adjust camera to fit all tracks
                const numSamples = data.samples?.length || 1;
                const trackSpacing = parseFloat(document.getElementById('track-spacing').value);
                const centerY = (numSamples * trackSpacing) / 2;
                const regionSize = data.length || (currentEnd - currentStart);

                // For nucleotide regions, zoom out more to see the wider track
                const cameraZ = regionSize <= 500
                    ? Math.max(regionSize * 0.5, 30)  // Further back for wide nucleotide views
                    : 40;

                controls.target.set(0, centerY, 0);
                camera.position.set(0, centerY + 15, cameraZ);
                controls.update();
            }} catch (err) {{
                console.error('Failed to load region:', err);
            }}
        }}

        function setupInteraction() {{
            const raycaster = new THREE.Raycaster();
            const mouse = new THREE.Vector2();
            let hoveredObject = null;

            renderer.domElement.addEventListener('mousemove', (event) => {{
                mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
                mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

                raycaster.setFromCamera(mouse, camera);
                const intersects = raycaster.intersectObjects([...snpGroup.children, ...gapGroup.children]);

                if (intersects.length > 0) {{
                    const obj = intersects[0].object;
                    if (hoveredObject !== obj) {{
                        if (hoveredObject) {{
                            hoveredObject.material.emissiveIntensity = 0.2;
                        }}
                        hoveredObject = obj;
                        obj.material.emissiveIntensity = 0.6;
                        renderer.domElement.style.cursor = 'pointer';
                    }}
                }} else {{
                    if (hoveredObject) {{
                        hoveredObject.material.emissiveIntensity = 0.2;
                        hoveredObject = null;
                    }}
                    renderer.domElement.style.cursor = 'default';
                }}
            }});

            renderer.domElement.addEventListener('click', (event) => {{
                mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
                mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

                raycaster.setFromCamera(mouse, camera);
                const intersects = raycaster.intersectObjects([...snpGroup.children, ...gapGroup.children]);

                if (intersects.length > 0) {{
                    showInfo(intersects[0].object.userData);
                }}
            }});
        }}

        function showInfo(data) {{
            const panel = document.getElementById('info-panel');
            const title = document.getElementById('info-title');
            const content = document.getElementById('info-content');

            panel.style.display = 'block';

            if (data.type === 'snp') {{
                const statusColors = {{
                    'OnlyPipelineA': '#7B1FA2',
                    'SnippyOnly': '#7B1FA2',
                    'OnlyPipelineB': '#0288D1',
                    'CfsanOnly': '#0288D1',
                    'Concordant': '#4CAF50',
                    'Both': '#4CAF50',
                    'Discordant': '#ff5722'
                }};
                const statusLabels = {{
                    'OnlyPipelineA': 'Snippy Only',
                    'SnippyOnly': 'Snippy Only',
                    'OnlyPipelineB': 'CFSAN Only',
                    'CfsanOnly': 'CFSAN Only',
                    'Concordant': 'Concordant',
                    'Both': 'Concordant',
                    'Discordant': 'Discordant'
                }};
                const statusColor = statusColors[data.status] || '#888';
                const statusLabel = statusLabels[data.status] || data.status;

                title.textContent = `SNP at position ${{data.position.toLocaleString()}}`;
                content.innerHTML = `
                    <div class="info-row"><span class="info-label">Sample</span><span class="info-value">${{data.sample}}</span></div>
                    <div class="info-row"><span class="info-label">Reference</span><span class="info-value">${{data.ref || 'N/A'}}</span></div>
                    <div class="info-row"><span class="info-label">Alt Allele</span><span class="info-value">${{data.alt || 'N/A'}}</span></div>
                    <div class="info-row"><span class="info-label">Status</span><span class="info-value" style="color: ${{statusColor}}">${{statusLabel}}</span></div>
                    <div class="info-row"><span class="info-label">Quality</span><span class="info-value">${{data.quality.toFixed(1)}}</span></div>
                    ${{data.snippy ? `<div class="info-row"><span class="info-label">Snippy DP</span><span class="info-value">${{data.snippy.depth || 'N/A'}}</span></div>` : ''}}
                    ${{data.cfsan ? `<div class="info-row"><span class="info-label">CFSAN DP</span><span class="info-value">${{data.cfsan.depth || 'N/A'}}</span></div>` : ''}}
                `;
            }} else if (data.type === 'gap') {{
                title.textContent = `Gap (${{data.gapType.toUpperCase()}})`;
                content.innerHTML = `
                    <div class="info-row"><span class="info-label">Sample</span><span class="info-value">${{data.sample}}</span></div>
                    <div class="info-row"><span class="info-label">Start</span><span class="info-value">${{data.start.toLocaleString()}}</span></div>
                    <div class="info-row"><span class="info-label">End</span><span class="info-value">${{data.end.toLocaleString()}}</span></div>
                    <div class="info-row"><span class="info-label">Size</span><span class="info-value">${{data.size.toLocaleString()}} bp</span></div>
                    <div class="info-row"><span class="info-label">Source</span><span class="info-value">${{data.gapType === 'mm2' ? 'minimap2' : data.gapType === 'snippy' ? 'Snippy (N positions)' : 'CFSAN (zero coverage)'}}</span></div>
                `;
            }}
        }}

        function closeInfoPanel() {{
            document.getElementById('info-panel').style.display = 'none';
        }}

        function resetCamera() {{
            const trackSpacing = parseFloat(document.getElementById('track-spacing').value);
            const numSamples = samples.length || 5;
            const centerY = (numSamples * trackSpacing) / 2;
            camera.position.set(0, centerY + 15, 40);
            controls.target.set(0, centerY, 0);
            controls.update();
        }}

        function toggleAutoRotate() {{
            autoRotate = !autoRotate;
            controls.autoRotate = autoRotate;
            controls.autoRotateSpeed = 1.0;
        }}

        // Jump to first region with gaps
        function jumpToGaps() {{
            // Use known gap regions from report
            const gapStarts = [11500, 15800, 23500, 57500, 95300, 161100];  // Known gap positions
            const randomGap = gapStarts[Math.floor(Math.random() * gapStarts.length)];
            document.getElementById('region-start').value = Math.max(1, randomGap - 500);
            loadRegion();
        }}

        // Navigate forward/backward
        function navigate(direction) {{
            const size = parseInt(document.getElementById('region-size').value) || 5000;
            const currentStartInput = document.getElementById('region-start');
            let newStart = parseInt(currentStartInput.value) + (direction * size);
            newStart = Math.max(1, Math.min(newStart, referenceLength - size));
            currentStartInput.value = newStart;
            loadRegion();
        }}

        // Load snake/spiral view of whole genome
        async function loadSnakeGenome() {{
            document.getElementById('loading').style.display = 'flex';
            document.getElementById('loading').querySelector('div:last-child').textContent = 'Loading full genome as snake...';

            // Clear existing
            while(trackGroup.children.length) trackGroup.remove(trackGroup.children[0]);
            while(snpGroup.children.length) snpGroup.remove(snpGroup.children[0]);
            while(gapGroup.children.length) gapGroup.remove(gapGroup.children[0]);

            try {{
                // Load ALL data
                const [reportRes, snpRes] = await Promise.all([
                    fetch('/report.json'),
                    fetch(`/api/region?start=1&end=${{referenceLength}}`)
                ]);
                const report = await reportRes.json();
                const data = await snpRes.json();

                // Snake layout parameters
                const bpPerRow = 50000;  // 50kb per row
                const rowLength = 80;    // 3D units per row
                const rowGap = 15;       // Vertical gap between snake rows (space for samples)
                const numRows = Math.ceil(referenceLength / bpPerRow);
                const scaleFactor = rowLength / bpPerRow;
                const trackSpacing = 2;  // Space between sample tracks
                const numSamples = data.samples?.length || 0;
                const snpScale = parseFloat(document.getElementById('snp-scale')?.value) || 1;
                const gapScale = parseFloat(document.getElementById('gap-scale')?.value) || 1;

                console.log(`Snake: ${{numRows}} rows, ${{bpPerRow}} bp/row, ${{numSamples}} samples`);

                // Helper: convert bp position to 3D coordinates
                // Rows go into Z (depth), samples stack in Y, position along X
                function bpTo3D(bp) {{
                    const row = Math.floor(bp / bpPerRow);
                    const posInRow = bp % bpPerRow;
                    const x = (posInRow * scaleFactor) - rowLength/2;  // Position along X
                    const z = -row * rowGap;  // Rows go back into Z (depth)
                    return {{ x, y: 0, z, row }};
                }}

                // 1. CREATE SAMPLE TRACKS FOR EACH ROW (rows go into Z depth)
                for (let row = 0; row < numRows; row++) {{
                    const rowZ = -row * rowGap;  // Rows go back into Z

                    // Reference track (thin, at y=-0.5)
                    if (filters['reference']) {{
                        const refGeom = new THREE.BoxGeometry(rowLength, 0.15, 1.5);
                        const refMat = new THREE.MeshPhongMaterial({{ color: COLORS.reference, shininess: 60 }});
                        const refMesh = new THREE.Mesh(refGeom, refMat);
                        refMesh.position.set(0, -0.5, rowZ);
                        trackGroup.add(refMesh);
                    }}

                    // Sample tracks (stacked in Y, each row at different Z)
                    data.samples?.forEach((sampleName, idx) => {{
                        const yPos = idx * trackSpacing;
                        const trackColor = COLORS.trackColors[idx % COLORS.trackColors.length];

                        const trackGeom = new THREE.BoxGeometry(rowLength, 0.2, 1.2);
                        const trackMat = new THREE.MeshPhongMaterial({{
                            color: trackColor,
                            transparent: true,
                            opacity: 0.85,
                            shininess: 30
                        }});
                        const track = new THREE.Mesh(trackGeom, trackMat);
                        track.position.set(0, yPos, rowZ);
                        trackGroup.add(track);
                    }});

                    // Position label at start of each row
                    const rowStartBp = row * bpPerRow;
                    const canvas = document.createElement('canvas');
                    const ctx = canvas.getContext('2d');
                    canvas.width = 256;
                    canvas.height = 64;
                    ctx.fillStyle = '#fff';
                    ctx.font = 'bold 24px Arial';
                    ctx.fillText(`${{(rowStartBp/1000000).toFixed(2)}}M`, 10, 40);
                    const texture = new THREE.CanvasTexture(canvas);
                    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture, transparent: true }}));
                    sprite.scale.set(6, 1.5, 1);
                    sprite.position.set(-rowLength/2 - 5, 0, rowZ);
                    trackGroup.add(sprite);
                }}

                // Sample name labels (at front, z=0)
                data.samples?.forEach((sampleName, idx) => {{
                    const canvas = document.createElement('canvas');
                    const ctx = canvas.getContext('2d');
                    canvas.width = 256;
                    canvas.height = 64;
                    ctx.fillStyle = '#' + COLORS.trackColors[idx % COLORS.trackColors.length].toString(16).padStart(6, '0');
                    ctx.font = 'bold 24px Arial';
                    ctx.fillText(sampleName, 10, 40);
                    const texture = new THREE.CanvasTexture(canvas);
                    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture, transparent: true }}));
                    sprite.scale.set(10, 2.5, 1);
                    sprite.position.set(-rowLength/2 - 12, idx * trackSpacing, 5);
                    trackGroup.add(sprite);
                }});

                // 2. ADD SNPs AS PILLARS
                document.getElementById('loading').querySelector('div:last-child').textContent =
                    `Rendering ${{(data.snps || []).length}} SNPs...`;

                (data.snps || []).forEach(snp => {{
                    const sampleIdx = data.samples?.indexOf(snp.sample);
                    if (sampleIdx === -1 || sampleIdx === undefined) return;

                    // Check filter
                    const status = snp.status || 'unknown';
                    let filterKey;
                    if (status === 'OnlyPipelineA' || status === 'SnippyOnly') filterKey = 'snippy-only';
                    else if (status === 'OnlyPipelineB' || status === 'CfsanOnly') filterKey = 'cfsan-only';
                    else if (status === 'Concordant' || status === 'Both') filterKey = 'concordant';
                    else filterKey = 'discordant';
                    if (!filters[filterKey]) return;

                    const coord = bpTo3D(snp.position);
                    const yBase = sampleIdx * trackSpacing;  // Y is sample position

                    // Color based on status
                    let snpColor, emissiveColor;
                    if (status === 'OnlyPipelineA' || status === 'SnippyOnly') {{
                        snpColor = COLORS.snpSnippyOnly;
                        emissiveColor = 0x4A148C;
                    }} else if (status === 'OnlyPipelineB' || status === 'CfsanOnly') {{
                        snpColor = COLORS.snpCfsanOnly;
                        emissiveColor = 0x01579B;
                    }} else if (status === 'Concordant' || status === 'Both') {{
                        snpColor = COLORS.snpConcordant;
                        emissiveColor = 0x1B5E20;
                    }} else {{
                        snpColor = COLORS.snpDiscordant;
                        emissiveColor = 0xBF360C;
                    }}

                    const quality = snp.snippy?.qual || snp.cfsan?.qual || 100;
                    const height = Math.min((quality / 500) * snpScale + 0.3, 3);

                    const snpGeom = new THREE.CylinderGeometry(0.05, 0.05, height, 6);
                    const snpMat = new THREE.MeshPhongMaterial({{
                        color: snpColor,
                        shininess: 80,
                        emissive: emissiveColor,
                        emissiveIntensity: 0.2
                    }});
                    const snpMesh = new THREE.Mesh(snpGeom, snpMat);
                    snpMesh.position.set(coord.x, yBase + height/2 + 0.2, coord.z);  // X along row, Y up, Z is row depth
                    snpGroup.add(snpMesh);
                }});

                // 3. ADD GAPS AS COLORED BLOCKS
                document.getElementById('loading').querySelector('div:last-child').textContent =
                    `Rendering gaps...`;

                (data.gaps || []).forEach(gap => {{
                    const sampleIdx = data.samples?.indexOf(gap.sample);
                    if (sampleIdx === -1 || sampleIdx === undefined) return;

                    const gapFilterKey = 'gap-' + (gap.type || 'mm2');
                    if (!filters[gapFilterKey]) return;

                    const startCoord = bpTo3D(gap.start);
                    const endCoord = bpTo3D(gap.end);

                    // Only render if gap is within same row
                    let gapWidth, gapX;
                    if (startCoord.row !== endCoord.row) {{
                        // Gap spans rows - render only in start row
                        const rowEndBp = (startCoord.row + 1) * bpPerRow;
                        const endInRow = bpTo3D(rowEndBp - 1);
                        gapWidth = Math.abs(endInRow.x - startCoord.x) * gapScale;
                        gapX = (startCoord.x + endInRow.x) / 2;
                    }} else {{
                        gapWidth = Math.abs(endCoord.x - startCoord.x) * gapScale;
                        gapX = (startCoord.x + endCoord.x) / 2;
                    }}

                    const yPos = sampleIdx * trackSpacing;  // Y is sample position

                    let color;
                    switch(gap.type) {{
                        case 'snippy': color = COLORS.gapSnippy; break;
                        case 'cfsan': color = COLORS.gapCfsan; break;
                        default: color = COLORS.gapMm2;
                    }}

                    const gapGeom = new THREE.BoxGeometry(Math.max(0.08, gapWidth), 0.3, 0.8);
                    const gapMat = new THREE.MeshPhongMaterial({{
                        color: color,
                        transparent: true,
                        opacity: 0.8
                    }});
                    const gapMesh = new THREE.Mesh(gapGeom, gapMat);
                    gapMesh.position.set(gapX, yPos - 0.2, startCoord.z);  // Z is row depth
                    gapGroup.add(gapMesh);
                }});

                // Adjust camera - look down at the plane from above/front
                const centerZ = -numRows * rowGap / 2;
                const maxY = numSamples * trackSpacing;
                camera.position.set(0, maxY + 30, centerZ + 80);  // Above and in front
                controls.target.set(0, maxY / 2, centerZ);
                controls.update();

                document.getElementById('loading').style.display = 'none';
                document.getElementById('objects').textContent =
                    `Snake: ${{numRows}} rows × ${{numSamples}} samples, ${{snpGroup.children.length}} SNPs, ${{gapGroup.children.length}} gaps`;

            }} catch(e) {{
                console.error('Snake load error:', e);
                document.getElementById('loading').querySelector('div:last-child').textContent = 'Error: ' + e.message;
            }}
        }}

        async function loadSnakeView() {{
            // Legacy - redirect to new snake
            loadSnakeGenome();
            return;

            const binsPerRow = 50;
            const rowHeight = 4;
            const binSize = Math.ceil(referenceLength / 1000);
            const totalBins = Math.ceil(referenceLength / binSize);
            const numRows = Math.ceil(totalBins / binsPerRow);
            const trackWidth = 80;
            const binWidth = trackWidth / binsPerRow;
            const tubeRadius = 0.4;

            // Fetch gap density data
            const response = await fetch('/report.json');
            const report = await response.json();

            // Build gap density map
            const gapDensity = new Array(totalBins).fill(0);
            (report.reference_gaps || []).forEach(rg => {{
                (rg.reference_uncovered || []).forEach(gap => {{
                    const startBin = Math.floor(gap.start / binSize);
                    const endBin = Math.floor(gap.end / binSize);
                    for (let b = startBin; b <= endBin && b < totalBins; b++) {{
                        gapDensity[b]++;
                    }}
                }});
            }});

            // Find max density for normalization
            const maxDensity = Math.max(...gapDensity, 1);

            // Create S-shaped snake path
            const points = [];
            for (let row = 0; row < numRows; row++) {{
                const yPos = -row * rowHeight;
                const direction = row % 2 === 0 ? 1 : -1;

                // Horizontal segment
                const startX = direction === 1 ? -trackWidth/2 : trackWidth/2;
                const endX = direction === 1 ? trackWidth/2 : -trackWidth/2;

                // Add points along horizontal
                for (let i = 0; i <= 10; i++) {{
                    const t = i / 10;
                    points.push(new THREE.Vector3(
                        startX + (endX - startX) * t,
                        yPos,
                        0
                    ));
                }}

                // Add curved connector to next row (if not last row)
                if (row < numRows - 1) {{
                    const nextY = -(row + 1) * rowHeight;
                    const curveX = endX;

                    // Create smooth U-turn curve
                    for (let i = 1; i <= 8; i++) {{
                        const t = i / 8;
                        const angle = Math.PI * t;
                        const curveRadius = rowHeight / 2;
                        points.push(new THREE.Vector3(
                            curveX + (direction === 1 ? 1 : -1) * curveRadius * (1 - Math.cos(angle)),
                            yPos - curveRadius * Math.sin(angle) - curveRadius,
                            0
                        ));
                    }}
                }}
            }}

            // Create tube geometry for the snake
            const curve = new THREE.CatmullRomCurve3(points);
            const tubeGeometry = new THREE.TubeGeometry(curve, points.length * 4, tubeRadius, 8, false);

            // Create gradient material - green to blue along the genome
            const tubeMaterial = new THREE.MeshPhongMaterial({{
                color: 0x4CAF50,
                transparent: true,
                opacity: 0.8,
                shininess: 30
            }});
            const tubeMesh = new THREE.Mesh(tubeGeometry, tubeMaterial);
            trackGroup.add(tubeMesh);

            // Add gap markers along the snake
            for (let row = 0; row < numRows; row++) {{
                const yPos = -row * rowHeight;
                const direction = row % 2 === 0 ? 1 : -1;

                for (let col = 0; col < binsPerRow; col++) {{
                    const binIndex = row * binsPerRow + (direction === 1 ? col : binsPerRow - 1 - col);
                    if (binIndex >= totalBins) continue;

                    const density = gapDensity[binIndex];
                    if (density > 0) {{
                        const intensity = density / maxDensity;
                        const xPos = (col - binsPerRow/2 + 0.5) * binWidth;

                        // Gap sphere on the snake
                        const gapGeometry = new THREE.SphereGeometry(tubeRadius + 0.1 + intensity * 0.3, 8, 8);
                        const gapMaterial = new THREE.MeshPhongMaterial({{
                            color: 0xff5722,
                            transparent: true,
                            opacity: 0.4 + intensity * 0.4,
                            emissive: 0xff5722,
                            emissiveIntensity: 0.2
                        }});
                        const gapMesh = new THREE.Mesh(gapGeometry, gapMaterial);
                        gapMesh.position.set(xPos, yPos, 0);
                        gapMesh.userData = {{
                            type: 'snakeBin',
                            binIndex: binIndex,
                            startPos: binIndex * binSize,
                            endPos: Math.min((binIndex + 1) * binSize, referenceLength),
                            gapCount: density
                        }};
                        gapGroup.add(gapMesh);
                    }}
                }}

                // Add position labels
                const rowStartPos = row * binsPerRow * binSize;
                const labelX = (row % 2 === 0) ? -trackWidth/2 - 5 : trackWidth/2 + 5;

                const canvas = document.createElement('canvas');
                const ctx = canvas.getContext('2d');
                canvas.width = 128;
                canvas.height = 32;
                ctx.fillStyle = '#aaa';
                ctx.font = 'bold 16px Arial';
                ctx.textAlign = row % 2 === 0 ? 'right' : 'left';
                ctx.fillText(`${{(rowStartPos/1000000).toFixed(2)}}M`, row % 2 === 0 ? 120 : 8, 22);

                const texture = new THREE.CanvasTexture(canvas);
                const spriteMaterial = new THREE.SpriteMaterial({{ map: texture, transparent: true }});
                const sprite = new THREE.Sprite(spriteMaterial);
                sprite.scale.set(5, 1.2, 1);
                sprite.position.set(labelX, yPos, 0);
                trackGroup.add(sprite);
            }}

            // Add start/end markers
            const startGeom = new THREE.ConeGeometry(0.8, 1.5, 8);
            const startMat = new THREE.MeshPhongMaterial({{ color: 0x4CAF50, emissive: 0x4CAF50, emissiveIntensity: 0.3 }});
            const startMarker = new THREE.Mesh(startGeom, startMat);
            startMarker.position.set(-trackWidth/2 - 2, 0, 0);
            startMarker.rotation.z = Math.PI / 2;
            trackGroup.add(startMarker);

            const endRow = numRows - 1;
            const endX = (endRow % 2 === 0) ? trackWidth/2 + 2 : -trackWidth/2 - 2;
            const endGeom = new THREE.OctahedronGeometry(0.8);
            const endMat = new THREE.MeshPhongMaterial({{ color: 0xf44336, emissive: 0xf44336, emissiveIntensity: 0.3 }});
            const endMarker = new THREE.Mesh(endGeom, endMat);
            endMarker.position.set(endX, -endRow * rowHeight, 0);
            trackGroup.add(endMarker);

            // Adjust camera for snake view
            const centerY = -numRows * rowHeight / 2;
            camera.position.set(0, centerY, 100);
            controls.target.set(0, centerY, 0);
            controls.update();

            document.getElementById('objects').textContent =
                `Snake View: ${{numRows}} rows, ${{gapGroup.children.length}} gap regions`;
        }}

        // ============ FULL GENOME 3D (GPU ACCELERATED) ============
        async function loadFullGenome3D() {{
            document.getElementById('loading').style.display = 'flex';
            document.getElementById('loading').querySelector('div:last-child').textContent = 'Loading full genome (2.9M bp)...';

            // Clear existing
            while(trackGroup.children.length) trackGroup.remove(trackGroup.children[0]);
            while(snpGroup.children.length) snpGroup.remove(snpGroup.children[0]);
            while(gapGroup.children.length) gapGroup.remove(gapGroup.children[0]);

            try {{
                // Load ALL data
                const [reportRes, snpRes] = await Promise.all([
                    fetch('/report.json'),
                    fetch(`/api/region?start=1&end=${{referenceLength}}`)
                ]);
                const report = await reportRes.json();
                const regionData = await snpRes.json();

                const snps = regionData.snps || [];
                console.log('Loaded SNPs:', snps.length);
                console.log('Sample SNP:', snps[0]);

                document.getElementById('loading').querySelector('div:last-child').textContent =
                    `Processing ${{snps.length}} SNPs...`;

                // Layout parameters - BIGGER for visibility
                const bpPerUnit = 500;   // 500 bp = 1 unit (more dense)
                const rowWidth = 60;     // Width of each row
                const rowSpacing = 12;   // More space between rows
                const bpPerRow = rowWidth * bpPerUnit;  // 30kb per row
                const numRows = Math.ceil(referenceLength / bpPerRow);

                console.log('Layout:', {{ bpPerUnit, rowWidth, rowSpacing, bpPerRow, numRows }});

                // Helper: convert bp position to 3D coordinates (snake layout)
                function bpTo3D(bp) {{
                    const row = Math.floor(bp / bpPerRow);
                    const posInRow = bp % bpPerRow;
                    const direction = row % 2 === 0 ? 1 : -1;
                    const x = direction === 1
                        ? (posInRow / bpPerUnit) - rowWidth/2
                        : rowWidth/2 - (posInRow / bpPerUnit);
                    const y = -row * rowSpacing;
                    return {{ x, y, z: 0, row, direction }};
                }}

                // 1. CREATE REFERENCE TRACK (thin line at top)
                const refMat = new THREE.MeshPhongMaterial({{ color: 0x1a3a1a, transparent: true, opacity: 0.5 }});
                for (let row = 0; row < numRows; row++) {{
                    const y = -row * rowSpacing;
                    const refGeom = new THREE.BoxGeometry(rowWidth, 0.3, 0.5);
                    const refMesh = new THREE.Mesh(refGeom, refMat);
                    refMesh.position.set(0, y + 3, 0);
                    trackGroup.add(refMesh);
                }}

                // Get sample names and gaps
                const sampleGaps = {{}};
                (report.reference_gaps || []).forEach(rg => {{
                    sampleGaps[rg.sample] = rg.reference_uncovered || [];
                }});
                const sampleNames = Object.keys(sampleGaps);
                console.log('Samples:', sampleNames);

                // Sample colors
                const sampleColors = [0x5C6BC0, 0x26A69A, 0xFFCA28, 0xEF5350, 0x7E57C2, 0x66BB6A];

                // 2. CREATE SAMPLE TRACKS WITH GAPS CUT OUT
                sampleNames.forEach((sampleName, sampleIdx) => {{
                    const gaps = sampleGaps[sampleName];
                    const sampleColor = sampleColors[sampleIdx % sampleColors.length];
                    const sampleMat = new THREE.MeshPhongMaterial({{
                        color: sampleColor,
                        emissive: sampleColor,
                        emissiveIntensity: 0.15,
                        transparent: true,
                        opacity: 0.9
                    }});

                    // Sort gaps by start position
                    const sortedGaps = [...gaps].sort((a, b) => a.start - b.start);

                    // Build covered segments (regions between gaps)
                    const segments = [];
                    let lastEnd = 0;
                    sortedGaps.forEach(gap => {{
                        if (gap.start > lastEnd) {{
                            segments.push({{ start: lastEnd, end: gap.start }});
                        }}
                        lastEnd = Math.max(lastEnd, gap.end);
                    }});
                    if (lastEnd < referenceLength) {{
                        segments.push({{ start: lastEnd, end: referenceLength }});
                    }}

                    console.log(`${{sampleName}}: ${{segments.length}} covered segments, ${{gaps.length}} gaps`);

                    // Create geometry for covered segments using InstancedMesh
                    const segGeom = new THREE.BoxGeometry(1, 1.2, 1.5);
                    const instancedSegs = new THREE.InstancedMesh(segGeom, sampleMat, segments.length);
                    const matrix = new THREE.Matrix4();
                    const tempPos = new THREE.Vector3();
                    const tempScale = new THREE.Vector3();
                    const tempQuat = new THREE.Quaternion();

                    segments.forEach((seg, i) => {{
                        const startCoord = bpTo3D(seg.start);
                        const endCoord = bpTo3D(seg.end);

                        // Handle segments within same row
                        if (startCoord.row === endCoord.row) {{
                            const width = Math.abs(endCoord.x - startCoord.x);
                            const centerX = (startCoord.x + endCoord.x) / 2;

                            tempPos.set(centerX, startCoord.y - sampleIdx * 2.5, sampleIdx * 0.3);
                            tempScale.set(Math.max(0.1, width), 1, 1);
                            matrix.compose(tempPos, tempQuat, tempScale);
                            instancedSegs.setMatrixAt(i, matrix);
                        }} else {{
                            // Segment spans multiple rows - just show first row portion
                            const rowEnd = (startCoord.row + 1) * bpPerRow;
                            const endInRow = bpTo3D(Math.min(seg.end, rowEnd));
                            const width = Math.abs(endInRow.x - startCoord.x);
                            const centerX = (startCoord.x + endInRow.x) / 2;

                            tempPos.set(centerX, startCoord.y - sampleIdx * 2.5, sampleIdx * 0.3);
                            tempScale.set(Math.max(0.1, width), 1, 1);
                            matrix.compose(tempPos, tempQuat, tempScale);
                            instancedSegs.setMatrixAt(i, matrix);
                        }}
                    }});

                    instancedSegs.instanceMatrix.needsUpdate = true;
                    trackGroup.add(instancedSegs);

                    // Add sample label
                    const labelCanvas = document.createElement('canvas');
                    const labelCtx = labelCanvas.getContext('2d');
                    labelCanvas.width = 256;
                    labelCanvas.height = 64;
                    labelCtx.fillStyle = '#' + sampleColor.toString(16).padStart(6, '0');
                    labelCtx.font = 'bold 28px Arial';
                    labelCtx.fillText(sampleName, 10, 40);
                    const labelTex = new THREE.CanvasTexture(labelCanvas);
                    const labelSprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: labelTex, transparent: true }}));
                    labelSprite.scale.set(10, 2.5, 1);
                    labelSprite.position.set(-rowWidth/2 - 12, -sampleIdx * 2.5, 0);
                    trackGroup.add(labelSprite);
                }});

                // Add position labels every 10 rows
                for (let row = 0; row < numRows; row += 10) {{
                    const pos = row * bpPerRow;
                    const canvas = document.createElement('canvas');
                    const ctx = canvas.getContext('2d');
                    canvas.width = 256;
                    canvas.height = 64;
                    ctx.fillStyle = '#888';
                    ctx.font = 'bold 28px Arial';
                    ctx.fillText(`${{(pos/1000000).toFixed(1)}}M`, 10, 45);
                    const texture = new THREE.CanvasTexture(canvas);
                    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture, transparent: true }}));
                    sprite.scale.set(10, 2.5, 1);
                    sprite.position.set(rowWidth/2 + 8, -row * rowSpacing, 0);
                    trackGroup.add(sprite);
                }}

                // 2. CREATE SNIPPY/CFSAN GAP MARKERS (mm2 gaps shown as track breaks above)
                const gapMarkerColors = {{
                    snippy: 0xffeb3b,   // Yellow for snippy N-regions
                    cfsan: 0x00bcd4    // Cyan for cfsan gaps
                }};

                let totalGaps = 0;
                const gapGeometry = new THREE.BoxGeometry(1, 0.8, 2);

                (report.reference_gaps || []).forEach((rg, sampleIdx) => {{
                    // Only render snippy and cfsan gap markers (mm2 gaps are track breaks)
                    ['snippy', 'cfsan'].forEach((gapType, typeIdx) => {{
                        const gaps = gapType === 'snippy' ? (rg.snippy_gaps || []) : (rg.cfsan_gaps || []);

                        if (gaps.length === 0) return;

                        const gapMaterial = new THREE.MeshPhongMaterial({{
                            color: gapMarkerColors[gapType],
                            emissive: gapMarkerColors[gapType],
                            emissiveIntensity: 0.4,
                            transparent: true,
                            opacity: 0.9
                        }});

                        // Use InstancedMesh for efficiency
                        const instancedGaps = new THREE.InstancedMesh(gapGeometry, gapMaterial, gaps.length);
                        const matrix = new THREE.Matrix4();
                        const tempPos = new THREE.Vector3();
                        const tempScale = new THREE.Vector3();
                        const tempQuat = new THREE.Quaternion();

                        gaps.forEach((gap, i) => {{
                            const coord = bpTo3D(gap.start);
                            const gapWidth = Math.max(0.3, Math.min((gap.end - gap.start) / bpPerUnit, 8));

                            // Position below the sample tracks
                            tempPos.set(
                                coord.x + (coord.direction === 1 ? gapWidth/2 : -gapWidth/2),
                                coord.y - sampleIdx * 2.5 - 1.2,  // Just below the sample track
                                2 + typeIdx * 2  // Stagger in Z
                            );
                            tempScale.set(gapWidth, 0.6, 1);

                            matrix.compose(tempPos, tempQuat, tempScale);
                            instancedGaps.setMatrixAt(i, matrix);
                        }});

                        instancedGaps.instanceMatrix.needsUpdate = true;
                        gapGroup.add(instancedGaps);
                        totalGaps += gaps.length;
                    }});
                }});

                console.log('Total gaps rendered:', totalGaps);

                // 3. CREATE SNPs using Points - BIGGER size
                const snpPositions = new Float32Array(snps.length * 3);
                const snpColors = new Float32Array(snps.length * 3);

                const colorMap = {{
                    'OnlyPipelineA': new THREE.Color(0xaa00ff),  // Bright purple
                    'SnippyOnly': new THREE.Color(0xaa00ff),
                    'OnlyPipelineB': new THREE.Color(0x00bcd4),  // Cyan
                    'CfsanOnly': new THREE.Color(0x00bcd4),
                    'Concordant': new THREE.Color(0x00ff00),     // Bright green
                    'Both': new THREE.Color(0x00ff00),
                    'Discordant': new THREE.Color(0xff0000),     // Bright red
                    'default': new THREE.Color(0xff0000)
                }};

                snps.forEach((snp, i) => {{
                    const coord = bpTo3D(snp.position);
                    snpPositions[i * 3] = coord.x;
                    snpPositions[i * 3 + 1] = coord.y + 2;  // Above track
                    snpPositions[i * 3 + 2] = 2;

                    const color = colorMap[snp.status] || colorMap['default'];
                    snpColors[i * 3] = color.r;
                    snpColors[i * 3 + 1] = color.g;
                    snpColors[i * 3 + 2] = color.b;
                }});

                const snpGeometry = new THREE.BufferGeometry();
                snpGeometry.setAttribute('position', new THREE.BufferAttribute(snpPositions, 3));
                snpGeometry.setAttribute('color', new THREE.BufferAttribute(snpColors, 3));

                const snpMaterial = new THREE.PointsMaterial({{
                    size: 2.5,  // BIGGER points
                    vertexColors: true,
                    transparent: true,
                    opacity: 0.95,
                    sizeAttenuation: true
                }});

                const snpPoints = new THREE.Points(snpGeometry, snpMaterial);
                snpGroup.add(snpPoints);

                console.log('SNPs rendered:', snps.length);

                // 4. Adjust camera to see everything - start closer
                const centerY = -numRows * rowSpacing / 2;
                camera.position.set(0, centerY, 80);  // Closer view
                controls.target.set(0, centerY, 0);
                controls.minDistance = 20;
                controls.maxDistance = 500;
                controls.update();

                document.getElementById('loading').style.display = 'none';
                document.getElementById('objects').textContent =
                    `Full 3D: ${{snps.length.toLocaleString()}} SNPs, ${{totalGaps.toLocaleString()}} gaps, ${{numRows}} rows (scroll to zoom)`;

            }} catch (e) {{
                console.error('Error loading full genome:', e);
                document.getElementById('loading').querySelector('div:last-child').textContent = 'Error: ' + e.message;
            }}
        }}

        // ============ FULL GENOME MAP ============
        let genomeMapData = null;

        async function openGenomeMap() {{
            document.getElementById('genome-map-modal').style.display = 'block';
            document.getElementById('map-status').textContent = 'Loading all genome data...';

            if (!genomeMapData) {{
                // Load all data
                try {{
                    const [reportRes, snpRes] = await Promise.all([
                        fetch('/report.json'),
                        fetch(`/api/region?start=1&end=${{referenceLength}}`)
                    ]);
                    const report = await reportRes.json();
                    const regionData = await snpRes.json();

                    genomeMapData = {{
                        gaps: {{}},
                        snps: regionData.snps || []
                    }};

                    // Process gaps by sample
                    (report.reference_gaps || []).forEach(rg => {{
                        const sampleName = rg.sample;
                        genomeMapData.gaps[sampleName] = {{
                            mm2: rg.reference_uncovered || [],
                            snippy: rg.snippy_gaps || [],
                            cfsan: rg.cfsan_gaps || []
                        }};
                    }});

                    document.getElementById('map-status').textContent =
                        `Loaded: ${{genomeMapData.snps.length}} SNPs, ${{Object.keys(genomeMapData.gaps).length}} samples`;
                }} catch (e) {{
                    document.getElementById('map-status').textContent = 'Error loading data: ' + e.message;
                    return;
                }}
            }}

            renderGenomeMap();

            // Setup canvas click handler
            const canvas = document.getElementById('genome-map-canvas');
            canvas.onclick = (e) => {{
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                const zoom = parseInt(document.getElementById('map-zoom').value);
                const bpPerRow = Math.min(10000, Math.ceil(referenceLength / 100)) * zoom;
                const row = Math.floor(y / 25);
                const col = x;
                const pos = row * bpPerRow + col * zoom;

                if (pos > 0 && pos <= referenceLength) {{
                    closeGenomeMap();
                    const size = parseInt(document.getElementById('region-size').value) || 5000;
                    document.getElementById('region-start').value = Math.max(1, pos - size/2);
                    loadRegion();
                }}
            }};

            // Setup tooltip
            canvas.onmousemove = (e) => {{
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                const zoom = parseInt(document.getElementById('map-zoom').value);
                const bpPerRow = Math.min(10000, Math.ceil(referenceLength / 100)) * zoom;
                const row = Math.floor(y / 25);
                const col = x;
                const pos = row * bpPerRow + col * zoom;

                if (pos > 0 && pos <= referenceLength) {{
                    const tooltip = document.getElementById('map-tooltip');
                    tooltip.style.display = 'block';
                    tooltip.style.left = (e.clientX + 15) + 'px';
                    tooltip.style.top = (e.clientY + 15) + 'px';
                    tooltip.innerHTML = `<b>Position:</b> ${{pos.toLocaleString()}} bp<br><small>Click to view in 3D</small>`;
                }}
            }};

            canvas.onmouseleave = () => {{
                document.getElementById('map-tooltip').style.display = 'none';
            }};
        }}

        function closeGenomeMap() {{
            document.getElementById('genome-map-modal').style.display = 'none';
        }}

        function renderGenomeMap() {{
            if (!genomeMapData) return;

            const canvas = document.getElementById('genome-map-canvas');
            const ctx = canvas.getContext('2d');
            const zoom = parseInt(document.getElementById('map-zoom').value);

            // Calculate dimensions
            const rowHeight = 25;
            const sampleHeight = 4;
            const sampleNames = Object.keys(genomeMapData.gaps);
            const numSamples = sampleNames.length;
            const bpPerRow = Math.min(10000, Math.ceil(referenceLength / 100)) * zoom;
            const canvasWidth = Math.ceil(bpPerRow / zoom);
            const numRows = Math.ceil(referenceLength / bpPerRow);

            canvas.width = canvasWidth;
            canvas.height = numRows * rowHeight;

            // Clear
            ctx.fillStyle = '#0a0a10';
            ctx.fillRect(0, 0, canvas.width, canvas.height);

            document.getElementById('map-status').textContent =
                `Rendering ${{referenceLength.toLocaleString()}} bp at ${{zoom}} bp/px (${{numRows}} rows x ${{canvasWidth}} px)...`;

            // Build gap lookup for fast access (pixel-level)
            const gapArrays = {{}};
            sampleNames.forEach((sample, sampleIdx) => {{
                const sampleGaps = genomeMapData.gaps[sample];
                gapArrays[sample] = {{
                    mm2: new Uint8Array(Math.ceil(referenceLength / zoom)),
                    snippy: new Uint8Array(Math.ceil(referenceLength / zoom)),
                    cfsan: new Uint8Array(Math.ceil(referenceLength / zoom))
                }};

                // Fill gap arrays
                (sampleGaps.mm2 || []).forEach(g => {{
                    for (let p = Math.floor(g.start/zoom); p <= Math.floor(g.end/zoom) && p < gapArrays[sample].mm2.length; p++) {{
                        gapArrays[sample].mm2[p] = 1;
                    }}
                }});
                (sampleGaps.snippy || []).forEach(g => {{
                    for (let p = Math.floor(g.start/zoom); p <= Math.floor(g.end/zoom) && p < gapArrays[sample].snippy.length; p++) {{
                        gapArrays[sample].snippy[p] = 1;
                    }}
                }});
                (sampleGaps.cfsan || []).forEach(g => {{
                    for (let p = Math.floor(g.start/zoom); p <= Math.floor(g.end/zoom) && p < gapArrays[sample].cfsan.length; p++) {{
                        gapArrays[sample].cfsan[p] = 1;
                    }}
                }});
            }});

            // Draw rows
            for (let row = 0; row < numRows; row++) {{
                const rowStart = row * bpPerRow;
                const rowY = row * rowHeight;

                // Draw base track (green = covered)
                ctx.fillStyle = '#1a3a1a';
                ctx.fillRect(0, rowY, canvasWidth, rowHeight - 2);

                // Draw position marker
                ctx.fillStyle = '#444';
                ctx.font = '9px monospace';
                ctx.fillText(`${{(rowStart/1000000).toFixed(2)}}M`, 2, rowY + 10);

                // Draw gaps for each sample
                sampleNames.forEach((sample, sampleIdx) => {{
                    const yOffset = rowY + 12 + sampleIdx * sampleHeight;

                    for (let px = 0; px < canvasWidth; px++) {{
                        const binIdx = Math.floor((rowStart + px * zoom) / zoom);
                        if (binIdx >= gapArrays[sample].mm2.length) continue;

                        // Priority: mm2 gaps (gray), then snippy (orange), then cfsan (blue)
                        if (gapArrays[sample].mm2[binIdx]) {{
                            ctx.fillStyle = '#9e9e9e';
                            ctx.fillRect(px, yOffset, 1, sampleHeight - 1);
                        }} else if (gapArrays[sample].snippy[binIdx]) {{
                            ctx.fillStyle = '#e65100';
                            ctx.fillRect(px, yOffset, 1, sampleHeight - 1);
                        }} else if (gapArrays[sample].cfsan[binIdx]) {{
                            ctx.fillStyle = '#0d47a1';
                            ctx.fillRect(px, yOffset, 1, sampleHeight - 1);
                        }} else {{
                            ctx.fillStyle = '#2a5a2a';
                            ctx.fillRect(px, yOffset, 1, sampleHeight - 1);
                        }}
                    }}
                }});
            }}

            // Draw SNPs on top
            ctx.globalAlpha = 0.9;
            genomeMapData.snps.forEach(snp => {{
                const pos = snp.position;
                const row = Math.floor(pos / bpPerRow);
                const px = Math.floor((pos % bpPerRow) / zoom);
                const rowY = row * rowHeight;

                // Color by status
                const status = snp.status || '';
                if (status === 'OnlyPipelineA' || status === 'SnippyOnly') {{
                    ctx.fillStyle = '#7B1FA2';
                }} else if (status === 'OnlyPipelineB' || status === 'CfsanOnly') {{
                    ctx.fillStyle = '#0288D1';
                }} else if (status === 'Concordant' || status === 'Both') {{
                    ctx.fillStyle = '#4CAF50';
                }} else {{
                    ctx.fillStyle = '#ff5722';
                }}

                // Draw SNP marker (small dot at top of row)
                ctx.beginPath();
                ctx.arc(px, rowY + 6, 2, 0, Math.PI * 2);
                ctx.fill();
            }});
            ctx.globalAlpha = 1.0;

            document.getElementById('map-status').textContent =
                `${{referenceLength.toLocaleString()}} bp | ${{genomeMapData.snps.length}} SNPs | ${{numSamples}} samples | ${{zoom}} bp/pixel | Click to navigate`;
        }}

        // Draw genome overview
        function drawOverview() {{
            const canvas = document.getElementById('overview-canvas');
            if (!canvas) return;
            const ctx = canvas.getContext('2d');
            const width = canvas.width;
            const height = canvas.height;

            // Clear
            ctx.fillStyle = '#2a2a34';
            ctx.fillRect(0, 0, width, height);

            // Draw reference bar
            ctx.fillStyle = '#4CAF50';
            ctx.fillRect(0, height/2 - 3, width, 6);

            // Draw current view indicator
            const viewStart = (currentStart / referenceLength) * width;
            const viewWidth = Math.max(2, ((currentEnd - currentStart) / referenceLength) * width);
            ctx.fillStyle = 'rgba(255, 255, 255, 0.3)';
            ctx.fillRect(viewStart, 0, viewWidth, height);
            ctx.strokeStyle = '#fff';
            ctx.lineWidth = 2;
            ctx.strokeRect(viewStart, 0, viewWidth, height);

            // Draw gap density markers (simplified - known gap regions)
            const gapRegions = [
                [11000, 25000], [57000, 58000], [95000, 96000], [161000, 162000],
                [232000, 234000], [280000, 290000], [400000, 420000], [560000, 570000],
                [700000, 760000], [880000, 900000], [1070000, 1090000], [1280000, 1310000],
                [1470000, 1510000], [1700000, 1720000], [1890000, 1920000], [2080000, 2100000],
                [2350000, 2380000], [2540000, 2560000], [2700000, 2720000], [2870000, 2900000]
            ];
            ctx.fillStyle = 'rgba(158, 158, 158, 0.8)';
            gapRegions.forEach(([start, end]) => {{
                const x = (start / referenceLength) * width;
                const w = Math.max(1, ((end - start) / referenceLength) * width);
                ctx.fillRect(x, height/2 - 8, w, 16);
            }});

            // Add click handler
            canvas.onclick = (e) => {{
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const clickPos = Math.floor((x / width) * referenceLength);
                const size = parseInt(document.getElementById('region-size').value) || 5000;
                document.getElementById('region-start').value = Math.max(1, clickPos - size/2);
                loadRegion();
            }};
        }}

        // Update slider display and re-render
        function updateSlider(sliderId) {{
            const slider = document.getElementById(sliderId);
            const valueSpan = document.getElementById(sliderId + '-value');
            if (valueSpan) {{
                valueSpan.textContent = parseFloat(slider.value).toFixed(1);
            }}
            // Re-render with cached data if available
            if (lastLoadedData) {{
                createTracks(lastLoadedData);
            }}
        }}

        // Toggle filter and update visualization
        function toggleFilter(filterName) {{
            filters[filterName] = !filters[filterName];

            // Update UI
            const legendItem = document.querySelector(`[data-filter="${{filterName}}"]`);
            const statusIcon = legendItem.querySelector('.filter-status');

            if (filters[filterName]) {{
                legendItem.classList.remove('filtered');
                statusIcon.textContent = '\u2713';  // checkmark
                statusIcon.classList.remove('inactive');
            }} else {{
                legendItem.classList.add('filtered');
                statusIcon.textContent = '\u2717';  // X mark
                statusIcon.classList.add('inactive');
            }}

            // Re-render with cached data
            if (lastLoadedData) {{
                createTracks(lastLoadedData);
            }}
        }}

        function onWindowResize() {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }}

        function animate() {{
            requestAnimationFrame(animate);

            controls.update();
            renderer.render(scene, camera);

            // FPS counter
            frameCount++;
            const now = performance.now();
            if (now - lastTime >= 1000) {{
                document.getElementById('fps').textContent = 'FPS: ' + frameCount;
                frameCount = 0;
                lastTime = now;
            }}
        }}

        // Start
        init();
    </script>
</body>
</html>
"##, reference_name = reference_name, reference_length = reference_length, report_json = report_json)
}

/// Load gaps from Snippy BAM using samtools depth
fn load_snippy_gaps(sample_name: &str, ref_len: usize, min_depth: usize) -> Vec<(usize, usize)> {
    let bam_path = format!(
        "/home/IZSNT/a.deruvo/snps-study/results/{}_snippy/{}_snippy.bam",
        sample_name, sample_name
    );
    load_gaps_from_bam(&bam_path, ref_len, min_depth)
}

/// Load gaps from CFSAN BAM using samtools depth
fn load_cfsan_gaps(sample_name: &str, ref_len: usize, min_depth: usize) -> Vec<(usize, usize)> {
    let bam_path = format!(
        "/home/IZSNT/a.deruvo/snps-study/results/cfsan/samples/{}/reads.sorted.bam",
        sample_name
    );
    load_gaps_from_bam(&bam_path, ref_len, min_depth)
}

/// Generic function to load gaps from any BAM file
/// Returns gaps as 1-based positions: [(start, end), ...] where both are 1-based inclusive
fn load_gaps_from_bam(bam_path: &str, ref_len: usize, min_depth: usize) -> Vec<(usize, usize)> {
    use std::process::Command;

    let output = Command::new("samtools")
        .args(&["depth", "-a", bam_path])
        .output();

    let mut gaps = Vec::new();

    if let Ok(output) = output {
        if output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            let mut in_gap = false;
            let mut gap_start = 0;

            for line in stdout.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 3 {
                    let pos: usize = fields[1].parse().unwrap_or(0); // 1-based from samtools
                    let depth: usize = fields[2].parse().unwrap_or(0);

                    if depth < min_depth {
                        if !in_gap {
                            gap_start = pos; // Keep 1-based
                            in_gap = true;
                        }
                    } else {
                        if in_gap {
                            gaps.push((gap_start, pos)); // Store as [start, end) 1-based
                            in_gap = false;
                        }
                    }
                }
            }

            if in_gap {
                gaps.push((gap_start, ref_len + 1)); // End is exclusive, so ref_len + 1
            }
        }
    }

    log::info!("Loaded {} gap regions from {}", gaps.len(), bam_path);
    gaps
}

/// Load Snippy SNPs from VCF file (decomposing MNPs and complex variants)
fn load_snippy_snps(sample_name: &str) -> Vec<(usize, String, String, f64, usize)> {
    // Try to find VCF file
    let vcf_path = format!(
        "/home/IZSNT/a.deruvo/snps-study/results/{}_snippy/{}_snippy.filt.vcf",
        sample_name, sample_name
    );

    let mut snps = Vec::new();

    if let Ok(content) = std::fs::read_to_string(&vcf_path) {
        for line in content.lines() {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 8 {
                continue;
            }

            let info = fields[7];

            // Skip indels (insertions/deletions) - only handle substitutions
            if info.contains("TYPE=ins") || info.contains("TYPE=del") {
                continue;
            }

            let pos: usize = fields[1].parse().unwrap_or(0);
            let ref_allele = fields[3];
            let alt_allele = fields[4];
            let qual: f64 = fields[5].parse().unwrap_or(0.0);

            // Extract depth from INFO field
            let depth: usize = info
                .split(';')
                .find(|s| s.starts_with("DP="))
                .and_then(|s| s[3..].parse().ok())
                .unwrap_or(0);

            // Decompose: compare ref and alt base by base
            let ref_bytes = ref_allele.as_bytes();
            let alt_bytes = alt_allele.as_bytes();
            let min_len = ref_bytes.len().min(alt_bytes.len());

            for i in 0..min_len {
                let ref_base = ref_bytes[i] as char;
                let alt_base = alt_bytes[i] as char;

                // Only add if bases differ (actual SNP)
                if ref_base != alt_base {
                    snps.push((
                        pos + i,
                        ref_base.to_string(),
                        alt_base.to_string(),
                        qual,
                        depth
                    ));
                }
            }
        }
    }

    snps
}

/// Load CFSAN SNPs from VCF file
fn load_cfsan_snps(sample_name: &str) -> Vec<(usize, String, String, f64, usize)> {
    let vcf_path = format!(
        "/home/IZSNT/a.deruvo/snps-study/results/cfsan/samples/{}/var.flt.vcf",
        sample_name
    );

    let mut snps = Vec::new();

    if let Ok(content) = std::fs::read_to_string(&vcf_path) {
        for line in content.lines() {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 10 {
                continue;
            }

            let pos: usize = fields[1].parse().unwrap_or(0);
            let ref_allele = fields[3];
            let alt_allele = fields[4];

            // Skip indels (different length ref/alt)
            if ref_allele.len() != alt_allele.len() {
                continue;
            }

            // Extract depth from INFO field (ADP=)
            let info = fields[7];
            let depth: usize = info
                .split(';')
                .find(|s| s.starts_with("ADP="))
                .and_then(|s| s[4..].parse().ok())
                .unwrap_or(0);

            // CFSAN doesn't have QUAL in the same way, use depth as proxy
            let qual: f64 = depth as f64;

            // Decompose MNPs if any (compare base by base)
            let ref_bytes = ref_allele.as_bytes();
            let alt_bytes = alt_allele.as_bytes();

            for i in 0..ref_bytes.len() {
                let ref_base = ref_bytes[i] as char;
                let alt_base = alt_bytes[i] as char;

                if ref_base != alt_base {
                    snps.push((
                        pos + i,
                        ref_base.to_string(),
                        alt_base.to_string(),
                        qual,
                        depth
                    ));
                }
            }
        }
    }

    snps
}

/// Generate simple nucleotide-level alignment viewer
fn generate_nucleotide_html(reference_fasta: &str, report: &Report, _start: usize, _chars_per_line: usize) -> String {
    // Parse reference sequence from FASTA
    let ref_seq: String = reference_fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<Vec<_>>()
        .join("")
        .to_uppercase();

    let ref_len = ref_seq.len();
    let sample_names: Vec<&str> = report.reference_gaps.iter()
        .map(|rg| rg.sample.as_str())
        .collect();

    // Build gap regions as JSON for client-side rendering
    let mut gaps_json = String::from("{");
    for (i, rg) in report.reference_gaps.iter().enumerate() {
        if i > 0 { gaps_json.push(','); }
        gaps_json.push_str(&format!("\"{}\":[", rg.sample));
        for (j, region) in rg.reference_uncovered.iter().enumerate() {
            if j > 0 { gaps_json.push(','); }
            gaps_json.push_str(&format!("[{},{}]", region.start, region.end));
        }
        gaps_json.push(']');
    }
    gaps_json.push('}');

    // Load Snippy gaps for each sample (from BAM)
    let mut snippy_gaps_json = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { snippy_gaps_json.push(','); }
        let gaps = load_snippy_gaps(sample, ref_len, 1); // min depth 1
        snippy_gaps_json.push_str(&format!("\"{}\":[", sample));
        for (j, (start, end)) in gaps.iter().enumerate() {
            if j > 0 { snippy_gaps_json.push(','); }
            snippy_gaps_json.push_str(&format!("[{},{}]", start, end));
        }
        snippy_gaps_json.push(']');
    }
    snippy_gaps_json.push('}');

    // Load CFSAN gaps for each sample (from BAM)
    let mut cfsan_gaps_json = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { cfsan_gaps_json.push(','); }
        let gaps = load_cfsan_gaps(sample, ref_len, 1); // min depth 1
        cfsan_gaps_json.push_str(&format!("\"{}\":[", sample));
        for (j, (start, end)) in gaps.iter().enumerate() {
            if j > 0 { cfsan_gaps_json.push(','); }
            cfsan_gaps_json.push_str(&format!("[{},{}]", start, end));
        }
        cfsan_gaps_json.push(']');
    }
    cfsan_gaps_json.push('}');

    // Load Snippy SNPs for each sample
    let mut snippy_snps_json = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { snippy_snps_json.push(','); }
        let snps = load_snippy_snps(sample);
        snippy_snps_json.push_str(&format!("\"{}\":{{", sample));
        for (j, (pos, ref_al, alt_al, qual, depth)) in snps.iter().enumerate() {
            if j > 0 { snippy_snps_json.push(','); }
            snippy_snps_json.push_str(&format!(
                "\"{}\":[\"{}\",\"{}\",{:.1},{}]",
                pos, ref_al, alt_al, qual, depth
            ));
        }
        snippy_snps_json.push_str("}");
        log::info!("Loaded {} Snippy SNPs for {}", snps.len(), sample);
    }
    snippy_snps_json.push('}');

    // Load CFSAN SNPs for each sample
    let mut cfsan_snps_json = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { cfsan_snps_json.push(','); }
        let snps = load_cfsan_snps(sample);
        cfsan_snps_json.push_str(&format!("\"{}\":{{", sample));
        for (j, (pos, ref_al, alt_al, qual, depth)) in snps.iter().enumerate() {
            if j > 0 { cfsan_snps_json.push(','); }
            cfsan_snps_json.push_str(&format!(
                "\"{}\":[\"{}\",\"{}\",{:.1},{}]",
                pos, ref_al, alt_al, qual, depth
            ));
        }
        cfsan_snps_json.push_str("}");
        log::info!("Loaded {} CFSAN SNPs for {}", snps.len(), sample);
    }
    cfsan_snps_json.push('}');

    let samples_json = serde_json::to_string(&sample_names).unwrap_or_default();

    // Generate HTML with client-side rendering
    format!(r##"<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>CoreGuard - Nucleotide View</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Courier New', monospace;
            background: #0a0a0a;
            color: #e0e0e0;
            font-size: 11px;
            line-height: 1.3;
        }}
        .header {{
            background: #1a1a2e;
            padding: 8px 15px;
            border-bottom: 2px solid #e94560;
            position: sticky;
            top: 0;
            z-index: 100;
            display: flex;
            align-items: center;
            gap: 15px;
            flex-wrap: wrap;
        }}
        .header h1 {{ color: #e94560; font-size: 1.1rem; }}
        .nav {{ display: flex; gap: 8px; align-items: center; }}
        .nav input {{
            width: 100px;
            padding: 4px 8px;
            background: #2a2a3e;
            border: 1px solid #444;
            color: #fff;
            border-radius: 3px;
            font-family: inherit;
        }}
        .nav button {{
            padding: 4px 12px;
            background: #e94560;
            border: none;
            color: white;
            cursor: pointer;
            border-radius: 3px;
            font-family: inherit;
        }}
        .nav button:hover {{ background: #ff6b6b; }}
        .filters {{ display: flex; gap: 5px; align-items: center; }}
        .filters label {{ color: #888; }}
        .filter-btn {{
            padding: 3px 8px;
            background: #2a2a3e;
            border: 1px solid #444;
            color: #888;
            cursor: pointer;
            border-radius: 3px;
            font-size: 0.75rem;
        }}
        .filter-btn:hover {{ border-color: #666; color: #ccc; }}
        .filter-btn.active {{ background: #4a4a6e; border-color: #e94560; color: #fff; }}
        .row-check {{ cursor: pointer; margin-right: 3px; opacity: 0.6; }}
        .row-check:hover {{ opacity: 1; }}
        .row.hidden .sample, .row.hidden .ref {{ display: none; }}
        .row.hidden {{ opacity: 0.4; }}
        .info {{ color: #888; font-size: 0.8rem; }}
        .kpis {{ display: flex; gap: 15px; }}
        .kpi {{ font-size: 0.75rem; color: #888; }}
        .kpi-val {{ font-weight: bold; color: #4fc3f7; font-size: 0.85rem; }}
        .chunk-header {{
            background: #1a1a2e;
            color: #e94560;
            padding: 4px 10px;
            margin-top: 10px;
            font-size: 0.8rem;
            border-left: 3px solid #e94560;
        }}
        .block.compact {{ margin-bottom: 20px; }}
        .nuc {{ cursor: default; }}
        .nuc:hover, .gap:hover, .dim:hover {{ background: #333; }}
        #alignment {{
            padding: 5px 10px;
            padding-bottom: 40px;
        }}
        .block {{ margin-bottom: 8px; }}
        .pos-marker {{ color: #555; font-size: 10px; }}
        .row {{ white-space: pre; }}
        .ref {{ color: #4fc3f7; }}
        .sample {{ color: #ccc; }}
        .gap {{ color: #e94560; background: #2a1a1a; }}
        .snp-snippy {{ background: #2d4a1a; color: #8bc34a; cursor: pointer; }}
        .snp-snippy:hover {{ background: #4a7a2a; }}
        .snp-cfsan {{ background: #1a2d4a; color: #64b5f6; cursor: pointer; }}
        .snp-cfsan:hover {{ background: #2a4a7a; }}
        .snp-both {{ background: #4a1a4a; color: #ce93d8; cursor: pointer; }}
        .snp-both:hover {{ background: #7a2a7a; }}
        .snp-consensus {{ background: #1a4a2a; color: #4caf50; font-weight: bold; }}
        .snp-uncertain {{ background: #4a4a1a; color: #ffeb3b; font-style: italic; }}
        .snp-discord {{ background: #4a1a1a; color: #ff5722; text-decoration: underline; }}
        .dim {{ color: #333; }}
        .lbl {{ color: #666; display: inline-block; width: 75px; text-align: right; margin-right: 8px; }}
        .lbl.sub {{ color: #555; font-size: 9px; font-style: italic; }}
        .popup {{
            position: fixed;
            background: #1a1a2e;
            border: 1px solid #e94560;
            border-radius: 5px;
            padding: 10px 15px;
            z-index: 1000;
            font-size: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.5);
            max-width: 300px;
        }}
        .popup h4 {{ color: #e94560; margin-bottom: 8px; }}
        .popup .row {{ margin: 4px 0; }}
        .popup .label {{ color: #888; }}
        .popup .value {{ color: #4fc3f7; }}
        .popup .close {{ position: absolute; top: 5px; right: 10px; cursor: pointer; color: #888; }}
        .legend {{
            padding: 6px 15px;
            background: #1a1a1a;
            border-top: 1px solid #333;
            position: fixed;
            bottom: 0;
            width: 100%;
            display: flex;
            gap: 20px;
            font-size: 0.75rem;
        }}
        .legend span {{ margin-right: 15px; }}
        .legend .gap {{ color: #e94560; }}
        .legend .match {{ color: #666; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>CoreGuard Nucleotide View</h1>
        <div class="nav">
            <label>Position:</label>
            <input type="number" id="pos-input" value="0" min="0" max="{max_pos}">
            <button onclick="goToPos()">Go</button>
            <button onclick="prevPage()">&larr;</button>
            <button onclick="nextPage()">&rarr;</button>
        </div>
        <div class="filters">
            <label>Filter:</label>
            <button id="filter-all" class="filter-btn active" onclick="setFilter('all')">All</button>
            <button id="filter-snippy" class="filter-btn" onclick="setFilter('snippy')">Snippy only</button>
            <button id="filter-cfsan" class="filter-btn" onclick="setFilter('cfsan')">CFSAN only</button>
            <button id="filter-both" class="filter-btn" onclick="setFilter('both')">Both</button>
            <button id="filter-consensus" class="filter-btn" onclick="setFilter('consensus')">Consensus</button>
            <button id="filter-discord" class="filter-btn" onclick="setFilter('discord')">Discordant</button>
            <button id="filter-ingaps" class="filter-btn" onclick="setFilter('ingaps')">In mm2 gaps</button>
        </div>
        <div class="kpis">
            <span class="kpi"><span class="kpi-val" id="kpi-snippy">-</span> Snippy SNPs</span>
            <span class="kpi"><span class="kpi-val" id="kpi-cfsan">-</span> CFSAN SNPs</span>
            <span class="kpi"><span class="kpi-val" id="kpi-gaps">-</span> bp gaps</span>
            <span class="kpi"><span class="kpi-val" id="kpi-consensus">-</span> consensus</span>
        </div>
        <div class="info" id="info">
            Reference: {ref_len} bp | Samples: {num_samples} | Loading...
        </div>
    </div>

    <div id="alignment"></div>

    <div class="legend">
        <span><span class="ref">ATCG</span> = Reference</span>
        <span><span class="gap">----</span> = Gap (no coverage)</span>
        <span><span class="snp-consensus">T</span> = Consensus SNP (both agree)</span>
        <span><span class="snp-uncertain">T</span> = Uncertain (only 1 pipeline)</span>
        <span><span class="snp-discord">?</span> = Discordant (pipelines disagree)</span>
        <span><span class="snp-snippy">A</span> = Snippy SNP</span>
        <span><span class="snp-cfsan">A</span> = CFSAN SNP</span>
    </div>
    <div id="popup" class="popup" style="display:none;"></div>

    <script>
        const refSeq = "{ref_seq}";
        const refLen = {ref_len};
        const samples = {samples_json};
        const gaps = {gaps_json};
        const snippyGaps = {snippy_gaps_json};
        const cfsanGaps = {cfsan_gaps_json};
        const snippySnps = {snippy_snps_json};
        const cfsanSnps = {cfsan_snps_json};

        let currentPos = 0;
        let charsPerLine = 100;
        const linesPerPage = 40;
        const labelWidth = 83; // pixels for label
        let activeFilters = new Set(); // multiple filters can be active
        let filteredPositions = []; // positions matching current filter
        let filterIndex = 0; // current position in filtered array
        // Visibility state: track hidden rows
        let hiddenRows = new Set(); // e.g., "TE15676", "TE15676-snippy", "TE15676-cfsan"

        function toggleRow(rowId) {{
            if (hiddenRows.has(rowId)) {{
                hiddenRows.delete(rowId);
            }} else {{
                hiddenRows.add(rowId);
            }}
            rerender();
        }}

        function isRowVisible(rowId) {{
            return !hiddenRows.has(rowId);
        }}

        function rerender() {{
            if (activeFilters.size === 0) {{
                render();
            }} else {{
                renderFiltered();
            }}
        }}

        // Cached total KPIs (computed once at startup)
        let totalKPIs = {{ snippy: 0, cfsan: 0, consensus: 0, gaps: 0 }};

        // Compute total KPIs once at startup
        function computeTotalKPIs() {{
            let snippyCount = 0;
            let cfsanCount = 0;
            let consensusCount = 0;
            let totalGapBp = 0;

            for (let pos = 1; pos <= refLen; pos++) {{
                let allSnippy = true;
                let allCfsan = true;
                for (const sample of samples) {{
                    if (!getSnippySnp(sample, pos)) allSnippy = false;
                    if (!getCfsanSnp(sample, pos)) allCfsan = false;
                }}
                if (allSnippy) snippyCount++;
                if (allCfsan) cfsanCount++;
                if (allSnippy && allCfsan) consensusCount++;
            }}

            for (const sample of samples) {{
                const sampleGaps = gaps[sample] || [];
                for (const [start, end] of sampleGaps) {{
                    totalGapBp += (end - start);
                }}
            }}

            totalKPIs = {{ snippy: snippyCount, cfsan: cfsanCount, consensus: consensusCount, gaps: totalGapBp }};
        }}

        // Update KPIs based on current filter
        function updateKPIs() {{
            if (activeFilters.size === 0) {{
                // No filter: show totals
                document.getElementById('kpi-snippy').textContent = totalKPIs.snippy.toLocaleString();
                document.getElementById('kpi-cfsan').textContent = totalKPIs.cfsan.toLocaleString();
                document.getElementById('kpi-gaps').textContent = totalKPIs.gaps.toLocaleString();
                document.getElementById('kpi-consensus').textContent = totalKPIs.consensus.toLocaleString();
            }} else {{
                // Filter active: count in filtered positions
                let snippyInFilter = 0;
                let cfsanInFilter = 0;
                let consensusInFilter = 0;
                let gapsInFilter = 0;

                for (const pos of filteredPositions) {{
                    const hasSnippy = hasSnippyConsensus(pos);
                    const hasCfsan = hasCfsanConsensus(pos);
                    if (hasSnippy) snippyInFilter++;
                    if (hasCfsan) cfsanInFilter++;
                    if (hasSnippy && hasCfsan) consensusInFilter++;

                    for (const sample of samples) {{
                        if (isGap(sample, pos - 1)) gapsInFilter++;
                    }}
                }}

                document.getElementById('kpi-snippy').textContent = `${{snippyInFilter}}/${{totalKPIs.snippy}}`;
                document.getElementById('kpi-cfsan').textContent = `${{cfsanInFilter}}/${{totalKPIs.cfsan}}`;
                document.getElementById('kpi-gaps').textContent = `${{gapsInFilter}} pos`;
                document.getElementById('kpi-consensus').textContent = `${{consensusInFilter}}/${{totalKPIs.consensus}}`;
            }}
        }}

        function computeKPIs() {{
            computeTotalKPIs();
            updateKPIs();
        }}

        // Get Snippy SNP at position for sample
        function getSnippySnp(sample, pos) {{
            const sampleSnps = snippySnps[sample];
            if (!sampleSnps) return null;
            return sampleSnps[pos] || null; // [ref, alt, qual, depth]
        }}

        // Get CFSAN SNP at position for sample
        function getCfsanSnp(sample, pos) {{
            const sampleSnps = cfsanSnps[sample];
            if (!sampleSnps) return null;
            return sampleSnps[pos] || null; // [ref, alt, qual, depth]
        }}

        // Check if position has consensus Snippy SNP (all samples)
        function hasSnippyConsensus(pos) {{
            for (const sample of samples) {{
                if (!getSnippySnp(sample, pos)) return false;
            }}
            return true;
        }}

        // Check if position has consensus CFSAN SNP (all samples)
        function hasCfsanConsensus(pos) {{
            for (const sample of samples) {{
                if (!getCfsanSnp(sample, pos)) return false;
            }}
            return true;
        }}

        // Pre-compute filtered positions
        function computeFilteredPositions() {{
            filteredPositions = [];
            for (let pos = 1; pos <= refLen; pos++) {{
                const hasSnippy = hasSnippyConsensus(pos);
                const hasCfsan = hasCfsanConsensus(pos);

                // Check if any sample has a gap at this position
                let hasGap = false;
                for (const sample of samples) {{
                    if (isGap(sample, pos - 1)) {{ hasGap = true; break; }}
                }}

                // AND logic: position must match ALL active filters
                let matchAll = true;
                if (activeFilters.has('snippy') && !hasSnippy) matchAll = false;
                if (activeFilters.has('cfsan') && !hasCfsan) matchAll = false;
                if (activeFilters.has('both') && !(hasSnippy && hasCfsan)) matchAll = false;
                if (activeFilters.has('consensus') && !(hasSnippy && hasCfsan)) matchAll = false;
                if (activeFilters.has('discord') && !(hasSnippy !== hasCfsan)) matchAll = false;
                if (activeFilters.has('ingaps') && !hasGap) matchAll = false;

                // Must have at least one SNP to be relevant
                const hasAnySNP = hasSnippy || hasCfsan;
                if (matchAll && hasAnySNP) filteredPositions.push(pos);
            }}
            filterIndex = 0;
            console.log(`Filters [${{Array.from(activeFilters).join(', ')}}]: ${{filteredPositions.length}} positions`);
        }}

        // Toggle filter (multiple can be active)
        function setFilter(filter) {{
            if (filter === 'all') {{
                // Reset: clear all filters, show normal view
                activeFilters.clear();
                document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
                document.getElementById('filter-all').classList.add('active');
                currentPos = 0;
                render();
            }} else {{
                // Toggle this filter
                const btn = document.getElementById('filter-' + filter);
                if (activeFilters.has(filter)) {{
                    activeFilters.delete(filter);
                    btn.classList.remove('active');
                }} else {{
                    activeFilters.add(filter);
                    btn.classList.add('active');
                }}
                // Remove 'all' active state if any filter is selected
                document.getElementById('filter-all').classList.remove('active');

                if (activeFilters.size === 0) {{
                    // No filters selected, show normal view
                    document.getElementById('filter-all').classList.add('active');
                    currentPos = 0;
                    render();
                    updateKPIs();
                }} else {{
                    computeFilteredPositions();
                    filterIndex = 0;
                    renderFiltered();
                    updateKPIs();
                }}
            }}
        }}

        // Calculate chars per line based on window width
        function calcCharsPerLine() {{
            const charWidth = 7.2; // approximate monospace char width at 11px
            const availableWidth = window.innerWidth - labelWidth - 30;
            return Math.floor(availableWidth / charWidth);
        }}

        // Check if position is in a mm2 gap for given sample
        function isGap(sample, pos) {{
            const sampleGaps = gaps[sample] || [];
            for (const [start, end] of sampleGaps) {{
                if (pos >= start && pos < end) return true;
            }}
            return false;
        }}

        // Check if position is in a Snippy gap for given sample
        function isSnippyGap(sample, pos) {{
            const sampleGaps = snippyGaps[sample] || [];
            for (const [start, end] of sampleGaps) {{
                if (pos >= start && pos < end) return true;
            }}
            return false;
        }}

        // Check if position is in a CFSAN gap for given sample
        function isCfsanGap(sample, pos) {{
            const sampleGaps = cfsanGaps[sample] || [];
            for (const [start, end] of sampleGaps) {{
                if (pos >= start && pos < end) return true;
            }}
            return false;
        }}

        // Render alignment
        function render() {{
            charsPerLine = calcCharsPerLine();
            const container = document.getElementById('alignment');
            let html = '';

            const totalLines = Math.min(linesPerPage, Math.ceil((refLen - currentPos) / charsPerLine));

            for (let line = 0; line < totalLines; line++) {{
                const lineStart = currentPos + line * charsPerLine;
                const lineEnd = Math.min(lineStart + charsPerLine, refLen);
                if (lineStart >= refLen) break;

                const refSlice = refSeq.substring(lineStart, lineEnd);

                html += '<div class="block">';
                html += `<span class="pos-marker">${{lineStart + 1}}</span>`;

                // REF row with checkbox
                const refHidden = hiddenRows.has('REF') ? ' hidden' : '';
                const refChecked = !hiddenRows.has('REF') ? 'checked' : '';
                html += `<div class="row${{refHidden}}"><input type="checkbox" class="row-check" ${{refChecked}} onchange="toggleRow('REF')"><span class="lbl">REF</span><span class="ref">${{refSlice}}</span></div>`;

                // Per-sample rows
                for (const sample of samples) {{
                    // mm2 alignment row - shows reconstructed sample sequence
                    const sampleHidden = hiddenRows.has(sample) ? ' hidden' : '';
                    const sampleChecked = !hiddenRows.has(sample) ? 'checked' : '';
                    let mm2Line = '';
                    for (let i = lineStart; i < lineEnd; i++) {{
                        const pos = i + 1;
                        const isG = isGap(sample, i);
                        if (isG) {{
                            mm2Line += '<span class="gap">-</span>';
                        }} else {{
                            // Check for consensus SNP (both pipelines agree)
                            const snippySnp = getSnippySnp(sample, pos);
                            const cfsanSnp = getCfsanSnp(sample, pos);
                            if (snippySnp && cfsanSnp && snippySnp[1] === cfsanSnp[1]) {{
                                // Consensus: both pipelines found same ALT
                                mm2Line += `<span class="snp-consensus" title="pos ${{pos}}: consensus SNP">${{snippySnp[1]}}</span>`;
                            }} else if (snippySnp && cfsanSnp && snippySnp[1] !== cfsanSnp[1]) {{
                                // Discordant: different ALT alleles
                                mm2Line += `<span class="snp-discord" title="pos ${{pos}}: discordant (Snippy=${{snippySnp[1]}}, CFSAN=${{cfsanSnp[1]}})">${{refSeq[i]}}</span>`;
                            }} else if (snippySnp || cfsanSnp) {{
                                // Only one pipeline found SNP - show with uncertainty marker
                                const snp = snippySnp || cfsanSnp;
                                const pipeline = snippySnp ? 'Snippy' : 'CFSAN';
                                mm2Line += `<span class="snp-uncertain" title="pos ${{pos}}: only ${{pipeline}} (${{snp[1]}})">${{snp[1]}}</span>`;
                            }} else {{
                                // No SNP - show reference
                                mm2Line += refSeq[i];
                            }}
                        }}
                    }}
                    html += `<div class="row${{sampleHidden}}"><input type="checkbox" class="row-check" ${{sampleChecked}} onchange="toggleRow('${{sample}}')"><span class="lbl">${{sample}}</span><span class="sample">${{mm2Line}}</span></div>`;

                    // Snippy SNPs row
                    const snippyId = `${{sample}}-snippy`;
                    const snippyHidden = hiddenRows.has(snippyId) ? ' hidden' : '';
                    const snippyChecked = !hiddenRows.has(snippyId) ? 'checked' : '';
                    let snippyLine = '';
                    for (let i = lineStart; i < lineEnd; i++) {{
                        const pos = i + 1;
                        const snpData = getSnippySnp(sample, pos);
                        if (snpData) {{
                            const [ref, alt, qual, depth] = snpData;
                            snippyLine += `<span class="snp-snippy" data-pos="${{pos}}" data-sample="${{sample}}" data-pipeline="Snippy" data-ref="${{ref}}" data-alt="${{alt}}" data-qual="${{qual}}" data-depth="${{depth}}">${{alt}}</span>`;
                        }} else if (isSnippyGap(sample, pos)) {{
                            snippyLine += '<span class="gap">-</span>';
                        }} else {{
                            snippyLine += '<span class="dim">.</span>';
                        }}
                    }}
                    html += `<div class="row${{snippyHidden}}"><input type="checkbox" class="row-check" ${{snippyChecked}} onchange="toggleRow('${{snippyId}}')"><span class="lbl sub">snippy</span><span class="sample">${{snippyLine}}</span></div>`;

                    // CFSAN SNPs row
                    const cfsanId = `${{sample}}-cfsan`;
                    const cfsanHidden = hiddenRows.has(cfsanId) ? ' hidden' : '';
                    const cfsanChecked = !hiddenRows.has(cfsanId) ? 'checked' : '';
                    let cfsanLine = '';
                    for (let i = lineStart; i < lineEnd; i++) {{
                        const pos = i + 1;
                        const snpData = getCfsanSnp(sample, pos);
                        if (snpData) {{
                            const [ref, alt, qual, depth] = snpData;
                            cfsanLine += `<span class="snp-cfsan" data-pos="${{pos}}" data-sample="${{sample}}" data-pipeline="CFSAN" data-ref="${{ref}}" data-alt="${{alt}}" data-qual="${{qual}}" data-depth="${{depth}}">${{alt}}</span>`;
                        }} else if (isCfsanGap(sample, pos)) {{
                            cfsanLine += '<span class="gap">-</span>';
                        }} else {{
                            cfsanLine += '<span class="dim">.</span>';
                        }}
                    }}
                    html += `<div class="row${{cfsanHidden}}"><input type="checkbox" class="row-check" ${{cfsanChecked}} onchange="toggleRow('${{cfsanId}}')"><span class="lbl sub">cfsan</span><span class="sample">${{cfsanLine}}</span></div>`;
                }}

                html += '</div>';
            }}

            container.innerHTML = html;
            document.getElementById('pos-input').value = currentPos;
            document.getElementById('info').textContent =
                `Reference: ${{refLen.toLocaleString()}} bp | Samples: ${{samples.length}} | ${{charsPerLine}} bp/line | Position: ${{currentPos.toLocaleString()}}-${{Math.min(currentPos + charsPerLine * linesPerPage, refLen).toLocaleString()}}`;

            // Add click handlers for SNPs
            document.querySelectorAll('.snp-snippy, .snp-cfsan').forEach(el => {{
                el.addEventListener('click', showSnpPopup);
            }});
        }}

        // Show SNP popup
        function showSnpPopup(e) {{
            const el = e.target;
            const popup = document.getElementById('popup');
            const pos = el.dataset.pos;
            const sample = el.dataset.sample;
            const pipeline = el.dataset.pipeline;
            const ref = el.dataset.ref;
            const alt = el.dataset.alt;
            const qual = parseFloat(el.dataset.qual).toFixed(1);
            const depth = el.dataset.depth;

            popup.innerHTML = `
                <span class="close" onclick="closePopup()">&times;</span>
                <h4>${{pipeline}} SNP</h4>
                <div class="row"><span class="label">Position:</span> <span class="value">${{pos}}</span></div>
                <div class="row"><span class="label">Sample:</span> <span class="value">${{sample}}</span></div>
                <div class="row"><span class="label">Change:</span> <span class="value">${{ref}} &rarr; ${{alt}}</span></div>
                <div class="row"><span class="label">Quality:</span> <span class="value">${{qual}}</span></div>
                <div class="row"><span class="label">Depth:</span> <span class="value">${{depth}}x</span></div>
            `;

            popup.style.display = 'block';
            popup.style.left = (e.pageX + 10) + 'px';
            popup.style.top = (e.pageY - 10) + 'px';
        }}

        function closePopup() {{
            document.getElementById('popup').style.display = 'none';
        }}

        // Close popup when clicking outside
        document.addEventListener('click', (e) => {{
            if (!e.target.classList.contains('snp') && !e.target.closest('.popup')) {{
                closePopup();
            }}
        }});

        // Render filtered view (only SNP positions)
        function renderFiltered() {{
            charsPerLine = calcCharsPerLine();
            const container = document.getElementById('alignment');
            let html = '';

            if (filteredPositions.length === 0) {{
                html = '<div class="block"><p style="color:#888;padding:20px;">No positions match this filter.</p></div>';
                container.innerHTML = html;
                document.getElementById('info').textContent = `Filters: ${{Array.from(activeFilters).join(' + ')}} | 0 positions found`;
                return;
            }}

            const posPerPage = charsPerLine * 10; // more positions per page in compact view
            const startIdx = filterIndex;
            const endIdx = Math.min(startIdx + posPerPage, filteredPositions.length);
            const positions = filteredPositions.slice(startIdx, endIdx);

            // Compact view: all filtered positions in sequence, hover for position
            html += '<div class="block compact">';

            // REF row with checkbox
            const refHidden = hiddenRows.has('REF') ? ' hidden' : '';
            const refChecked = !hiddenRows.has('REF') ? 'checked' : '';
            let refLine = '';
            for (const pos of positions) {{
                refLine += `<span class="nuc" title="pos ${{pos}}">${{refSeq[pos - 1]}}</span>`;
            }}
            html += `<div class="row${{refHidden}}"><input type="checkbox" class="row-check" ${{refChecked}} onchange="toggleRow('REF')"><span class="lbl">REF</span><span class="ref">${{refLine}}</span></div>`;

            // Per-sample rows
            for (const sample of samples) {{
                // mm2 row - shows reconstructed sample sequence
                const sampleHidden = hiddenRows.has(sample) ? ' hidden' : '';
                const sampleChecked = !hiddenRows.has(sample) ? 'checked' : '';
                let sampleLine = '';
                for (const pos of positions) {{
                    const idx = pos - 1;
                    if (isGap(sample, idx)) {{
                        sampleLine += `<span class="gap" title="pos ${{pos}}">-</span>`;
                    }} else {{
                        // Check for consensus SNP
                        const snippySnp = getSnippySnp(sample, pos);
                        const cfsanSnp = getCfsanSnp(sample, pos);
                        if (snippySnp && cfsanSnp && snippySnp[1] === cfsanSnp[1]) {{
                            sampleLine += `<span class="snp-consensus" title="pos ${{pos}}: consensus">${{snippySnp[1]}}</span>`;
                        }} else if (snippySnp && cfsanSnp && snippySnp[1] !== cfsanSnp[1]) {{
                            sampleLine += `<span class="snp-discord" title="pos ${{pos}}: discordant">${{refSeq[idx]}}</span>`;
                        }} else if (snippySnp || cfsanSnp) {{
                            const snp = snippySnp || cfsanSnp;
                            sampleLine += `<span class="snp-uncertain" title="pos ${{pos}}: uncertain">${{snp[1]}}</span>`;
                        }} else {{
                            sampleLine += `<span class="nuc" title="pos ${{pos}}">${{refSeq[idx]}}</span>`;
                        }}
                    }}
                }}
                html += `<div class="row${{sampleHidden}}"><input type="checkbox" class="row-check" ${{sampleChecked}} onchange="toggleRow('${{sample}}')"><span class="lbl">${{sample}}</span><span class="sample">${{sampleLine}}</span></div>`;

                // Snippy row
                const snippyId = `${{sample}}-snippy`;
                const snippyHidden = hiddenRows.has(snippyId) ? ' hidden' : '';
                const snippyChecked = !hiddenRows.has(snippyId) ? 'checked' : '';
                let snippyLine = '';
                for (const pos of positions) {{
                    const snpData = getSnippySnp(sample, pos);
                    if (snpData) {{
                        const [ref, alt, qual, depth] = snpData;
                        snippyLine += `<span class="snp-snippy" title="pos ${{pos}}" data-pos="${{pos}}" data-sample="${{sample}}" data-pipeline="Snippy" data-ref="${{ref}}" data-alt="${{alt}}" data-qual="${{qual}}" data-depth="${{depth}}">${{alt}}</span>`;
                    }} else if (isSnippyGap(sample, pos)) {{
                        snippyLine += `<span class="gap" title="pos ${{pos}}">-</span>`;
                    }} else {{
                        snippyLine += `<span class="dim" title="pos ${{pos}}">.</span>`;
                    }}
                }}
                html += `<div class="row${{snippyHidden}}"><input type="checkbox" class="row-check" ${{snippyChecked}} onchange="toggleRow('${{snippyId}}')"><span class="lbl sub">snippy</span><span class="sample">${{snippyLine}}</span></div>`;

                // CFSAN row
                const cfsanId = `${{sample}}-cfsan`;
                const cfsanHidden = hiddenRows.has(cfsanId) ? ' hidden' : '';
                const cfsanChecked = !hiddenRows.has(cfsanId) ? 'checked' : '';
                let cfsanLine = '';
                for (const pos of positions) {{
                    const snpData = getCfsanSnp(sample, pos);
                    if (snpData) {{
                        const [ref, alt, qual, depth] = snpData;
                        cfsanLine += `<span class="snp-cfsan" title="pos ${{pos}}" data-pos="${{pos}}" data-sample="${{sample}}" data-pipeline="CFSAN" data-ref="${{ref}}" data-alt="${{alt}}" data-qual="${{qual}}" data-depth="${{depth}}">${{alt}}</span>`;
                    }} else if (isCfsanGap(sample, pos)) {{
                        cfsanLine += `<span class="gap" title="pos ${{pos}}">-</span>`;
                    }} else {{
                        cfsanLine += `<span class="dim" title="pos ${{pos}}">.</span>`;
                    }}
                }}
                html += `<div class="row${{cfsanHidden}}"><input type="checkbox" class="row-check" ${{cfsanChecked}} onchange="toggleRow('${{cfsanId}}')"><span class="lbl sub">cfsan</span><span class="sample">${{cfsanLine}}</span></div>`;
            }}

            html += '</div>';

            container.innerHTML = html;
            document.getElementById('pos-input').value = positions[0];
            document.getElementById('info').textContent =
                `Filters: ${{Array.from(activeFilters).join(' + ')}} | Showing ${{startIdx + 1}}-${{endIdx}} of ${{filteredPositions.length}} positions (hover for pos)`;

            // Add click handlers
            document.querySelectorAll('.snp-snippy, .snp-cfsan').forEach(el => {{
                el.addEventListener('click', showSnpPopup);
            }});
        }}

        function goToPos() {{
            const pos = parseInt(document.getElementById('pos-input').value) || 0;
            if (currentFilter === 'all') {{
                currentPos = Math.max(0, Math.min(pos, refLen - 1));
                render();
            }} else {{
                // Find nearest filtered position
                const idx = filteredPositions.findIndex(p => p >= pos);
                filterIndex = idx >= 0 ? idx : 0;
                renderFiltered();
            }}
        }}

        function prevPage() {{
            if (currentFilter === 'all') {{
                currentPos = Math.max(0, currentPos - charsPerLine * linesPerPage);
                render();
            }} else {{
                filterIndex = Math.max(0, filterIndex - charsPerLine * linesPerPage);
                renderFiltered();
            }}
        }}

        function nextPage() {{
            if (currentFilter === 'all') {{
                currentPos = Math.min(refLen - charsPerLine, currentPos + charsPerLine * linesPerPage);
                render();
            }} else {{
                filterIndex = Math.min(filteredPositions.length - 1, filterIndex + charsPerLine * linesPerPage);
                renderFiltered();
            }}
        }}

        document.getElementById('pos-input').addEventListener('keypress', e => {{
            if (e.key === 'Enter') goToPos();
        }});

        window.addEventListener('resize', () => {{
            if (activeFilters.size === 0) render();
            else renderFiltered();
        }});
        computeKPIs();
        render();
    </script>
</body>
</html>
"##,
        max_pos = ref_len.saturating_sub(1),
        ref_len = ref_len,
        num_samples = sample_names.len(),
        ref_seq = ref_seq,
        samples_json = samples_json,
        gaps_json = gaps_json,
        snippy_gaps_json = snippy_gaps_json,
        cfsan_gaps_json = cfsan_gaps_json,
        snippy_snps_json = snippy_snps_json,
        cfsan_snps_json = cfsan_snps_json
    )
}

/// Generate WASM-powered nucleotide viewer HTML
fn generate_nucleotide_wasm_html(reference_fasta: &str, report: &Report) -> String {
    // Extract reference sequence
    let ref_seq: String = reference_fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<Vec<_>>()
        .join("");

    let ref_len = ref_seq.len();
    let sample_names: Vec<&str> = report.reference_gaps.iter()
        .map(|rg| rg.sample.as_str())
        .collect();

    // Build JSON data for initialization
    let samples_json = format!("[{}]",
        sample_names.iter()
            .map(|s| format!("\"{}\"", s))
            .collect::<Vec<_>>()
            .join(",")
    );

    // Build mm2 gaps JSON
    let mut mm2_gaps_data = String::from("{");
    for (i, rg) in report.reference_gaps.iter().enumerate() {
        if i > 0 { mm2_gaps_data.push(','); }
        mm2_gaps_data.push_str(&format!("\"{}\":[", rg.sample));
        for (j, region) in rg.reference_uncovered.iter().enumerate() {
            if j > 0 { mm2_gaps_data.push(','); }
            mm2_gaps_data.push_str(&format!("[{},{}]", region.start, region.end));
        }
        mm2_gaps_data.push(']');
    }
    mm2_gaps_data.push('}');

    // Load Snippy gaps
    let mut snippy_gaps_data = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { snippy_gaps_data.push(','); }
        let gaps = load_snippy_gaps(sample, ref_len, 1);
        snippy_gaps_data.push_str(&format!("\"{}\":[", sample));
        for (j, (start, end)) in gaps.iter().enumerate() {
            if j > 0 { snippy_gaps_data.push(','); }
            snippy_gaps_data.push_str(&format!("[{},{}]", start, end));
        }
        snippy_gaps_data.push(']');
    }
    snippy_gaps_data.push('}');

    // Load CFSAN gaps
    let mut cfsan_gaps_data = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { cfsan_gaps_data.push(','); }
        let gaps = load_cfsan_gaps(sample, ref_len, 1);
        cfsan_gaps_data.push_str(&format!("\"{}\":[", sample));
        for (j, (start, end)) in gaps.iter().enumerate() {
            if j > 0 { cfsan_gaps_data.push(','); }
            cfsan_gaps_data.push_str(&format!("[{},{}]", start, end));
        }
        cfsan_gaps_data.push(']');
    }
    cfsan_gaps_data.push('}');

    // Load Snippy SNPs
    let mut snippy_snps_data = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { snippy_snps_data.push(','); }
        let snps = load_snippy_snps(sample);
        snippy_snps_data.push_str(&format!("\"{}\":[", sample));
        for (j, (pos, ref_al, alt_al, qual, depth)) in snps.iter().enumerate() {
            if j > 0 { snippy_snps_data.push(','); }
            snippy_snps_data.push_str(&format!(
                "[{},\"{}\",\"{}\",{:.1},{}]",
                pos, ref_al, alt_al, qual, depth
            ));
        }
        snippy_snps_data.push(']');
        log::info!("WASM: Loaded {} Snippy SNPs for {}", snps.len(), sample);
    }
    snippy_snps_data.push('}');

    // Load CFSAN SNPs
    let mut cfsan_snps_data = String::from("{");
    for (i, sample) in sample_names.iter().enumerate() {
        if i > 0 { cfsan_snps_data.push(','); }
        let snps = load_cfsan_snps(sample);
        cfsan_snps_data.push_str(&format!("\"{}\":[", sample));
        for (j, (pos, ref_al, alt_al, qual, depth)) in snps.iter().enumerate() {
            if j > 0 { cfsan_snps_data.push(','); }
            cfsan_snps_data.push_str(&format!(
                "[{},\"{}\",\"{}\",{:.1},{}]",
                pos, ref_al, alt_al, qual, depth
            ));
        }
        cfsan_snps_data.push(']');
        log::info!("WASM: Loaded {} CFSAN SNPs for {}", snps.len(), sample);
    }
    cfsan_snps_data.push('}');

    format!(r##"<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>CoreGuard WASM Nucleotide Viewer</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Courier New', monospace;
            background: #0a0a0f;
            color: #e0e0e0;
            font-size: 12px;
            padding: 10px;
        }}
        .header {{
            display: flex;
            align-items: center;
            gap: 15px;
            margin-bottom: 10px;
            padding: 10px;
            background: #1a1a2e;
            border-radius: 5px;
            flex-wrap: wrap;
        }}
        .header h2 {{ color: #e94560; margin-right: 20px; }}
        .header input {{
            background: #0a0a0f;
            border: 1px solid #333;
            color: #fff;
            padding: 5px 10px;
            border-radius: 3px;
            width: 120px;
        }}
        .header button {{
            background: #e94560;
            border: none;
            color: #fff;
            padding: 5px 15px;
            border-radius: 3px;
            cursor: pointer;
        }}
        .header button:hover {{ background: #ff6b6b; }}
        .kpi {{
            background: #1a2a1a;
            padding: 5px 10px;
            border-radius: 3px;
            font-size: 11px;
        }}
        .kpi .val {{ color: #4caf50; font-weight: bold; }}
        #status {{
            color: #888;
            font-style: italic;
        }}
        #alignment {{
            background: #0f0f1a;
            padding: 10px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        .block {{ margin-bottom: 5px; }}
        .row {{ white-space: nowrap; line-height: 1.4; }}
        .lbl {{
            color: #666;
            display: inline-block;
            width: 80px;
            text-align: right;
            margin-right: 8px;
        }}
        .lbl.sub {{ color: #555; font-size: 9px; font-style: italic; }}
        .ref {{ color: #888; }}
        .sample {{ color: #ccc; }}
        .gap {{ color: #e94560; background: #2a1a1a; }}
        .snp-consensus {{ background: #1a4a2a; color: #4caf50; font-weight: bold; }}
        .snp-uncertain {{ background: #4a4a1a; color: #ffeb3b; font-style: italic; }}
        .snp-discord {{ background: #4a1a1a; color: #ff5722; }}
        .snp-snippy {{ background: #2d4a1a; color: #8bc34a; }}
        .snp-cfsan {{ background: #1a2d4a; color: #64b5f6; }}
        .dim {{ color: #333; }}
        .legend {{
            margin-top: 10px;
            padding: 10px;
            background: #1a1a2e;
            border-radius: 5px;
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            font-size: 11px;
        }}
        .perf {{
            color: #888;
            font-size: 10px;
            margin-left: auto;
        }}
        .filters {{
            display: flex;
            gap: 5px;
            flex-wrap: wrap;
        }}
        .filter-btn {{
            background: #2a2a3a;
            border: 1px solid #444;
            color: #888;
            padding: 4px 10px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 10px;
        }}
        .filter-btn:hover {{ border-color: #666; color: #ccc; }}
        .filter-btn.active {{ background: #4a1a4a; border-color: #e94560; color: #fff; }}
        .filter-btn.active.snippy {{ background: #2d4a1a; border-color: #8bc34a; }}
        .filter-btn.active.cfsan {{ background: #1a2d4a; border-color: #64b5f6; }}
        .filter-btn.active.both {{ background: #4a4a1a; border-color: #ffeb3b; }}
        .filter-btn.active.consensus {{ background: #1a4a2a; border-color: #4caf50; }}
        .filter-btn.active.discord {{ background: #4a1a1a; border-color: #ff5722; }}
        .filter-btn.active.ingaps {{ background: #3a1a2a; border-color: #e94560; }}
        .nuc {{ cursor: default; }}
        .compact .nuc:hover, .compact .gap:hover, .compact .snp-snippy:hover, .compact .snp-cfsan:hover,
        .compact .snp-consensus:hover, .compact .snp-uncertain:hover {{
            outline: 1px solid #fff;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h2>CoreGuard WASM</h2>
        <span>Position:</span>
        <input type="number" id="pos-input" value="0" min="0" max="{max_pos}">
        <button onclick="goToPos()">Go</button>
        <button onclick="prevPage()">← Prev</button>
        <button onclick="nextPage()">Next →</button>
        <div class="filters">
            <button class="filter-btn snippy" data-filter="snippy" onclick="toggleFilter('snippy')">Snippy only</button>
            <button class="filter-btn cfsan" data-filter="cfsan" onclick="toggleFilter('cfsan')">CFSAN only</button>
            <button class="filter-btn both" data-filter="both" onclick="toggleFilter('both')">Both</button>
            <button class="filter-btn consensus" data-filter="consensus" onclick="toggleFilter('consensus')">Consensus</button>
            <button class="filter-btn discord" data-filter="discord" onclick="toggleFilter('discord')">Discordant</button>
            <button class="filter-btn ingaps" data-filter="ingaps" onclick="toggleFilter('ingaps')">In mm2 gaps</button>
            <button class="filter-btn" onclick="clearFilters()">Clear</button>
        </div>
        <div class="kpi">Snippy: <span class="val" id="kpi-snippy">-</span></div>
        <div class="kpi">CFSAN: <span class="val" id="kpi-cfsan">-</span></div>
        <div class="kpi">Gaps: <span class="val" id="kpi-gaps">-</span> bp</div>
        <div class="kpi">Consensus: <span class="val" id="kpi-consensus">-</span></div>
        <span id="status">Loading WASM...</span>
        <span class="perf" id="perf"></span>
    </div>

    <div id="alignment">Initializing WASM module...</div>

    <div class="legend">
        <span><span class="ref">ATCG</span> = Reference</span>
        <span><span class="gap">-</span> = Gap</span>
        <span><span class="snp-consensus">T</span> = Consensus SNP</span>
        <span><span class="snp-uncertain">T</span> = Uncertain</span>
        <span><span class="snp-snippy">A</span> = Snippy</span>
        <span><span class="snp-cfsan">A</span> = CFSAN</span>
    </div>

    <script type="module">
        import init, {{ GenomeData }} from '/wasm/coreguard_wasm.js';

        const refSeq = "{ref_seq}";
        const refLen = {ref_len};
        const samples = {samples_json};
        const mm2Gaps = {mm2_gaps_data};
        const snippyGaps = {snippy_gaps_data};
        const cfsanGaps = {cfsan_gaps_data};
        const snippySnps = {snippy_snps_data};
        const cfsanSnps = {cfsan_snps_data};

        let genomeData = null;
        let currentPos = 0;
        let charsPerLine = 100;
        const linesPerPage = 30;
        let activeFilters = new Set();
        let filteredPositions = [];
        let filterIndex = 0;

        async function initWasm() {{
            const status = document.getElementById('status');
            const t0 = performance.now();

            try {{
                status.textContent = 'Loading WASM module...';
                await init();

                status.textContent = 'Creating data structures...';
                genomeData = new GenomeData();
                genomeData.set_reference(refSeq);

                // Add samples
                for (const sample of samples) {{
                    genomeData.add_sample(sample);
                }}

                status.textContent = 'Loading mm2 gaps...';
                for (const [sample, gaps] of Object.entries(mm2Gaps)) {{
                    for (const [start, end] of gaps) {{
                        genomeData.add_mm2_gap(sample, start, end);
                    }}
                }}

                status.textContent = 'Loading Snippy gaps...';
                for (const [sample, gaps] of Object.entries(snippyGaps)) {{
                    for (const [start, end] of gaps) {{
                        genomeData.add_snippy_gap(sample, start, end);
                    }}
                }}

                status.textContent = 'Loading CFSAN gaps...';
                for (const [sample, gaps] of Object.entries(cfsanGaps)) {{
                    for (const [start, end] of gaps) {{
                        genomeData.add_cfsan_gap(sample, start, end);
                    }}
                }}

                status.textContent = 'Loading Snippy SNPs...';
                for (const [sample, snps] of Object.entries(snippySnps)) {{
                    for (const [pos, ref, alt, qual, depth] of snps) {{
                        genomeData.add_snippy_snp(sample, pos, ref.charCodeAt(0), alt.charCodeAt(0), qual, depth);
                    }}
                }}

                status.textContent = 'Loading CFSAN SNPs...';
                for (const [sample, snps] of Object.entries(cfsanSnps)) {{
                    for (const [pos, ref, alt, qual, depth] of snps) {{
                        genomeData.add_cfsan_snp(sample, pos, ref.charCodeAt(0), alt.charCodeAt(0), qual, depth);
                    }}
                }}

                status.textContent = 'Finalizing...';
                genomeData.finalize();

                const t1 = performance.now();
                status.textContent = `Ready (init: ${{(t1-t0).toFixed(0)}}ms)`;

                updateKPIs();
                render();

            }} catch (e) {{
                status.textContent = 'Error: ' + e;
                console.error(e);
            }}
        }}

        function updateKPIs() {{
            if (!genomeData) return;
            document.getElementById('kpi-snippy').textContent = genomeData.get_total_snippy_snps().toLocaleString();
            document.getElementById('kpi-cfsan').textContent = genomeData.get_total_cfsan_snps().toLocaleString();
            document.getElementById('kpi-gaps').textContent = genomeData.get_total_gap_bases().toLocaleString();
            // Consensus count requires filtering - simplified for now
            document.getElementById('kpi-consensus').textContent = '-';
        }}

        function render() {{
            if (!genomeData) return;

            if (activeFilters.size > 0) {{
                renderFiltered();
                return;
            }}

            const t0 = performance.now();
            const container = document.getElementById('alignment');

            // Calculate visible region
            const end = Math.min(currentPos + charsPerLine * linesPerPage, refLen);

            // Use WASM to render HTML
            const html = genomeData.render_region(JSON.stringify(samples), currentPos, end);

            container.innerHTML = html;
            document.getElementById('pos-input').value = currentPos;

            const t1 = performance.now();
            document.getElementById('perf').textContent = `render: ${{(t1-t0).toFixed(1)}}ms`;
            document.getElementById('status').textContent = `Pos ${{currentPos.toLocaleString()}}-${{end.toLocaleString()}} of ${{refLen.toLocaleString()}}`;
        }}

        function renderFiltered() {{
            const t0 = performance.now();
            const container = document.getElementById('alignment');

            // Get filtered positions from WASM
            const filtersStr = Array.from(activeFilters).join(',');
            const posJson = genomeData.get_filtered_positions(JSON.stringify(samples), filtersStr);
            filteredPositions = JSON.parse(posJson);

            if (filteredPositions.length === 0) {{
                container.innerHTML = '<div class="block"><p style="color:#888;padding:20px;">No positions match the selected filters.</p></div>';
                document.getElementById('status').textContent = `Filters: ${{filtersStr}} | 0 positions`;
                return;
            }}

            const posPerPage = charsPerLine * 10;
            const html = genomeData.render_filtered(JSON.stringify(samples), posJson, filterIndex, posPerPage);

            container.innerHTML = html;
            document.getElementById('pos-input').value = filteredPositions[filterIndex] || 0;

            const t1 = performance.now();
            document.getElementById('perf').textContent = `render: ${{(t1-t0).toFixed(1)}}ms`;
            document.getElementById('status').textContent = `Filters: ${{filtersStr}} | ${{filteredPositions.length.toLocaleString()}} positions`;
        }}

        window.toggleFilter = function(filter) {{
            const btn = document.querySelector(`[data-filter="${{filter}}"]`);
            if (activeFilters.has(filter)) {{
                activeFilters.delete(filter);
                btn.classList.remove('active');
            }} else {{
                activeFilters.add(filter);
                btn.classList.add('active');
            }}
            filterIndex = 0;
            render();
        }};

        window.clearFilters = function() {{
            activeFilters.clear();
            document.querySelectorAll('.filter-btn').forEach(btn => btn.classList.remove('active'));
            filterIndex = 0;
            render();
        }};

        window.goToPos = function() {{
            const pos = parseInt(document.getElementById('pos-input').value) || 0;
            if (activeFilters.size > 0) {{
                // Find nearest filtered position
                const idx = filteredPositions.findIndex(p => p >= pos);
                filterIndex = idx >= 0 ? idx : 0;
                renderFiltered();
            }} else {{
                currentPos = Math.max(0, Math.min(pos, refLen - 1));
                render();
            }}
        }};

        window.prevPage = function() {{
            if (activeFilters.size > 0) {{
                filterIndex = Math.max(0, filterIndex - charsPerLine * 10);
                renderFiltered();
            }} else {{
                currentPos = Math.max(0, currentPos - charsPerLine * linesPerPage);
                render();
            }}
        }};

        window.nextPage = function() {{
            if (activeFilters.size > 0) {{
                filterIndex = Math.min(filteredPositions.length - 1, filterIndex + charsPerLine * 10);
                renderFiltered();
            }} else {{
                currentPos = Math.min(refLen - charsPerLine, currentPos + charsPerLine * linesPerPage);
                render();
            }}
        }};

        // Initialize
        initWasm();
    </script>
</body>
</html>
"##,
        max_pos = ref_len.saturating_sub(1),
        ref_len = ref_len,
        ref_seq = ref_seq,
        samples_json = samples_json,
        mm2_gaps_data = mm2_gaps_data,
        snippy_gaps_data = snippy_gaps_data,
        cfsan_gaps_data = cfsan_gaps_data,
        snippy_snps_data = snippy_snps_data,
        cfsan_snps_data = cfsan_snps_data
    )
}
