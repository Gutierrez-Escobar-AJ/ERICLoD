#!/usr/bin/env python3
import re
import os
import json
import random
import shutil
import argparse
import logging
import subprocess
import gzip
from typing import List, Dict, Union
from collections import defaultdict
from pyfaidx import Fasta

##########################################
# ARGUMENT PARSER
##########################################

def parse_args():
    parser = argparse.ArgumentParser(description="Synthetic TP53 Validation Kit Generator")
    parser.add_argument("--vcf", type=str, default="clinvar_tp53_full.vcf", help="Path to ClinVar VCF file")
    parser.add_argument("--ref", type=str, default="GRCh38.fa", help="Path to reference FASTA")
    parser.add_argument("--vaf", type=float, default=1.0, help="Variant allele frequency (0-1)")
    parser.add_argument("--paired", action="store_true", help="Simulate paired-end reads")
    parser.add_argument("--rules", type=str, default="missense_variant=3", help="Selection rules (e.g. missense_variant=3,frameshift_variant=2)")
    parser.add_argument("--output", type=str, default="tp53_sim_output", help="Output directory")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument("--coverage", type=int, default=1000, help="Total coverage per region")
    parser.add_argument("--read-length", type=int, default=150, help="Read length")
    parser.add_argument("--mean-fraglen", type=int, default=300, help="Mean fragment length for paired-end")
    parser.add_argument("--fraglen-std", type=int, default=10, help="Standard deviation of fragment length")
    parser.add_argument("--art-system", type=str, default="HS25", help="ART platform system")
    parser.add_argument("--flank", type=int, default=250, help="Flanking region around variants")
    parser.add_argument("--check-ref-match", action="store_true", help="Enable checking and reporting of reference mismatches")
    return parser.parse_args()

def parse_rules(rule_str: str) -> Dict[str, int]:
    rules = {}
    for item in rule_str.split(','):
        try:
            key, value = item.split('=')
            rules[key.strip()] = int(value)
        except ValueError as e:
            logging.error(f"Invalid rule format for '{item}': {e}")
    return rules

##########################################
# 1. LOAD CLINVAR VARIANTS
##########################################

def parse_info_field(info_str: str) -> Dict[str, Union[str, List[str]]]:
    info_dict = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            if key == "MC":
                info_dict[key] = value.split(',')
            else:
                info_dict[key] = value.split('|') if '|' in value else value
    return info_dict

def extract_consequences(mc_list: Union[str, List[str]]) -> List[str]:
    """Returns list of consequence terms, logs problematic entries"""
    consequences = []
    if not mc_list:
        return consequences
    
    mc_iter = mc_list if isinstance(mc_list, list) else [mc_list]
    
    for mc in mc_iter:
        if not isinstance(mc, str):
            logging.warning(f"Non-string MC entry: {mc}")
            continue
        if '|' not in mc:
            logging.debug(f"Malformed MC entry (missing '|'): {mc}")
            continue
        parts = mc.split('|')
        if len(parts) < 2:
            logging.debug(f"Incomplete MC entry: {mc}")
            continue
        consequences.append(parts[1])
    
    return consequences

def load_clinvar_variants(vcf_path: str) -> List[Dict]:
    """Load and validate TP53 variants from ClinVar VCF"""
    variants = []
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    
    try:
        with open_fn(vcf_path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith("#"):
                    continue
                
                try:
                    fields = line.strip().split("\t")
                    if len(fields) < 8:
                        logging.warning(f"Skipping malformed line {line_num}: insufficient fields")
                        continue

                    chrom, pos, _, ref, alt, _, _, info_str = fields[:8]
                    info = parse_info_field(info_str)
                    
                    # Validate gene information
                    gene_info = info.get("GENEINFO", "")
                    if not gene_info:
                        logging.debug(f"Skipping line {line_num}: No GENEINFO")
                        continue
                        
                    if isinstance(gene_info, list):
                        gene_match = any("TP53" in g for g in gene_info)
                    else:
                        gene_match = "TP53" in gene_info
                        
                    if not gene_match:
                        logging.debug(f"Skipping line {line_num}: Not TP53 ({gene_info})")
                        continue

                    # Extract consequences
                    consequences = extract_consequences(info.get("MC", []))
                    if not consequences:
                        logging.info(f"Skipping line {line_num}: No valid consequences")
                        continue
                    if all("UTR" in c for c in consequences):
                        logging.info(f"Skipping line {line_num}: UTR-only consequences")
                        continue

                    variants.append({
                        "chrom": chrom.replace("chr", ""),
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                        "hgvs_g": info.get("CLNHGVS", [""])[0],
                        "pathogenicity": info.get("CLNSIG", ["Uncertain"])[0],
                        "consequence": consequences,
                        "gene": "TP53"
                    })
                    
                except Exception as e:
                    logging.warning(f"Error processing line {line_num}: {str(e)}")
                    continue

    except Exception as e:
        logging.error(f"File reading failed: {str(e)}")
        raise

    logging.info(f"Loaded {len(variants)} TP53 variants with valid consequences")
    return variants
    
##########################################
# 2. SELECT VARIANTS & MERGE REGIONS
##########################################

def select_variants(variants: List[Dict], rules: Dict[str, int]) -> List[Dict]:
    selected = []
    pool = variants.copy()
    for cons, count in rules.items():
        candidates = [v for v in pool if any(cons in c for c in v["consequence"])]
        if len(candidates) < count:
            logging.warning(f"Requested {count} variants for consequence '{cons}' but only found {len(candidates)}.")
        chosen = random.sample(candidates, min(count, len(candidates)))
        selected.extend(chosen)
        pool = [v for v in pool if v not in chosen]
    return selected

def merge_regions(variants: List[Dict], flank: int) -> List[Dict]:
    regions = defaultdict(list)
    for var in variants:
        chrom = var["chrom"]
        start = max(var["pos"] - flank, 1)
        end = var["pos"] + flank
        regions[chrom].append((start, end, var))
    
    merged = []
    for chrom in regions:
        intervals = sorted(regions[chrom], key=lambda x: x[0])
        current_start, current_end = intervals[0][0], intervals[0][1]
        current_vars = [intervals[0][2]]
        for start, end, var in intervals[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
                current_vars.append(var)
            else:
                merged.append({"chrom": chrom, "start": current_start, "end": current_end, "variants": current_vars})
                current_start, current_end = start, end
                current_vars = [var]
        merged.append({"chrom": chrom, "start": current_start, "end": current_end, "variants": current_vars})
    return merged

##########################################
# 3. EXTRACT & MUTATE SEQUENCES
##########################################

def extract_sequence(ref_fasta: str, chrom: str, start: int, end: int) -> str:
    try:
        with Fasta(ref_fasta) as fasta:
            chrom_key = chrom if chrom in fasta else f"chr{chrom}" if f"chr{chrom}" in fasta else chrom
            return str(fasta[chrom_key][start-1:end]).upper()
    except Exception as e:
        logging.error(f"Failed to extract {chrom}:{start}-{end}: {e}")
        raise

def apply_mutations(seq: str, variants: List[Dict], region_start: int, mismatch_report: list = None) -> str:
    seq_list = list(seq)
    for var in variants:
        rel_pos = var["pos"] - region_start
        ref_len = len(var["ref"])
        # Check if variant is within bounds
        if rel_pos < 0 or rel_pos + ref_len > len(seq_list):
            logging.error(f"Variant {var['chrom']}:{var['pos']} out of bounds for region starting at {region_start}")
            if mismatch_report is not None:
                mismatch_report.append({
                    "chrom": var["chrom"],
                    "pos": var["pos"],
                    "expected": var["ref"],
                    "found": "NA (out of bounds)",
                    "region_start": region_start,
                    "error": "out_of_bounds"
                })
            continue
        if seq_list[rel_pos:rel_pos+ref_len] != list(var["ref"]):
            found_bases = ''.join(seq_list[rel_pos:rel_pos+ref_len])
            logging.error(f"Ref mismatch at {var['chrom']}:{var['pos']}. Expected {var['ref']}, found {found_bases}")
            if mismatch_report is not None:
                mismatch_report.append({
                    "chrom": var["chrom"],
                    "pos": var["pos"],
                    "expected": var["ref"],
                    "found": found_bases,
                    "region_start": region_start,
                    "error": "ref_mismatch"
                })
            continue
        seq_list[rel_pos:rel_pos+ref_len] = list(var["alt"])
    return ''.join(seq_list)

##########################################
# 4. SIMULATE READS
##########################################

def simulate_reads_art(fasta_path: str, out_prefix: str, coverage: int, args: argparse.Namespace):
    # Check if ART Illumina tool is available
    if not shutil.which("art_illumina"):
        logging.error("ART Illumina tool ('art_illumina') is not available in your PATH. Please install or add it to the PATH.")
        raise FileNotFoundError("ART Illumina tool not found")
    
    cmd = [
        "art_illumina",
        "-ss", args.art_system,
        "-sam",
        "-i", fasta_path,
        "-l", str(args.read_length),
        "-f", str(coverage),
        "-o", out_prefix
    ]
    if args.paired:
        cmd += ["-p", "-m", str(args.mean_fraglen), "-s", str(args.fraglen_std)]
    logging.info(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Simulation failed: {e}")
        raise

##########################################
# 5. OUTPUT & PACKAGING
##########################################

def generate_truth_table(variants: List[Dict], output_path: str):
    with open(output_path, 'w') as f:
        f.write("chrom\tpos\tref\talt\thgvs_g\tpathogenicity\tconsequence\n")
        for v in variants:
            f.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['hgvs_g']}\t{v['pathogenicity']}\t{','.join(v['consequence'])}\n")

def package_output(output_dir: str):
    shutil.make_archive(output_dir, 'zip', output_dir)

def write_simulation_summary(output_dir: str, summary_stats: Dict):
    summary_file = os.path.join(output_dir, "simulation_summary.tsv")
    with open(summary_file, 'w') as f:
        f.write("metric\tvalue\n")
        for key, value in summary_stats.items():
            f.write(f"{key}\t{value}\n")
    logging.info(f"Simulation summary written to {summary_file}")

def write_mismatch_report(output_dir: str, mismatch_report: List[Dict]):
    report_file = os.path.join(output_dir, "ref_mismatch_summary.tsv")
    with open(report_file, 'w') as f:
        f.write("chrom\tpos\texpected\tfound\tregion_start\terror\n")
        for record in mismatch_report:
            f.write(f"{record['chrom']}\t{record['pos']}\t{record['expected']}\t{record['found']}\t{record['region_start']}\t{record['error']}\n")
    logging.info(f"Reference mismatch report written to {report_file}")

##########################################
# MAIN
##########################################

def main():
    args = parse_args()
    if args.seed is not None:
        random.seed(args.seed)
    if not (0 <= args.vaf <= 1):
        raise ValueError("VAF must be between 0 and 1")
    
    os.makedirs(args.output, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

    # 1. Load variants
    variants = load_clinvar_variants(args.vcf)
    
    # 2. Select variants and merge regions
    selected = select_variants(variants, parse_rules(args.rules))
    merged_regions = merge_regions(selected, args.flank)
    
    # Initialize mismatch records if check-ref-match is enabled
    mismatch_records = [] if args.check_ref_match else None
    regions_simulated = 0

    # 3. Process each merged region
    for i, region in enumerate(merged_regions):
        chrom, start, end = region["chrom"], region["start"], region["end"]
        try:
            # Extract reference sequence
            ref_seq = extract_sequence(args.ref, chrom, start, end)
            # Apply mutations to the reference (with mismatch reporting if enabled)
            mut_seq = apply_mutations(ref_seq, region["variants"], start, mismatch_records)
            
            # Write FASTA files for original and mutated sequences
            region_dir = os.path.join(args.output, f"region_{i+1}")
            os.makedirs(region_dir, exist_ok=True)
            with open(os.path.join(region_dir, "original.fa"), 'w') as f:
                f.write(f">original\n{ref_seq}\n")
            with open(os.path.join(region_dir, "mutated.fa"), 'w') as f:
                f.write(f">mutated\n{mut_seq}\n")
            
            # Simulate reads for wild-type and mutated sequences
            wt_coverage = args.coverage * (1 - args.vaf)
            simulate_reads_art(os.path.join(region_dir, "original.fa"), os.path.join(region_dir, "wt"), wt_coverage, args)
            simulate_reads_art(os.path.join(region_dir, "mutated.fa"), os.path.join(region_dir, "mut"), args.coverage * args.vaf, args)
            
            regions_simulated += 1
        
        except Exception as e:
            logging.error(f"Skipping region {chrom}:{start}-{end}: {e}")
            continue
    
    # Generate truth table and package output
    truth_table_path = os.path.join(args.output, "truth_table.tsv")
    generate_truth_table(selected, truth_table_path)
    package_output(args.output)
    logging.info("✅ Simulation complete")
    
    # Write simulation summary table
    summary_stats = {
        "total_variants_loaded": len(variants),
        "total_selected_variants": len(selected),
        "total_merged_regions": len(merged_regions),
        "total_regions_simulated": regions_simulated,
        "total_ref_mismatches": len(mismatch_records) if mismatch_records is not None else 0
    }
    write_simulation_summary(args.output, summary_stats)
    
    # Write reference mismatch report if enabled and mismatches exist
    if args.check_ref_match and mismatch_records and len(mismatch_records) > 0:
        write_mismatch_report(args.output, mismatch_records)

if __name__ == "__main__":
    main()
    
# Example of how to run:
# python3 ERICLoD.py --vcf clinvar_tp53_full.vcf --ref GRCh38.fa --vaf 0.1 --paired --rules missense_variant=5,frameshift_variant=2 --seed 42 --check-ref-match --output my_simulation

