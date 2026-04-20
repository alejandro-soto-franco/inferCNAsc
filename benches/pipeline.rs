//! Criterion benches for the inferCNAsc Rust core.
//!
//! `cargo bench --no-default-features` to run. Each bench seeds a synthetic
//! expression matrix drawn from a simple Normal(0, 1) per-gene distribution
//! with a planted chr1 loss in a random subset of cells, then times the full
//! pipeline and each stage in isolation.

use criterion::{Criterion, criterion_group, criterion_main};
use infercnasc::{assign_cnas_to_cells, find_cnas, smooth_expression};
use ndarray::Array2;
use std::hint::black_box;

fn synth_matrix(n_cells: usize, n_genes: usize, n_chroms: usize) -> Array2<f64> {
    let mut data = Vec::with_capacity(n_cells * n_genes);
    let mut seed: u64 = 0x1234_5678_ABCD_EF01;
    for c in 0..n_cells {
        for g in 0..n_genes {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u = ((seed >> 33) as f64) / ((1u64 << 31) as f64) - 1.0;
            // Plant a chr1 loss in half the cells: chr1 = gene 0..n_genes/n_chroms.
            let chrom_idx = g * n_chroms / n_genes;
            let planted = if chrom_idx == 0 && c % 2 == 0 {
                -2.0
            } else {
                0.0
            };
            data.push(u + planted);
        }
    }
    Array2::from_shape_vec((n_cells, n_genes), data).unwrap()
}

fn chroms(n_genes: usize, n_chroms: usize) -> Vec<String> {
    (0..n_genes)
        .map(|g| format!("{}", 1 + g * n_chroms / n_genes))
        .collect()
}

fn coords(n_genes: usize) -> (Vec<i64>, Vec<i64>, Vec<String>) {
    let mut starts = Vec::with_capacity(n_genes);
    let mut ends = Vec::with_capacity(n_genes);
    let mut names = Vec::with_capacity(n_genes);
    for g in 0..n_genes {
        starts.push((g as i64) * 1000);
        ends.push((g as i64) * 1000 + 999);
        names.push(format!("G{g}"));
    }
    (starts, ends, names)
}

fn bench_pipeline(c: &mut Criterion) {
    for (n_cells, n_genes, n_chroms) in [(500usize, 2000usize, 23usize), (2000, 10000, 23)] {
        let label = format!("{n_cells}x{n_genes}");
        let expr = synth_matrix(n_cells, n_genes, n_chroms);
        let chrom_strs = chroms(n_genes, n_chroms);
        let chrom_refs: Vec<&str> = chrom_strs.iter().map(String::as_str).collect();
        let (starts, ends, names) = coords(n_genes);
        let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();

        let mut group = c.benchmark_group(format!("pipeline_{label}"));
        group.sample_size(10);

        group.bench_function("smooth", |b| {
            b.iter(|| {
                let smoothed =
                    smooth_expression(black_box(&expr), black_box(&chrom_refs), 50).unwrap();
                black_box(smoothed);
            });
        });

        let smoothed = smooth_expression(&expr, &chrom_refs, 50).unwrap();

        group.bench_function("find_cnas", |b| {
            b.iter(|| {
                let (g, l) = find_cnas(black_box(&smoothed), 1.0);
                black_box((g, l));
            });
        });

        let (gains, losses) = find_cnas(&smoothed, 1.0);

        group.bench_function("assign", |b| {
            b.iter(|| {
                let recs = assign_cnas_to_cells(
                    black_box(&gains),
                    black_box(&losses),
                    black_box(&chrom_refs),
                    black_box(&starts),
                    black_box(&ends),
                    black_box(&name_refs),
                    3,
                );
                black_box(recs);
            });
        });

        group.bench_function("end_to_end", |b| {
            b.iter(|| {
                let s = smooth_expression(&expr, &chrom_refs, 50).unwrap();
                let (g, l) = find_cnas(&s, 1.0);
                let recs = assign_cnas_to_cells(&g, &l, &chrom_refs, &starts, &ends, &name_refs, 3);
                black_box(recs);
            });
        });

        group.finish();
    }
}

criterion_group!(benches, bench_pipeline);
criterion_main!(benches);
