use crate::types::InsertKey;
use crate::utils::{BASES, base_index, ratio, split_base_change};
use rayon::prelude::*;
use std::collections::HashMap;

/// Lower and upper quantiles of the reported (95%) percentile interval, in thousandths so the
/// quantile index is integer arithmetic.
const CI_LOW_PER_MILLE: usize = 25;
const CI_HIGH_PER_MILLE: usize = 975;

/// A block's `(key id, count)` pairs.
type SparseBlock = Vec<(u32, u32)>;
/// A block's genomic position: `(tid, start)`.
type BlockPosition = (i32, i64);

/// Per-block mismatch counts kept for block-bootstrap resampling.
///
/// Each block is one processing chunk of the genome. Counts are stored sparsely against an
/// interned key id, since a genome-wide run holds thousands of blocks and a string-keyed
/// map per block would cost several GB.
#[derive(Default)]
pub struct BlockCounts {
    index: HashMap<InsertKey, u32>,
    keys: Vec<InsertKey>,
    /// Each block's genomic position `(tid, start)` and its sparse counts.
    blocks: Vec<(BlockPosition, SparseBlock)>,
}

impl BlockCounts {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record the counts of the block starting at `position` (`(tid, start)`). Empty blocks
    /// are dropped: resampling is over blocks that contributed data, so empty stretches
    /// (off-target regions, gaps) don't dilute it. Blocks may arrive in any order; they are
    /// resampled in genomic order so a seed always draws the same blocks.
    pub fn add_block(&mut self, position: BlockPosition, counts: &HashMap<InsertKey, usize>) {
        if counts.is_empty() {
            return;
        }
        let block = counts
            .iter()
            .map(|(key, &count)| {
                let id = *self.index.entry(key.clone()).or_insert_with(|| {
                    self.keys.push(key.clone());
                    (self.keys.len() - 1) as u32
                });
                let count = u32::try_from(count).expect("per-block count exceeds u32");
                (id, count)
            })
            .collect();
        self.blocks.push((position, block));
    }

    /// Number of non-empty blocks recorded.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Blocks in genomic order, whatever order they were added in.
    fn in_order(&self) -> Vec<&(BlockPosition, SparseBlock)> {
        let mut ordered: Vec<&(BlockPosition, SparseBlock)> = self.blocks.iter().collect();
        ordered.sort_by_key(|(position, _)| *position);
        ordered
    }
}

/// `(ref, read)` base indices of an A/C/G/T `base_change`, or `None` for anything else.
fn base_change_indices(base_change: &str) -> Option<(usize, usize)> {
    let (ref_base, read_base) = split_base_change(base_change)?;
    Some((base_index(ref_base)?, base_index(read_base)?))
}

/// Linear-interpolation quantile (R's default, type 7) at `per_mille` thousandths. Reorders
/// `values` but needs only two order statistics, so it selects rather than sorts. Kept identical to
/// tractslip's (`src/bootstrap.rs`).
fn quantile(values: &mut [f32], per_mille: usize) -> f64 {
    let h = (values.len() - 1) * per_mille;
    let (lo, frac) = (h / 1000, h % 1000);
    let (_, &mut a, above) = values.select_nth_unstable_by(lo, f32::total_cmp);
    let b = above.iter().copied().min_by(f32::total_cmp).unwrap_or(a);
    let t = f64::from(u16::try_from(frac).expect("below 1000")) / 1000.0;
    (f64::from(b) - f64::from(a)).mul_add(t, f64::from(a))
}

/// One generator per replicate, forked in order from one seeded with `seed`, so replicate `r`
/// draws the same blocks whatever thread runs it. tractslip (`src/bootstrap.rs`) draws the same
/// way; keep the two in step.
fn replicate_rngs(seed: u64, replicates: usize) -> Vec<fastrand::Rng> {
    let mut root = fastrand::Rng::with_seed(seed);
    (0..replicates).map(|_| root.fork()).collect()
}

/// How many times each of `n_blocks` blocks is drawn in one replicate.
fn multiplicities(n_blocks: usize, rng: &mut fastrand::Rng) -> Vec<u64> {
    let mut multiplicity = vec![0u64; n_blocks];
    for _ in 0..n_blocks {
        multiplicity[rng.usize(0..n_blocks)] += 1;
    }
    multiplicity
}

/// Draws bootstrap replicates of the per-key counts. Shared by every interval type, so for a
/// given seed each statistic sees exactly the same resampled blocks.
struct Resampler<'a> {
    blocks: Vec<&'a SparseBlock>,
    discount: Vec<u64>,
    n_keys: usize,
}

impl<'a> Resampler<'a> {
    fn new(blocks: &'a BlockCounts, discounts: &HashMap<InsertKey, usize>) -> Self {
        Self {
            blocks: blocks
                .in_order()
                .into_iter()
                .map(|(_, block)| block)
                .collect(),
            discount: blocks
                .keys
                .iter()
                .map(|key| discounts.get(key).copied().unwrap_or(0) as u64)
                .collect(),
            n_keys: blocks.keys.len(),
        }
    }

    /// Per-key counts with each block weighted by `multiplicity`, minus discounts
    /// (saturating at zero).
    fn counts(&self, multiplicity: &[u64]) -> Vec<u64> {
        let mut counts = vec![0u64; self.n_keys];
        for (&block, &m) in self.blocks.iter().zip(multiplicity) {
            if m == 0 {
                continue;
            }
            for &(id, count) in block {
                counts[id as usize] += m * count as u64;
            }
        }
        for (count, &d) in counts.iter_mut().zip(&self.discount) {
            *count = count.saturating_sub(d);
        }
        counts
    }

    /// Counts for one replicate: `blocks.len()` blocks drawn with replacement.
    fn replicate(&self, rng: &mut fastrand::Rng) -> Vec<u64> {
        self.counts(&multiplicities(self.blocks.len(), rng))
    }

    /// Counts from the full data: every block once.
    fn full(&self) -> Vec<u64> {
        self.counts(&vec![1; self.blocks.len()])
    }
}

/// The (2.5%, 97.5%) interval of statistic `k` across replicates.
fn percentile_interval(replicates: &[Vec<f32>], k: usize) -> (f64, f64) {
    let mut values: Vec<f32> = replicates.iter().map(|rep| rep[k]).collect();
    (
        quantile(&mut values, CI_LOW_PER_MILLE),
        quantile(&mut values, CI_HIGH_PER_MILLE),
    )
}

/// Mismatch-rate statistics reported per (read_num, position) by `bootstrap_intervals`:
/// the 12 A/C/G/T mismatch classes, the 6 strand-folded classes, and all mismatches.
const N_RAW_CLASSES: usize = 12;
const FOLDED_CLASSES: [(usize, usize); 6] = [(1, 0), (1, 2), (1, 3), (3, 0), (3, 1), (3, 2)];
const N_RATE_STATS: usize = N_RAW_CLASSES + FOLDED_CLASSES.len() + 1;

/// Index of the reverse complement of an A/C/G/T index.
fn complement_index(b: usize) -> usize {
    3 - b
}

/// Raw mismatch classes in order: for each reference base, its three read bases.
fn raw_classes() -> impl Iterator<Item = (usize, usize)> {
    (0..4).flat_map(|r| (0..4).filter(move |&b| b != r).map(move |b| (r, b)))
}

/// A mismatch rate at one position, with its 95% block-bootstrap interval.
#[derive(Debug, Clone, PartialEq)]
pub struct RateInterval {
    pub read_num: u8,
    pub position: usize,
    /// A raw class (`G>A`), a strand-folded class (`C>T/G>A`), or `N>N` for all mismatches.
    pub class: String,
    pub count: u64,
    /// Bases read against this class's reference base(s): `rate`'s denominator.
    pub ref_total: u64,
    pub rate: f64,
    pub ci_low: f64,
    pub ci_high: f64,
    /// Every A/C/G/T base observed at this position: `rate_all`'s denominator, so each class's
    /// `rate_all` is its share of the `N>N` rate (the raw classes, like the folded ones, sum to it).
    pub all_bases: u64,
    pub rate_all: f64,
    pub rate_all_ci_low: f64,
    pub rate_all_ci_high: f64,
}

/// Every statistic's (numerator, denominator) for one (read_num, position) group, from its
/// `[ref][read]` base counts pooled over reference orders.
fn rate_stats(bases: &[[u64; 4]; 4]) -> [(u64, u64); N_RATE_STATS] {
    let ref_total = |r: usize| bases[r].iter().sum::<u64>();
    let mut stats = [(0u64, 0u64); N_RATE_STATS];
    for (i, (r, b)) in raw_classes().enumerate() {
        stats[i] = (bases[r][b], ref_total(r));
    }
    for (i, &(r, b)) in FOLDED_CLASSES.iter().enumerate() {
        let (cr, cb) = (complement_index(r), complement_index(b));
        stats[N_RAW_CLASSES + i] = (bases[r][b] + bases[cr][cb], ref_total(r) + ref_total(cr));
    }
    let mismatches = (0..4).map(|r| ref_total(r) - bases[r][r]).sum::<u64>();
    stats[N_RATE_STATS - 1] = (mismatches, (0..4).map(ref_total).sum());
    stats
}

fn rate_class_label(stat: usize) -> String {
    let base = |i: usize| BASES[i];
    if stat < N_RAW_CLASSES {
        let (r, b) = raw_classes().nth(stat).expect("raw class index in range");
        format!("{}>{}", base(r), base(b))
    } else if stat < N_RATE_STATS - 1 {
        let (r, b) = FOLDED_CLASSES[stat - N_RAW_CLASSES];
        let (cr, cb) = (complement_index(r), complement_index(b));
        format!("{}>{}/{}>{}", base(r), base(b), base(cr), base(cb))
    } else {
        "N>N".to_string()
    }
}

/// Frequency statistics: each key's count over its `(read_num, position, ref_base)` group
/// total, exactly as `normalize_mismatch_counts` normalizes.
struct FrequencyStats {
    /// Group of each key, or `None` for keys without a parseable reference base.
    key_group: Vec<Option<usize>>,
    n_groups: usize,
}

impl FrequencyStats {
    fn new(keys: &[InsertKey]) -> Self {
        let mut group_ids: HashMap<(u8, usize, char), usize> = HashMap::new();
        let key_group = keys
            .iter()
            .map(|key| {
                let (ref_base, _) = split_base_change(&key.base_change)?;
                let next = group_ids.len();
                Some(
                    *group_ids
                        .entry((key.read_num, key.base_position, ref_base))
                        .or_insert(next),
                )
            })
            .collect();
        Self {
            key_group,
            n_groups: group_ids.len(),
        }
    }

    /// One value per key.
    fn values(&self, counts: &[u64]) -> Vec<f32> {
        let mut group_totals = vec![0u64; self.n_groups];
        for (count, group) in counts.iter().zip(&self.key_group) {
            if let Some(g) = group {
                group_totals[*g] += count;
            }
        }
        counts
            .iter()
            .zip(&self.key_group)
            .map(|(&count, group)| group.map_or(0.0, |g| ratio(count, group_totals[g]) as f32))
            .collect()
    }

    fn intervals(
        &self,
        keys: &[InsertKey],
        replicates: &[Vec<f32>],
    ) -> HashMap<InsertKey, (f64, f64)> {
        (0..keys.len())
            .into_par_iter()
            .filter(|&k| self.key_group[k].is_some())
            .map(|k| (keys[k].clone(), percentile_interval(replicates, k)))
            .collect()
    }
}

/// Every rate statistic's (numerator, reference denominator) for one (read_num, position)
/// group, and all bases observed there.
type GroupStats = ([(u64, u64); N_RATE_STATS], u64);

/// Rate statistics per (read_num, position), pooled over reference orders.
struct RateStats {
    /// `(group, ref index, read index)` of each A/C/G/T key.
    key_slot: Vec<Option<(usize, usize, usize)>>,
    /// Each group's `(read_num, position)`, indexed by group.
    groups: Vec<(u8, usize)>,
}

impl RateStats {
    fn new(keys: &[InsertKey]) -> Self {
        let mut group_ids: HashMap<(u8, usize), usize> = HashMap::new();
        let mut groups = Vec::new();
        let key_slot = keys
            .iter()
            .map(|key| {
                let (r, b) = base_change_indices(&key.base_change)?;
                let g = *group_ids
                    .entry((key.read_num, key.base_position))
                    .or_insert_with(|| {
                        groups.push((key.read_num, key.base_position));
                        groups.len() - 1
                    });
                Some((g, r, b))
            })
            .collect();
        Self { key_slot, groups }
    }

    fn group_stats(&self, counts: &[u64]) -> Vec<GroupStats> {
        let mut bases = vec![[[0u64; 4]; 4]; self.groups.len()];
        for (&count, slot) in counts.iter().zip(&self.key_slot) {
            if let Some((g, r, b)) = slot {
                bases[*g][*r][*b] += count;
            }
        }
        bases
            .iter()
            .map(|b| (rate_stats(b), b.iter().flatten().sum()))
            .collect()
    }

    /// For each group and statistic, `[rate, rate_all]`.
    fn values(&self, counts: &[u64]) -> Vec<f32> {
        self.group_stats(counts)
            .iter()
            .flat_map(|(stats, all)| {
                stats
                    .iter()
                    .flat_map(move |&(num, den)| [ratio(num, den) as f32, ratio(num, *all) as f32])
            })
            .collect()
    }

    fn intervals(&self, full_counts: &[u64], replicates: &[Vec<f32>]) -> Vec<RateInterval> {
        let full = self.group_stats(full_counts);
        let mut order: Vec<usize> = (0..self.groups.len()).collect();
        order.sort_by_key(|&g| self.groups[g]);
        order
            .into_iter()
            .flat_map(|g| {
                let (read_num, position) = self.groups[g];
                let (stats, all_bases) = full[g];
                (0..N_RATE_STATS).filter_map(move |stat| {
                    let (count, ref_total) = stats[stat];
                    if ref_total == 0 {
                        return None;
                    }
                    let index = 2 * (g * N_RATE_STATS + stat);
                    let (ci_low, ci_high) = percentile_interval(replicates, index);
                    let (rate_all_ci_low, rate_all_ci_high) =
                        percentile_interval(replicates, index + 1);
                    Some(RateInterval {
                        read_num,
                        position,
                        class: rate_class_label(stat),
                        count,
                        ref_total,
                        rate: ratio(count, ref_total),
                        ci_low,
                        ci_high,
                        all_bases,
                        rate_all: ratio(count, all_bases),
                        rate_all_ci_low,
                        rate_all_ci_high,
                    })
                })
            })
            .collect()
    }
}

/// Block-bootstrap 95% intervals, all from one set of replicates.
#[derive(Debug, Default)]
pub struct BootstrapIntervals {
    /// Percentile interval for each key's normalized frequency. Keys without a parseable
    /// reference base get none.
    pub frequencies: HashMap<InsertKey, (f64, f64)>,
    /// Per-position rate intervals (see `RateInterval`), when requested.
    pub rates: Option<Vec<RateInterval>>,
}

/// Block-bootstrap 95% percentile intervals for each key's normalized frequency and,
/// with `with_rates`, for per-position mismatch rates.
///
/// Each replicate draws `blocks.len()` blocks with replacement, sums their counts, and
/// subtracts `discounts` (the per-key amounts `apply_external_discounts` removed from the
/// full-data counts, saturating at zero). Both statistics are computed from the same
/// replicate counts, so resampling runs once.
///
/// Frequencies are normalized exactly as `normalize_mismatch_counts` does: within
/// `(read_num, position, ref_base)` groups.
///
/// Rates are per (read_num, position), pooled over reference orders:
/// - each raw A/C/G/T class: `count / ref_total` of its reference base
/// - each strand-folded class: `C>T/G>A` = `(C>T + G>A) / (C + G)`
/// - all mismatches: `N>N`
///
/// Each rate also gets `rate_all`, its count over every base observed at the position, so
/// the classes' `rate_all` values add up to the `N>N` rate. Positions where a statistic's
/// denominator is zero are omitted for that statistic.
pub fn bootstrap_intervals(
    blocks: &BlockCounts,
    replicates: usize,
    seed: u64,
    discounts: &HashMap<InsertKey, usize>,
    with_rates: bool,
) -> BootstrapIntervals {
    if blocks.is_empty() || replicates == 0 {
        return BootstrapIntervals {
            rates: with_rates.then(Vec::new),
            ..Default::default()
        };
    }
    let resampler = Resampler::new(blocks, discounts);
    let frequency_stats = FrequencyStats::new(&blocks.keys);
    let rate_stats = with_rates.then(|| RateStats::new(&blocks.keys));

    let (frequency_replicates, rate_replicates): (Vec<Vec<f32>>, Vec<Vec<f32>>) =
        replicate_rngs(seed, replicates)
            .into_par_iter()
            .map(|mut rng| {
                let counts = resampler.replicate(&mut rng);
                let rates = rate_stats
                    .as_ref()
                    .map_or_else(Vec::new, |r| r.values(&counts));
                (frequency_stats.values(&counts), rates)
            })
            .unzip();

    BootstrapIntervals {
        frequencies: frequency_stats.intervals(&blocks.keys, &frequency_replicates),
        rates: rate_stats.map(|r| r.intervals(&resampler.full(), &rate_replicates)),
    }
}

/// Write `BootstrapIntervals::rates` rows as TSV.
pub fn write_rate_intervals(
    rows: &[RateInterval],
    position_label: &str,
    path: &str,
) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "read_num\t{position_label}\tclass\tcount\tref_total\trate\tci_low\tci_high\t\
         all_bases\trate_all\trate_all_ci_low\trate_all_ci_high"
    )?;
    for row in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.6e}\t{:.6e}\t{:.6e}\t{}\t{:.6e}\t{:.6e}\t{:.6e}",
            row.read_num,
            row.position,
            row.class,
            row.count,
            row.ref_total,
            row.rate,
            row.ci_low,
            row.ci_high,
            row.all_bases,
            row.rate_all,
            row.rate_all_ci_low,
            row.rate_all_ci_high
        )?;
    }
    w.flush()
}

/// Blocks whose expected mismatch count for a class is below this are left out of that
/// class's residual histogram and dispersion: Pearson residuals of small counts are too
/// skewed to compare against N(0, 1).
const MIN_EXPECTED_COUNT: f64 = 5.0;

/// Interior residual-histogram edges; the first and last bins are open-ended tails.
const RESIDUAL_BIN_WIDTH: f64 = 0.5;
const RESIDUAL_BIN_LIMIT: f64 = 4.0;
const N_RESIDUAL_BINS: usize = (2.0 * RESIDUAL_BIN_LIMIT / RESIDUAL_BIN_WIDTH) as usize + 2;

/// Between-block heterogeneity of one mismatch class (one read, pooled over positions).
#[derive(Debug, Clone)]
pub struct ClassDiagnostics {
    pub read_num: u8,
    pub base_change: String,
    /// Pooled rate over every block: the null each block is compared against.
    pub pooled_rate: f64,
    /// Blocks with enough expected counts to contribute a residual.
    pub blocks_used: usize,
    /// Blocks left out for having fewer than `MIN_EXPECTED_COUNT` expected mismatches.
    pub blocks_sparse: usize,
    /// Pearson dispersion `sum(z^2) / (blocks_used - 1)`: about 1 when blocks differ only by
    /// counting noise, well above 1 with real between-block variation or outlier blocks.
    pub dispersion: f64,
    /// Residual counts per bin: `(-inf, -4)`, `[-4, -3.5)`, ..., `[3.5, 4)`, `[4, inf)`.
    pub histogram: [usize; N_RESIDUAL_BINS],
    /// Block with the largest `|z|`, and its residual.
    pub most_extreme: Option<(BlockPosition, f64)>,
}

/// Lower edge of residual bin `i` (`-inf` for the first).
fn residual_bin_low(i: usize) -> f64 {
    if i == 0 {
        f64::NEG_INFINITY
    } else {
        -RESIDUAL_BIN_LIMIT + (i - 1) as f64 * RESIDUAL_BIN_WIDTH
    }
}

/// Upper edge of residual bin `i` (`inf` for the last).
fn residual_bin_high(i: usize) -> f64 {
    if i + 1 == N_RESIDUAL_BINS {
        f64::INFINITY
    } else {
        -RESIDUAL_BIN_LIMIT + i as f64 * RESIDUAL_BIN_WIDTH
    }
}

fn residual_bin(z: f64) -> usize {
    if z < -RESIDUAL_BIN_LIMIT {
        0
    } else if z >= RESIDUAL_BIN_LIMIT {
        N_RESIDUAL_BINS - 1
    } else {
        1 + ((z + RESIDUAL_BIN_LIMIT) / RESIDUAL_BIN_WIDTH) as usize
    }
}

/// Standard normal CDF, via the Abramowitz & Stegun 7.1.26 erf approximation (|error| < 1.5e-7).
fn normal_cdf(x: f64) -> f64 {
    if x.is_infinite() {
        return if x > 0.0 { 1.0 } else { 0.0 };
    }
    let t_arg = x.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.327_591_1 * t_arg);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let erf = 1.0 - poly * (-t_arg * t_arg).exp();
    0.5 * (1.0 + erf.copysign(x))
}

/// Per-class Pearson residuals of each block against the pooled rate.
///
/// For a block with `n` bases read against a class's reference base and `k` of that
/// mismatch, `z = (k - n p) / sqrt(n p (1 - p))`, with `p` the pooled rate. If blocks differ
/// only by counting noise, `z` is roughly N(0, 1); a wide histogram or dispersion well above
/// 1 means some blocks really differ, which the bootstrap intervals then reflect. Counts are
/// before `--discount-table`. Only A/C/G/T mismatch classes on reads 1 and 2 are reported.
pub fn block_residual_diagnostics(blocks: &BlockCounts) -> Vec<ClassDiagnostics> {
    // Per key: (read index, ref index, read-base index), or None for anything unreported.
    let key_class: Vec<Option<(usize, usize, usize)>> = blocks
        .keys
        .iter()
        .map(|key| {
            let read = match key.read_num {
                1 => 0,
                2 => 1,
                _ => return None,
            };
            let (ref_idx, read_idx) = base_change_indices(&key.base_change)?;
            Some((read, ref_idx, read_idx))
        })
        .collect();

    // Per block: counts[read][ref][read_base]; the ref total is the sum over read bases.
    let per_block: Vec<(BlockPosition, [[[u64; 4]; 4]; 2])> = blocks
        .in_order()
        .iter()
        .map(|(position, block)| {
            let mut counts = [[[0u64; 4]; 4]; 2];
            for &(id, count) in block {
                if let Some((read, r, b)) = key_class[id as usize] {
                    counts[read][r][b] += count as u64;
                }
            }
            (*position, counts)
        })
        .collect();

    let mut diagnostics = Vec::new();
    for read in 0..2 {
        for r in 0..4 {
            for b in (0..4).filter(|&b| b != r) {
                let (mut k_all, mut n_all) = (0u64, 0u64);
                for (_, counts) in &per_block {
                    k_all += counts[read][r][b];
                    n_all += counts[read][r].iter().sum::<u64>();
                }
                if k_all == 0 || k_all == n_all {
                    continue;
                }
                let p = k_all as f64 / n_all as f64;

                let mut diag = ClassDiagnostics {
                    read_num: read as u8 + 1,
                    base_change: format!("{}>{}", BASES[r], BASES[b]),
                    pooled_rate: p,
                    blocks_used: 0,
                    blocks_sparse: 0,
                    dispersion: f64::NAN,
                    histogram: [0; N_RESIDUAL_BINS],
                    most_extreme: None,
                };
                let mut sum_sq = 0.0;
                for (position, counts) in &per_block {
                    let n = counts[read][r].iter().sum::<u64>() as f64;
                    let expected = n * p;
                    if expected < MIN_EXPECTED_COUNT {
                        if n > 0.0 {
                            diag.blocks_sparse += 1;
                        }
                        continue;
                    }
                    let z = (counts[read][r][b] as f64 - expected) / (expected * (1.0 - p)).sqrt();
                    diag.blocks_used += 1;
                    diag.histogram[residual_bin(z)] += 1;
                    sum_sq += z * z;
                    if diag
                        .most_extreme
                        .is_none_or(|(_, worst)| z.abs() > worst.abs())
                    {
                        diag.most_extreme = Some((*position, z));
                    }
                }
                if diag.blocks_used > 1 {
                    diag.dispersion = sum_sq / (diag.blocks_used - 1) as f64;
                }
                diagnostics.push(diag);
            }
        }
    }
    diagnostics
}

/// Write each class's residual histogram beside the counts N(0, 1) would give.
pub fn write_block_diagnostics(
    diagnostics: &[ClassDiagnostics],
    path: &str,
) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "read_num\tbase_change\tpooled_rate\tblocks_used\tdispersion\tbin_low\tbin_high\tblocks\texpected_blocks"
    )?;
    for diag in diagnostics {
        for (i, &observed) in diag.histogram.iter().enumerate() {
            let (low, high) = (residual_bin_low(i), residual_bin_high(i));
            let expected = diag.blocks_used as f64 * (normal_cdf(high) - normal_cdf(low));
            writeln!(
                w,
                "{}\t{}\t{:.6e}\t{}\t{:.3}\t{}\t{}\t{}\t{:.2}",
                diag.read_num,
                diag.base_change,
                diag.pooled_rate,
                diag.blocks_used,
                diag.dispersion,
                low,
                high,
                observed,
                expected
            )?;
        }
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ReferenceOrder;

    fn key(base_change: &str, base_position: usize) -> InsertKey {
        InsertKey {
            base_change: base_change.to_string(),
            read_num: 1,
            base_position,
            reference_order: ReferenceOrder::First,
        }
    }

    fn block(ct: usize, cc: usize) -> HashMap<InsertKey, usize> {
        HashMap::from([(key("C>T", 1), ct), (key("C>C", 1), cc)])
    }

    fn frequency_intervals(
        blocks: &BlockCounts,
        replicates: usize,
        seed: u64,
        discounts: &HashMap<InsertKey, usize>,
    ) -> HashMap<InsertKey, (f64, f64)> {
        bootstrap_intervals(blocks, replicates, seed, discounts, false).frequencies
    }

    fn rate_intervals(
        blocks: &BlockCounts,
        replicates: usize,
        seed: u64,
        discounts: &HashMap<InsertKey, usize>,
    ) -> Vec<RateInterval> {
        bootstrap_intervals(blocks, replicates, seed, discounts, true)
            .rates
            .expect("rates were requested")
    }

    #[test]
    fn requesting_rates_leaves_frequency_intervals_unchanged() {
        let mut blocks = BlockCounts::new();
        for i in 0..10 {
            blocks.add_block((0, i as i64), &block(2 * i, 100 - 2 * i));
        }
        let both = bootstrap_intervals(&blocks, 200, 5, &HashMap::new(), true);
        assert_eq!(
            both.frequencies,
            frequency_intervals(&blocks, 200, 5, &HashMap::new())
        );
        assert!(!both.rates.expect("rates were requested").is_empty());
    }

    #[test]
    fn quantile_matches_sorted_interpolation() {
        let mut shuffled = [3.0f32, 10.0, 0.0, 2.0, 1.0];
        // Type 7: h = 4 * 0.975 = 3.9, between 3.0 and 10.0.
        assert!((quantile(&mut shuffled, 975) - (3.0 + 0.9 * 7.0)).abs() < 1e-9);
        assert!((quantile(&mut shuffled, 500) - 2.0).abs() < 1e-9);
        assert_eq!(quantile(&mut shuffled, 1000), 10.0);
        assert_eq!(quantile(&mut [4.0f32], 25), 4.0);
    }

    #[test]
    fn replicate_draws_are_pinned() {
        // The same values tractslip pins: seed 1, 10 blocks, replicates 0 and 7.
        let mut rngs = replicate_rngs(1, 8);
        assert_eq!(
            multiplicities(10, &mut rngs[0]),
            [0, 1, 1, 1, 1, 3, 0, 1, 2, 0]
        );
        assert_eq!(
            multiplicities(10, &mut rngs[7]),
            [1, 2, 1, 0, 0, 1, 0, 3, 2, 0]
        );
    }

    #[test]
    fn identical_blocks_give_a_zero_width_interval() {
        let mut blocks = BlockCounts::new();
        for i in 0..10 {
            blocks.add_block((0, i), &block(1, 9));
        }
        let intervals = frequency_intervals(&blocks, 200, 1, &HashMap::new());
        let (low, high) = intervals[&key("C>T", 1)];
        assert!((low - 0.1).abs() < 1e-6 && (high - 0.1).abs() < 1e-6);
    }

    #[test]
    fn interval_brackets_the_pooled_frequency_and_is_reproducible() {
        let mut blocks = BlockCounts::new();
        // Pooled C>T frequency: (0 + 2 + 4 + ... + 18) / (10 * 100) = 0.09.
        for i in 0..10 {
            blocks.add_block((0, i as i64), &block(2 * i, 100 - 2 * i));
        }
        let intervals = frequency_intervals(&blocks, 500, 7, &HashMap::new());
        let (low, high) = intervals[&key("C>T", 1)];
        assert!(
            low < 0.09 && 0.09 < high,
            "({low}, {high}) should contain 0.09"
        );
        assert!(low > 0.0 && high < 0.18);
        assert_eq!(
            intervals,
            frequency_intervals(&blocks, 500, 7, &HashMap::new())
        );
    }

    #[test]
    fn empty_blocks_are_not_resampled() {
        let mut blocks = BlockCounts::new();
        blocks.add_block((0, 0), &HashMap::new());
        blocks.add_block((0, 1), &block(1, 9));
        assert_eq!(blocks.len(), 1);
    }

    #[test]
    fn discounts_are_subtracted_before_normalizing() {
        let mut blocks = BlockCounts::new();
        for i in 0..5 {
            blocks.add_block((0, i), &block(2, 8));
        }
        // Totals: C>T 10, C>C 40. Discounting 5 C>T leaves 5 / 45.
        let discounts = HashMap::from([(key("C>T", 1), 5)]);
        let (low, high) = frequency_intervals(&blocks, 100, 1, &discounts)[&key("C>T", 1)];
        assert!((low - 5.0 / 45.0).abs() < 1e-6 && (high - 5.0 / 45.0).abs() < 1e-6);
    }

    fn class<'a>(diags: &'a [ClassDiagnostics], base_change: &str) -> &'a ClassDiagnostics {
        diags
            .iter()
            .find(|d| d.read_num == 1 && d.base_change == base_change)
            .unwrap_or_else(|| panic!("no diagnostics for {base_change}"))
    }

    #[test]
    fn homogeneous_blocks_have_zero_residuals() {
        let mut blocks = BlockCounts::new();
        for i in 0..10 {
            blocks.add_block((0, i), &block(10, 90));
        }
        let diags = block_residual_diagnostics(&blocks);
        let ct = class(&diags, "C>T");
        assert!((ct.pooled_rate - 0.1).abs() < 1e-12);
        assert_eq!(ct.blocks_used, 10);
        assert!(ct.dispersion.abs() < 1e-12);
        // Every residual is exactly 0, which falls in the [0, 0.5) bin.
        assert_eq!(ct.histogram[residual_bin(0.0)], 10);
        assert_eq!(residual_bin_low(residual_bin(0.0)), 0.0);
    }

    #[test]
    fn outlier_block_inflates_dispersion_and_is_named() {
        let mut blocks = BlockCounts::new();
        for i in 0..20 {
            blocks.add_block((0, i * 1000), &block(10, 990));
        }
        blocks.add_block((3, 5000), &block(200, 800));
        let diags = block_residual_diagnostics(&blocks);
        let ct = class(&diags, "C>T");
        assert!(ct.dispersion > 10.0, "dispersion {}", ct.dispersion);
        let (position, z) = ct.most_extreme.unwrap();
        assert_eq!(position, (3, 5000));
        assert!(z > 4.0);
        assert_eq!(ct.histogram[N_RESIDUAL_BINS - 1], 1);
    }

    #[test]
    fn sparse_blocks_are_left_out() {
        let mut blocks = BlockCounts::new();
        for i in 0..5 {
            blocks.add_block((0, i), &block(10, 90));
        }
        // Expected C>T here is 10 * 0.1 = 1, below the minimum.
        blocks.add_block((0, 99), &block(1, 9));
        let diags = block_residual_diagnostics(&blocks);
        let ct = class(&diags, "C>T");
        assert_eq!((ct.blocks_used, ct.blocks_sparse), (5, 1));
    }

    #[test]
    fn residual_bins_cover_the_line() {
        assert_eq!(residual_bin(-10.0), 0);
        assert_eq!(residual_bin(-4.0), 1);
        assert_eq!(residual_bin(3.99), N_RESIDUAL_BINS - 2);
        assert_eq!(residual_bin(4.0), N_RESIDUAL_BINS - 1);
        let total: f64 = (0..N_RESIDUAL_BINS)
            .map(|i| normal_cdf(residual_bin_high(i)) - normal_cdf(residual_bin_low(i)))
            .sum();
        assert!((total - 1.0).abs() < 1e-9);
        assert!((normal_cdf(1.96) - 0.975).abs() < 1e-4);
    }

    fn rates_for<'a>(rows: &'a [RateInterval], class: &str) -> &'a RateInterval {
        rows.iter()
            .find(|r| r.read_num == 1 && r.position == 1 && r.class == class)
            .unwrap_or_else(|| panic!("no {class} rate"))
    }

    #[test]
    fn rate_intervals_pool_reference_orders_and_fold_strands() {
        let mut second_order = key("C>T", 1);
        second_order.reference_order = ReferenceOrder::Second;
        let mut blocks = BlockCounts::new();
        for i in 0..8 {
            blocks.add_block(
                (0, i),
                &HashMap::from([
                    (key("C>T", 1), 3),
                    (second_order.clone(), 1),
                    (key("C>C", 1), 96),
                    (key("G>A", 1), 2),
                    (key("G>G", 1), 98),
                ]),
            );
        }
        let rows = rate_intervals(&blocks, 100, 1, &HashMap::new());

        // Both reference orders pool into one C>T rate: (3 + 1) / 100 per block.
        let ct = rates_for(&rows, "C>T");
        assert_eq!((ct.count, ct.ref_total), (32, 800));
        assert!((ct.rate - 0.04).abs() < 1e-12);
        // Identical blocks: every replicate matches the full-data rate.
        assert!((ct.ci_low - 0.04).abs() < 1e-6 && (ct.ci_high - 0.04).abs() < 1e-6);

        let folded = rates_for(&rows, "C>T/G>A");
        assert_eq!((folded.count, folded.ref_total), (48, 1600));

        let all = rates_for(&rows, "N>N");
        assert_eq!((all.count, all.ref_total), (48, 1600));

        // Over all bases, the classes decompose N>N: raw C>T + G>A = folded C>T/G>A = N>N here.
        assert_eq!((ct.all_bases, folded.all_bases), (1600, 1600));
        assert!((ct.rate_all - 32.0 / 1600.0).abs() < 1e-12);
        let ga = rates_for(&rows, "G>A");
        assert!((ct.rate_all + ga.rate_all - folded.rate_all).abs() < 1e-12);
        assert!((folded.rate_all - all.rate_all).abs() < 1e-12);
        assert!((all.rate_all - all.rate).abs() < 1e-12);
        assert!((ct.rate_all_ci_low - ct.rate_all).abs() < 1e-6);

        // A and T were never observed, so their classes are omitted.
        assert!(
            rows.iter()
                .all(|r| !r.class.starts_with('A') && !r.class.starts_with('T'))
        );
    }

    #[test]
    fn rate_interval_brackets_the_pooled_rate() {
        let mut blocks = BlockCounts::new();
        for i in 0..20 {
            blocks.add_block((0, i as i64), &block(i, 100 - i));
        }
        let rows = rate_intervals(&blocks, 500, 3, &HashMap::new());
        let ct = rates_for(&rows, "C>T");
        assert!(ct.ci_low < ct.rate && ct.rate < ct.ci_high, "{ct:?}");
    }
}
