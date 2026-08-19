# Remaining code review items

## Clippy — remaining warnings

`cargo clippy --workspace --all-targets` is clean (aside from an upstream `proc-macro-error2` future-incompat notice unrelated to this crate).

**Type complexity** — done via `RescalingMatrix`/`TidNameMap`/`RegionList`/`MdMismatches`/`MdMatches` type aliases in `types.rs`:
- [x] `io.rs` — `load_rescaling_matrix`, `load_rescaling_matrix_from_reader`, `write_rescaling_matrix_output` return `RescalingMatrix`
- [x] `processing.rs` — `build_tid_map_and_regions` returns `(TidNameMap, RegionList)`
- [x] `utils.rs` — `parse_md_tag` returns `(MdMismatches, MdMatches)`

**Manual flatten** (not auto-fixable):
- [x] `tasmanian-rescale-quality.rs:128` — `for record_result in bam.records() { if let Ok(mut record) = ...` → `for mut record in bam.records().flatten()`

**Collapsible if** (auto-fixed via `cargo clippy --fix`):
- [x] `io.rs`, `processing.rs` (7 sites) — collapsed nested `if let`/`if` into `if let ... && ...` let-chains

**Functions with too many arguments** — these lower-level processing functions still need the same `ProcessingConfig`/`ProcessingContext` treatment applied to the main hot path:
- [ ] `create_mismatch_key` (12 args, `processing.rs:69`)
- [ ] `compare_and_count` (15 args, `processing.rs:144`)
- [ ] `count_softclip_mismatches` (16 args, `processing.rs:334`)
- [ ] `process_overlap_region` (18 args, `processing.rs:484`)
- [ ] `process_record` (17 args, `processing.rs:636`)
- [ ] `process_single_record` (11 args, `processing.rs:813`)
- [ ] `base_position_for_mode` (8 args, `processing.rs:1214`)

## Design / refactoring

- [ ] `process_paired_reads_with_overlap` (194 lines) — three distinct concerns in one pass:
  1. build `read1/read2_overlap_map` for inconsistency detection
  2. count mismatches in the overlap (first read only)
  3. count mismatches outside the overlap (both reads)
  Split into functions

- [ ] `compare_record_to_reference` — soft-clip post-loop is a near-duplicate of the aligned-base loop body; extract a shared `emit_insert_key` helper called by both.

- [x] `write_output` and `write_normalized_output` — identical file-vs-stdout dispatch boilerplate; extract a shared writer helper.

- [ ] `frequencies_to_rescaling_matrix` — returns `scaling_factor = 1.0` for every entry (documented placeholder). Either implement the actual conversion or replace the body with `unimplemented!()` so callers know it isn't ready.

- [ ] `merge_reads_into_insert_position_mode` — confirmed unused (not called by any binary, not exercised by any test) as of 2026-08-19. Kept intentionally rather than deleted; revisit if it stays dead.
- [ ] `print_read_pair_inconsistency_table` (`io.rs:127`) — same situation: `pub`, re-exported from `lib.rs`, zero callers, zero tests, as of 2026-08-19.

## Test coverage gaps (audit 2026-08-19)

- [ ] `--plot` (`launch_visualization` in `io.rs`, invoked from `tasmanian-mismatch.rs`) shells out to `scripts/visualize.py` and has no integration test — a broken plotting script or Python/bokeh env would go undetected by `cargo test`. Full unit/integration coverage otherwise confirmed solid: 68/68 tests pass, assertions check real output values (not just exit codes), no tautological or disabled assertions found.

## Module structure

- [x] Move `ProcessingConfig` and `RegionSpec` (now `GenomicRegion`) from `processing.rs` to `types.rs`. (`ProcessingContext` intentionally kept in `processing.rs` as an implementation detail.)
- [ ] Move `print_position_table` and `print_main_output` out of `utils.rs` into `io.rs` or a new `output.rs`.
- [ ] Split position-mode math out of `processing.rs` — `insert_mode_read_position`, `read_mode_read_position`, `base_position_for_mode`, `correct_read_len_with_mode` into `positions.rs`.
- [ ] Move `src/workflow.md` and `src/control_flow.md` to `docs/`.
- [x] Replace `num_cpus::get()` in `tasmanian-rescale-quality.rs` with `std::thread::available_parallelism()` and drop the `num_cpus` dependency.
