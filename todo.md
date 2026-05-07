# Remaining code review items

## Clippy — remaining warnings

**Type complexity** — add a `RescalingMatrix` type alias to eliminate the repeated `HashMap<(u8, u16, char, char), f32>`:
- [ ] `io.rs:223,236,472` — complex return type in `load_rescaling_matrix`, `load_rescaling_matrix_from_reader`, `write_rescaling_matrix_output`
- [ ] `processing.rs:1249` — complex return type in `build_tid_map_and_regions`
- [ ] `utils.rs:166` — complex return type in `parse_md_tag`

**Manual flatten** (not auto-fixable):
- [x] `tasmanian-rescale-quality.rs:128` — `for record_result in bam.records() { if let Ok(mut record) = ...` → `for mut record in bam.records().flatten()`

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

- [ ] `merge_reads_into_insert_position_mode` — verify whether this is called anywhere in production paths; if not, delete it.

## Module structure

- [x] Move `ProcessingConfig` and `RegionSpec` (now `GenomicRegion`) from `processing.rs` to `types.rs`. (`ProcessingContext` intentionally kept in `processing.rs` as an implementation detail.)
- [ ] Move `print_position_table` and `print_main_output` out of `utils.rs` into `io.rs` or a new `output.rs`.
- [ ] Split position-mode math out of `processing.rs` — `insert_mode_read_position`, `read_mode_read_position`, `base_position_for_mode`, `correct_read_len_with_mode` into `positions.rs`.
- [ ] Move `src/workflow.md` and `src/control_flow.md` to `docs/`.
- [x] Replace `num_cpus::get()` in `tasmanian-rescale-quality.rs` with `std::thread::available_parallelism()` and drop the `num_cpus` dependency.
