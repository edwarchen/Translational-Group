# Unified TCR/BCR Basic Repertoire Pipeline Design

## Goal

Extend the existing TCR pipeline so that it can also process BCR multiplex
amplicon libraries while preserving current TCR production behavior.

The first release covers basic immune repertoire statistics only:

- CDR3 clonotype export
- CDR3 length filtering and distribution
- V/J usage statistics
- VDJTools diversity, segment usage, and spectratype analysis

The first release does not include immunoglobulin isotype analysis, somatic
hypermutation analysis, clonal lineage inference, or lineage trees.

## Experimental Model

TCR and BCR libraries are prepared and delivered as separate FASTQ inputs.
They share one analysis implementation but run in distinct receptor modes.

The current BCR multiplex panel contains `IGH`, `IGK`, and `IGL` primers in one
library. The design leaves room for additional TCR loci such as `TRA`, `TRD`,
and `TRG`, but the first release does not claim support for loci without
production primers and validation data.

## Implementation Approach

Use centralized shell condition branches for the first release.

This is intentionally narrower than a YAML-driven panel registry. It minimizes
the initial production change while keeping receptor-specific decisions in one
place. If the number of supported panels grows, these branches can later move
to a configuration file without changing the data flow.

Do not create separate TCR and BCR pipeline copies. Both modes use the same
quality control, MiXCR, export, and chain-level statistics stages.

## Command-Line Interface

Add a receptor option to the batch entry point:

```text
--receptor TCR | BCR
--chain <chain selection>
--adapter-pattern <info.csv adapter marker>
```

Allowed combinations:

| Receptor | Allowed chain selections |
| --- | --- |
| `TCR` | `TRA`, `TRB`, `BOTH` |
| `BCR` | `IGH`, `IGK`, `IGL`, `BCR_ALL` |

Expansion rules:

```text
BOTH    -> TRA TRB
BCR_ALL -> IGH IGK IGL
```

Compatibility requirements:

- `--receptor` is optional and defaults to `TCR`.
- Existing commands that only pass `--chain TRB`, `--chain TRA`, or
  `--chain BOTH` continue to work.
- BCR runs must pass `--receptor BCR` explicitly.
- Existing TCR runs retain the current default `--adapter-pattern Adapter_TCR_`.
- BCR runs must pass `--adapter-pattern` explicitly until the upstream BCR
  adapter marker is standardized.
- Invalid receptor values and invalid receptor/chain combinations fail before
  sample download or analysis begins.

Examples:

```bash
# Existing TCR behavior remains valid.
bash tcr_process_v3.sh ... --chain TRB

# BCR multiplex library: align once, export three chain-level outputs.
bash tcr_process_v3.sh ... \
  --receptor BCR \
  --chain BCR_ALL \
  --adapter-pattern "$BCR_ADAPTER_PATTERN"
```

## Primer Selection

Retain the existing pre-MiXCR V/J primer filtering behavior.

The single-sample script selects primer files by receptor:

```bash
TCR_V_PRIMERS=/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/V_primer.txt
TCR_J_RC_PRIMERS=/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_rc_primer.txt
BCR_V_PRIMERS=""
BCR_J_RC_PRIMERS=""
```

The implementation keeps the current TCR server paths as defaults. The BCR
paths are deliberately empty deployment inputs: they must be replaced with
production server absolute paths before the first BCR run. The script validates
the selected primer files before processing reads and fails with a specific
error if either value is empty or either file is missing. There is no silent
fallback from BCR primers to TCR primers.

The `Panel` column currently retained from `info.csv` is not used for routing
because its values are not yet stable. A future release may map stable business
panel names to analysis settings.

## Data Flow

The batch entry point passes `receptor` and `chain` through the task generator
to the single-sample pipeline.

```text
FASTQ
  -> fastp quality trimming and adapter removal
  -> FLASH paired-end merge
  -> receptor-specific multiplex V-primer filter
  -> receptor-specific multiplex J-primer filter
  -> MiXCR align
  -> MiXCR assemble
  -> loop over expanded chains
       -> MiXCR exportClones --chains <chain>
       -> CDR3 length filter
       -> CDR3 length distribution
       -> VDJTools format conversion
       -> V/J usage statistics
       -> VDJTools diversity, segment usage, and spectratype
```

For `--receptor BCR --chain BCR_ALL`, MiXCR alignment and assembly run once.
The export and downstream chain-level analysis loop runs for `IGH`, `IGK`, and
`IGL` separately. The three chains are never merged into one repertoire
statistic.

## Script Changes

### `scripts/tcr_process_v3.sh`

- Add `--receptor`, defaulting to `TCR`.
- Normalize receptor and chain values to uppercase.
- Validate the receptor/chain matrix.
- Preserve `Adapter_TCR_` as the TCR adapter pattern default and require an
  explicit `--adapter-pattern` for BCR until the upstream value stabilizes.
- Keep the existing `Panel` value in `filtered_info.csv`.
- Add the normalized receptor value to generated sample metadata so the task
  generator can pass it to the single-sample script.
- Include receptor information in dry-run output.

### `scripts/tcr_get_shell_fixed_primers_v3.py`

- Read the receptor value produced by the batch entry point.
- Default missing receptor values to `TCR` for compatibility with existing
  matched sample tables.
- Validate receptor/chain combinations before generating commands.
- Pass receptor and chain to the single-sample script.

### `scripts/TCR_analysis_pipeline.v7.sh`

- Accept receptor as an additional optional argument, defaulting to `TCR`.
- Centralize primer selection in one receptor `case` statement.
- Expand chain selections in one receptor-aware chain `case` statement.
- Validate selected primer files before the filtering steps.
- Preserve the existing MiXCR align and assemble behavior for this first
  release.
- Continue running the current chain-level downstream analysis for each
  expanded chain.
- For BCR runs, skip TCR-only final summary helpers with an explicit warning if
  they are not BCR-compatible. Chain-level output generation must still
  complete successfully.

### Batch QC helpers

- Keep the current TCR batch QC behavior unchanged.
- Skip TCR-only batch QC helpers for BCR mode with an explicit warning until
  they are verified or extended for `IGH`, `IGK`, and `IGL`.
- Continue generating generic fastp QC for BCR batches.

### Documentation

- Update the operation manual with TCR and BCR invocation examples.
- Document the supported receptor/chain matrix.
- State the first-release BCR scope and the excluded advanced analyses.
- Document the required deployment step for BCR primer server paths.

## Output Contract

Keep the existing chain-qualified naming convention. For example, an `IGH`
result includes:

```text
Map_Clone_Analysis/<sample>.clonotypes.IGH.raw.txt
Map_Clone_Analysis/<sample>.clonotypes.IGH.txt
Map_Clone_Analysis/convert.<sample>.clonotypes.IGH.txt
Stat_Picture/<sample>.clonotypes.IGH.cdr3.stat.txt
Stat_Picture/<sample>.IGH.VJ.stat.txt
Stat_Picture/vdjtools/<sample>.IGH.*
```

`IGK` and `IGL` use the same pattern. Existing TCR output names do not change.

## Error Handling

Fail early with clear messages for:

- unsupported receptor values
- unsupported receptor/chain combinations
- BCR runs without an explicit adapter pattern
- missing BCR primer paths or missing selected primer files
- missing single-sample pipeline scripts

Examples:

```text
ERROR: receptor must be TCR or BCR, got '...'
ERROR: chain 'IGH' is not valid for receptor 'TCR'
ERROR: --adapter-pattern is required for BCR runs
ERROR: BCR V primer file not found: ...
```

External TCR-specific summary helpers must not cause a completed set of BCR
chain-level results to be discarded. If they are not BCR-compatible, emit a
warning and skip those helpers for BCR mode.

## Test Strategy

### Batch CLI tests

Extend the existing CLI tests to cover:

- legacy TCR invocation without `--receptor`
- explicit `--receptor TCR --chain TRB`
- `--receptor BCR --chain BCR_ALL`
- `--receptor BCR --chain IGH`
- BCR invocation without `--adapter-pattern`
- invalid receptor values
- invalid combinations such as `--receptor TCR --chain IGH`

### Task generator tests

Verify generated commands:

- retain existing TCR behavior when receptor metadata is absent
- include `BCR BCR_ALL` for BCR multiplex samples
- reject invalid receptor/chain combinations

### Single-sample smoke tests

Use minimal command stubs to verify:

- TCR selects TCR primer files
- BCR selects BCR primer files
- `BCR_ALL` exports `IGH`, `IGK`, and `IGL` exactly once each
- missing BCR primer files fail before read processing
- BCR mode does not fail solely because a TCR-only final summary helper is
  skipped

## Rollout

1. Add the shell branches, validations, tests, and manual updates.
2. Replace the BCR primer path deployment inputs with production server paths.
3. Run the BCR smoke tests.
4. Validate one known BCR sample and inspect `IGH`, `IGK`, and `IGL` outputs
   independently.
5. Compare basic metrics against an expected or manually reviewed result set.
6. Enable routine BCR batch runs after validation.

## Future Extensions

Possible follow-up work:

- add stable business `Panel` aliases when the upstream field stabilizes
- migrate receptor and panel branches to a YAML registry as panels grow
- add validated `TRA`, `TRD`, and `TRG` support when primers and test data are
  available
- add optional BCR isotype, SHM, and lineage analysis as separate downstream
  modules
