# skater-cli

`skater-cli` provides the `skater` and `tt-skater` command-line programs for preparing, running, and collecting SKaTER-seq and TT-SKaTER optimization analyses respectively.

## Installation

Clone this repository and install its two commands (`skater` and `tt-skater`) as uv-managed tools:

```bash
git clone git@github.com:TheRealGambleLab/skater-cli.git
cd skater-cli
uv tool install --no-sources .
uv tool update-shell
```

`--no-sources` uses the project's resolved dependency sources instead of local development checkouts. 
`uv tool update-shell` adds uv's tool directory to your shell `PATH`.

### Alternative: editable local component repositories

For development across SKaTER components, check out the component repositories
beside this repository and install the CLI in editable mode without
`--no-sources`. The project's
`[tool.uv.sources]` entries will use those local checkouts as editable
dependencies:

```text
parent-directory/
├── skater-cli/
├── skater-formats/
├── skater-load/
├── skater-annotate/
├── skater-equations/
├── skater-optimize/
├── tt-optimize/
├── skater-output/
└── skater-slurm/
```

From `skater-cli/`, install this layout with:

```bash
uv tool install --editable .
uv tool update-shell
```

## Configure an analysis

Start from the supplied template that matches your experiment:

```bash
cp configs/skater_config.yml my_skater.yml
# or
cp configs/ttseq_config.yml my_ttseq.yml
```

Edit the copied file to set paths to the annotation and sequencing files, output directory, experimental values, and (if applicable) cluster settings.
Comments in each template describe the available fields. 
Keep the sample order in the configuration consistent with the listed time points.

## Run commands

After installation, run the commands directly from any directory:

```bash
skater --help
tt-skater --help
```

Typical SKaTER workflow:

```bash
# Build an annotation file from the values in the configuration.
skater annotate --config my_skater.yml

# Create the genome database used by the analysis.
skater unpack --config my_skater.yml

# Load coverage for a single gene
skater load --config my_skater.yml --gene GENE_NAME

# Compile event functions required to optimize the gene
skater compile --config my_skater.yml --event EVENT_NAME

# Run optimization for gene
skater run --config my_skater.yml --gene GENE_NAME --run 1

# Write the collected results.
skater output --config my_skater.yml --out results.tsv
```

The TT-seq interface uses the same core commands, plus `pro` for generating
PRO-seq coverage:

```bash
tt-skater annotate --config my_ttseq.yml
tt-skater unpack --config my_ttseq.yml
tt-skater pro --config my_ttseq.yml --out pro_coverage.bw
tt-skater load --config my_ttseq.yml --gene GENE_NAME
tt-skater compile --config my_ttseq.yml --event EVENT_NAM
tt-skater run --config my_ttseq.yml --gene GENE_NAME --run 1
tt-skater output --config my_ttseq.yml --out results.tsv
```

## Slurm execution

Add `--slurm` to submit an individual command with `sbatch`. The wrapper
accepts resource options such as `--partition`, `--time_limit`, `--mem`,
`--cpu`, and `--log` before the subcommand:

```bash
skater --slurm --partition quick --time_limit 02:00:00 --mem 6gb --cpu 1 load --config my_skater.yml --gene GENE_NAME
```

For cluster-managed runs, use `watcher --config ...`. The `cluster` section of
the YAML file controls its user, partitions, resource limits, and number of
runs. `watcher --method cleanup` performs the corresponding cleanup mode.

```bash
skater watcher --config my_skater.yml
tt-skater watcher --config my_ttseq.yml
```

Append `--debug` to show debug logging, or `--resource_monitor` to include
process memory and CPU usage in log messages.
