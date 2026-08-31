# Public-data placement

Primary matrices are not committed. Download them from SRA/GEO and place the reconstructed files exactly as shown below.

## SRP361778 / Liao cohort

Create `data/liao2022/` with sparse gene-by-cell R matrices:

`onfh1.rds`, `onfh2.rds`, `onfh3.rds`, `onfh4.rds`, `onfh5.rds`, `onfh6.rds`, `hoa1.rds`, `hoa2.rds`, `hoa3.rds`, `fnf1.rds`, and `fnf2.rds`.

Each RDS must have gene symbols as row names and original cell barcodes as column names. The virtual-knockout barcode manifest expects the same original barcodes before the analysis prefix is attached.

The repository starts from these processed matrices and does not claim to reproduce the original study's Cell Ranger/SoupX preprocessing from FASTQ files. For the two matrices used by the official-R virtual-knockout release, exact byte sizes and SHA256 values are frozen in `virtual_knockout/vko_input_checksums.tsv`. Verify them with:

```powershell
python qa\check_repository.py --check-vko-data
```

An exact full-atlas rebuild requires the corresponding processed RDS files from the source-data release or authors; substituting newly reprocessed FASTQs may change cell barcodes, retained genes, and downstream values.

## GSE169396

Create `data/gse169396/` containing the GEO supplementary 10x files for GSM5201883–GSM5201886, for example:

`GSM5201883_S1_matrix.mtx.gz`, `GSM5201883_S1_barcodes.tsv.gz`, and `GSM5201883_S1_features.tsv.gz`, with the same pattern for S2–S4.

## GSE290411

Create `data/gse290411/` containing the supplementary 10x files for GSM8812280–GSM8812283. The deposited names are consumed exactly by `analysis/prepare_matrices.R`, including the absence of an underscore before `matrix`, `barcodes`, and `features`.

## GSE123568

Place these files directly under `data/`:

- `GSE123568_series_matrix.txt.gz`
- `GSE123568_family.soft.gz`

## Curated interaction/regulon resources

Place the release-matched exports used by the code directly under `data/`:

- `cellchatdb_interactions.csv`
- `cellchatdb_complex_named.csv`
- `dorothea_ABC.tsv`

`data/datasets.tsv` is the machine-readable source/role manifest. File availability can be checked before analysis with `python qa/check_repository.py --check-data`.
