import csv
import gzip
import os
from pathlib import Path


repo = Path(os.environ.get("FAP_REPO_ROOT", Path.cwd())).resolve()
data_root = Path(os.environ.get("FAP_DATA_ROOT", repo / "data")).resolve()
series = data_root / "GSE39582" / "GSE39582_series_matrix.txt.gz"
out = repo / "config" / "gse39582_sample_metadata.txt"
keys = (
    "!Sample_geo_accession",
    "!Sample_title",
    "!Sample_source_name_ch1",
    "!Sample_characteristics_ch1",
)

with gzip.open(series, "rt", encoding="utf-8", errors="replace") as src, out.open("w", encoding="utf-8") as dst:
    sample_ids = []
    datasets = []
    for line in src:
        if line.startswith(keys):
            dst.write(line)
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1"):
            values = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if values and all(x.lower().startswith("dataset:") for x in values):
                datasets = [x.split(":", 1)[1].strip() for x in values]
        if line.startswith("!series_matrix_table_begin"):
            break

if len(sample_ids) != len(datasets):
    raise RuntimeError(f"sample/dataset length mismatch: {len(sample_ids)} vs {len(datasets)}")

map_path = out.with_name("gse39582_sample_dataset.csv")
with map_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["sample", "dataset"])
    writer.writerows(zip(sample_ids, datasets))

print(out)
print(map_path)
