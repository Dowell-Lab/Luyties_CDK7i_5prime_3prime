Scripts for three separate pipelines are in this directory.

1) 5' and 3 metagene creation
  - `metagene_pipeline.sh` calls `metagene_pipeline_annotation_windows.sbatch`
  - `metagene_pipeline.sh` then calls `metagene_pipeline_5prime.sbatch` and `metagene_pipeline_3prime.sbatch`
  - Those latter two scripts both call `metagene_pipeline_geneends.py`

2) Pause ratio scripts
  - `pause_ratio.sbatch` calls `featurecounts_pauseratio.r` for counting, then `pause_ratio.py`

3) 'Processivity,' or relative 3' elongation, was calculated as in `metagene_pipeline_processivity.sbatch` and `metagene_pipeline_processivity.py`, the first of which is also called in `metagene_pipeline.sh`

Any custom gene lists/gene windows associated with these pipelines are present in the `annotation` directory
