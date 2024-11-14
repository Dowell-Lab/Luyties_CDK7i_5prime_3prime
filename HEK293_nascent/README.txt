Scripts for the metagene pipeline are in this directory.

- `metagene_pipeline.sh` calls `metagene_pipeline_annotation_windows.sbatch`
- `metagene_pipeline.sh` then calls `metagene_pipeline_5prime.sbatch` and `metagene_pipeline_3prime.sbatch`
- Those latter two scripts both call `metagene_pipeline_geneends.py` for the TTseq data or `metagene_pipeline_geneends_NETseq.py` for the NETseq data

Any custom gene lists/gene windows associated with this pipelines are present in the `annotation` directory
