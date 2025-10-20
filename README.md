# bm4322 promoter analysis

student: 210657g
genome: gca_900637025.1

## quick start

```bash
just install
just all
```

## assignment tasks

### task 1: ppm construction

- extract upstream regions (-15 to -5) from 1100 genes
- manually select 100 promoters with wawwwt pattern
- build position probability matrix with pseudocounts for c/g

### task 2: statistical alignment

- score remaining 1000 regions using ppm
- detect promoter presence/absence

### task 3: cross-validation

- test ppm on other students' 1000 samples
- compare detection rates across genomes

## project structure

```
data/                           # genome data
results/
  task1_ppm/                    # task 1 outputs
  task2_alignment/              # task 2 outputs
  task3_cross_validation/       # task 3 outputs
  figures/                      # visualizations
src/                            # source code
logs/                           # execution logs
```

## commands

- `just task-1` - run task 1 (ppm construction)
- `just task-2` - run task 2 (statistical alignment)
- `just task-3` - run task 3 (cross-validation)
- `just all` - complete analysis (all tasks)
- `just clean` - remove generated files
