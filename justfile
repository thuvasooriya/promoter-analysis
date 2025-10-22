default:
    just --list

# install dependencies
setup:
    uv sync
    mkdir -p logs

# run task 1 (ppm construction)
task-1: setup
    uv run python src/main.py --task 1

# run task 2 (statistical alignment with empirical threshold)
task-2: setup
    uv run python src/main.py --task 2 --threshold-method empirical

# run task 2 with consensus threshold
task-2-consensus: setup
    uv run python src/main.py --task 2 --threshold-method consensus

# compare both threshold methods
compare-thresholds: task-2 task-2-consensus
    uv run python src/threshold_comparison.py

# run task 3 (cross-validation)
task-3: setup
    uv run python src/main.py --task 3

# run complete analysis (all tasks)
all: setup
    uv run python src/main.py --task all

# generate visualizations
visualize:
    uv run python src/visualizations.py

report:
    cd report && typst compile gsp.typ 210657G-a1.pdf

report-w:
    cd report && typst watch gsp.typ 210657G-a1.pdf

# clean up generated files
clean:
    rm -rf results/*
    rm -rf logs/*
    rm -rf src/__pycache__/*
    rm -rf .pytest_cache

# check data structure
check-data:
    @echo "checking data structure..."
    @ls -la data/
    @echo ""
    @echo "checking student data directories:"
    @find data -maxdepth 2 -type d -name "*_GCA_*"
