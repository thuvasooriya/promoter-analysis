default:
    just --list

# install dependencies
setup:
    uv sync
    mkdir -p logs

# run task 1 (ppm construction)
task-1: setup
    uv run python src/main.py --task 1

# run task 2 (statistical alignment)
task-2: setup
    uv run python src/main.py --task 2

# run task 3 (cross-validation)
task-3: setup
    uv run python src/main.py --task 3

# run complete analysis (all tasks)
all: setup
    uv run python src/main.py --task all

# generate visualizations
visualize:
    uv run python src/visualizations.py

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
