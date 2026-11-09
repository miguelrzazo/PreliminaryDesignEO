#!/bin/bash
# Install the JCC wrapper and physical data without changing the system Python.
set -euo pipefail

task_root="${TFG_OREKIT_ROOT:-/private/tmp/tfg-orekit}"
if [[ "$(uname -s)" != Darwin || "$(uname -m)" != arm64 ]]; then
    echo "This bootstrap targets macOS Apple Silicon. Use environment-orekit.yml on other platforms." >&2
    exit 1
fi
mkdir -p "$task_root/bootstrap" "$task_root/data"
export MAMBA_ROOT_PREFIX="$task_root/mamba"
export CONDA_PKGS_DIRS="$task_root/packages"
micromamba="$task_root/bootstrap/bin/micromamba"
if [[ ! -x "$micromamba" ]]; then
    curl --fail --location --retry 2 --connect-timeout 20 \
        'https://micro.mamba.pm/api/micromamba/osx-arm64/latest' \
        --output "$task_root/bootstrap/micromamba.tar.bz2.part"
    mv "$task_root/bootstrap/micromamba.tar.bz2.part" "$task_root/bootstrap/micromamba.tar.bz2"
    tar -xjf "$task_root/bootstrap/micromamba.tar.bz2" -C "$task_root/bootstrap" bin/micromamba
fi

action=create
if [[ -f "$task_root/env/conda-meta/history" ]]; then
    action=install
fi
"$micromamba" "$action" --yes --prefix "$task_root/env" --override-channels \
    --channel conda-forge 'python=3.11' 'orekit=13.1' \
    'numpy>=1.26' 'scipy>=1.11' 'pandas>=2.1' 'matplotlib>=3.8' \
    'pyyaml>=6' 'pytest>=7.4' pytest-cov plotly pip

if [[ ! -s "$task_root/data/orekit-data.zip" ]]; then
    (
        cd "$task_root/data"
        "$task_root/env/bin/python" -c \
            'import orekit; orekit.initVM(); from orekit.pyhelpers import download_orekit_data_curdir; download_orekit_data_curdir()'
    )
fi
export OREKIT_DATA_PATH="$task_root/data/orekit-data.zip"
"$task_root/env/bin/python" - <<'PY'
import os
import orekit
orekit.initVM()
from orekit.pyhelpers import setup_orekit_curdir
from org.orekit.time import TimeScalesFactory
setup_orekit_curdir(os.environ["OREKIT_DATA_PATH"])
print("Orekit JVM and UTC data ready:", TimeScalesFactory.getUTC())
PY
echo
echo "Installation complete. No shell configuration was changed."
echo "Python: $task_root/env/bin/python"
echo "Data:   $OREKIT_DATA_PATH"
