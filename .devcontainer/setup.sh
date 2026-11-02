#!/usr/bin/env bash
# Runs ONCE when the Codespace is created. Installs micromamba and every tool
# the workshop needs into a conda environment called "pfdr".
# Participants never need to run this themselves.
set -euo pipefail

MM_VERSION=2.10.0
MM_ROOT="$HOME/micromamba"
MM_BIN="$HOME/.local/bin/micromamba"

echo ">>> Installing micromamba ${MM_VERSION}"
mkdir -p "$HOME/.local/bin"
curl -sSL "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MM_VERSION}-0.tar.bz2" \
  | tar -xj -C "$HOME/.local" bin/micromamba

export MAMBA_ROOT_PREFIX="$MM_ROOT"
echo ">>> Creating the 'pfdr' environment (takes ~3-5 min)"
"$MM_BIN" create -y -q -f environment.yml
"$MM_BIN" clean -a -y -q

# Activate the environment automatically in every new terminal
if ! grep -q "pfdr workshop" "$HOME/.bashrc"; then
cat >> "$HOME/.bashrc" <<'RC'
# >>> pfdr workshop >>>
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$HOME/.local/bin/micromamba" shell hook -s bash)"
micromamba activate pfdr
# <<< pfdr workshop <<<
RC
fi

chmod +x steps/*.sh tests/*.sh
echo ">>> Setup complete. Open a NEW terminal and start with: bash steps/00_reference.sh"
