#!/usr/bin/env bash
set -Eeuo pipefail

# PEPATAC ONE-BUTTON INSTALLER — ROBUST ANALYSIS VERSION
# WSL Ubuntu / Linux
#
# Purpose:
#   Prepare this computer to run the PEPATAC ATAC-seq pipeline.
#
# Installs/checks:
#   - Basic Linux packages
#   - Miniconda, if missing
#   - Conda environment: FetchPA
#   - Alignment + peak-calling tools: bowtie2/bowtie2-inspect, samtools, macs3, bedtools,
#     samblaster, skewer, preseq, HOMER, ucsc-wigtobigwig, ucsc-bedtobigbed
#   - QC tools: fastqc, deeptools (bamCoverage, computeMatrix, plotHeatmap)
#   - Trimming tools: trim-galore, cutadapt
#   - PEPATAC repo (pipeline code)
#   - looper (pipeline runner)
#   - PEPATACr (R reporting package)
#   - R packages: GenomicDistributions, DiffBind, DESeq2, csaw, ChIPseeker,
#     GenomicFeatures, clusterProfiler, assembly-matched TxDb packages,
#     species-matched OrgDb packages, SQLite annotation support, ggrepel, pheatmap, etc.
#   - Rust + gtars (fragmentation scoring)
#   - Refgenie (sequence/alignment asset manager — init only; assets pulled by runner)
#   - Environment export files for reproducibility
#
# Does NOT do yet:
#   - Download/build Refgenie genome assets
#   - Download HOMER genome packages
#   - Configure per-sample looper project files
#   - Run alignment, peak-calling, or reporting
#
# Those belong in the PEPATAC runner so genome and sample choices are
# confirmed only after the user reviews a run configuration.

# USER SETTINGS

INSTALLER_VERSION="1.28-annotation-integrity-audit-deps"
ENV_NAME="FetchPA"
PEPATAC_DIR="$HOME/pepatac"
REFGENIE_DIR="$HOME/refgenie"
REFGENIE_CONFIG="$REFGENIE_DIR/refgenie.yaml"

INSTALL_LOG="$HOME/pepatac_install_$(date +%Y%m%d_%H%M%S).log"

# Detect whether we're attached to a real interactive terminal. Used only
# to decide whether to emit color codes.
IS_TTY=0
if [[ -t 1 ]]; then
    IS_TTY=1
fi

# Mirror all installer output to a log file while still showing it live.
exec > >(tee -a "$INSTALL_LOG") 2>&1

LAST_SECTION="startup"
trap 'echo ""; echo "  [ERROR] Installer failed during: $LAST_SECTION"; echo "  [ERROR] See log: $INSTALL_LOG"' ERR

# Helpers

# Minimal color helpers. No-ops when not attached to a terminal.
if [[ "$IS_TTY" -eq 1 ]]; then
    C_RESET=$'\033[0m'
    C_DIM=$'\033[2m'
    C_BOLD=$'\033[1m'
    C_CYAN=$'\033[36m'
    C_GREEN=$'\033[32m'
    C_YELLOW=$'\033[33m'
    C_MAGENTA=$'\033[35m'
else
    C_RESET=""; C_DIM=""; C_BOLD=""; C_CYAN=""; C_GREEN=""; C_YELLOW=""; C_MAGENTA=""
fi

# section TITLE
# Prints a plain banner.
section() {
    LAST_SECTION="$1"
    local title="$1"
    echo ""
    echo "${C_CYAN}============================================================${C_RESET}"
    echo "${C_BOLD}  $title${C_RESET}"
    echo "${C_CYAN}============================================================${C_RESET}"
}

# run_with_spinner LABEL -- CMD [ARGS...]
# Runs a command, printing plain [INFO]/[OK]/[WARN] status lines. On
# failure, shows the last bit of the command's output to help diagnose
# the problem.
run_with_spinner() {
    local label="$1"; shift
    [[ "${1:-}" == "--" ]] && shift

    info "$label"

    local out_file
    out_file="$(mktemp)"

    local status=0
    if "$@" >"$out_file" 2>&1; then
        status=0
    else
        status=$?
    fi

    if [[ $status -eq 0 ]]; then
        ok "$label"
    else
        warn "$label (failed — see details below)"
        echo "  ---- last output ----"
        tail -n 25 "$out_file" || true
        echo "  ----------------------"
    fi

    rm -f "$out_file"
    return $status
}

info()    { echo "  [INFO]  $*"; }
ok()      { echo "  [OK]    $*"; }
warn()    { echo "  [WARN]  $*"; }
err()     { echo "  [ERROR] $*"; }
die()     { echo "  [ERROR] $*"; exit 1; }

command_exists() { command -v "$1" >/dev/null 2>&1; }

# pip_install_self_healing LABEL -- CMD [ARGS...]
# Like run_with_spinner, but specifically for `pip install` commands. Some
# environments end up with a package's .dist-info directory present but its
# METADATA file missing — usually from an interrupted previous install, a
# full disk, or antivirus/indexing software touching files mid-write (all
# more likely on WSL, where the conda env can live on an NTFS-backed mount).
# pip then raises a raw OSError ("No such file or directory: .../METADATA")
# instead of a normal, retryable error. If that exact signature is detected,
# force-reinstall just the broken package(s) and retry the original install
# once before giving up.
pip_install_self_healing() {
    local label="$1"; shift
    [[ "${1:-}" == "--" ]] && shift

    info "$label"

    local out_file
    out_file="$(mktemp)"

    local status=0
    if "$@" >"$out_file" 2>&1; then
        status=0
    else
        status=$?
    fi

    if [[ $status -eq 0 ]]; then
        ok "$label"
        rm -f "$out_file"
        return 0
    fi

    # Extract package name(s) from either of two symptoms of the same
    # underlying issue (an incompletely-written dist-info directory):
    #   1. "No such file or directory: .../<pkg>-<version>.dist-info/METADATA"
    #   2. "no RECORD file was found for <pkg>" (surfaces if something
    #      upstream already attempted --force-reinstall and hit the
    #      uninstall step first)
    local broken_pkgs
    broken_pkgs="$(
        {
            grep -oE '[A-Za-z0-9_.+-]+-[0-9][A-Za-z0-9_.+-]*\.dist-info/METADATA' "$out_file" \
                | sed -E 's/-[0-9][A-Za-z0-9_.+-]*\.dist-info\/METADATA$//'
            grep -oE 'no RECORD file was found for [A-Za-z0-9_.+-]+' "$out_file" \
                | awk '{print $NF}'
        } | sort -u || true
    )"

    if [[ -n "$broken_pkgs" ]]; then
        warn "$label (failed — detected corrupted package metadata, attempting repair)"
        local pkg repaired_all=1

        # Purge pip's wheel/build cache first: a stale cached artifact could
        # itself be the corrupted one (e.g. from an earlier interrupted
        # download), and reinstalling from it would just recreate the
        # problem.
        in_env python -m pip cache purge >/dev/null 2>&1 || true

        # Resolve the env's site-packages directory so we can physically
        # remove the broken dist-info folder. --ignore-installed alone is
        # not enough: it skips pip's uninstall step (which needs a RECORD
        # file we don't have) but leaves the old, broken dist-info directory
        # sitting on disk. pip then finds two dist-info dirs for the same
        # package on the next run and can pick the broken one again,
        # reproducing the exact same failure even after a "successful" repair.
        local site_packages_dir
        site_packages_dir="$(in_env python -c \
            'import sysconfig; print(sysconfig.get_paths()["purelib"])' 2>/dev/null || true)"

        for pkg in $broken_pkgs; do
            info "Repairing corrupted package: $pkg"

            if [[ -n "$site_packages_dir" && -d "$site_packages_dir" ]]; then
                local stale_dir
                for stale_dir in "$site_packages_dir/${pkg}"-*.dist-info \
                                  "$site_packages_dir/${pkg//-/_}"-*.dist-info; do
                    [[ -d "$stale_dir" ]] || continue
                    info "Removing stale metadata directory: $stale_dir"
                    rm -rf "$stale_dir"
                done
            else
                warn "Could not resolve site-packages for '$ENV_NAME' — skipping stale directory cleanup for $pkg."
            fi

            local repair_out
            repair_out="$(mktemp)"
            if in_env python -m pip install --no-deps --no-cache-dir "$pkg" >"$repair_out" 2>&1; then
                ok "Repaired: $pkg"
            else
                warn "Could not repair $pkg automatically."
                echo "  ---- repair attempt output for $pkg ----"
                tail -n 20 "$repair_out" || true
                echo "  ------------------------------------------"
                repaired_all=0
            fi
            rm -f "$repair_out"
        done

        info "Retrying: $label"
        local retry_out
        retry_out="$(mktemp)"
        if "$@" >"$retry_out" 2>&1; then
            ok "$label (succeeded after repair)"
            rm -f "$out_file" "$retry_out"
            return 0
        else
            status=$?
            warn "$label (still failing after repair attempt)"
            echo "  ---- last output ----"
            tail -n 25 "$retry_out" || true
            echo "  ----------------------"
            if [[ "$repaired_all" -eq 0 ]]; then
                warn "One or more packages could not be repaired automatically."
            fi
            warn "If this persists, the '$ENV_NAME' environment may need a full rebuild:"
            warn "    conda env remove -n $ENV_NAME   # then re-run this installer"
            warn "Also worth checking free disk space: df -h \"\$HOME\""
            rm -f "$retry_out"
        fi
    else
        warn "$label (failed — see details below)"
        echo "  ---- last output ----"
        tail -n 25 "$out_file" || true
        echo "  ----------------------"
    fi

    rm -f "$out_file"
    return $status
}

# ─────────────────────────────────────────────────────────────
# Path normalisation  (shared with run.sh / diff_analysis.sh / explore.sh)
# ─────────────────────────────────────────────────────────────

# normalize_path VALUE
# Accepts paths pasted from either Linux/WSL or Windows File Explorer and
# returns a Linux path suitable for tools running inside WSL.
# Supported examples:
#   /home/dustin/project
#   ~/project
#   C:\\Users\\Dustin\\project      -> /mnt/c/Users/Dustin/project
#   F:/ATACseq/project             -> /mnt/f/ATACseq/project
#   \\\\wsl.localhost\\Ubuntu\\home\\dustin\\project -> /home/dustin/project
#   \\\\wsl$\\Ubuntu\\home\\dustin\\project         -> /home/dustin/project
# Surrounding quotes from Explorer's "Copy as path" are removed.
# Unrecognised Windows-style paths (e.g. UNC \\server\share) are returned
# unchanged and a warning is printed to stderr so the caller can surface it.
normalize_path() {
    local input="${1:-}"
    local slashified rest converted drive

    # Remove a pasted Windows carriage return, if present.
    input="${input//$'\r'/}"

    # Strip one matching pair of surrounding single or double quotes.
    if [[ ${#input} -ge 2 ]]; then
        if [[ "${input:0:1}" == '"' && "${input: -1}" == '"' ]]; then
            input="${input:1:${#input}-2}"
        elif [[ "${input:0:1}" == "'" && "${input: -1}" == "'" ]]; then
            input="${input:1:${#input}-2}"
        fi
    fi

    # Expand a Linux home shortcut.
    case "$input" in
        "~")  input="$HOME" ;;
        "~/"*) input="$HOME/${input:2}" ;;
    esac

    # Convert backslashes only for recognizing Windows-style paths.
    slashified="${input//\\//}"

    # Windows Explorer view of a WSL distribution:
    # //wsl.localhost/Ubuntu/home/dustin/... or //wsl$/Ubuntu/home/dustin/...
    case "${slashified,,}" in
        //wsl.localhost/*|//wsl\$/*)
            rest="${slashified#//*/}"  # remove //wsl.../
            if [[ "$rest" == */* ]]; then
                rest="${rest#*/}"      # remove distribution name
                input="/${rest}"
            else
                input="/"
            fi
            ;;
        *)
            # Windows drive path. Prefer WSL's converter, with a portable
            # /mnt/<drive> fallback for environments where wslpath is absent.
            if [[ "$slashified" =~ ^([A-Za-z]):/(.*)$ ]]; then
                converted=""
                if command -v wslpath >/dev/null 2>&1; then
                    converted="$(wslpath -u "$input" 2>/dev/null || true)"
                fi
                if [[ -n "$converted" ]]; then
                    input="$converted"
                else
                    drive="${BASH_REMATCH[1],,}"
                    input="/mnt/${drive}/${BASH_REMATCH[2]}"
                fi
            # Looks like a Windows-style path (starts with // or \\ after
            # backslash conversion) but didn't match any known pattern.
            # UNC shares (\\server\share), unusual prefixes, etc. land here.
            # Return the input unchanged and warn so the user isn't left
            # wondering why a later "directory not found" error appeared.
            elif [[ "$slashified" == //* ]]; then
                echo "  ⚠  normalize_path: unrecognised Windows-style path — returning as-is: $input" >&2
                echo "     Expected formats: C:\\\\path, F:/path, \\\\\\\\wsl.localhost\\\\Distro\\\\path" >&2
            fi
            ;;
    esac

    printf '%s\n' "$input"
}

# normalize_path_var VARIABLE_NAME
# In-place wrapper used immediately after interactive path prompts.
normalize_path_var() {
    local -n _path_ref="$1"
    _path_ref="$(normalize_path "$_path_ref")"
}

safe_append_bashrc() {
    local line="$1"
    grep -Fxq "$line" "$HOME/.bashrc" 2>/dev/null || echo "$line" >> "$HOME/.bashrc"
}

# Run a command inside the conda env without needing conda activate.
in_env() {
    conda run --no-capture-output -n "$ENV_NAME" "$@"
}

# Same idea, but with inherited R settings removed. This protects Bioconductor
# installs from weird host/session variables like R_ARCH.
in_env_clean() {
    conda run --no-capture-output -n "$ENV_NAME" \
        env -u R_ARCH -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$@"
}

# Install packages into the env from the base conda process.
# Important: do NOT run `conda install` through `conda run`; that creates
# fragile nested-conda behavior during package post-link scripts.
conda_install_env() {
    conda install -n "$ENV_NAME" \
        --override-channels \
        -c conda-forge \
        -c bioconda \
        "$@" \
        -y
}

verify_r_packages() {
    local label="$1"
    shift
    local pkg_list=("$@")
    in_env_clean Rscript --vanilla - "${pkg_list[@]}" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
missing <- args[!vapply(args, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
}
cat("    R packages OK:", paste(args, collapse = ", "), "\n")
REOF
}

# ------------------------------------------------------------
# Pinned tool versions
# ------------------------------------------------------------
# Known-good versions for the standalone CLI tools that this installer/
# pipeline has been built and tested against. These are intentionally
# pinned individually since each is a standalone tool with no
# cross-package version coupling.
#
# Bioconductor packages (DiffBind, DESeq2, clusterProfiler, etc.) are generally
# NOT pinned individually here: Bioconductor releases ship as a bundle of
# interdependent package versions tied to a specific R version, so arbitrary
# one-package pins can create solver conflicts. csaw is the deliberate exception
# below because DiffBind's background-bin normalization requires it at runtime,
# and this pipeline is fixed to R 4.3 / Bioconductor 3.18. Pinning the matching
# csaw release prevents a future installer run from silently selecting a build
# for a newer R/Bioconductor generation.
#
# trim-galore is deliberately pinned to the last Perl-based release
# (0.6.11) rather than the newer Rust rewrite (v2.x). The Rust rewrite is
# billed as a drop-in replacement, but it's recent and this installer
# hasn't been validated against it yet.
declare -A PINNED_VERSIONS
PINNED_VERSIONS["bowtie2"]="2.5.4"
PINNED_VERSIONS["samtools"]="1.21"     # 1.21 satisfies libzlib >=1.3.1 and is compatible with preseq
PINNED_VERSIONS["bedtools"]="2.31.1"
PINNED_VERSIONS["fastqc"]="0.12.1"
PINNED_VERSIONS["macs3"]="3.0.2"
PINNED_VERSIONS["deeptools"]="3.5.6"
# preseq is intentionally NOT pinned — its zlib requirement conflicts with
# pinned samtools >=1.23. Letting conda resolve preseq freely alongside the
# samtools pin avoids the LibMambaUnsatisfiableError seen with 1.23.1+preseq=3.2.0.
PINNED_VERSIONS["samblaster"]="0.1.26"
PINNED_VERSIONS["skewer"]="0.2.2"
PINNED_VERSIONS["trim-galore"]="0.6.11"
PINNED_VERSIONS["cutadapt"]="5.2"

# DiffBind calls csaw when dba.normalize(background=TRUE) generates genome-wide
# background bins. R 4.3 maps to Bioconductor 3.18, whose csaw release is 1.36.0.
# Keep these two values synchronized when deliberately moving the R/Bioc stack.
CSAW_VERSION="1.36.0"
CSAW_BIOC_VERSION="3.18"

# VERSION MODE — Ask once before anything else runs
# pinned     → use every pin as-is, no conda search, no prompts (fastest, safest)
# interactive → check conda for each tool and ask permission to upgrade (original behaviour)
#
# Non-interactive sessions (piped, CI, no TTY) always behave as pinned.

VERSION_MODE="pinned"   # default; overridden below when interactive

if [[ "$IS_TTY" -eq 1 ]]; then
    echo ""
    echo "${C_CYAN}============================================================${C_RESET}"
    echo "${C_BOLD}  Tool Version Mode${C_RESET}"
    echo "${C_CYAN}============================================================${C_RESET}"
    echo ""
    echo "  How should tool versions be chosen?"
    echo ""
    echo "    ${C_BOLD}1.  Use pinned versions${C_RESET}  (default — recommended)"
    echo "        Installs pinned versions of the core CLI tools this pipeline"
    echo "        was built and tested against. A few tools (preseq, HOMER,"
    echo "        refgenie) and the R/Bioconductor stack are intentionally left"
    echo "        unpinned -- see PINNED_VERSIONS comments above. Fastest, no"
    echo "        extra network lookups, no prompts."
    echo ""
    echo "    ${C_BOLD}2.  Let me decide${C_RESET}"
    echo "        Checks conda for the latest available version of each tool."
    echo "        If a newer version exists you will be asked whether to use it."
    echo "        Useful if a pin is known to be broken or you want cutting-edge"
    echo "        tools — but newer versions are not tested against this pipeline."
    echo ""
    echo "  ${C_DIM}Pinned versions:${C_RESET}"
    for _tool in bowtie2 samtools bedtools fastqc macs3 deeptools samblaster skewer trim-galore cutadapt; do
        printf "    ${C_DIM}%-20s %s${C_RESET}\n" "$_tool" "${PINNED_VERSIONS[$_tool]}"
    done
    printf "    ${C_DIM}%-20s %s${C_RESET}\n" "preseq" "(unpinned — conda resolves)"
    printf "    ${C_DIM}%-20s %s${C_RESET}\n" "csaw (R/Bioc)" "${CSAW_VERSION} (Bioconductor ${CSAW_BIOC_VERSION})"
    echo ""

    while true; do
        read -r -p "  Choice [1/2, default: 1]: " _ver_choice
        _ver_choice="${_ver_choice:-1}"
        case "$_ver_choice" in
            1)
                VERSION_MODE="pinned"
                ok "Using pinned versions."
                break
                ;;
            2)
                VERSION_MODE="interactive"
                ok "Will check for newer versions and ask per tool."
                break
                ;;
            *)
                warn "Enter 1 or 2."
                ;;
        esac
    done
    unset _tool _ver_choice
else
    info "Non-interactive session — using pinned versions silently."
fi

# resolve_tool_version TOOL
# In pinned mode:      returns the pin immediately, no network call, no prompt.
# In interactive mode: checks conda for the latest, asks if a newer version
#                      exists, remembers the answer for the rest of this run.
#                      Falls back to pin silently in non-TTY sessions.
declare -A RESOLVED_VERSIONS
resolve_tool_version() {
    local tool="$1"
    local pinned="${PINNED_VERSIONS[$tool]:-}"

    # Cache hit — already resolved this session.
    if [[ -n "${RESOLVED_VERSIONS[$tool]:-}" ]]; then
        echo "${RESOLVED_VERSIONS[$tool]}"
        return 0
    fi

    # Pinned mode: return immediately, no conda search needed.
    if [[ "$VERSION_MODE" == "pinned" ]]; then
        RESOLVED_VERSIONS["$tool"]="$pinned"
        echo "$pinned"
        return 0
    fi

    # Interactive mode: look up the latest available version.
    local latest=""
    latest="$(conda search --override-channels -c conda-forge -c bioconda "$tool" 2>/dev/null \
        | awk -v t="$tool" '$1==t {print $2}' \
        | sort -V \
        | tail -n 1)"

    if [[ -z "$latest" || "$latest" == "$pinned" ]]; then
        RESOLVED_VERSIONS["$tool"]="$pinned"
        echo "$pinned"
        return 0
    fi

    # A newer version exists — ask the user.
    warn "$tool: newer version available  (installed pin: ${pinned}  →  latest: ${latest})" >&2
    local reply
    read -r -p "  Use ${latest} instead of the pinned ${pinned}? [y/N]: " reply
    if [[ "${reply,,}" == "y" ]]; then
        ok "  $tool → $latest" >&2
        RESOLVED_VERSIONS["$tool"]="$latest"
        echo "$latest"
    else
        ok "  $tool → $pinned (keeping pin)" >&2
        RESOLVED_VERSIONS["$tool"]="$pinned"
        echo "$pinned"
    fi
}

# STEP 0 — Disk space check

section "STEP 0 — Disk space check"

# Rough budget for a brand-new machine:
#   Miniconda itself                ~0.5 GB
#   FetchPA conda env                ~5-8 GB (many bioinformatics tools + R/Bioc)
#   PEPATAC repo + gtars build      ~0.5 GB
#   Project folders + QC output     ~1 GB (grows a lot once real samples are run)
# We require 10 GB free on $HOME as a conservative floor for *this installer*.
# Refgenie genome assets and HOMER genome packages (downloaded later by the
# runner) need considerably more and are NOT covered by this check.
MIN_FREE_GB=10

AVAILABLE_KB=$(df -Pk "$HOME" | awk 'NR==2 {print $4}')
AVAILABLE_GB=$(( AVAILABLE_KB / 1024 / 1024 ))

info "Free space on $HOME: ${AVAILABLE_GB} GB (need at least ${MIN_FREE_GB} GB for this installer)"

if (( AVAILABLE_GB < MIN_FREE_GB )); then
    die "Not enough free disk space. Found ${AVAILABLE_GB} GB free on $HOME, need at least ${MIN_FREE_GB} GB. Free up space and re-run this script. (Note: Refgenie genome assets and HOMER genomes downloaded later will need substantially more space than this installer alone.)"
else
    ok "Enough free disk space to continue."
fi

# STEP 1 — Basic system packages

section "STEP 1 — System packages"

if command_exists apt-get; then
    info "About to install/check these system packages via apt-get (requires sudo):"
    cat <<'PKGLIST'
    build-essential   git              wget             curl
    unzip             ca-certificates  bzip2            gzip
    fontconfig        perl             pkg-config       libxml2-dev
    libuv1-dev        libssl-dev       libcurl4-openssl-dev  zlib1g-dev
    libbz2-dev        liblzma-dev
PKGLIST
    info "These are standard build tools and libraries needed to compile R/Bioconductor packages."

    # Refresh sudo credentials BEFORE backgrounding anything, so a password
    # prompt (if needed) appears normally instead of hiding behind the spinner.
    sudo -v

    run_with_spinner "Updating apt package lists" -- sudo apt-get update -qq
    run_with_spinner "Installing system packages" -- sudo apt-get install -y \
        build-essential \
        coreutils \
        util-linux \
        git \
        wget \
        curl \
        unzip \
        ca-certificates \
        bzip2 \
        gzip \
        fontconfig \
        perl \
        pkg-config \
        libxml2-dev \
        libuv1-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev
else
    warn "apt-get not found — assuming system packages are already present."
fi

ok "System package step complete."

# STEP 2 — Miniconda

section "STEP 2 — Miniconda"

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"

# Pinned to the exact installer build this installer has been run against,
# instead of "Miniconda3-latest-Linux-x86_64.sh". That filename is a moving
# target -- repo.anaconda.com overwrites its contents (and thus its
# checksum) every time a new Miniconda point release ships, so a checksum
# pinned against "-latest" breaks on the very next release. Anaconda keeps
# every dated build permanently available under its own versioned
# filename, so pinning to one of those stays downloadable long-term.
# Verified directly against https://repo.anaconda.com/miniconda/ 2026-07-30.
# Update both of these together (from that same index page) if this needs
# to move to a newer build.
MINICONDA_INSTALLER="Miniconda3-py314_26.5.3-1-Linux-x86_64.sh"
MINICONDA_SHA256="42cfece170da342364a78d629e06b94dfd81b0f2717d7655729100d888d606b4"

if [[ -d "$HOME/miniconda3" ]]; then
    ok "Conda already appears to be installed at $HOME/miniconda3."
elif command_exists conda; then
    # A conda is on PATH, but not at the one location this installer (and
    # every other PEPATAC script -- run.sh, diff_analysis.sh, explorer.sh)
    # hardcodes. Letting this fall through as "already installed" was the
    # old behavior; it then died confusingly a few lines below at the
    # conda.sh sourcing step instead of here, with no explanation of why.
    # Failing fast with the actual base path and a concrete fix is much
    # less confusing for a first-time user with an existing Miniforge/
    # Anaconda/system Conda install than a bare "file not found" later.
    EXISTING_CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    die "Found a 'conda' command on PATH (base: ${EXISTING_CONDA_BASE:-unknown}), but this installer and every other PEPATAC script here only ever look for Conda at \$HOME/miniconda3 -- they will not find or use a Conda installed elsewhere, even though this check alone would have let you proceed. Fix with one of: (1) symlink your existing install so the pipeline can find it: ln -s \"${EXISTING_CONDA_BASE:-/path/to/your/conda}\" \"$HOME/miniconda3\" -- then re-run this script; or (2) remove/rename the existing 'conda' from PATH and re-run so this installer sets up its own Miniconda at $HOME/miniconda3 instead."
else
    ARCH="$(uname -m)"
    if [[ "$ARCH" != "x86_64" ]]; then
        die "Detected CPU architecture '$ARCH', but this installer downloads the x86_64 build of Miniconda. If you're on Windows-on-ARM (e.g. Surface Pro X/11), this installer won't work as written. Get the matching Miniconda installer for $ARCH from https://repo.anaconda.com/miniconda/ and install it manually to $HOME/miniconda3, then re-run this script."
    fi
    info "Installing Miniconda ($MINICONDA_INSTALLER)..."
    cd "$HOME"
    run_with_spinner "Downloading Miniconda installer" -- \
        wget -q -O "$MINICONDA_INSTALLER" \
        "https://repo.anaconda.com/miniconda/$MINICONDA_INSTALLER"

    ACTUAL_SHA256="$(sha256sum "$MINICONDA_INSTALLER" | awk '{print $1}')"
    if [[ "$ACTUAL_SHA256" != "$MINICONDA_SHA256" ]]; then
        rm -f "$MINICONDA_INSTALLER"
        die "Miniconda installer checksum mismatch. Expected $MINICONDA_SHA256, got $ACTUAL_SHA256. Refusing to run an unverified installer -- re-download, or check https://repo.anaconda.com/miniconda/ for a current build and update MINICONDA_INSTALLER/MINICONDA_SHA256 above."
    fi
    ok "Checksum verified."

    run_with_spinner "Running Miniconda installer" -- \
        bash "$MINICONDA_INSTALLER" -b -p "$HOME/miniconda3"
    "$HOME/miniconda3/bin/conda" init bash
    ok "Miniconda installed."
fi

if [[ -f "$CONDA_SH" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_SH"
else
    die "Cannot find Conda startup script at: $CONDA_SH"
fi

unset R_ARCH R_LIBS R_LIBS_USER R_LIBS_SITE || true

# STEP 3 — Conda Terms of Service / config

section "STEP 3 — Conda Terms of Service and config"

accept_conda_tos_channel() {
    local channel="$1"

    if conda tos accept --override-channels --channel "$channel" >/dev/null 2>&1; then
        ok "Accepted Conda Terms of Service for: $channel"
    else
        warn "Could not auto-accept Conda Terms of Service for: $channel"
        warn "If Conda fails later with CondaToSNonInteractiveError, run:"
        warn "  conda tos accept --override-channels --channel $channel"
    fi
}

if conda tos accept --help >/dev/null 2>&1; then
    accept_conda_tos_channel "https://repo.anaconda.com/pkgs/main"
    accept_conda_tos_channel "https://repo.anaconda.com/pkgs/r"
else
    ok "This Conda version does not require/use 'conda tos accept'."
fi

conda config --set channel_priority strict >/dev/null 2>&1 || true
conda config --set show_channel_urls true >/dev/null 2>&1 || true

ok "Conda configured."

# STEP 4 — Conda environment

section "STEP 4 — Conda environment: $ENV_NAME"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    ok "Environment '$ENV_NAME' already exists."
else
    info "Creating environment '$ENV_NAME'..."
    run_with_spinner "Creating conda environment '$ENV_NAME'" -- \
        conda create -n "$ENV_NAME" \
        --override-channels \
        -c conda-forge \
        -c bioconda \
        python=3.10 \
        pip \
        r-base=4.3 \
        -y
    ok "Environment created."
fi

# STEP 5 — PEPATAC bioinformatics CLI tools

section "STEP 5 — PEPATAC bioinformatics CLI tools"

info "About to install these tools into the '$ENV_NAME' conda environment:"
cat <<'TOOLLIST'
    bowtie2          - short-read aligner
    samtools         - BAM filtering, sorting, indexing
    bedtools         - genome arithmetic (peaks, coverage)
    fastqc           - pre/post-trim quality control
    deeptools        - bamCoverage, computeMatrix, plotHeatmap
    macs3            - ATAC peak calling
    preseq           - library complexity estimation
    samblaster       - duplicate marking in stream
    skewer           - adapter trimming (ATAC-specific)
    trim-galore      - trimming wrapper
    cutadapt         - trimming backend used by Trim Galore
    ucsc-wigtobigwig - wig → bigWig conversion
    ucsc-bedtobigbed - bed → bigBed conversion
    homer            - motif enrichment analysis
    refgenie         - reference genome manager
    refgenconf       - refgenie config library
    pigz             - faster gzip-compatible compression/decompression
    wget, curl, unzip, tree, dos2unix - general utilities
    r-base=4.3       - R language runtime (pinned version)
    pip              - Python package installer
TOOLLIST

if [[ "$VERSION_MODE" == "pinned" ]]; then
    info "Using pinned versions (selected at startup)."
else
    info "Checking latest available versions against pins (interactive mode)..."
fi

BOWTIE2_VER="$(resolve_tool_version bowtie2)";       info "  bowtie2: $BOWTIE2_VER"
SAMTOOLS_VER="$(resolve_tool_version samtools)";     info "  samtools: $SAMTOOLS_VER"
BEDTOOLS_VER="$(resolve_tool_version bedtools)";     info "  bedtools: $BEDTOOLS_VER"
FASTQC_VER="$(resolve_tool_version fastqc)";         info "  fastqc: $FASTQC_VER"
MACS3_VER="$(resolve_tool_version macs3)";           info "  macs3: $MACS3_VER"
DEEPTOOLS_VER="$(resolve_tool_version deeptools)";   info "  deeptools: $DEEPTOOLS_VER"
SAMBLASTER_VER="$(resolve_tool_version samblaster)"; info "  samblaster: $SAMBLASTER_VER"
info "  preseq: (unpinned — conda will resolve compatible version)"
SKEWER_VER="$(resolve_tool_version skewer)";         info "  skewer: $SKEWER_VER"
TRIMGALORE_VER="$(resolve_tool_version trim-galore)"; info "  trim-galore: $TRIMGALORE_VER"
CUTADAPT_VER="$(resolve_tool_version cutadapt)";     info "  cutadapt: $CUTADAPT_VER"

if ! run_with_spinner "Installing PEPATAC bioinformatics CLI tools" -- conda_install_env \
    bowtie2="$BOWTIE2_VER" \
    samtools="$SAMTOOLS_VER" \
    bedtools="$BEDTOOLS_VER" \
    fastqc="$FASTQC_VER" \
    macs3="$MACS3_VER" \
    deeptools="$DEEPTOOLS_VER" \
    preseq \
    samblaster="$SAMBLASTER_VER" \
    skewer="$SKEWER_VER" \
    trim-galore="$TRIMGALORE_VER" \
    cutadapt="$CUTADAPT_VER" \
    ucsc-wigtobigwig \
    ucsc-bedtobigbed \
    homer \
    refgenie \
    refgenconf \
    pigz \
    wget \
    curl \
    unzip \
    tree \
    dos2unix \
    r-base=4.3 \
    pip
then
    warn "Installing PEPATAC bioinformatics CLI tools failed."
    warn "The most common cause is a zlib/libzlib version conflict between"
    warn "samtools and preseq. This installer unpins preseq to avoid that,"
    warn "but if the error persists check the conda output above for the"
    warn "conflicting package and report it with the install log:"
    warn "  $INSTALL_LOG"
    die "Stopping: PEPATAC bioinformatics CLI tools were not fully installed."
fi

ok "PEPATAC bioinformatics CLI tools installed."

# STEP 6 — R / Bioconductor packages (core PEPATAC reporting)

section "STEP 6 — R packages: core PEPATAC reporting"

# Install binary conda R packages first (faster, avoids source compilation).
# Bioconductor data packages (GenomicDistributions, GenomicDistributionsData)
# are installed separately because their post-link scripts have historically
# failed on some machines; we use a BiocManager fallback for those.

info "Installing binary R helper packages via conda..."
info "About to install these R packages (used by PEPATACr reports):"
cat <<'CRANLIST'
    r-remotes             - install R packages from GitHub/local
    r-biocmanager         - installs Bioconductor packages from R
    r-matrix, r-fs        - common utilities used by Bioconductor
    r-xml, r-xml2         - XML parsing (used by many Bioc packages)
    r-bslib, r-rmarkdown  - report rendering
    r-htmlwidgets         - interactive HTML output
    r-shiny, r-dt         - interactive tables in PEPATAC reports
    r-gplots              - classic R plotting (used by PEPATAC heatmaps)
    r-pheatmap            - correlation and clustered heatmaps used by explorer tools
CRANLIST
run_with_spinner "Installing binary R helper packages" -- conda_install_env \
    r-remotes \
    r-gplots \
    r-pheatmap \
    r-biocmanager \
    r-matrix \
    r-fs \
    r-xml \
    r-xml2 \
    r-bslib \
    r-rmarkdown \
    r-htmlwidgets \
    r-shiny \
    r-dt \
    libxml2 \
    pkg-config \
    libuv

# Verify the binary R dependencies landed before proceeding.
verify_r_packages "conda R dependency install" \
    Matrix fs XML xml2 BiocManager remotes pheatmap

# CRAN helpers used by PEPATAC/PEPATACr (not available as conda binaries).
# optigrab has no CRAN release, so it's installed from GitHub -- pinned to
# a specific commit (the repo's HEAD as of 2026-07-30) rather than always
# tracking whatever the default branch currently points to, so re-running
# this installer later doesn't silently pick up unreviewed upstream changes.
OPTIGRAB_COMMIT="5166b8bad103ec4ef95114eb07101a511c39a1a4"
info "Installing CRAN helpers (pepr, optigrab)..."
in_env_clean Rscript --vanilla - "$OPTIGRAB_COMMIT" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
optigrab_commit <- args[1]
options(repos=c(CRAN="https://cloud.r-project.org"))
for (pkg in c("pepr")) {
    if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
if (!requireNamespace("optigrab", quietly=TRUE)) {
    remotes::install_github(paste0("decisionpatterns/optigrab@", optigrab_commit))
}
REOF
verify_r_packages "CRAN helper install" pepr optigrab

# GenomicDistributions packages are needed for PEPATAC QC reports. Try conda
# first; fall back to BiocManager if post-link scripts fail (a known issue on
# some machines for GenomicInfoDbData).
if verify_r_packages "GenomicDistributions precheck" \
    GenomicDistributions GenomicDistributionsData >/dev/null 2>&1; then
    ok "GenomicDistributions packages already present."
else
    info "Installing GenomicDistributions packages..."
    if run_with_spinner "Installing GenomicDistributions via conda" -- \
        conda_install_env bioconductor-genomicdistributions bioconductor-genomicdistributionsdata; then
        ok "GenomicDistributions packages installed with conda."
    else
        warn "Conda install of GenomicDistributions packages failed; trying BiocManager fallback."
        in_env_clean Rscript --vanilla -e '
options(repos=BiocManager::repositories())
for (pkg in c("GenomicDistributions", "GenomicDistributionsData")) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        BiocManager::install(pkg, ask=FALSE, update=FALSE)
    }
}
'
    fi
fi
verify_r_packages "GenomicDistributions final check" \
    GenomicDistributions GenomicDistributionsData

ok "Core PEPATAC R packages installed and verified."

# STEP 6b — DiffBind + differential analysis R packages
# These packages are required by PEPATAC_diff_analysis.sh.
# They are kept separate from the PEPATAC core packages above
# because they are only needed for differential accessibility
# analysis, not for the alignment/peak-calling pipeline.
# edgeR and limma are explicitly required by DiffBind's RLE normalization
# (DBA_NORM_RLE). csaw is required when background=TRUE generates genome-wide
# background bins. None of these dependencies is optional for the default path.

section "STEP 6b — DiffBind + DESeq2 differential analysis packages"

# Try conda first for the heaviest Bioconductor packages.
# This avoids source compilation and is much faster.
DIFFBIND_CONDA_PKGS=(
    bioconductor-diffbind
    bioconductor-deseq2
    bioconductor-edger
    bioconductor-limma
    "bioconductor-csaw=${CSAW_VERSION}"
    bioconductor-biocparallel
    bioconductor-genomicranges
    bioconductor-rsamtools
    bioconductor-summarizedexperiment
    bioconductor-annotationdbi
    bioconductor-genomicfeatures
    bioconductor-genomeinfodb
    bioconductor-iranges
    r-rsqlite
    bioconductor-chipseeker
    bioconductor-clusterprofiler
    bioconductor-enrichplot
    bioconductor-go.db
    bioconductor-org.hs.eg.db
    bioconductor-org.mm.eg.db
    bioconductor-org.rn.eg.db
    bioconductor-org.dm.eg.db
    bioconductor-org.dr.eg.db
    r-ggrepel
)

info "Installing DiffBind and dependencies via conda (this may take several minutes)..."
if run_with_spinner "Installing DiffBind packages via conda" -- \
    conda_install_env "${DIFFBIND_CONDA_PKGS[@]}"; then
    ok "DiffBind packages installed via conda."
else
    warn "Conda install of DiffBind packages failed; falling back to BiocManager."
    in_env_clean Rscript --vanilla -e '
options(repos=c(CRAN="https://cloud.r-project.org"))
bioc_pkgs <- c(
    "DiffBind", "DESeq2", "edgeR", "limma", "csaw",
    "BiocParallel", "GenomicRanges", "Rsamtools",
    "SummarizedExperiment", "AnnotationDbi", "GenomicFeatures", "GenomeInfoDb", "IRanges",
    "ChIPseeker", "clusterProfiler", "enrichplot", "GO.db",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "TxDb.Rnorvegicus.UCSC.rn7.refGene",
    "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
    "TxDb.Drerio.UCSC.danRer11.refGene"
)
cran_pkgs <- c("ggrepel", "dplyr", "tidyr", "RSQLite")
for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        cat("  Installing", pkg, "from CRAN\n")
        install.packages(pkg, quiet=TRUE)
    }
}
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        cat("  Installing", pkg, "from Bioconductor\n")
        BiocManager::install(pkg, ask=FALSE, update=FALSE)
    }
}
'
fi

# Enforce the exact csaw release even when the broad conda transaction fell
# back to BiocManager. Presence alone is insufficient: a csaw build from a
# different Bioconductor generation can be ABI/API-incompatible with R 4.3 and
# the rest of this environment.
info "Ensuring pinned csaw ${CSAW_VERSION} (Bioconductor ${CSAW_BIOC_VERSION})..."
in_env_clean Rscript --vanilla - "$CSAW_VERSION" "$CSAW_BIOC_VERSION" <<'RCSAW'
args <- commandArgs(trailingOnly = TRUE)
expected_version <- args[1]
bioc_version <- args[2]

installed_version <- if (requireNamespace("csaw", quietly = TRUE)) {
    as.character(packageVersion("csaw"))
} else {
    NA_character_
}

if (is.na(installed_version) || !identical(installed_version, expected_version)) {
    cat(sprintf("  Installing pinned csaw %s from Bioconductor %s (currently: %s)\n",
                expected_version, bioc_version,
                ifelse(is.na(installed_version), "not installed", installed_version)))
    options(repos = BiocManager::repositories(version = bioc_version))
    BiocManager::install(
        "csaw",
        version = bioc_version,
        ask = FALSE,
        update = FALSE,
        force = TRUE
    )
}

if (!requireNamespace("csaw", quietly = TRUE)) {
    stop("csaw did not install successfully")
}
actual_version <- as.character(packageVersion("csaw"))
if (!identical(actual_version, expected_version)) {
    stop(sprintf("csaw version mismatch after install: expected %s, found %s",
                 expected_version, actual_version))
}
cat(sprintf("  csaw version verified: %s\n", actual_version))
RCSAW
ok "Pinned csaw ${CSAW_VERSION} installed and verified."

# TxDb packages are Bioconductor annotation-data packages. Install them with
# BiocManager even when the conda core-package transaction succeeds, because
# not every Bioconductor release exposes every TxDb through the same conda label.
in_env_clean Rscript --vanilla - <<'RTXDB'
annotation_pkgs <- c(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "TxDb.Rnorvegicus.UCSC.rn7.refGene",
    "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
    "TxDb.Drerio.UCSC.danRer11.refGene"
)
missing <- annotation_pkgs[!vapply(annotation_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) {
    cat("  Installing assembly-matched TxDb packages:", paste(missing, collapse=", "), "\n")
    BiocManager::install(missing, ask=FALSE, update=FALSE)
}
RTXDB

verify_r_packages "DiffBind dependency install" \
    DiffBind DESeq2 edgeR limma csaw BiocParallel GenomicRanges \
    Rsamtools SummarizedExperiment AnnotationDbi GenomicFeatures GenomeInfoDb IRanges RSQLite ChIPseeker \
    ggrepel dplyr tidyr clusterProfiler enrichplot GO.db \
    org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db org.Dm.eg.db org.Dr.eg.db \
    TxDb.Hsapiens.UCSC.hg38.knownGene \
    TxDb.Mmusculus.UCSC.mm10.knownGene \
    TxDb.Rnorvegicus.UCSC.rn7.refGene \
    TxDb.Dmelanogaster.UCSC.dm6.ensGene \
    TxDb.Drerio.UCSC.danRer11.refGene

in_env_clean Rscript --vanilla - "$CSAW_VERSION" <<'RCSAWVERIFY'
expected <- commandArgs(trailingOnly = TRUE)[1]
actual <- as.character(packageVersion("csaw"))
if (!identical(actual, expected)) {
    stop(sprintf("csaw version mismatch: expected %s, found %s", expected, actual))
}
cat(sprintf("    csaw exact version OK: %s\n", actual))
RCSAWVERIFY

ok "DiffBind and differential analysis packages installed and verified."

# STEP 7 — Clone / update PEPATAC repo

section "STEP 7 — PEPATAC repo"

# Pinned to a specific commit (the repo's HEAD as of 2026-07-30) rather
# than always tracking the default branch. This previously pulled whatever
# was newest on every run, which means the exact pipeline code in use
# could silently change between installs with no record of what changed.
# Update PEPATAC_COMMIT deliberately -- after validating the new commit
# against this installer and PEPATACr -- rather than removing the pin.
PEPATAC_COMMIT="69c3fad7226152a5f4e040657cd1cb22a4aef6c5"

if [[ -d "$PEPATAC_DIR/.git" ]]; then
    ok "Repo exists at $PEPATAC_DIR — fetching and checking out the pinned commit."
    run_with_spinner "Fetching PEPATAC updates" -- git -C "$PEPATAC_DIR" fetch origin
    run_with_spinner "Checking out pinned PEPATAC commit ($PEPATAC_COMMIT)" -- \
        git -C "$PEPATAC_DIR" checkout "$PEPATAC_COMMIT"
else
    run_with_spinner "Cloning PEPATAC repo" -- \
        git clone https://github.com/databio/pepatac.git "$PEPATAC_DIR"
    run_with_spinner "Checking out pinned PEPATAC commit ($PEPATAC_COMMIT)" -- \
        git -C "$PEPATAC_DIR" checkout "$PEPATAC_COMMIT"
    ok "Cloned PEPATAC to $PEPATAC_DIR"
fi

# STEP 8 — Python requirements for PEPATAC

section "STEP 8 — PEPATAC Python requirements"

# Some conda-created environments do not expose a bare `pip` executable.
# Use `python -m pip`, and install pip first if it is missing.
if ! in_env python -m pip --version >/dev/null 2>&1; then
    warn "pip is missing inside the '$ENV_NAME' environment. Installing pip..."
    conda_install_env pip
fi

# pip and looper pinned to their current PyPI releases (checked 2026-07-30)
# rather than always installing whatever is newest. looper is the actual
# pipeline-execution tool PEPATAC runs through, so an unpinned version can
# change run behavior between installs with no record of what changed.
# pip's own version doesn't affect any scientific output -- it's pinned
# here mainly for install-to-install consistency; if a future package
# needs a newer pip resolver feature, bump PIP_VERSION deliberately rather
# than dropping the pin back to "--upgrade" with no version.
PIP_VERSION="26.2"
LOOPER_VERSION="2.1.1"

pip_install_self_healing "Upgrading pip" -- in_env python -m pip install "pip==$PIP_VERSION"
pip_install_self_healing "Installing PEPATAC Python requirements" -- \
    in_env python -m pip install -r "$PEPATAC_DIR/requirements.txt"
pip_install_self_healing "Installing looper" -- in_env python -m pip install "looper==$LOOPER_VERSION"

ok "Python requirements installed."

# STEP 9 — PEPATACr

section "STEP 9 — PEPATACr"

run_with_spinner "Installing PEPATACr from local repo" -- \
    in_env_clean Rscript --vanilla -e \
    "remotes::install_local('$PEPATAC_DIR/PEPATACr', dependencies=FALSE)"

in_env_clean Rscript --vanilla -e \
    'if (!requireNamespace("PEPATACr", quietly=TRUE)) stop("PEPATACr did not install/load"); cat("    PEPATACr loaded OK\n")'

ok "PEPATACr installed and verified."

# STEP 10 — pepatac.py shortcut

section "STEP 10 — pepatac.py shortcut"

mkdir -p "$HOME/bin"
ln -sf "$PEPATAC_DIR/pipelines/pepatac.py" "$HOME/bin/pepatac.py"
chmod +x "$PEPATAC_DIR/pipelines/pepatac.py"
safe_append_bashrc 'export PATH="$HOME/bin:$PATH"'
export PATH="$HOME/bin:$PATH"

ok "pepatac.py linked to ~/bin/pepatac.py"

# ============================================================
# ============================================================
# STEP 11 — Rust + gtars

section "STEP 11 — Rust + gtars"

# Pinned gtars git tag — 0.8.0 has a known bug where gtars uniwig silently
# fails to process BAM chromosomes, producing no per-chromosome bigwig files,
# then floods the log with hundreds of "Error opening file: No such file or
# directory" errors when it tries to merge bigwigs that were never created.
#
# v0.9.0 is the version this installer has actually been building and
# running against. A previous version of this comment claimed a bare
# "v0.9.2" release was needed instead -- that tag does not exist for the
# gtars CLI. Checked directly against the upstream repo's tags: the only
# v0.9.2-labeled tags are gtars-python-v0.9.2 and gtars-r-v0.9.2, which
# version the separate Python/R language bindings, not the CLI this
# installer builds from source. `git checkout v0.9.2` here would simply
# fail with "did not match any file(s) known to git". v0.9.0 is a real,
# valid tag and is what gets built below -- do not change GTARS_VERSION
# to "0.9.2" based on that stale claim.
#
# gtars-cli is not published to crates.io; it must be built from the repo
# at the matching git tag. If a newer core CLI tag (a bare vX.Y.Z, not a
# language-binding-specific one) is validated against PEPATAC in the
# future, update GTARS_VERSION then.
GTARS_VERSION="0.9.0"
GTARS_TAG="v${GTARS_VERSION}"
GTARS_DIR="$HOME/gtars"

# Install Rust via rustup (the official method — cargo is not on conda).
safe_append_bashrc 'export PATH="$HOME/.cargo/bin:$PATH"'
export PATH="$HOME/.cargo/bin:$PATH"

if command_exists rustc; then
    ok "Rust already installed: $(rustc --version)"
else
    info "Installing Rust via rustup..."
    run_with_spinner "Installing Rust" -- \
        bash -c 'curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path'
    ok "Rust installed."
fi

# Reload cargo into this session.
if [[ -f "$HOME/.cargo/env" ]]; then
    source "$HOME/.cargo/env"
fi

# Check whether the correct version is already installed.
GTARS_INSTALLED_VERSION="$(gtars --version 2>/dev/null | awk '{print $2}' || echo '')"
if [[ "$GTARS_INSTALLED_VERSION" == "$GTARS_VERSION" ]]; then
    ok "gtars $GTARS_VERSION already installed."
else
    if [[ -n "$GTARS_INSTALLED_VERSION" ]]; then
        warn "gtars $GTARS_INSTALLED_VERSION is installed but the required version is $GTARS_VERSION — reinstalling."
    else
        info "Installing gtars $GTARS_VERSION from git tag $GTARS_TAG..."
    fi

    # Clone or update the repo, then checkout the pinned tag before building.
    # gtars-cli is not published to crates.io — it must be built from source.
    if [[ -d "$GTARS_DIR/.git" ]]; then
        info "gtars repo exists — fetching tags."
        run_with_spinner "Fetching gtars tags" -- git -C "$GTARS_DIR" fetch --tags
    else
        run_with_spinner "Cloning gtars repo" -- \
            git clone https://github.com/databio/gtars.git "$GTARS_DIR"
    fi

    run_with_spinner "Checking out gtars $GTARS_TAG" -- \
        git -C "$GTARS_DIR" checkout "$GTARS_TAG"

    run_with_spinner "Building gtars $GTARS_VERSION (cargo install — this may take a few minutes)" -- \
        cargo install \
        --path "$GTARS_DIR/gtars-cli" \
        --features "uniwig" \
        --force

    GTARS_INSTALLED_VERSION="$(gtars --version 2>/dev/null | awk '{print $2}' || echo '')"
    if [[ "$GTARS_INSTALLED_VERSION" == "$GTARS_VERSION" ]]; then
        ok "gtars $GTARS_VERSION installed successfully."
    else
        warn "gtars version after install: '$GTARS_INSTALLED_VERSION' (expected $GTARS_VERSION)."
        warn "This may cause bigwig generation failures in PEPATAC."
    fi
fi

# STEP 12 — HOMER (genome downloads deferred to runner)

section "STEP 12 — HOMER genome downloads deferred to runner"

# HOMER software is installed via conda in STEP 5. HOMER genome packages are
# intentionally NOT installed here because the runner checks and downloads the
# matching HOMER genome only after the user chooses a genome and enables motif
# analysis. This prevents installer-time downloads for genomes that may never
# be used.

HOMER_CONFIGURE="$(in_env find "$HOME/miniconda3/envs/$ENV_NAME" \
    -name configureHomer.pl 2>/dev/null | head -1)"

if [[ -z "$HOMER_CONFIGURE" ]]; then
    warn "configureHomer.pl not found. HOMER software may not be installed correctly."
else
    ok "HOMER software found: $HOMER_CONFIGURE"
    ok "HOMER genome downloads will be handled by the runner when needed."
fi

# STEP 13 — Refgenie init (genome assets deferred to runner)

section "STEP 13 — Refgenie init"

mkdir -p "$REFGENIE_DIR"
export REFGENIE="$REFGENIE_CONFIG"
safe_append_bashrc 'export REFGENIE="$HOME/refgenie/refgenie.yaml"'

if [[ -f "$REFGENIE_CONFIG" ]]; then
    ok "Refgenie config already exists: $REFGENIE_CONFIG"
else
    run_with_spinner "Initializing Refgenie config" -- \
        in_env refgenie init \
        -c "$REFGENIE_CONFIG" \
        -s http://refgenomes.databio.org
    ok "Refgenie initialized."
fi

ok "Refgenie genome assets will be downloaded by the runner after genome selection."

# ── Repair broken chrom sizes symlinks ───────────────────────
# If a previous install or run left a dead symlink pointing at a path
# that no longer exists (e.g. a file that was on a Windows drive that
# got unmounted), regenerate the actual file from the FASTA that refgenie
# already has. This is the most common cause of the preflight error:
#   "Chrom sizes symlink is broken or unreadable: ... -> /home/.../genomes/hg38.chrom.sizes"

section "STEP 13b — Repair broken refgenie chrom sizes symlinks"

REFGENIE_ALIAS_DIR="$REFGENIE_DIR/alias"
REPAIRED=0
SKIPPED=0

if [[ -d "$REFGENIE_ALIAS_DIR" ]]; then
    # Walk every genome alias refgenie knows about.
    for genome_dir in "$REFGENIE_ALIAS_DIR"/*/; do
        genome="$(basename "$genome_dir")"
        chrom_sizes_link="$genome_dir/fasta/default/${genome}.chrom.sizes"
        fasta_file="$genome_dir/fasta/default/${genome}.fa"

        # Only act when the symlink exists but is broken (points nowhere).
        if [[ -L "$chrom_sizes_link" ]] && [[ ! -e "$chrom_sizes_link" ]]; then
            warn "Broken chrom sizes symlink detected for $genome:"
            warn "  $chrom_sizes_link -> $(readlink "$chrom_sizes_link")"

            if [[ ! -f "$fasta_file" ]]; then
                warn "  FASTA not found at $fasta_file — cannot repair automatically."
                warn "  Re-run the runner to re-download genome assets for $genome."
                (( SKIPPED++ )) || true
                continue
            fi

            # The target path the symlink was supposed to point at.
            # Recreate it so the symlink resolves again.
            TARGET="$(readlink "$chrom_sizes_link")"
            TARGET_DIR="$(dirname "$TARGET")"
            mkdir -p "$TARGET_DIR"

            info "Regenerating chrom sizes for $genome from FASTA..."
            if in_env samtools faidx "$fasta_file" 2>/dev/null && \
               cut -f1,2 "${fasta_file}.fai" > "$TARGET"; then
                ok "Repaired: $TARGET"
                (( REPAIRED++ )) || true
            else
                warn "samtools faidx failed for $genome — trying fallback via refgenie pull."
                # Ask refgenie to re-pull just the chrom sizes asset.
                in_env refgenie pull -c "$REFGENIE_CONFIG" "${genome}/fasta:chrom_sizes" 2>/dev/null || \
                    warn "refgenie pull also failed — runner will re-download assets at next run."
                (( SKIPPED++ )) || true
            fi
        fi
    done

    if [[ $REPAIRED -gt 0 ]]; then
        ok "Repaired $REPAIRED broken chrom sizes symlink(s)."
    fi
    if [[ $SKIPPED -gt 0 ]]; then
        warn "$SKIPPED genome(s) could not be auto-repaired — runner will handle asset download."
    fi
    if [[ $REPAIRED -eq 0 && $SKIPPED -eq 0 ]]; then
        ok "No broken chrom sizes symlinks found."
    fi
else
    ok "No refgenie alias directory yet — nothing to repair."
fi

# STEP 14 — Convenience environment checker

section "STEP 14 — Convenience checker"

mkdir -p "$HOME/bin"
safe_append_bashrc 'export PATH="$HOME/bin:$PATH"'

cat > "$HOME/bin/pepatac_check.sh" <<'CHECKEOF'
#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="FetchPA"
CSAW_VERSION="1.36.0"
GTARS_VERSION="0.9.0"
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"

if [[ -f "$CONDA_SH" ]]; then
    source "$CONDA_SH"
else
    echo "[ERROR] Conda startup script not found: $CONDA_SH"
    exit 1
fi

echo "PEPATAC environment check"
echo "=========================="
echo "Environment: $ENV_NAME"
echo ""

missing=0
for tool in pepatac.py looper bowtie2 bowtie2-inspect samtools macs3 bedtools fastqc bamCoverage \
            computeMatrix plotHeatmap skewer cutadapt trim_galore samblaster preseq \
            refgenie findMotifsGenome.pl Rscript; do
    if conda run --no-capture-output -n "$ENV_NAME" bash -c "command -v $tool" >/dev/null 2>&1; then
        printf "  OK      %s\n" "$tool"
    else
        printf "  MISSING %s\n" "$tool"
        missing=1
    fi
done

# gtars gets its own check rather than folding into the loop above:
# existence on PATH isn't enough -- v0.8.0 has a known bug (see the
# GTARS_VERSION comment in PEPATAC_install.sh STEP 11), so a stale or
# wrong-version binary still on PATH needs to fail this check, not pass it.
# Checked bare (not wrapped in conda run), matching how it's installed --
# gtars is a cargo binary on the login PATH via ~/.cargo/env, not part of
# the FetchPA conda environment; wrapping it in conda run risks a false
# MISSING if ~/.cargo/bin isn't inherited into that subprocess's PATH.
#
# The `|| echo ''` matters under set -e + pipefail (both active in this
# script): a genuinely missing gtars makes `gtars --version` exit non-zero
# before awk ever runs, and pipefail propagates that failure out of the
# whole pipeline -- without the fallback, this bare assignment would
# terminate the script right here instead of falling through to the
# MISSING message below, which is the one case this check most needs to
# handle gracefully. Matches the same protective pattern STEP 11's own
# gtars checks already use.
gtars_ver_found="$(gtars --version 2>/dev/null | awk '{print $2}' || echo '')"
if [[ "$gtars_ver_found" == "$GTARS_VERSION" ]]; then
    printf "  OK      gtars (%s)\n" "$gtars_ver_found"
else
    printf "  MISSING gtars (found '%s', need %s)\n" "${gtars_ver_found:-none}" "$GTARS_VERSION"
    missing=1
fi

echo ""
echo "R package check"
echo "==============="

set +e
conda run --no-capture-output -n "$ENV_NAME" \
    env -u R_ARCH -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
    Rscript --vanilla - "$CSAW_VERSION" <<'RCHECK'
expected_csaw <- commandArgs(trailingOnly = TRUE)[1]
required_pkgs <- c(
  "PEPATACr", "pepr", "optigrab",
  "pheatmap",
  "GenomicDistributions", "GenomicDistributionsData",
  "DiffBind", "DESeq2", "edgeR", "limma", "csaw",
  "BiocParallel", "GenomicRanges", "Rsamtools",
  "SummarizedExperiment", "AnnotationDbi", "GenomicFeatures", "GenomeInfoDb", "IRanges", "RSQLite",
  "ChIPseeker", "clusterProfiler", "enrichplot", "GO.db",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "TxDb.Rnorvegicus.UCSC.rn7.refGene",
  "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
  "TxDb.Drerio.UCSC.danRer11.refGene",
  "ggrepel", "dplyr", "tidyr",
  "Matrix", "fs", "XML", "xml2", "BiocManager", "remotes"
)
status <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
csaw_version_ok <- status[["csaw"]] &&
  identical(as.character(packageVersion("csaw")), expected_csaw)
for (pkg in required_pkgs) {
  if (!status[[pkg]]) {
    cat(sprintf("  %-16s %s\n", "MISSING", pkg))
  } else if (pkg == "csaw" && !csaw_version_ok) {
    cat(sprintf("  %-16s %s (installed %s; expected %s)\n",
                "VERSION_MISMATCH", pkg,
                as.character(packageVersion("csaw")), expected_csaw))
  } else {
    cat(sprintf("  %-16s %s\n", "OK", pkg))
  }
}
if (!all(status) || !csaw_version_ok) quit(status = 1)
RCHECK
r_status=$?
set -e

echo ""
echo "Refgenie config check"
echo "====================="
REFGENIE_CONFIG="$HOME/refgenie/refgenie.yaml"
if [[ -f "$REFGENIE_CONFIG" ]]; then
    echo "  OK      $REFGENIE_CONFIG"
else
    echo "  MISSING $REFGENIE_CONFIG"
    missing=1
fi

if [[ "$r_status" -ne 0 ]]; then
    missing=1
fi

echo ""
if [[ "$missing" -eq 0 ]]; then
    echo "All PEPATAC tools, R packages, and config are present."
else
    echo "Some tools, R packages, or config are missing. Re-run PEPATAC_install.sh and check the install log."
    exit 1
fi
CHECKEOF

chmod +x "$HOME/bin/pepatac_check.sh"
ok "Wrote checker: $HOME/bin/pepatac_check.sh"

# STEP 15 — Verification

section "STEP 15 — Verification"

MISSING_TOOLS=()

echo ""
echo "  CLI tools:"
for tool in pepatac.py looper bowtie2 bowtie2-inspect samtools macs3 bedtools fastqc bamCoverage \
            computeMatrix plotHeatmap skewer cutadapt trim_galore samblaster preseq \
            refgenie findMotifsGenome.pl Rscript; do
    if in_env bash -c "command -v $tool" >/dev/null 2>&1; then
        echo "    OK      $tool"
    else
        echo "    MISSING $tool"
        MISSING_TOOLS+=("$tool")
    fi
done

for tool in sha256sum flock; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "    OK      $tool (system)"
    else
        echo "    MISSING $tool (system)"
        MISSING_TOOLS+=("$tool")
    fi
done

# gtars gets its own check rather than folding into the loop above:
# existence on PATH isn't enough -- v0.8.0 has a known bug (see the
# GTARS_VERSION comment in STEP 11 above), so a stale or wrong-version
# binary left on PATH (e.g. from a system-wide install predating this
# script, or a reinstall that didn't actually replace what's resolved
# first on PATH) needs to fail this check, not silently pass it. Checked
# bare (not wrapped in in_env/conda run), matching STEP 11's own gtars
# checks -- gtars is a cargo binary on the login PATH via ~/.cargo/env,
# not part of the FetchPA conda environment.
#
# The `|| echo ''` matters under set -e + pipefail (both active in this
# script): a genuinely missing gtars makes this bare assignment terminate
# the script right here instead of falling through to the MISSING message
# below -- the one case this check most needs to handle gracefully.
GTARS_VERSION_FOUND="$(gtars --version 2>/dev/null | awk '{print $2}' || echo '')"
if [[ "$GTARS_VERSION_FOUND" == "$GTARS_VERSION" ]]; then
    echo "    OK      gtars ($GTARS_VERSION_FOUND)"
else
    echo "    MISSING gtars (found '${GTARS_VERSION_FOUND:-none}', need $GTARS_VERSION)"
    MISSING_TOOLS+=("gtars_version_mismatch")
fi

echo ""
echo "  R packages:"

set +e
R_PACKAGE_CHECK=$(in_env_clean Rscript --vanilla - "$CSAW_VERSION" <<'RCHECK'
expected_csaw <- commandArgs(trailingOnly = TRUE)[1]
required_pkgs <- c(
  "PEPATACr", "pepr", "optigrab",
  "pheatmap",
  "GenomicDistributions", "GenomicDistributionsData",
  "DiffBind", "DESeq2", "edgeR", "limma", "csaw",
  "BiocParallel", "GenomicRanges", "Rsamtools",
  "SummarizedExperiment", "AnnotationDbi", "GenomicFeatures", "GenomeInfoDb", "IRanges", "RSQLite",
  "ChIPseeker", "clusterProfiler", "enrichplot", "GO.db",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "TxDb.Rnorvegicus.UCSC.rn7.refGene",
  "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
  "TxDb.Drerio.UCSC.danRer11.refGene",
  "ggrepel", "dplyr", "tidyr",
  "Matrix", "fs", "XML", "xml2", "BiocManager", "remotes"
)
status <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
csaw_version_ok <- status[["csaw"]] &&
  identical(as.character(packageVersion("csaw")), expected_csaw)
for (pkg in required_pkgs) {
  if (!status[[pkg]]) {
    cat(sprintf("MISSING\t%s\n", pkg))
  } else if (pkg == "csaw" && !csaw_version_ok) {
    cat(sprintf("VERSION_MISMATCH\t%s (installed %s; expected %s)\n",
                pkg, as.character(packageVersion("csaw")), expected_csaw))
  } else {
    cat(sprintf("OK\t%s\n", pkg))
  }
}
if (!all(status) || !csaw_version_ok) quit(status = 1)
RCHECK
)
R_PACKAGE_CHECK_STATUS=$?
set -e

if [[ -n "$R_PACKAGE_CHECK" ]]; then
    echo "$R_PACKAGE_CHECK" | awk -F '\t' '{printf "    %-8s %s\n", $1, $2}'
else
    echo "    [ERROR] R package check produced no output (Rscript exit $R_PACKAGE_CHECK_STATUS)."
    echo "    The R environment itself may be broken -- check the FetchPA conda env directly."
fi

# Gated on Rscript's own exit status, not on text-matching the captured
# output. The R heredoc above calls quit(status=1) whenever any package is
# missing/mismatched, so R_PACKAGE_CHECK_STATUS covers that case -- but it
# also covers the case a pure text-match against MISSING/VERSION_MISMATCH
# cannot: R crashing before any cat() line runs at all (e.g. requireNamespace
# itself erroring unexpectedly, or Rscript failing to start), which would
# leave R_PACKAGE_CHECK empty and let a broken R environment through as if
# every package were fine. The previous `) || true` here discarded the exit
# status entirely, so only the pattern-match (which that empty-output case
# defeats) stood between a crash and a false "COMPLETE."
if [[ "$R_PACKAGE_CHECK_STATUS" -ne 0 ]]; then
    MISSING_TOOLS+=("one_or_more_required_R_packages")
fi

echo ""
echo "  Refgenie config:"
if [[ -f "$REFGENIE_CONFIG" ]]; then
    echo "    OK      $REFGENIE_CONFIG"
else
    echo "    MISSING $REFGENIE_CONFIG"
    MISSING_TOOLS+=("refgenie_config")
fi

echo ""
echo "  Versions:"
in_env bowtie2 --version 2>/dev/null | head -n 1 | sed 's/^/    /' || true
in_env samtools --version 2>/dev/null | head -n 1 | sed 's/^/    /' || true
in_env fastqc --version 2>/dev/null | sed 's/^/    /' || true
in_env macs3 --version 2>/dev/null | sed 's/^/    /' || true
in_env bamCoverage --version 2>/dev/null | sed 's/^/    /' || true
in_env Rscript --version 2>&1 | sed 's/^/    /' || true
in_env_clean Rscript --vanilla -e 'cat("    PEPATACr: ", as.character(packageVersion("PEPATACr")), "\n", sep="")' 2>/dev/null || true
in_env_clean Rscript --vanilla -e 'cat("    DiffBind: ", as.character(packageVersion("DiffBind")), "\n", sep="")' 2>/dev/null || true
in_env_clean Rscript --vanilla -e 'cat("    DESeq2: ", as.character(packageVersion("DESeq2")), "\n", sep="")' 2>/dev/null || true
in_env_clean Rscript --vanilla -e 'cat("    csaw: ", as.character(packageVersion("csaw")), " (expected 1.36.0)\n", sep="")' 2>/dev/null || true
gtars --version 2>/dev/null | sed 's/^/    gtars: /' || true

if [[ ${#MISSING_TOOLS[@]} -eq 0 ]]; then
    ok "All PEPATAC tool, R package, and config checks passed."
else
    warn "Some tools, R packages, or config are missing — review the output above."
fi

# STEP 16 — Export environment

section "STEP 16 — Export conda environment"

conda env export -n "$ENV_NAME" > "$HOME/pepatac_environment_full.yml"
conda env export -n "$ENV_NAME" --from-history > "$HOME/pepatac_environment_minimal.yml"

ok "Saved full environment: $HOME/pepatac_environment_full.yml"
ok "Saved minimal environment: $HOME/pepatac_environment_minimal.yml"

# DONE

# A nonempty MISSING_TOOLS must stop this from claiming success. Without
# this, a partially-broken install still prints "INSTALL COMPLETE" and
# exits 0, so a wrapper script or a person skimming the output has no
# signal that anything needs attention.
if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    section "INSTALL INCOMPLETE"
    err "The following did not pass verification:"
    for item in "${MISSING_TOOLS[@]}"; do
        err "  - $item"
    done
    err ""
    err "Environment files were still exported above (STEP 16) for debugging,"
    err "but this install is not ready to run PEPATAC. This script is safe to"
    err "re-run -- it skips work that's already done -- or install the"
    err "missing piece(s) manually. Then re-run $HOME/bin/pepatac_check.sh to"
    err "confirm before using the pipeline."
    exit 1
fi

section "INSTALL COMPLETE"

cat <<DONE

  PEPATAC tools, R packages, and pipeline code are installed and ready.

  Installer log:
    $INSTALL_LOG

  Main paths:
    PEPATAC pipeline:    $PEPATAC_DIR
    Refgenie config:     $REFGENIE_CONFIG
    Refgenie data:       $REFGENIE_DIR

  To check the install later:
    pepatac_check.sh

  Genome assets:
    Refgenie manages only FASTA/chrom-size/Bowtie2 assets after genome selection.
    Bioconductor TxDb/OrgDb packages provide peak/gene annotation; HOMER genomes are motif-only.

  To activate manually:
    conda activate $ENV_NAME

  Next steps:
    Run PEPATAC_run_flexible_paths.sh to process samples through the
    pipeline. Once you have differential-accessibility results,
    PEPATAC_diff_analysis.sh and PEPATAC_explore_organized_outputs.sh
    handle the DiffBind/DESeq2 analysis and interactive results
    exploration. All three are interactive and prompt for what they need.

DONE
