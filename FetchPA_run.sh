#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# PEPATAC RUNNER
# Interactive sample processor for local machines.
#
# Flow:
#   1. Collect all run configuration first (or restore from a
#      previous run's resume state).
#   2. Show a final run summary.
#   3. Only after confirmation: download/install missing assets.
#   4. Run preflight checks.
#   5. Take off and process samples (skipping any that already
#      PASSED in a previous run).
# ============================================================

ENV_NAME="FetchPA"
PEPATAC_DIR="$HOME/pepatac"
REFGENIE_CONFIG="$HOME/refgenie/refgenie.yaml"
PIPELINE="$PEPATAC_DIR/pipelines/pepatac.py"
BLACKLIST_DIR="$HOME/pepatac_blacklists"

RUNNER_VERSION="1.34-all-prompts-explicit"
SCRIPT_VERSION="$RUNNER_VERSION"
# Schema 3: custom-profile-only snapshot (CUSTOM_* fields).
# Schema 4: adds genome-agnostic ANNOTATION_* fields, written for every
# genome (built-in or custom) -- see write_reference_snapshot(). A schema-3
# file on disk is still fully readable by this version (all its CUSTOM_*
# fields are unchanged); this version just also writes the newer fields
# going forward, and diff_analysis.sh reads them when present.
REFERENCE_SNAPSHOT_SCHEMA=4
RUN_STARTED="$(date '+%Y-%m-%d %H:%M:%S')"
RUN_ID="$(date '+%Y%m%d_%H%M%S')"



SUPPORTED_GENOMES=("mm10" "hg38" "rn7" "dm6" "danRer11")

# Effective (mappable) genome size for MACS peak calling, passed to pepatac.py
# via -gs/--genome-size. Without this flag pepatac.py defaults to 2.7e9 (the
# human value) for EVERY genome -- silently wrong for anything else, since
# MACS's significance calls depend on this number.
#
# hg38, mm10, dm6 use the effective genome sizes documented in current
# MACS3 (3.0.4) callpeak docs -- https://macs3-project.github.io/MACS/docs/callpeak.html
# (2,913,022,398 / 2,652,783,500 / 142,573,017 respectively, verified
# 2026-07-28). Note this differs from the older, still widely-cited MACS2-era
# shortcuts (hs=2.7e9, mm=1.87e9, dm=1.2e8) -- both have been "official" at
# different points; these are the current documented values. rn7 and
# danRer11 have no official MACS shortcut and no consensus published value
# -- rather than hardcode a guess, their size is
# computed at runtime from the actual chrom.sizes file for the assembly in
# use (sum of primary contigs x ~0.90, the "mappable fraction" approximation
# MACS's own docs cite for genomes without a precompiled value). Override
# either genome's value by exporting PEPATAC_GENOME_SIZE before running this
# script if you have a more precise mappability-based estimate.
declare -A GENOME_SIZE_MAP=(
    [hg38]="2913022398"
    [mm10]="2652783500"
    [dm6]="142573017"
)

# resolve_effective_genome_size GENOME CHROM_SIZES_FILE
# Sets RUN_GENOME_SIZE and GENOME_SIZE_METHOD. The method is persisted in
# the immutable reference snapshot so downstream analyses can distinguish a
# built-in MACS3 value from a user override or chrom-sizes estimate.
resolve_effective_genome_size() {
    local genome="$1"
    local chrom_sizes_file="$2"
    local computed=""

    if [[ -n "${PEPATAC_GENOME_SIZE:-}" ]]; then
        if [[ "$PEPATAC_GENOME_SIZE" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] && \
           awk "BEGIN{exit !($PEPATAC_GENOME_SIZE > 0)}" 2>/dev/null; then
            RUN_GENOME_SIZE="$PEPATAC_GENOME_SIZE"
            GENOME_SIZE_METHOD="environment_override"
            return 0
        fi
        die "PEPATAC_GENOME_SIZE is set but is not a positive number: '$PEPATAC_GENOME_SIZE'"
    fi

    if [[ -n "${GENOME_SIZE_MAP[$genome]+x}" ]]; then
        RUN_GENOME_SIZE="${GENOME_SIZE_MAP[$genome]}"
        GENOME_SIZE_METHOD="built_in_macs3"
        return 0
    fi

    if [[ -n "${CUSTOM_GENOME_SIZE[$genome]:-}" ]]; then
        RUN_GENOME_SIZE="${CUSTOM_GENOME_SIZE[$genome]}"
        GENOME_SIZE_METHOD="user_supplied"
        return 0
    fi

    if [[ -f "$chrom_sizes_file" ]]; then
        if is_supported_genome "$genome"; then
            computed=$(awk '
                $1 !~ /_/ && tolower($1) !~ /(random|unplaced|unlocalized|alt|hap|fix|decoy)/ {
                    total += $2
                }
                END { printf "%.0f", total * 0.90 }
            ' "$chrom_sizes_file")
        else
            # Do not exclude underscores for arbitrary assemblies: scaffold_1,
            # contig_00001, etc. are often the primary sequence names.
            computed=$(awk '
                tolower($1) !~ /(random|unplaced|unlocalized|alt|hap|fix|decoy)/ {
                    total += $2
                }
                END { printf "%.0f", total * 0.90 }
            ' "$chrom_sizes_file")
        fi

        if [[ "$computed" =~ ^[0-9]+$ ]] && [[ "$computed" -gt 0 ]]; then
            RUN_GENOME_SIZE="$computed"
            GENOME_SIZE_METHOD="chrom_sizes_90_percent_estimate"
            return 0
        fi

        die "Computed effective genome size for '$genome' from $chrom_sizes_file was '${computed:-<empty>}' (not a positive number)."$'\n'"       Set PEPATAC_GENOME_SIZE to a known value and re-run, or check chrom.sizes for unusual naming."
    fi

    die "Could not determine an effective genome size for '$genome' (chrom.sizes not found: $chrom_sizes_file)."$'\n'"       Set PEPATAC_GENOME_SIZE to a known value and re-run."
}

# Backward-compatible value-only wrapper for callers/tests that only need the
# number. The main run path calls resolve_effective_genome_size directly so
# GENOME_SIZE_METHOD is retained in the current shell.
genome_size_for() {
    resolve_effective_genome_size "$1" "$2"
    printf '%s\n' "$RUN_GENOME_SIZE"
}

# estimate_memory_per_sample_mb GENOME_INDEX_PREFIX
# Rough per-sample peak memory estimate in MB, used only for the
# oversubscription warning below -- not exact. Bowtie2's actual resident
# memory usage roughly tracks the on-disk size of its index files, so this
# measures the real index for the genome actually in use rather than
# guessing a fixed number. Adds ~20% + 1GB headroom for auxiliary processes
# that run alongside it per sample (samtools, pigz, FastQC). Falls back to
# 4096 MB (close to what's actually been observed for hg38 on this pipeline)
# if the index files can't be found/measured -- e.g. before genome assets
# have been downloaded yet.
estimate_memory_per_sample_mb() {
    local index_prefix="$1"
    local index_dir total_bytes=0
    index_dir=$(dirname "$index_prefix" 2>/dev/null)

    if [[ -d "$index_dir" ]]; then
        total_bytes=$(find "$index_dir" -maxdepth 1 -name "$(basename "$index_prefix")*.bt2*" \
            -exec stat -c%s {} + 2>/dev/null | awk '{sum+=$1} END{print sum+0}')
    fi

    if [[ "$total_bytes" -gt 0 ]]; then
        echo $(( total_bytes * 12 / 10 / 1024 / 1024 + 1024 ))
    else
        echo 4096
    fi
}

# Known blacklist URLs for common genomes.
# Downloaded only after the final run confirmation.
#
# Pinned to a specific commit rather than the mutable `master` branch, so a
# future push to that repo can't silently change which regions get filtered
# out from under a published analysis. Content verified via BLACKLIST_SHA256
# below the same way MINICONDA_SHA256 is verified in PEPATAC_install.sh.
# Pinned commit is master's HEAD as of 2026-08-16 (Blacklist v2 itself was
# last touched in 2019; this repo has seen no substantive changes since).
# Update both the commit and the matching SHA256 together if this ever
# needs to move to a newer blacklist release.
BLACKLIST_COMMIT="61a04d2c5e49341d76735d485c61f0d1177d08a8"
declare -A BLACKLIST_URLS
BLACKLIST_URLS["mm10"]="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/${BLACKLIST_COMMIT}/lists/mm10-blacklist.v2.bed.gz"
BLACKLIST_URLS["hg38"]="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/${BLACKLIST_COMMIT}/lists/hg38-blacklist.v2.bed.gz"
BLACKLIST_URLS["dm6"]="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/${BLACKLIST_COMMIT}/lists/dm6-blacklist.v2.bed.gz"
declare -A BLACKLIST_SHA256
BLACKLIST_SHA256["mm10"]="febafb843c6df492f9a9fc418f8796762ee899d9864330fb509ae2d38ddc0b46"
BLACKLIST_SHA256["hg38"]="c92e763af17271446194991e71917ac220593a5a3d40a06667be24178ef08cf2"
BLACKLIST_SHA256["dm6"]="6174aa0304859a2ee836de7f1eaabc4e7ab7fc4c75ffd4b7587172b1e234026e"
# rn7 and danRer11 do not have published blacklists.

# ─────────────────────────────────────────────────────────────
# Genome metadata — Refgenie availability and UCSC sequence-asset URLs.
# ─────────────────────────────────────────────────────────────
# GENOME_ON_REFGENIE: whether the genome is on the standard Refgenie server.
# UCSC_FASTA_URL:     UCSC bigZips FASTA (.fa.gz) for local builds.
# UCSC_CHROMSIZES_URL: UCSC chrom sizes file for local builds.
# When adding a new genome, populate Refgenie availability plus FASTA and chrom-size URLs.
declare -A GENOME_ON_REFGENIE
GENOME_ON_REFGENIE["mm10"]="yes"
GENOME_ON_REFGENIE["hg38"]="yes"
GENOME_ON_REFGENIE["rn7"]="no"
GENOME_ON_REFGENIE["dm6"]="yes"
GENOME_ON_REFGENIE["danRer11"]="no"

declare -A UCSC_FASTA_URL
UCSC_FASTA_URL["mm10"]="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
UCSC_FASTA_URL["hg38"]="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
UCSC_FASTA_URL["rn7"]="https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz"
UCSC_FASTA_URL["dm6"]="https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
UCSC_FASTA_URL["danRer11"]="https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz"

declare -A UCSC_CHROMSIZES_URL
UCSC_CHROMSIZES_URL["mm10"]="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
UCSC_CHROMSIZES_URL["hg38"]="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
UCSC_CHROMSIZES_URL["rn7"]="https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.chrom.sizes"
UCSC_CHROMSIZES_URL["dm6"]="https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"
UCSC_CHROMSIZES_URL["danRer11"]="https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.chrom.sizes"

# ─────────────────────────────────────────────────────────────
# User-defined genome profiles — assemblies outside the five built-ins.
#
# A profile is accepted only when the user supplies the same scientific
# components hard-coded for the built-in genomes:
#   • an official UCSC assembly FASTA and chromosome-size file derived from the assembly ID
#   • an assembly-matched TxDb package
#   • a species-matched OrgDb package
#   • an optional KEGG organism code
#
# The TxDb and OrgDb are validated after final confirmation, then frozen as
# immutable SQLite databases. Arbitrary GTF/GFF annotation is intentionally
# not supported by this mode.
# ─────────────────────────────────────────────────────────────
declare -A CUSTOM_FASTA_URL
declare -A CUSTOM_CHROMSIZES_URL
declare -A CUSTOM_UCSC_ASSEMBLY
declare -A CUSTOM_TXDB_PKG
declare -A CUSTOM_TXDB_VERSION
declare -A CUSTOM_ORGDB_PKG
declare -A CUSTOM_ORGDB_VERSION
declare -A CUSTOM_KEGG_ORG
declare -A CUSTOM_TXDB_SQLITE
declare -A CUSTOM_TXDB_SQLITE_SHA256
declare -A CUSTOM_ORGDB_SQLITE
declare -A CUSTOM_ORGDB_SQLITE_SHA256
declare -A CUSTOM_TSS_BED
declare -A CUSTOM_TSS_BED_SHA256
declare -A CUSTOM_FEATURE_BED
declare -A CUSTOM_FEATURE_BED_SHA256
declare -A CUSTOM_GENES_BED
declare -A CUSTOM_GENES_BED_SHA256
declare -A CUSTOM_EXONS_BED
declare -A CUSTOM_EXONS_BED_SHA256
declare -A CUSTOM_INTRONS_BED
declare -A CUSTOM_INTRONS_BED_SHA256
declare -A CUSTOM_PROMOTERS_BED
declare -A CUSTOM_PROMOTERS_BED_SHA256
declare -A CUSTOM_ANNOTATION_FINGERPRINT
declare -A CUSTOM_PROFILE_DIR
declare -A CUSTOM_PROFILE_BUILD_CONF
declare -A CUSTOM_PROFILE_BUILD_CONF_SHA256
declare -A CUSTOM_FASTA_SHA256
declare -A CUSTOM_CHROM_SIZES_SHA256
declare -A CUSTOM_BT2_INDEX_SHA256
# Blank means estimate 90% of the assembly length from chrom.sizes.
declare -A CUSTOM_GENOME_SIZE

CUSTOM_GENOME_REGISTRY_ROOT="$HOME/pepatac_custom_genomes"
CUSTOM_GENOME_STAGING_ROOT="$CUSTOM_GENOME_REGISTRY_ROOT/.staging"
CUSTOM_GENOME_LOCK_ROOT="$CUSTOM_GENOME_REGISTRY_ROOT/.locks"

# ENA Portal API endpoint for SRA/ENA accession resolution.
ENA_FILEREPORT_BASE="https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS="run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,library_strategy,instrument_platform,sample_title,experiment_title,read_count,base_count"

# ─────────────────────────────────────────────────────────────
# Persistent user defaults
# ─────────────────────────────────────────────────────────────
# Appending " -newdefault" after a path at any prompt saves it
# as the bracketed default for every future run, without editing
# this script. Source this before any DEFAULT_* value is used
# in a prompt, so saved overrides take effect immediately.
USER_DEFAULTS_FILE="$HOME/.pepatac_run_defaults.sh"
if [[ -f "$USER_DEFAULTS_FILE" ]]; then
    # shellcheck source=/dev/null
    source "$USER_DEFAULTS_FILE"
fi

DEFAULT_OUTPUT_BASE="${DEFAULT_OUTPUT_BASE:-$HOME/pepatac_output}"
DEFAULT_BLACKLIST_DIR="${DEFAULT_BLACKLIST_DIR:-$HOME/pepatac_blacklists}"
DEFAULT_SRA_DOWNLOAD_DIR="${DEFAULT_SRA_DOWNLOAD_DIR:-$DEFAULT_OUTPUT_BASE/sra_downloads}"

# strip_newdefault_marker VALUE_VAR
# If the user's raw input ends in " -newdefault", strips the marker,
# trims trailing whitespace, and returns 0 (marker present) or 1 (absent).
# Modifies VALUE_VAR in place via nameref.
strip_newdefault_marker() {
    local -n _val_ref="$1"
    if [[ "$_val_ref" =~ ^(.*[^[:space:]])[[:space:]]+-newdefault[[:space:]]*$ ]]; then
        _val_ref="${BASH_REMATCH[1]}"
        return 0
    fi
    return 1
}

# save_user_default KEY_NAME VALUE
# Persists VALUE under KEY_NAME into USER_DEFAULTS_FILE, replacing any
# prior line for that key. Safe to call repeatedly.
save_user_default() {
    local key_name="$1"
    local value="$2"
    local tmp_file
    tmp_file="$(mktemp)"

    {
        echo "# Auto-generated by PEPATAC_run.sh - persistent path defaults."
        echo "# Saved via '-newdefault' at a path prompt. Safe to edit or delete."
        if [[ -f "$USER_DEFAULTS_FILE" ]]; then
            grep -v "^${key_name}=" "$USER_DEFAULTS_FILE" 2>/dev/null | grep -v '^#' || true
        fi
        printf '%s=%q\n' "$key_name" "$value"
    } > "$tmp_file"

    mv "$tmp_file" "$USER_DEFAULTS_FILE"
    ok "Saved as new default for future runs: $key_name=$value"
    echo -e "  ${DIM}(stored in $USER_DEFAULTS_FILE — delete that file anytime to reset)${RESET}"
}

# ─────────────────────────────────────────────────────────────
# Terminal colors / formatting
# ─────────────────────────────────────────────────────────────

BOLD="\e[1m"
DIM="\e[2m"
CYAN="\e[36m"
GREEN="\e[32m"
YELLOW="\e[33m"
RED="\e[31m"
BLUE="\e[34m"
MAGENTA="\e[35m"
RESET="\e[0m"

header() {
    echo ""
    echo -e "${BOLD}${CYAN}══════════════════════════════════════════════════${RESET}"
    echo -e "${BOLD}${CYAN}  $1${RESET}"
    echo -e "${BOLD}${CYAN}══════════════════════════════════════════════════${RESET}"
    echo ""
}

label()   { echo -e "  ${BOLD}${BLUE}▸${RESET} $*"; }
ok()      { echo -e "  ${GREEN}✔${RESET}  $*"; }
warn()    { echo -e "  ${YELLOW}⚠${RESET}  $*"; }
err()     { echo -e "  ${RED}✘${RESET}  $*"; }
die()     { err "$*"; exit 1; }
blank()   { echo ""; }

command_exists() { command -v "$1" >/dev/null 2>&1; }


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
# Also refuses paths containing whitespace: Pypiper/PEPATAC split some
# file paths on spaces internally, independent of how carefully this
# script quotes the value, so a space here causes a confusing failure
# much later instead of a clear one now.
normalize_path_var() {
    local -n _path_ref="$1"
    _path_ref="$(normalize_path "$_path_ref")"
    if [[ "$_path_ref" == *[[:space:]]* ]]; then
        die "Path contains whitespace, which Pypiper/PEPATAC cannot handle reliably: $_path_ref"$'\n'"       Please move/rename to a path with no spaces (e.g. ${_path_ref// /_})."
    fi
}

in_env() {
    conda run --no-capture-output -n "$ENV_NAME" "$@"
}

# in_env_clean: like in_env, but strips R_ARCH/R_LIBS/R_LIBS_USER/R_LIBS_SITE
# before running the command. Without this, a user's own ambient R library
# environment variables (very commonly set in a shell profile, e.g.
# R_LIBS_USER pointing at a personal library) could shadow the pipeline's
# own isolated conda R library -- silently loading a different, possibly
# incompatible or version-mismatched TxDb/OrgDb/Bioconductor package than
# the one this environment actually has installed and pinned. Used for
# every Rscript invocation that touches annotation packages specifically
# (ensure_profile_annotation_packages, build_profile_annotation_assets,
# validate_frozen_annotation_profile, and the built-in package-version
# resolution in ensure_builtin_annotation_assets), matching the identical
# definition already present and used in diff_analysis.sh/explore.sh/
# install.sh -- this was the one place in the four scripts where it had
# been called but never actually defined.
in_env_clean() {
    conda run --no-capture-output -n "$ENV_NAME" \
        env -u R_ARCH -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$@"
}

asset_ok() {
    in_env refgenie seek -c "$REFGENIE_CONFIG" "$1" >/dev/null 2>&1
}

# ─────────────────────────────────────────────────────────────
# Reference file validation
# ─────────────────────────────────────────────────────────────

validate_readable_reference_file() {
    local path="$1"
    local label="$2"
    local target

    [[ -n "$path" ]] || die "$label path is empty."

    if [[ -L "$path" ]]; then
        target="$(readlink -f "$path" 2>/dev/null || true)"
        if [[ -z "$target" || ! -e "$target" || ! -r "$target" ]]; then
            # Dead symlink — try to repair before giving up.
            # There are two distinct failure modes:
            #
            # Mode A — chrom sizes target missing, but FASTA data is present:
            #   hg38.chrom.sizes -> /home/dustin/genomes/hg38.chrom.sizes  (dead)
            #   hg38.fa -> ../../../../data/.../hg38.fa  (resolves OK)
            #   Fix: samtools faidx + cut to regenerate the chrom sizes file.
            #
            # Mode B — FASTA data directory itself is missing:
            #   hg38.fa -> ../../../../data/.../hg38.fa  (also dead)
            #   This means refgenie has the alias skeleton but the actual
            #   downloaded data was never written or got deleted.
            #   Fix: refgenie pull hg38/fasta (re-downloads the full FASTA),
            #   then regenerate chrom sizes from the freshly downloaded file.

            local dead_target
            dead_target="$(readlink "$path" 2>/dev/null || echo UNKNOWN)"
            warn "$label symlink is broken: $path -> $dead_target"
            warn "Attempting automatic repair..."

            # Derive genome name and FASTA symlink path from the broken path.
            local link_dir genome_name fasta_link fasta_real
            link_dir="$(dirname "$path")"
            genome_name="$(basename "$path" .chrom.sizes)"
            fasta_link="$link_dir/${genome_name}.fa"

            # Resolve the FASTA symlink fully — it may itself be a symlink
            # chain into refgenie's data directory.
            fasta_real="$(readlink -f "$fasta_link" 2>/dev/null || true)"

            local repaired=false

            # Locate samtools directly in the conda env bin — avoids conda
            # activation issues that make `in_env samtools` fail silently
            # inside functions called before full env initialisation.
            local samtools_bin
            samtools_bin="$HOME/miniconda3/envs/${ENV_NAME}/bin/samtools"
            if [[ ! -x "$samtools_bin" ]]; then
                samtools_bin="$(command -v samtools 2>/dev/null || true)"
            fi

            # ── Mode B: FASTA data missing — pull the full fasta asset ──────
            if [[ -z "$fasta_real" || ! -f "$fasta_real" ]]; then
                warn "FASTA data is missing (refgenie alias exists but data directory is gone)."
                warn "Pulling hg38/fasta from refgenie — this will take a few minutes..."
                if in_env refgenie pull -c "$REFGENIE_CONFIG" "${genome_name}/fasta" 2>&1 | \
                       grep -v "^$" | sed 's/^/  /'; then
                    # After pull, re-resolve the FASTA symlink.
                    fasta_real="$(readlink -f "$fasta_link" 2>/dev/null || true)"
                    if [[ -n "$fasta_real" && -f "$fasta_real" ]]; then
                        ok "FASTA restored via refgenie pull: $fasta_real"
                    else
                        fasta_real=""
                    fi
                fi
            fi

            # ── Mode A: FASTA present, chrom sizes target missing ────────────
            if [[ -n "$fasta_real" && -f "$fasta_real" && -x "$samtools_bin" ]]; then
                warn "FASTA present — regenerating chrom sizes..."
                local target_dir
                target_dir="$(dirname "$dead_target")"
                mkdir -p "$target_dir" 2>/dev/null || true
                "$samtools_bin" faidx "$fasta_real" 2>/dev/null || true
                if [[ -f "${fasta_real}.fai" ]] && \
                   cut -f1,2 "${fasta_real}.fai" > "$dead_target" 2>/dev/null && \
                   [[ -s "$dead_target" ]]; then
                    ok "Chrom sizes regenerated: $dead_target"
                    repaired=true
                else
                    warn "samtools faidx ran but chrom sizes could not be written."
                fi
            elif [[ -n "$fasta_real" && -f "$fasta_real" && ! -x "$samtools_bin" ]]; then
                warn "FASTA present but samtools not found at $samtools_bin"
            fi

            # ── Final fallback: chrom_sizes-only pull ────────────────────────
            if ! $repaired; then
                warn "Trying refgenie pull for chrom_sizes asset only..."
                if in_env refgenie pull -c "$REFGENIE_CONFIG" \
                       "${genome_name}/fasta:chrom_sizes" 2>/dev/null; then
                    local new_target
                    new_target="$(readlink -f "$path" 2>/dev/null || true)"
                    if [[ -n "$new_target" && -e "$new_target" && -s "$new_target" ]]; then
                        ok "Chrom sizes restored via refgenie pull."
                        repaired=true
                    fi
                fi
            fi

            if ! $repaired; then
                err "Could not repair broken chrom sizes symlink automatically."
                err "  Symlink : $path"
                err "  Target  : $dead_target"
                err "  FASTA   : ${fasta_real:-not found}"
                err ""
                err "  To fix manually, run:"
                err "    conda run -n FetchPA refgenie pull -c $REFGENIE_CONFIG ${genome_name}/fasta"
                err "    mkdir -p $target_dir"
                err "    samtools faidx $fasta_link"
                err "    cut -f1,2 ${fasta_link}.fai > $dead_target"
                die "$label symlink is broken and could not be repaired."
            fi

            # Re-resolve after repair.
            target="$(readlink -f "$path" 2>/dev/null || true)"
        fi
    fi

    [[ -e "$path" ]] || die "$label does not exist: $path"
    [[ -f "$path" ]] || die "$label is not a regular file: $path"
    [[ -r "$path" ]] || die "$label is not readable: $path"
    [[ -s "$path" ]] || die "$label is empty: $path"
}


validate_chrom_sizes_file() {
    local chrom_sizes="$1"
    local genome_size

    validate_readable_reference_file "$chrom_sizes" "Chrom sizes"

    awk '
        NF < 2 {bad=1; exit}
        $2 !~ /^[0-9]+$/ {bad=1; exit}
        {sum += $2; n++}
        END {if (bad || n == 0 || sum <= 0) exit 1}
    ' "$chrom_sizes" || die "Chrom sizes file is malformed: $chrom_sizes"

    genome_size=$(awk '{sum+=$2} END {printf "%.0f", sum}' "$chrom_sizes")
    ok "Chrom sizes verified: $chrom_sizes (genome size: $genome_size)"
}

validate_bowtie2_index_prefix() {
    local genome_index="$1"
    [[ -n "$genome_index" ]] || die "Bowtie2 index prefix is empty."
    local index_files count
    if ! index_files="$(bowtie2_index_files "$genome_index")"; then
        die "Bowtie2 index validation failed for prefix: $genome_index"
    fi
    count="$(grep -c . <<< "$index_files")"
    [[ "$count" -eq 6 ]] || die "Bowtie2 index has $count components; expected exactly 6."
    ok "Bowtie2 index verified: $genome_index ($count complete components)"
}

validate_genome_reference_integrity() {
    local genome_index="$1"
    local chrom_sizes="$2"

    label "Validating sequence/alignment reference integrity for $GENOME..."
    validate_bowtie2_index_prefix "$genome_index"
    validate_chrom_sizes_file "$chrom_sizes"
}

# ─────────────────────────────────────────────────────────────
# Abort / Ctrl+C cleanup
# ─────────────────────────────────────────────────────────────

CURRENT_PGID=""
CURRENT_CHILD_PID=""
CURRENT_SAMPLE=""
ABORTING=false

process_group_alive() {
    local pgid="${1:-}"
    [[ -n "$pgid" ]] || return 1
    kill -0 -- "-$pgid" >/dev/null 2>&1
}

terminate_process_group() {
    local pgid="${1:-}"
    [[ -n "$pgid" ]] || return 0

    if process_group_alive "$pgid"; then
        warn "Sending TERM to active process group: $pgid"
        kill -TERM -- "-$pgid" >/dev/null 2>&1 || true

        local n
        for n in {1..10}; do
            process_group_alive "$pgid" || return 0
            sleep 1
        done

        if process_group_alive "$pgid"; then
            warn "Process group $pgid is still running. Sending KILL."
            kill -KILL -- "-$pgid" >/dev/null 2>&1 || true
        fi
    fi
}

terminate_direct_child_tree() {
    local child="${1:-}"
    [[ -n "$child" ]] || return 0

    if kill -0 "$child" >/dev/null 2>&1; then
        warn "Stopping child process tree rooted at PID: $child"
        pkill -TERM -P "$child" >/dev/null 2>&1 || true
        kill -TERM "$child" >/dev/null 2>&1 || true
        sleep 3
        pkill -KILL -P "$child" >/dev/null 2>&1 || true
        kill -KILL "$child" >/dev/null 2>&1 || true
    fi
}

terminate_output_dir_processes() {
    local out="${OUTPUT_DIR:-}"
    [[ -n "$out" ]] || return 0

    local pid args killed_any=false

    while read -r pid args; do
        [[ -n "$pid" ]] || continue
        [[ "$pid" == "$$" ]] && continue
        [[ "${args:-}" == *"$out"* ]] || continue

        killed_any=true
        warn "Stopping run-related orphan PID $pid"
        kill -TERM "$pid" >/dev/null 2>&1 || true
    done < <(ps -eo pid=,args= 2>/dev/null || true)

    if $killed_any; then
        sleep 3
        while read -r pid args; do
            [[ -n "$pid" ]] || continue
            [[ "$pid" == "$$" ]] && continue
            [[ "${args:-}" == *"$out"* ]] || continue

            warn "Force-stopping stubborn run-related PID $pid"
            kill -KILL "$pid" >/dev/null 2>&1 || true
        done < <(ps -eo pid=,args= 2>/dev/null || true)
    fi
}

abort_cleanup() {
    local signal="${1:-INT}"

    if $ABORTING; then
        exit 130
    fi
    ABORTING=true

    trap - INT TERM HUP

    blank
    err "Abort requested ($signal). Stopping this run before exiting..."

    if [[ -n "${CURRENT_SAMPLE:-}" ]]; then
        warn "Active sample at abort: $CURRENT_SAMPLE"
    fi

    # Under concurrent sample processing, each background subshell tracks its
    # own child via a per-sample PID file (CURRENT_CHILD_PID/CURRENT_PGID only
    # ever cover a single foreground child and are subshell-local otherwise).
    local pid_file pid
    if [[ -n "${LOG_DIR:-}" && -d "$LOG_DIR" ]]; then
        for pid_file in "$LOG_DIR"/*.pid; do
            [[ -f "$pid_file" ]] || continue
            pid=$(cat "$pid_file" 2>/dev/null || true)
            [[ -n "$pid" ]] || continue
            warn "Stopping tracked sample process: PID $pid ($(basename "$pid_file" .pid))"
            terminate_process_group "$pid"
            terminate_direct_child_tree "$pid"
        done
    fi

    terminate_process_group "${CURRENT_PGID:-}"
    terminate_direct_child_tree "${CURRENT_CHILD_PID:-}"
    terminate_output_dir_processes
    cleanup_custom_reference_staging
    _release_custom_genome_build_lock

    blank
    err "Run aborted. Active child processes from this run were stopped."
    exit 130
}

trap 'abort_cleanup INT' INT
trap 'abort_cleanup TERM' TERM
trap 'abort_cleanup HUP' HUP
trap 'cleanup_all_staging' EXIT

# run_sample_command_with_cleanup LOG_FILE [--quiet] -- CMD [ARGS...]
# --quiet suppresses the live filtered terminal echo (used under concurrent
# sample processing, where interleaved output from multiple samples at once
# would be unreadable) -- the full command output is still written to
# LOG_FILE either way. Also tracks the child's PID in a per-sample .pid file
# next to LOG_FILE so abort_cleanup() can find and stop it even when this
# runs inside a backgrounded subshell (CURRENT_CHILD_PID/CURRENT_PGID are
# subshell-local in that case and never reach the parent's signal trap).
run_sample_command_with_cleanup() {
    local log_file="$1"
    shift

    local quiet=false
    if [[ "${1:-}" == "--quiet" ]]; then
        quiet=true
        shift
    fi

    local pid_file="${log_file%.log}.pid"
    local exit_code

    if $quiet; then
        if command_exists setsid; then
            setsid "$@" > "$log_file" 2>&1 &
            CURRENT_CHILD_PID="$!"
            CURRENT_PGID="$CURRENT_CHILD_PID"
        else
            "$@" > "$log_file" 2>&1 &
            CURRENT_CHILD_PID="$!"
            CURRENT_PGID=""
        fi
    elif command_exists setsid; then
        setsid "$@" > >(tee "$log_file" | \
            grep --line-buffered -E \
                "INFO|WARNING|ERROR|Trimming|Aligning|Calling peaks|Motif|done|FAIL|PASS|error" | \
            while IFS= read -r line; do
                echo -e "    ${DIM}$line${RESET}"
            done) 2>&1 &
        CURRENT_CHILD_PID="$!"
        CURRENT_PGID="$CURRENT_CHILD_PID"
    else
        warn "setsid was not found. Ctrl+C cleanup will use a less robust fallback."
        "$@" > >(tee "$log_file" | \
            grep --line-buffered -E \
                "INFO|WARNING|ERROR|Trimming|Aligning|Calling peaks|Motif|done|FAIL|PASS|error" | \
            while IFS= read -r line; do
                echo -e "    ${DIM}$line${RESET}"
            done) 2>&1 &
        CURRENT_CHILD_PID="$!"
        CURRENT_PGID=""
    fi

    echo "$CURRENT_CHILD_PID" > "$pid_file" 2>/dev/null || true

    # Capturing wait's result via if/else -- not a bare statement -- is what
    # exempts it from errexit. A nonzero wait as a plain top-level statement
    # under `set -e` would abort this entire function (and everything that
    # called it) the instant the tracked command failed, before exit_code
    # was ever assigned.
    if wait "$CURRENT_CHILD_PID"; then
        exit_code=0
    else
        exit_code="$?"
    fi

    rm -f "$pid_file"
    CURRENT_PGID=""
    CURRENT_CHILD_PID=""

    return "$exit_code"
}

# ─────────────────────────────────────────────────────────────
# General helpers
# ─────────────────────────────────────────────────────────────

is_supported_genome() {
    local query="$1"
    local g
    for g in "${SUPPORTED_GENOMES[@]}"; do
        [[ "$g" == "$query" ]] && return 0
    done
    return 1
}

# is_valid_custom_genome_name NAME
# User-defined genome profile names become directory names (registry, local FASTA/index
# cache) and are echoed into shell/R scripts elsewhere in the pipeline, so
# the charset is restricted the same way group/contrast labels are
# elsewhere in this codebase.
is_valid_custom_genome_name() {
    local name="$1"
    [[ -n "$name" ]] || return 1
    [[ "$name" =~ ^[A-Za-z0-9_]+$ ]] || return 1
    return 0
}

# ─────────────────────────────────────────────────────────────
# User-defined genome profile registry and immutable profile store
#
# Mutable state is limited to:
#   <root>/<assembly>/genome.conf
#   <root>/<assembly>/current -> assemblies/<fasta_sha256>
#
# Immutable data live under:
#   <root>/<assembly>/assemblies/<fasta_sha256>/
#       genome.fa
#       genome.fa.fai
#       genome.chrom.sizes
#       bowtie2/<assembly>.*
#       annotations/<annotation_fingerprint>/
#           txdb.sqlite
#           orgdb.sqlite
#           tss.bed
#           genes.bed
#           exons.bed
#           introns.bed
#           promoters.bed
#           features.bed
#           profile_build.conf
#
# The annotation fingerprint is derived from the frozen TxDb/OrgDb bytes and
# every generated BED hash. Existing hash-specific data are never overwritten.
# ─────────────────────────────────────────────────────────────

CUSTOM_REFERENCE_REGISTRY_SCHEMA=3
ACTIVE_CUSTOM_STAGE=""
CUSTOM_SELECTED_REGISTRY_SHA256=""
CUSTOM_REFERENCE_FROM_SNAPSHOT=false

custom_genome_registry_file() { printf '%s/%s/genome.conf\n' "$CUSTOM_GENOME_REGISTRY_ROOT" "$1"; }
custom_genome_current_link() { printf '%s/%s/current\n' "$CUSTOM_GENOME_REGISTRY_ROOT" "$1"; }
custom_genome_assembly_dir() { printf '%s/%s/assemblies/%s\n' "$CUSTOM_GENOME_REGISTRY_ROOT" "$1" "$2"; }

ucsc_profile_fasta_url() { printf 'https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.fa.gz\n' "$1" "$1"; }
ucsc_profile_chrom_sizes_url() { printf 'https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes\n' "$1" "$1"; }

sha256_of() {
    local f="$1"
    [[ -f "$f" ]] || { printf '\n'; return 0; }
    command_exists sha256sum || die "sha256sum is required for custom-reference integrity checks."
    sha256sum "$f" | awk '{print $1}'
}
require_sha256() {
    local f="$1" label_text="$2" digest
    digest="$(sha256_of "$f")"
    [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || die "Could not compute SHA-256 for $label_text: $f"
    printf '%s\n' "$digest"
}
verify_file_sha256() {
    local f="$1" expected="$2" label_text="$3" actual
    [[ -n "$expected" ]] || die "No registered SHA-256 is available for $label_text: $f"
    [[ -f "$f" ]] || die "$label_text is missing: $f"
    actual="$(require_sha256 "$f" "$label_text")"
    [[ "$actual" == "$expected" ]] || die "$label_text hash does not match the registry."$'\n'"       File: $f"$'\n'"       Registered: $expected"$'\n'"       Actual:     $actual"$'\n'"       Treat this as a different reference version; do not repair it in place."
}

bowtie2_index_files() {
    local prefix="$1"
    local -a small=("${prefix}.1.bt2" "${prefix}.2.bt2" "${prefix}.3.bt2" "${prefix}.4.bt2" "${prefix}.rev.1.bt2" "${prefix}.rev.2.bt2")
    local -a large=("${prefix}.1.bt2l" "${prefix}.2.bt2l" "${prefix}.3.bt2l" "${prefix}.4.bt2l" "${prefix}.rev.1.bt2l" "${prefix}.rev.2.bt2l")
    local f small_ok=true large_ok=true any_small=false any_large=false
    for f in "${small[@]}"; do [[ -e "$f" ]] && any_small=true; [[ -s "$f" ]] || small_ok=false; done
    for f in "${large[@]}"; do [[ -e "$f" ]] && any_large=true; [[ -s "$f" ]] || large_ok=false; done
    if $small_ok && ! $any_large; then printf '%s\n' "${small[@]}"
    elif $large_ok && ! $any_small; then printf '%s\n' "${large[@]}"
    elif $small_ok && $large_ok; then die "Bowtie2 index has both complete .bt2 and .bt2l sets: $prefix"
    else
        shopt -s nullglob; local existing=("${prefix}".*.bt2 "${prefix}".*.bt2l); shopt -u nullglob
        die "Bowtie2 index is partial or corrupted for prefix: $prefix (found ${#existing[@]} component files; expected one complete six-file set)."
    fi
}
bowtie2_index_sha256() {
    local prefix="$1" tmp f index_files digest
    tmp="$(mktemp)"
    index_files="$(bowtie2_index_files "$prefix")"
    while IFS= read -r f; do [[ -n "$f" ]] && printf '%s  %s\n' "$(require_sha256 "$f" "Bowtie2 index component")" "$(basename "$f")"; done <<< "$index_files" | sort -k2,2 > "$tmp"
    digest="$(require_sha256 "$tmp" "Bowtie2 index digest manifest")"; rm -f "$tmp"; printf '%s\n' "$digest"
}
verify_bowtie2_index_sha256() {
    local prefix="$1" expected="$2" actual
    [[ -n "$expected" ]] || die "No registered Bowtie2 index SHA-256 is available for: $prefix"
    actual="$(bowtie2_index_sha256 "$prefix")"
    [[ "$actual" == "$expected" ]] || die "Cached Bowtie2 index hash does not match the registry."$'\n'"       Prefix: $prefix"$'\n'"       Registered: $expected"$'\n'"       Actual:     $actual"
}

_acquire_custom_genome_build_lock() {
    local genome="$1" lock_file
    mkdir -p "$CUSTOM_GENOME_LOCK_ROOT"; lock_file="$CUSTOM_GENOME_LOCK_ROOT/${genome}.lock"
    command_exists flock || die "flock is required but was not found. Run install.sh first."
    eval "exec 200>\"$lock_file\""
    flock -w 7200 200 || die "Timed out after 2 hours waiting for the custom-reference lock: $lock_file"
}
_release_custom_genome_build_lock() { flock -u 200 2>/dev/null || true; exec 200>&- 2>/dev/null || true; }
cleanup_custom_reference_staging() { [[ -n "${ACTIVE_CUSTOM_STAGE:-}" && -d "$ACTIVE_CUSTOM_STAGE" ]] && rm -rf "$ACTIVE_CUSTOM_STAGE"; ACTIVE_CUSTOM_STAGE=""; }

# ACTIVE_BUILTIN_ANNOTATION_STAGE mirrors ACTIVE_CUSTOM_STAGE above -- set
# right before staging a built-in genome's annotation build, cleared right
# after it's either removed (failure/race loss) or moved into place
# (success). cleanup_all_staging (registered on the trap below in place of
# cleanup_custom_reference_staging alone) removes whichever of the two is
# still set if the script is interrupted mid-build, so Ctrl+C doesn't leave
# an abandoned annotation_<genome>_* directory under .staging/.
ACTIVE_BUILTIN_ANNOTATION_STAGE=""
cleanup_builtin_annotation_staging() { [[ -n "${ACTIVE_BUILTIN_ANNOTATION_STAGE:-}" && -d "$ACTIVE_BUILTIN_ANNOTATION_STAGE" ]] && rm -rf "$ACTIVE_BUILTIN_ANNOTATION_STAGE"; ACTIVE_BUILTIN_ANNOTATION_STAGE=""; }
cleanup_all_staging() { cleanup_custom_reference_staging; cleanup_builtin_annotation_staging; }

assert_custom_registry_unchanged_since_selection() {
    local genome="$1" reg_file current="__ABSENT__" expected="${CUSTOM_SELECTED_REGISTRY_SHA256:-__ABSENT__}"
    reg_file="$(custom_genome_registry_file "$genome")"
    if [[ -f "$reg_file" ]]; then current="$(require_sha256 "$reg_file" "current genome-profile registry")"; fi
    [[ "$current" == "$expected" ]] || die "The user-defined genome profile registry changed after you reviewed it."$'\n'"       Profile: $genome"$'\n'"       Expected state: $expected"$'\n'"       Current state:  $current"$'\n'"       Re-run and review the current profile before publishing or reusing it."
}

load_custom_genome_registry() {
    local genome="$1" reg_file line key val
    reg_file="$(custom_genome_registry_file "$genome")"; [[ -f "$reg_file" ]] || return 1
    local -A rg=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        key="${line%%=*}"; val="${line#*=}"; [[ "$key" =~ ^[A-Za-z0-9_]+$ ]] && rg["$key"]="$val"
    done < "$reg_file"
    CUSTOM_REG_SCHEMA="${rg[CUSTOM_REFERENCE_REGISTRY_SCHEMA]:-1}"
    [[ "$CUSTOM_REG_SCHEMA" =~ ^[0-9]+$ ]] || die "Custom registry has an invalid schema value: $CUSTOM_REG_SCHEMA"
    (( CUSTOM_REG_SCHEMA <= CUSTOM_REFERENCE_REGISTRY_SCHEMA )) || die "Custom registry schema $CUSTOM_REG_SCHEMA is newer than this runner understands (maximum $CUSTOM_REFERENCE_REGISTRY_SCHEMA)."
    CUSTOM_REG_UCSC_ASSEMBLY="${rg[CUSTOM_UCSC_ASSEMBLY]:-${rg[CUSTOM_GENOME_NAME]:-$genome}}"
    CUSTOM_REG_FASTA_SOURCE="${rg[CUSTOM_FASTA_SOURCE]:-${rg[CUSTOM_FASTA]:-}}"
    CUSTOM_REG_CHROM_SOURCE="${rg[CUSTOM_CHROM_SIZES_SOURCE]:-}"
    CUSTOM_REG_FASTA_CACHE="${rg[CUSTOM_FASTA_CACHE]:-${rg[CUSTOM_FASTA]:-}}"
    CUSTOM_REG_FASTA_SHA256="${rg[CUSTOM_FASTA_SHA256]:-}"
    CUSTOM_REG_CHROM_SIZES="${rg[CUSTOM_CHROM_SIZES]:-}"
    CUSTOM_REG_CHROM_SIZES_SHA256="${rg[CUSTOM_CHROM_SIZES_SHA256]:-}"
    CUSTOM_REG_BT2_INDEX="${rg[CUSTOM_BT2_INDEX]:-}"
    CUSTOM_REG_BT2_INDEX_SHA256="${rg[CUSTOM_BT2_INDEX_SHA256]:-}"
    CUSTOM_REG_TXDB_PKG="${rg[CUSTOM_TXDB_PACKAGE]:-}"
    CUSTOM_REG_TXDB_VERSION="${rg[CUSTOM_TXDB_PACKAGE_VERSION]:-}"
    CUSTOM_REG_ORGDB_PKG="${rg[CUSTOM_ORGDB_PACKAGE]:-${rg[CUSTOM_ORGDB_PKG]:-}}"
    CUSTOM_REG_ORGDB_VERSION="${rg[CUSTOM_ORGDB_PACKAGE_VERSION]:-}"
    CUSTOM_REG_KEGG_ORG="${rg[CUSTOM_KEGG_ORG]:-}"
    CUSTOM_REG_TXDB_SQLITE="${rg[CUSTOM_TXDB_SQLITE]:-}"
    CUSTOM_REG_TXDB_SQLITE_SHA256="${rg[CUSTOM_TXDB_SQLITE_SHA256]:-}"
    CUSTOM_REG_ORGDB_SQLITE="${rg[CUSTOM_ORGDB_SQLITE]:-}"
    CUSTOM_REG_ORGDB_SQLITE_SHA256="${rg[CUSTOM_ORGDB_SQLITE_SHA256]:-}"
    CUSTOM_REG_TSS_BED="${rg[CUSTOM_TSS_BED]:-}"; CUSTOM_REG_TSS_BED_SHA256="${rg[CUSTOM_TSS_BED_SHA256]:-}"
    CUSTOM_REG_FEATURE_BED="${rg[CUSTOM_FEATURE_BED]:-}"; CUSTOM_REG_FEATURE_BED_SHA256="${rg[CUSTOM_FEATURE_BED_SHA256]:-}"
    CUSTOM_REG_GENES_BED="${rg[CUSTOM_GENES_BED]:-}"; CUSTOM_REG_GENES_BED_SHA256="${rg[CUSTOM_GENES_BED_SHA256]:-}"
    CUSTOM_REG_EXONS_BED="${rg[CUSTOM_EXONS_BED]:-}"; CUSTOM_REG_EXONS_BED_SHA256="${rg[CUSTOM_EXONS_BED_SHA256]:-}"
    CUSTOM_REG_INTRONS_BED="${rg[CUSTOM_INTRONS_BED]:-}"; CUSTOM_REG_INTRONS_BED_SHA256="${rg[CUSTOM_INTRONS_BED_SHA256]:-}"
    CUSTOM_REG_PROMOTERS_BED="${rg[CUSTOM_PROMOTERS_BED]:-}"; CUSTOM_REG_PROMOTERS_BED_SHA256="${rg[CUSTOM_PROMOTERS_BED_SHA256]:-}"
    CUSTOM_REG_ANNOTATION_FINGERPRINT="${rg[CUSTOM_ANNOTATION_FINGERPRINT]:-}"
    CUSTOM_REG_PROFILE_DIR="${rg[CUSTOM_PROFILE_DIR]:-}"
    CUSTOM_REG_PROFILE_BUILD_CONF="${rg[CUSTOM_PROFILE_BUILD_CONF]:-}"; CUSTOM_REG_PROFILE_BUILD_CONF_SHA256="${rg[CUSTOM_PROFILE_BUILD_CONF_SHA256]:-}"
    CUSTOM_REG_GENOME_SIZE="${rg[CUSTOM_GENOME_SIZE]:-}"
    CUSTOM_REG_EFFECTIVE_GENOME_SIZE="${rg[EFFECTIVE_GENOME_SIZE]:-}"
    CUSTOM_REG_GENOME_SIZE_METHOD="${rg[GENOME_SIZE_METHOD]:-}"
    CUSTOM_REG_ASSEMBLY_DIR="${rg[CUSTOM_ASSEMBLY_DIR]:-}"
    CUSTOM_REG_SOURCE_FILE="$reg_file"
    CUSTOM_REGISTRY_SHA256="$(require_sha256 "$reg_file" "user-defined genome profile registry")"
}

restore_registered_file_from_source() {
    local cache="$1" source="$2" expected_sha="$3" label_text="$4" restore_dir tmp actual
    [[ ! -e "$cache" ]] || die "Internal error: restore requested for an existing file: $cache"
    [[ -n "$source" ]] || die "$label_text is missing and no source was recorded: $cache"
    restore_dir="$(dirname "$cache")"; mkdir -p "$restore_dir"; tmp="${restore_dir}/.$(basename "$cache").restore.$$"; rm -f "$tmp"
    _fetch_sequence_asset "$source" "$tmp" "$label_text restoration"
    actual="$(require_sha256 "$tmp" "$label_text restoration")"
    [[ "$actual" == "$expected_sha" ]] || { rm -f "$tmp"; die "The recorded source for $label_text no longer matches the registered hash."$'\n'"       Source: $source"$'\n'"       Registered: $expected_sha"$'\n'"       Downloaded: $actual"$'\n'"       Register it as a new reference version instead of repairing the old version."; }
    mv "$tmp" "$cache"; ok "Restored missing $label_text with the registered hash: $cache"
}
verify_or_restore_registered_file() { local cache="$1" source="$2" expected="$3" label="$4"; [[ -f "$cache" ]] && verify_file_sha256 "$cache" "$expected" "$label" || restore_registered_file_from_source "$cache" "$source" "$expected" "$label"; }

_write_profile_fields() {
    local genome="$1"
    printf 'CUSTOM_UCSC_ASSEMBLY=%s\n' "${CUSTOM_UCSC_ASSEMBLY[$genome]:-$genome}"
    printf 'CUSTOM_TXDB_PACKAGE=%s\n' "${CUSTOM_TXDB_PKG[$genome]:-}"
    printf 'CUSTOM_TXDB_PACKAGE_VERSION=%s\n' "${CUSTOM_TXDB_VERSION[$genome]:-}"
    printf 'CUSTOM_ORGDB_PACKAGE=%s\n' "${CUSTOM_ORGDB_PKG[$genome]:-}"
    printf 'CUSTOM_ORGDB_PACKAGE_VERSION=%s\n' "${CUSTOM_ORGDB_VERSION[$genome]:-}"
    printf 'CUSTOM_KEGG_ORG=%s\n' "${CUSTOM_KEGG_ORG[$genome]:-}"
    printf 'CUSTOM_PROFILE_DIR=%s\n' "${CUSTOM_PROFILE_DIR[$genome]:-}"
    printf 'CUSTOM_ANNOTATION_FINGERPRINT=%s\n' "${CUSTOM_ANNOTATION_FINGERPRINT[$genome]:-}"
    local key
    for key in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
        local path_var="CUSTOM_${key}[$genome]" sha_var="CUSTOM_${key}_SHA256[$genome]"
        printf 'CUSTOM_%s=%s\n' "$key" "${!path_var:-}"
        printf 'CUSTOM_%s_SHA256=%s\n' "$key" "${!sha_var:-}"
    done
}

write_custom_genome_registry() {
    local genome="$1"
    local reg_dir="$CUSTOM_GENOME_REGISTRY_ROOT/$genome" reg_file tmp_file current_link tmp_link
    reg_file="$(custom_genome_registry_file "$genome")"; tmp_file="${reg_file}.tmp.$$"; mkdir -p "$reg_dir"
    {
        echo "# Auto-generated by PEPATAC_run.sh -- do not hand-edit."
        echo "# Mutable pointer to immutable validated genome-profile assets."
        printf 'CUSTOM_REFERENCE_REGISTRY_SCHEMA=%s\n' "$CUSTOM_REFERENCE_REGISTRY_SCHEMA"
        printf 'CUSTOM_GENOME_NAME=%s\n' "$genome"
        printf 'CUSTOM_ASSEMBLY_DIR=%s\n' "${CUSTOM_ASSEMBLY_DIR:-}"
        printf 'CUSTOM_FASTA_SOURCE=%s\n' "${CUSTOM_FASTA_URL[$genome]:-}"
        printf 'CUSTOM_CHROM_SIZES_SOURCE=%s\n' "${CUSTOM_CHROMSIZES_URL[$genome]:-}"
        printf 'CUSTOM_FASTA_CACHE=%s\n' "${LOCAL_CUSTOM_FASTA:-}"
        printf 'CUSTOM_FASTA_SHA256=%s\n' "${LOCAL_CUSTOM_FASTA_SHA256:-}"
        printf 'CUSTOM_CHROM_SIZES=%s\n' "${LOCAL_CHROM_SIZES:-}"
        printf 'CUSTOM_CHROM_SIZES_SHA256=%s\n' "${LOCAL_CHROM_SIZES_SHA256:-}"
        printf 'CUSTOM_BT2_INDEX=%s\n' "${LOCAL_BT2_INDEX:-}"
        printf 'CUSTOM_BT2_INDEX_SHA256=%s\n' "${LOCAL_BT2_INDEX_SHA256:-}"
        _write_profile_fields "$genome"
        printf 'CUSTOM_GENOME_SIZE=%s\n' "${CUSTOM_GENOME_SIZE[$genome]:-}"
        printf 'EFFECTIVE_GENOME_SIZE=%s\n' "${RUN_GENOME_SIZE:-}"
        printf 'GENOME_SIZE_METHOD=%s\n' "${GENOME_SIZE_METHOD:-}"
        printf 'REGISTERED=%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    } > "$tmp_file"
    mv "$tmp_file" "$reg_file"
    current_link="$(custom_genome_current_link "$genome")"; tmp_link="${current_link}.tmp.$$"; rm -f "$tmp_link"
    ln -s "assemblies/${LOCAL_CUSTOM_FASTA_SHA256}" "$tmp_link"; mv -Tf "$tmp_link" "$current_link"
    ok "Validated genome profile registry updated atomically: $reg_file"
}

write_reference_snapshot() {
    local genome="$1" snapshot_dir="$OUTPUT_DIR/reference_snapshot" snapshot_file tmp_file
    snapshot_file="$snapshot_dir/genome.conf"; tmp_file="${snapshot_file}.tmp.$$"; mkdir -p "$snapshot_dir"
    {
        echo "# Immutable validated genome-profile snapshot for this PEPATAC run."
        printf 'REFERENCE_SNAPSHOT_SCHEMA=%s\n' "$REFERENCE_SNAPSHOT_SCHEMA"
        printf 'CUSTOM_GENOME_NAME=%s\n' "$genome"
        printf 'CUSTOM_ASSEMBLY_DIR=%s\n' "${CUSTOM_ASSEMBLY_DIR:-}"
        printf 'CUSTOM_FASTA_SOURCE=%s\n' "${CUSTOM_FASTA_URL[$genome]:-}"
        printf 'CUSTOM_CHROM_SIZES_SOURCE=%s\n' "${CUSTOM_CHROMSIZES_URL[$genome]:-}"
        printf 'CUSTOM_FASTA_CACHE=%s\n' "${LOCAL_CUSTOM_FASTA:-}"
        printf 'CUSTOM_FASTA_SHA256=%s\n' "${LOCAL_CUSTOM_FASTA_SHA256:-}"
        printf 'CUSTOM_CHROM_SIZES=%s\n' "${LOCAL_CHROM_SIZES:-}"
        printf 'CUSTOM_CHROM_SIZES_SHA256=%s\n' "${LOCAL_CHROM_SIZES_SHA256:-}"
        printf 'CUSTOM_BT2_INDEX=%s\n' "${LOCAL_BT2_INDEX:-}"
        printf 'CUSTOM_BT2_INDEX_SHA256=%s\n' "${LOCAL_BT2_INDEX_SHA256:-}"
        _write_profile_fields "$genome"
        printf 'EFFECTIVE_GENOME_SIZE=%s\n' "${RUN_GENOME_SIZE:-}"
        printf 'GENOME_SIZE_METHOD=%s\n' "${GENOME_SIZE_METHOD:-}"
        # Schema 4: genome-agnostic annotation fields, written for every
        # genome (built-in or custom) from the already-resolved
        # RUN_ANNOTATION_* scalars -- see resolve_run_annotation_assets().
        # Deliberately not named CUSTOM_* like the block above: a built-in
        # genome's frozen TxDb/OrgDb has nothing to do with the
        # custom-profile system, and reusing that prefix here would be
        # exactly the misleading naming this schema bump exists to avoid.
        printf 'ANNOTATION_STATUS=%s\n' "${RUN_ANNOTATION_STATUS:-}"
        printf 'ANNOTATION_SOURCE=%s\n' "${RUN_ANNOTATION_SOURCE:-}"
        printf 'ANNOTATION_FINGERPRINT=%s\n' "${RUN_ANNOTATION_FINGERPRINT:-}"
        printf 'ANNOTATION_BUILDER_SCHEMA=%s\n' "${ANNOTATION_BUILDER_SCHEMA:-}"
        printf 'ANNOTATION_QC_FAILURE_REASON=%s\n' "${RUN_ANNOTATION_QC_FAILURE_REASON:-}"
        printf 'ANNOTATION_TXDB_PACKAGE=%s\n' "${RUN_ANNOTATION_TXDB_PACKAGE:-}"
        printf 'ANNOTATION_TXDB_PACKAGE_VERSION=%s\n' "${RUN_ANNOTATION_TXDB_VERSION:-}"
        printf 'ANNOTATION_TXDB_SQLITE=%s\n' "${RUN_ANNOTATION_TXDB_SQLITE:-}"
        printf 'ANNOTATION_TXDB_SQLITE_SHA256=%s\n' "${RUN_ANNOTATION_TXDB_SQLITE_SHA256:-}"
        printf 'ANNOTATION_ORGDB_PACKAGE=%s\n' "${RUN_ANNOTATION_ORGDB_PACKAGE:-}"
        printf 'ANNOTATION_ORGDB_PACKAGE_VERSION=%s\n' "${RUN_ANNOTATION_ORGDB_VERSION:-}"
        printf 'ANNOTATION_ORGDB_SQLITE=%s\n' "${RUN_ANNOTATION_ORGDB_SQLITE:-}"
        printf 'ANNOTATION_ORGDB_SQLITE_SHA256=%s\n' "${RUN_ANNOTATION_ORGDB_SQLITE_SHA256:-}"
        printf 'ANNOTATION_TSS_BED=%s\n' "${RUN_ANNOTATION_TSS_BED:-}"
        printf 'ANNOTATION_TSS_BED_SHA256=%s\n' "${RUN_ANNOTATION_TSS_BED_SHA256:-}"
        printf 'ANNOTATION_FEATURE_BED=%s\n' "${RUN_ANNOTATION_FEATURE_BED:-}"
        printf 'ANNOTATION_FEATURE_BED_SHA256=%s\n' "${RUN_ANNOTATION_FEATURE_BED_SHA256:-}"
        printf 'ANNOTATION_CHROM_SIZES=%s\n' "${RUN_ANNOTATION_CHROM_SIZES:-}"
        printf 'ANNOTATION_CHROM_SIZES_SHA256=%s\n' "${RUN_ANNOTATION_CHROM_SIZES_SHA256:-}"
        printf 'ANNOTATION_PROMOTER_UPSTREAM_BP=%s\n' "${ANNOTATION_PROMOTER_UPSTREAM_BP:-}"
        printf 'ANNOTATION_PROMOTER_DOWNSTREAM_BP=%s\n' "${ANNOTATION_PROMOTER_DOWNSTREAM_BP:-}"
        printf 'RUNNER_VERSION=%s\n' "$RUNNER_VERSION"
        printf 'SNAPSHOT_WRITTEN=%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    } > "$tmp_file"
    mv "$tmp_file" "$snapshot_file"; ok "Immutable reference/profile snapshot written: $snapshot_file"
}

_load_profile_fields_from_assoc() {
    local genome="$1" assoc_name="$2"; local -n a="$assoc_name"
    CUSTOM_UCSC_ASSEMBLY["$genome"]="${a[CUSTOM_UCSC_ASSEMBLY]:-$genome}"
    CUSTOM_TXDB_PKG["$genome"]="${a[CUSTOM_TXDB_PACKAGE]:-}"; CUSTOM_TXDB_VERSION["$genome"]="${a[CUSTOM_TXDB_PACKAGE_VERSION]:-}"
    CUSTOM_ORGDB_PKG["$genome"]="${a[CUSTOM_ORGDB_PACKAGE]:-${a[CUSTOM_ORGDB_PKG]:-}}"; CUSTOM_ORGDB_VERSION["$genome"]="${a[CUSTOM_ORGDB_PACKAGE_VERSION]:-}"
    CUSTOM_KEGG_ORG["$genome"]="${a[CUSTOM_KEGG_ORG]:-}"; CUSTOM_PROFILE_DIR["$genome"]="${a[CUSTOM_PROFILE_DIR]:-}"; CUSTOM_ANNOTATION_FINGERPRINT["$genome"]="${a[CUSTOM_ANNOTATION_FINGERPRINT]:-}"
    local key
    for key in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
        local path_var="CUSTOM_${key}[$genome]" sha_var="CUSTOM_${key}_SHA256[$genome]"
        printf -v "$path_var" '%s' "${a[CUSTOM_${key}]:-}"
        printf -v "$sha_var" '%s' "${a[CUSTOM_${key}_SHA256]:-}"
    done
}

load_reference_snapshot_into_runner() {
    local do_verify="${1:-true}" snapshot_file="$OUTPUT_DIR/reference_snapshot/genome.conf" line key val
    [[ -f "$snapshot_file" ]] || return 1
    local -A ss=()
    while IFS= read -r line || [[ -n "$line" ]]; do [[ -z "$line" || "$line" == \#* ]] && continue; key="${line%%=*}"; val="${line#*=}"; [[ "$key" =~ ^[A-Za-z0-9_]+$ ]] && ss["$key"]="$val"; done < "$snapshot_file"
    local schema="${ss[REFERENCE_SNAPSHOT_SCHEMA]:-0}"
    [[ "$schema" =~ ^[0-9]+$ ]] || die "Reference snapshot has an invalid schema value: $schema"
    (( schema <= REFERENCE_SNAPSHOT_SCHEMA )) || die "Reference snapshot schema $schema is newer than this runner understands (maximum $REFERENCE_SNAPSHOT_SCHEMA)."
    # Explicit enumeration rather than a loose ">= 3": a future schema bump
    # that changes field *meaning* rather than just adding fields should
    # have to be deliberately added here, not silently accepted because it
    # happens to be numerically higher than 3. This function only ever
    # reads CUSTOM_* fields (unchanged across schema 3 -> 4), so both are
    # equally valid inputs for it specifically.
    case "$schema" in
        3|4) ;;
        *) return 1 ;;
    esac
    LOCAL_CUSTOM_FASTA="${ss[CUSTOM_FASTA_CACHE]:-}"; LOCAL_CUSTOM_FASTA_SHA256="${ss[CUSTOM_FASTA_SHA256]:-}"
    LOCAL_CHROM_SIZES="${ss[CUSTOM_CHROM_SIZES]:-}"; LOCAL_CHROM_SIZES_SHA256="${ss[CUSTOM_CHROM_SIZES_SHA256]:-}"
    LOCAL_BT2_INDEX="${ss[CUSTOM_BT2_INDEX]:-}"; LOCAL_BT2_INDEX_SHA256="${ss[CUSTOM_BT2_INDEX_SHA256]:-}"
    CUSTOM_FASTA_URL["$GENOME"]="${ss[CUSTOM_FASTA_SOURCE]:-}"; CUSTOM_CHROMSIZES_URL["$GENOME"]="${ss[CUSTOM_CHROM_SIZES_SOURCE]:-}"
    _load_profile_fields_from_assoc "$GENOME" ss
    local snapshot_assembly="${CUSTOM_UCSC_ASSEMBLY[$GENOME]:-$GENOME}"
    [[ "${CUSTOM_FASTA_URL[$GENOME]}" == "$(ucsc_profile_fasta_url "$snapshot_assembly")" && \
       "${CUSTOM_CHROMSIZES_URL[$GENOME]}" == "$(ucsc_profile_chrom_sizes_url "$snapshot_assembly")" ]] \
        || die "This run snapshot is not tied to the official UCSC sequence sources for '$snapshot_assembly'."
    [[ "${CUSTOM_TXDB_PKG[$GENOME]:-}" == *".UCSC.${snapshot_assembly}."* ]] \
        || die "This run snapshot's TxDb package does not identify UCSC assembly '$snapshot_assembly'."
    CUSTOM_ASSEMBLY_DIR="${ss[CUSTOM_ASSEMBLY_DIR]:-}"; RUN_GENOME_SIZE="${ss[EFFECTIVE_GENOME_SIZE]:-}"; GENOME_SIZE_METHOD="${ss[GENOME_SIZE_METHOD]:-}"
    if [[ "$do_verify" == "true" ]]; then verify_custom_profile_assets "$GENOME"; fi
    GENOME_MODE="local_explicit"; GENOME_BUILD_METHOD="local_build"; CUSTOM_REFERENCE_FROM_SNAPSHOT=true
    ok "Loaded and verified this run's immutable validated genome profile snapshot."
}

validate_fai_matches_chrom_sizes() {
    local fasta="$1" chrom_sizes="$2"
    local fai="${fasta}.fai" left right diff_preview
    [[ -s "$fai" ]] || die "FASTA index is missing or empty: $fai"
    left="$(mktemp)"; right="$(mktemp)"; cut -f1,2 "$fai" > "$left"; cut -f1,2 "$chrom_sizes" > "$right"
    if ! cmp -s "$left" "$right"; then diff_preview="$(diff -u "$left" "$right" 2>/dev/null | head -20 || true)"; rm -f "$left" "$right"; die "FASTA .fai names/lengths do not exactly match chrom.sizes."$'\n'"$diff_preview"; fi
    rm -f "$left" "$right"; ok "FASTA .fai exactly matches chrom.sizes."
}

validate_bowtie2_matches_fasta() {
    local prefix="$1" fasta="$2" summary parsed expected diff_preview
    bowtie2_index_files "$prefix" >/dev/null
    command_exists bowtie2-inspect || die "bowtie2-inspect is required to validate a custom index."
    summary="$(mktemp)"; parsed="$(mktemp)"; expected="$(mktemp)"
    bowtie2-inspect --summary "$prefix" > "$summary" 2>&1 || { local tail_text="$(tail -20 "$summary" || true)"; rm -f "$summary" "$parsed" "$expected"; die "bowtie2-inspect failed for index: $prefix"$'\n'"$tail_text"; }
    awk '$1 ~ /^Sequence-[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $2 "\t" $3}' "$summary" > "$parsed"; cut -f1,2 "${fasta}.fai" > "$expected"
    [[ -s "$parsed" ]] || { rm -f "$summary" "$parsed" "$expected"; die "Could not parse sequence names and lengths from bowtie2-inspect --summary for: $prefix"; }
    if ! cmp -s "$parsed" "$expected"; then diff_preview="$(diff -u "$expected" "$parsed" 2>/dev/null | head -20 || true)"; rm -f "$summary" "$parsed" "$expected"; die "Bowtie2 index sequences/lengths do not match the FASTA .fai."$'\n'"$diff_preview"; fi
    rm -f "$summary" "$parsed" "$expected"; ok "Bowtie2 index sequence names and lengths match the FASTA."
}

validate_profile_bed() {
    local bed="$1" chrom_sizes="$2" label_text="$3" allow_empty="${4:-false}" report
    [[ -e "$bed" ]] || die "$label_text is missing: $bed"
    if [[ ! -s "$bed" ]]; then
        [[ "$allow_empty" == "true" ]] && return 0
        die "$label_text is empty: $bed"
    fi
    report="$(mktemp)"
    if ! awk -F'\t' -v report="$report" 'NR==FNR{if(NF>=2&&$2~/^[0-9]+$/)len[$1]=$2;next} /^#/||NF==0{next} {if(NF<3||!($1 in len)||$2!~/^[0-9]+$/||$3!~/^[0-9]+$/||$2<0||$3<=$2||$3>len[$1]){print "line " FNR ": " $0 > report; bad=1; exit} n++} END{if(bad||n==0)exit 1}' "$chrom_sizes" "$bed"; then
        local why="$(cat "$report" 2>/dev/null || true)"; rm -f "$report"
        die "$label_text failed BED coordinate validation: $why"
    fi
    rm -f "$report"
}

ensure_profile_annotation_packages() {
    local txdb_pkg="$1" orgdb_pkg="$2"
    in_env_clean Rscript --vanilla - "$txdb_pkg" "$orgdb_pkg" <<'RPKG'
args <- commandArgs(trailingOnly=TRUE); pkgs <- args[nzchar(args)]
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat(sprintf("  Installing required annotation package: %s\n", pkg))
    BiocManager::install(pkg, ask=FALSE, update=FALSE)
  }
  if (!requireNamespace(pkg, quietly=TRUE)) stop("Required annotation package could not be installed: ", pkg)
}
RPKG
}

build_profile_annotation_assets() {
    local assembly="$1" txdb_pkg="$2" orgdb_pkg="$3" chrom_sizes="$4" out_dir="$5"
    mkdir -p "$out_dir"
    in_env_clean Rscript --vanilla - "$assembly" "$txdb_pkg" "$orgdb_pkg" "$chrom_sizes" "$out_dir" <<'RPROFILE'
args <- commandArgs(trailingOnly=TRUE)
assembly <- args[1]; txdb_pkg <- args[2]; orgdb_pkg <- args[3]; chrom_file <- args[4]; out <- args[5]
suppressPackageStartupMessages({library(AnnotationDbi); library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(IRanges)})
load_pkg_object <- function(pkg, cls) {
  suppressPackageStartupMessages(library(pkg, character.only=TRUE))
  obj <- tryCatch(getExportedValue(pkg, pkg), error=function(e) NULL)
  if (!is.null(obj) && inherits(obj, cls)) return(obj)
  ns <- asNamespace(pkg)
  for (nm in ls(ns, all.names=TRUE)) {
    z <- tryCatch(get(nm, envir=ns), error=function(e) NULL)
    if (!is.null(z) && inherits(z, cls)) return(z)
  }
  stop(sprintf("Package %s does not expose a %s object.", pkg, cls))
}
txdb <- load_pkg_object(txdb_pkg, "TxDb"); orgdb <- load_pkg_object(orgdb_pkg, "OrgDb")
meta <- tryCatch(AnnotationDbi::metadata(txdb), error=function(e) data.frame())
reported <- NA_character_
if (nrow(meta)) { idx <- which(tolower(as.character(meta$name)) == "genome"); if (length(idx)) reported <- as.character(meta$value[idx[1]]) }
if (is.na(reported) || !nzchar(reported)) stop("The TxDb does not report a UCSC genome assembly in its metadata.")
if (!identical(tolower(reported), tolower(assembly))) stop(sprintf("TxDb assembly mismatch: selected UCSC assembly '%s', but %s reports '%s'.", assembly, txdb_pkg, reported))
chrom <- read.delim(chrom_file, header=FALSE, stringsAsFactors=FALSE, col.names=c("seqname","length"))
if (!nrow(chrom) || anyDuplicated(chrom$seqname)) stop("chrom.sizes is empty or contains duplicate sequence names.")
si <- GenomeInfoDb::seqinfo(txdb); txseq <- as.character(GenomeInfoDb::seqlevels(si)); txlen <- GenomeInfoDb::seqlengths(si)
missing <- setdiff(txseq, chrom$seqname)
if (length(missing)) {
  # A full UCSC TxDb (e.g. TxDb.Hsapiens.UCSC.hg38.knownGene) legitimately
  # models alt/patch/unlocalized contigs (chr1_GL383518v1_alt etc.) as part
  # of the complete reference assembly. A primary-assembly alignment index
  # -- the normal, often-recommended choice for ATAC-seq specifically to
  # avoid multi-mapping against alt contigs -- won't include those
  # sequences. That's not a mismatch between the wrong genome and the
  # wrong TxDb (the assembly-metadata check above already caught that
  # case); it's the completely standard situation of a full gene model
  # paired with a reduced alignment reference. Genes on those contigs can
  # never have aligned reads in this build anyway, so they're dropped from
  # the TxDb rather than treated as a fatal error.
  keep <- intersect(txseq, chrom$seqname)
  if (!length(keep)) stop(sprintf("TxDb and chrom.sizes share NO sequence names at all (checked %d TxDb sequences). This is a genuine mismatch, not just extra alt/patch contigs.", length(txseq)))
  cat(sprintf("  Note: dropping %d TxDb sequence(s) not present in this alignment reference (e.g. alt/patch contigs): %s%s\n",
              length(missing), paste(head(missing,5), collapse=", "), if (length(missing)>5) ", ..." else ""))
  txdb <- GenomeInfoDb::keepSeqlevels(txdb, keep, pruning.mode="coarse")
  si <- GenomeInfoDb::seqinfo(txdb); txseq <- as.character(GenomeInfoDb::seqlevels(si)); txlen <- GenomeInfoDb::seqlengths(si)
}
idx <- match(txseq, chrom$seqname)
if (any(is.na(txlen))) {
  bad <- txseq[is.na(txlen)]
  stop(sprintf("TxDb does not provide chromosome lengths for: %s", paste(head(bad,10), collapse=", ")))
}
if (any(txlen != chrom$length[idx])) { bad <- which(txlen != chrom$length[idx])[1]; stop(sprintf("TxDb chromosome length mismatch for %s: TxDb=%s, assembly=%s", txseq[bad], txlen[bad], chrom$length[idx[bad]])) }
normalize_ids <- function(ids, kt) { ids <- as.character(ids); if (kt %in% c("ENSEMBL","REFSEQ")) ids <- sub("\\.[0-9]+$", "", ids); ids }
tx_ids <- unique(as.character(AnnotationDbi::keys(txdb, keytype="GENEID"))); tx_ids <- tx_ids[!is.na(tx_ids)&nzchar(tx_ids)]
if (!length(tx_ids)) stop("The selected TxDb contains no GENEID keys.")
avail <- AnnotationDbi::keytypes(orgdb)
preferred <- if (grepl("knownGene",txdb_pkg,fixed=TRUE)) c("ENTREZID","ENSEMBL","REFSEQ","SYMBOL","FLYBASE","ZFIN") else if (grepl("ensGene",txdb_pkg,fixed=TRUE)) c("ENSEMBL","FLYBASE","ENTREZID","SYMBOL","REFSEQ","ZFIN") else c("ENTREZID","REFSEQ","ENSEMBL","SYMBOL","FLYBASE","ZFIN")
candidates <- intersect(preferred, avail); if (!length(candidates)) stop("TxDb and OrgDb share no supported identifier keytype.")
rates <- sapply(candidates, function(kt) { ok <- tryCatch(as.character(AnnotationDbi::keys(orgdb,keytype=kt)), error=function(e) character()); max(mean(tx_ids %in% ok), mean(normalize_ids(tx_ids,kt) %in% ok)) })
best <- candidates[which.max(rates)]; rate <- max(rates)
if (!is.finite(rate) || rate < 0.70) stop(sprintf("TxDb/OrgDb compatibility is too low: %.1f%% (minimum 70%%).",100*rate))
cat(sprintf("  Validated profile: %s; TxDb sequences in assembly 100.0%%; TxDb→OrgDb %s mapping %.1f%%.\n", assembly,best,100*rate))
AnnotationDbi::saveDb(txdb, file.path(out,"txdb.sqlite")); AnnotationDbi::saveDb(orgdb, file.path(out,"orgdb.sqlite"))
chrom_order <- setNames(seq_len(nrow(chrom)), chrom$seqname)
validate_core_ranges <- function(gr, label) {
  if (!length(gr)) return(invisible(TRUE))
  seqs <- as.character(seqnames(gr))
  idx <- match(seqs, chrom$seqname)
  if (any(is.na(idx))) stop(label, " contains sequence names absent from chrom.sizes.")
  bad <- which(start(gr) < 1L | end(gr) < start(gr) | end(gr) > chrom$length[idx])
  if (length(bad)) {
    i <- bad[1]
    stop(sprintf("%s contains an invalid/out-of-bounds interval: %s:%s-%s (chromosome length %s).",
                 label, seqs[i], start(gr)[i], end(gr)[i], chrom$length[idx[i]]))
  }
  invisible(TRUE)
}
clip_gr <- function(gr) {
  keep <- as.character(seqnames(gr)) %in% chrom$seqname
  gr <- gr[keep]
  if (!length(gr)) return(gr)
  lim <- chrom$length[match(as.character(seqnames(gr)),chrom$seqname)]
  new_start <- pmax(1L, start(gr))
  new_end <- pmin(end(gr), lim)
  valid <- new_start <= new_end
  gr <- gr[valid]
  if (!length(gr)) return(gr)
  ranges(gr) <- IRanges::IRanges(start=new_start[valid], end=new_end[valid])
  gr
}
name_col <- function(x, candidates, fallback) {
  for (n in candidates) if (n %in% colnames(mcols(x))) {
    v <- as.character(mcols(x)[[n]])
    v[is.na(v)|!nzchar(v)] <- fallback[is.na(v)|!nzchar(v)]
    return(v)
  }
  fallback
}
write_bed <- function(gr, path, fallback_prefix, allow_empty=FALSE) {
  gr <- clip_gr(gr)
  if (!length(gr)) {
    file.create(path)
    if (!allow_empty) stop("Generated empty BED: ", basename(path))
    return(invisible(data.frame()))
  }
  nm <- names(gr)
  if (is.null(nm)) nm <- paste0(fallback_prefix, "_", seq_along(gr))
  bad <- is.na(nm) | !nzchar(nm)
  nm[bad] <- paste0(fallback_prefix, "_", which(bad))
  ord <- order(chrom_order[as.character(seqnames(gr))], start(gr), end(gr), nm)
  gr <- gr[ord]; nm <- nm[ord]
  d <- unique(data.frame(as.character(seqnames(gr)), start(gr)-1L, end(gr), nm, 0,
                         as.character(strand(gr)), stringsAsFactors=FALSE))
  write.table(d,path,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  invisible(d)
}
category_bed <- function(gr, label) {
  gr <- clip_gr(gr)
  if (!length(gr)) return(NULL)
  gr <- GenomicRanges::reduce(gr, ignore.strand=TRUE)
  data.frame(as.character(seqnames(gr)), start(gr)-1L, end(gr), label, 0, ".",
             stringsAsFactors=FALSE)
}
subtract_ranges <- function(gr, masks) {
  gr <- clip_gr(gr)
  if (!length(gr)) return(gr)
  gr <- GenomicRanges::reduce(gr, ignore.strand=TRUE)
  masks <- clip_gr(masks)
  if (length(masks)) gr <- GenomicRanges::setdiff(gr, GenomicRanges::reduce(masks, ignore.strand=TRUE), ignore.strand=TRUE)
  gr
}

genes_gr <- suppressWarnings(GenomicFeatures::genes(txdb))
gene_names <- names(genes_gr); if (is.null(gene_names)) gene_names <- paste0("gene_",seq_along(genes_gr)); names(genes_gr) <- gene_names
tx <- GenomicFeatures::transcripts(txdb, columns=c("tx_name","gene_id"))
tx_names <- name_col(tx,c("tx_name","tx_id"),paste0("tx_",seq_along(tx))); names(tx) <- tx_names
ex <- GenomicFeatures::exons(txdb, columns=c("exon_name"))
ex_names <- name_col(ex,c("exon_name","exon_id"),paste0("exon_",seq_along(ex))); names(ex) <- ex_names
intr <- unlist(GenomicFeatures::intronsByTranscript(txdb,use.names=TRUE),use.names=TRUE)
if (length(intr)) {
  intr_names <- names(intr); if (is.null(intr_names)) intr_names <- paste0("intron_",seq_along(intr)); names(intr) <- intr_names
}
validate_core_ranges(genes_gr, "TxDb genes")
validate_core_ranges(tx, "TxDb transcripts")
validate_core_ranges(ex, "TxDb exons")
validate_core_ranges(intr, "TxDb introns")
tx_strand <- as.character(strand(tx))
if (any(!tx_strand %in% c("+", "-"))) stop("TxDb contains transcripts without a defined + or - strand; TSS/promoter generation is ambiguous.")
prom <- GenomicRanges::promoters(tx,upstream=2000,downstream=200); names(prom) <- tx_names
tss <- tx; plus <- tx_strand == "+"; end(tss)[plus] <- start(tss)[plus]; start(tss)[!plus] <- end(tss)[!plus]; names(tss) <- tx_names

write_bed(genes_gr,file.path(out,"genes.bed"),"gene")
write_bed(ex,file.path(out,"exons.bed"),"exon")
write_bed(intr,file.path(out,"introns.bed"),"intron",allow_empty=TRUE)
write_bed(prom,file.path(out,"promoters.bed"),"promoter")
write_bed(tss,file.path(out,"tss.bed"),"tss")

# PEPATAC's --anno-name input is a six-column BED with a feature class in
# column 4. Make a precedence-based, non-overlapping partition so one base is
# not counted simultaneously as promoter, exon, intron, and whole gene.
prom_part <- GenomicRanges::reduce(clip_gr(prom), ignore.strand=TRUE)
ex_part <- subtract_ranges(ex, prom_part)
intr_part <- subtract_ranges(intr, c(prom_part, ex_part))
gene_part <- subtract_ranges(genes_gr, c(prom_part, ex_part, intr_part))
feature_parts <- Filter(Negate(is.null), list(category_bed(prom_part,"Promoter"),
                                              category_bed(ex_part,"Exon"),
                                              category_bed(intr_part,"Intron"),
                                              category_bed(gene_part,"Gene")))
if (!length(feature_parts)) stop("Generated empty PEPATAC feature annotation BED.")
feature_rows <- unique(do.call(rbind, feature_parts))
feature_rows <- feature_rows[order(chrom_order[feature_rows[[1]]],feature_rows[[2]],feature_rows[[3]],feature_rows[[4]]),]
write.table(feature_rows,file.path(out,"features.bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
conf <- c(
    "PROFILE_VALIDATION_SCHEMA=1",
    sprintf("CUSTOM_UCSC_ASSEMBLY=%s", assembly),
    sprintf("CUSTOM_TXDB_PACKAGE=%s", txdb_pkg),
    sprintf("CUSTOM_TXDB_PACKAGE_VERSION=%s", as.character(packageVersion(txdb_pkg))),
    sprintf("CUSTOM_ORGDB_PACKAGE=%s", orgdb_pkg),
    sprintf("CUSTOM_ORGDB_PACKAGE_VERSION=%s", as.character(packageVersion(orgdb_pkg))),
    sprintf("CUSTOM_TXDB_REPORTED_ASSEMBLY=%s", reported),
    sprintf("CUSTOM_TXDB_ORGDB_KEYTYPE=%s", best),
    sprintf("CUSTOM_TXDB_ORGDB_MAPPING_RATE=%.8f", rate),
    sprintf("R_VERSION=%s", R.version.string),
    sprintf("BIOCONDUCTOR_VERSION=%s", if (requireNamespace("BiocManager", quietly=TRUE)) as.character(BiocManager::version()) else "unknown"),
    sprintf("ANNOTATIONDBI_VERSION=%s", as.character(packageVersion("AnnotationDbi"))),
    sprintf("GENOMICFEATURES_VERSION=%s", as.character(packageVersion("GenomicFeatures"))),
    sprintf("GENOMEINFODB_VERSION=%s", as.character(packageVersion("GenomeInfoDb"))),
    "PROMOTER_UPSTREAM_BP=2000",
    "PROMOTER_DOWNSTREAM_BP=200"
)
writeLines(conf, file.path(out, "profile_build.conf"))
RPROFILE
}

validate_frozen_annotation_profile() {
    local assembly="$1" txdb_sqlite="$2" orgdb_sqlite="$3" chrom_sizes="$4"
    in_env_clean Rscript --vanilla - "$assembly" "$txdb_sqlite" "$orgdb_sqlite" "$chrom_sizes" <<'RVALID'
args <- commandArgs(trailingOnly=TRUE); assembly <- args[1]
suppressPackageStartupMessages({library(AnnotationDbi);library(GenomicFeatures);library(GenomeInfoDb)})
txdb <- AnnotationDbi::loadDb(args[2]); orgdb <- AnnotationDbi::loadDb(args[3])
if (!inherits(txdb,"TxDb")) stop("Frozen TxDb SQLite did not load as a TxDb.")
if (!inherits(orgdb,"OrgDb")) stop("Frozen OrgDb SQLite did not load as an OrgDb.")
meta <- AnnotationDbi::metadata(txdb); idx <- which(tolower(as.character(meta$name))=="genome"); reported <- if(length(idx)) as.character(meta$value[idx[1]]) else NA_character_
if (is.na(reported)||tolower(reported)!=tolower(assembly)) stop("Frozen TxDb assembly metadata does not match profile assembly.")
chrom <- read.delim(args[4],header=FALSE,stringsAsFactors=FALSE,col.names=c("seqname","length")); si <- GenomeInfoDb::seqinfo(txdb); sn <- as.character(GenomeInfoDb::seqlevels(si)); sl <- GenomeInfoDb::seqlengths(si)
# saveDb()/loadDb() persist the TxDb's underlying data, not the
# session-level "active seqlevels" restriction build_profile_annotation_
# assets() applied via keepSeqlevels() before deriving genes/exons/BEDs --
# that restriction is a view filter, not part of the SQLite file, so a
# freshly reloaded TxDb always reports its FULL original sequence set
# again (alt/patch contigs included). The derived BED/GRanges output is
# unaffected by this (it was already correctly filtered at the point it
# was written), so this only needs to check that sequences present in
# BOTH the TxDb and chrom.sizes actually agree in length, not that the
# TxDb is a strict subset of chrom.sizes.
overlap <- intersect(sn,chrom$seqname)
if (!length(overlap)) stop("Frozen TxDb and chrom.sizes share NO sequence names at all. This is a genuine mismatch, not just extra alt/patch contigs.")
i <- match(overlap,chrom$seqname); ti <- match(overlap,sn)
if(any(is.na(sl[ti]))) stop("Frozen TxDb is missing one or more chromosome lengths for sequences shared with chrom.sizes.")
if(any(sl[ti]!=chrom$length[i])) stop("Frozen TxDb chromosome lengths do not match chrom.sizes for one or more shared sequences.")
RVALID
}

# ─────────────────────────────────────────────────────────────
# Built-in genome annotation QC assets
#
# --TSS-name/--anno-name are PEPATAC's own upstream QC inputs (TSS
# enrichment, feature-distribution plots) -- previously wired up only for
# user-defined custom genome profiles, leaving the five built-in genomes
# without this QC entirely. This closes that gap by reusing the exact same
# validated builder (build_profile_annotation_assets) that custom profiles
# already use, fed with each built-in genome's hardcoded TxDb/OrgDb
# package names -- the same ones diff_analysis.sh has always used.
#
# Storage is fingerprinted, not a fixed per-genome path: install.sh does
# not pin Bioconductor package versions, so a future BiocManager::install()
# could silently change what TxDb.Hsapiens.UCSC.hg38.knownGene resolves to
# under the same name. A fixed ~/pepatac_genomes/hg38/annotations/ path
# would then either silently keep serving BEDs built from the old package
# version, or overwrite them and change what later runs see -- exactly the
# stale-cache problem already fixed once for custom profiles. The
# fingerprint captures package name+version, the chrom.sizes this was
# validated against, the annotation-builder schema, and the promoter
# window, so any of those changing produces a new, non-colliding
# directory instead.
# ─────────────────────────────────────────────────────────────

# Promoter window used by build_profile_annotation_assets()'s R body
# (GenomicRanges::promoters(tx, upstream=2000, downstream=200)). Not
# passed as a parameter -- build_profile_annotation_assets is left
# unmodified and reused as-is for both custom and built-in genomes -- so
# these exist only to make that fixed window part of the fingerprint. If
# the hardcoded values inside build_profile_annotation_assets() ever
# change, update these to match.
ANNOTATION_BUILDER_SCHEMA="1"
ANNOTATION_PROMOTER_UPSTREAM_BP=2000
ANNOTATION_PROMOTER_DOWNSTREAM_BP=200

declare -A BUILTIN_TXDB_PKG=(
    [hg38]="TxDb.Hsapiens.UCSC.hg38.knownGene"
    [mm10]="TxDb.Mmusculus.UCSC.mm10.knownGene"
    [rn7]="TxDb.Rnorvegicus.UCSC.rn7.refGene"
    [dm6]="TxDb.Dmelanogaster.UCSC.dm6.ensGene"
    [danRer11]="TxDb.Drerio.UCSC.danRer11.refGene"
)
declare -A BUILTIN_ORGDB_PKG=(
    [hg38]="org.Hs.eg.db"
    [mm10]="org.Mm.eg.db"
    [rn7]="org.Rn.eg.db"
    [dm6]="org.Dm.eg.db"
    [danRer11]="org.Dr.eg.db"
)

# Generic, genome-agnostic run-level annotation state. Resolved exactly
# once per run (see resolve_run_annotation_assets below), then only ever
# read -- by the pepatac.py command builder and by write_reference_snapshot.
# Deliberately NOT folded into the CUSTOM_* associative arrays: those are
# woven through the custom-profile registry, immutable snapshot, resume
# state, and publish/hash-verification logic, and pulling built-in genomes
# into that same family would mean either teaching all of that machinery
# to handle a "this isn't really custom" case, or duplicating it. A flat
# set of scalars resolved once from either source touches none of the
# existing custom-profile code.
RUN_ANNOTATION_STATUS=""                 # "validated" | "disabled_by_override"
RUN_ANNOTATION_SOURCE=""                 # "custom_profile" | "builtin_txdb"
RUN_ANNOTATION_FINGERPRINT=""
RUN_ANNOTATION_QC_FAILURE_REASON=""
RUN_ANNOTATION_TXDB_PACKAGE=""
RUN_ANNOTATION_TXDB_VERSION=""
RUN_ANNOTATION_TXDB_SQLITE=""
RUN_ANNOTATION_TXDB_SQLITE_SHA256=""
RUN_ANNOTATION_ORGDB_PACKAGE=""
RUN_ANNOTATION_ORGDB_VERSION=""
RUN_ANNOTATION_ORGDB_SQLITE=""
RUN_ANNOTATION_ORGDB_SQLITE_SHA256=""
RUN_ANNOTATION_TSS_BED=""
RUN_ANNOTATION_TSS_BED_SHA256=""
RUN_ANNOTATION_FEATURE_BED=""
RUN_ANNOTATION_FEATURE_BED_SHA256=""
# Lets diff_analysis.sh independently repeat the TxDb-vs-assembly-length
# check for a built-in genome using a frozen snapshot, the same way it
# already can for a custom profile via CUSTOM_CHROM_SIZES -- previously
# only recorded for custom, so a built-in snapshot user had to trust that
# run.sh's build-time check (validate_frozen_annotation_profile) was
# sufficient with no way to re-verify it downstream.
RUN_ANNOTATION_CHROM_SIZES=""
RUN_ANNOTATION_CHROM_SIZES_SHA256=""

# compute_annotation_fingerprint ASSEMBLY CHROM_SHA TXDB_PKG TXDB_VERSION ORGDB_PKG ORGDB_VERSION
# Hashes a canonical, fixed-field-order input file -- deliberately NOT an
# iterated associative array, whose key order is not guaranteed stable
# across bash builds/versions and would make the "same" logical
# fingerprint hash differently on different machines.
compute_annotation_fingerprint() {
    local assembly="$1" chrom_sha="$2" txdb_pkg="$3" txdb_version="$4" orgdb_pkg="$5" orgdb_version="$6"
    local fp_input digest
    fp_input="$(mktemp)"
    {
        printf 'ASSEMBLY=%s\n' "$assembly"
        printf 'CHROM_SIZES_SHA256=%s\n' "$chrom_sha"
        printf 'TXDB_PACKAGE=%s\n' "$txdb_pkg"
        printf 'TXDB_VERSION=%s\n' "$txdb_version"
        printf 'ORGDB_PACKAGE=%s\n' "$orgdb_pkg"
        printf 'ORGDB_VERSION=%s\n' "$orgdb_version"
        printf 'ANNOTATION_BUILDER_SCHEMA=%s\n' "$ANNOTATION_BUILDER_SCHEMA"
        printf 'PROMOTER_UPSTREAM_BP=%s\n' "$ANNOTATION_PROMOTER_UPSTREAM_BP"
        printf 'PROMOTER_DOWNSTREAM_BP=%s\n' "$ANNOTATION_PROMOTER_DOWNSTREAM_BP"
    } > "$fp_input"
    digest="$(_soft_sha256_of "$fp_input")"
    rm -f "$fp_input"
    [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || return 1
    printf '%s\n' "$digest"
}

# _soft_sha256_of / _soft_verify_file_sha256 / _soft_validate_bed
# Non-dying counterparts of sha256_of/verify_file_sha256/validate_profile_bed,
# used ONLY by the built-in annotation path below. A failure there must be
# catchable so PEPATAC_ALLOW_MISSING_ANNOTATION_QC can offer a deliberate,
# recorded degradation instead of stopping the run outright. Custom genome
# profiles never call these -- they keep using the existing die-on-failure
# versions completely unchanged, which is exactly what makes "a custom
# profile's own validation always hard-stops, no override" true without an
# explicit custom-vs-built-in branch inside the validation logic itself.
_soft_sha256_of() {
    local f="$1"
    [[ -f "$f" ]] || return 1
    command_exists sha256sum || return 1
    local digest
    digest="$(sha256sum "$f" 2>/dev/null | awk '{print $1}')"
    [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || return 1
    printf '%s\n' "$digest"
}
_soft_verify_file_sha256() {
    local f="$1" expected="$2" actual
    [[ -n "$expected" && -f "$f" ]] || return 1
    actual="$(_soft_sha256_of "$f")" || return 1
    [[ "$actual" == "$expected" ]]
}
# Mirrors validate_profile_bed's coordinate-validation logic exactly (same
# awk pattern) but returns nonzero instead of dying. If that function's
# validation rules ever change, update this to match.
_soft_validate_bed() {
    local bed="$1" chrom_sizes="$2" allow_empty="${3:-false}"
    [[ -e "$bed" ]] || return 1
    if [[ ! -s "$bed" ]]; then
        [[ "$allow_empty" == "true" ]]
        return $?
    fi
    awk -F'\t' 'NR==FNR{if(NF>=2&&$2~/^[0-9]+$/)len[$1]=$2;next} /^#/||NF==0{next} {if(NF<3||!($1 in len)||$2!~/^[0-9]+$/||$3!~/^[0-9]+$/||$2<0||$3<=$2||$3>len[$1]){bad=1; exit} n++} END{if(bad||n==0)exit 1}' "$chrom_sizes" "$bed"
}

# _fail_builtin_annotation GENOME REASON
# Central failure point for ensure_builtin_annotation_assets(). Strict by
# default -- dies, matching every other reference-asset failure in this
# script (a failed bowtie2-build isn't allowed to gracefully degrade
# either, and this runs at the same pre-processing stage: before any FASTQ
# has been touched, so failing here costs nothing and the run is fully
# resumable). PEPATAC_ALLOW_MISSING_ANNOTATION_QC=1 is a deliberate,
# narrow escape hatch that disables PEPATAC's own TSS/feature QC for this
# run only -- recorded explicitly in the manifest and snapshot, never
# silent. Never consulted for custom genome profiles: those keep using the
# existing die-on-failure validation functions unchanged, so a custom
# profile's own validation always hard-stops regardless of this variable.
_fail_builtin_annotation() {
    local genome="$1" reason="$2"
    if [[ "${PEPATAC_ALLOW_MISSING_ANNOTATION_QC:-0}" == "1" ]]; then
        warn "Annotation QC assets could not be built or validated for '$genome': $reason"
        warn "PEPATAC_ALLOW_MISSING_ANNOTATION_QC=1 is set -- continuing WITHOUT --TSS-name/--anno-name for this run."
        warn "This run's PEPATAC QC report will be missing TSS enrichment and feature-distribution plots."
        RUN_ANNOTATION_STATUS="disabled_by_override"
        RUN_ANNOTATION_SOURCE="builtin_txdb"
        RUN_ANNOTATION_QC_FAILURE_REASON="$reason"
        return 0
    fi
    die "Could not build or validate annotation QC assets for $genome: $reason"$'\n'"       Set PEPATAC_ALLOW_MISSING_ANNOTATION_QC=1 to run without --TSS-name/--anno-name instead (recorded in the manifest)."
}

# _bind_builtin_annotation_assets DIR
# Populates RUN_ANNOTATION_* from a fingerprinted annotation directory
# (freshly published or reused), re-hashing and re-validating every file
# rather than trusting that a directory with the right name still has the
# right content. Reads the _BA_* scoped globals set by
# ensure_builtin_annotation_assets() just before calling this. Returns
# nonzero on any failure instead of dying, so the caller can route through
# _fail_builtin_annotation and respect the override.
_bind_builtin_annotation_assets() {
    local dir="$1"
    local manifest="$dir/asset_manifest.sha256"
    local txdb_sqlite="$dir/txdb.sqlite" orgdb_sqlite="$dir/orgdb.sqlite"
    local tss_bed="$dir/tss.bed" feature_bed="$dir/features.bed"

    # Verify every published file against the permanent manifest written at
    # publish time -- NOT freshly-computed hashes of whatever is currently
    # on disk. Re-hashing without a recorded expected value would bless any
    # directory with the right filenames regardless of content (reported
    # and reproduced: plain text files named txdb.sqlite/orgdb.sqlite would
    # otherwise pass).
    [[ -f "$manifest" ]] || return 1
    ( cd "$dir" && sha256sum --check --strict --status asset_manifest.sha256 ) || return 1

    # Re-run the same TxDb/OrgDb content validation the build path used --
    # confirms these actually load as TxDb/OrgDb objects with the right
    # assembly metadata and matching chrom.sizes, not just that their bytes
    # match what the manifest recorded (a manifest only proves the bytes
    # haven't changed since publish, not that publish itself was sound).
    validate_frozen_annotation_profile "$_BA_GENOME" "$txdb_sqlite" "$orgdb_sqlite" "$_BA_CHROM_SIZES" || return 1

    # Validate every generated BED, not only the two PEPATAC consumes.
    local _bf _bl _be
    for _bf in tss.bed:false features.bed:false genes.bed:false exons.bed:false introns.bed:true promoters.bed:false; do
        _bl="${_bf%%:*}"; _be="${_bf##*:}"
        _soft_validate_bed "$dir/$_bl" "$_BA_CHROM_SIZES" "$_be" || return 1
    done

    local txdb_sha orgdb_sha tss_sha feature_sha
    txdb_sha="$(_soft_sha256_of "$txdb_sqlite")" || return 1
    orgdb_sha="$(_soft_sha256_of "$orgdb_sqlite")" || return 1
    tss_sha="$(_soft_sha256_of "$tss_bed")" || return 1
    feature_sha="$(_soft_sha256_of "$feature_bed")" || return 1

    RUN_ANNOTATION_STATUS="validated"
    RUN_ANNOTATION_SOURCE="builtin_txdb"
    RUN_ANNOTATION_FINGERPRINT="$_BA_FINGERPRINT"
    RUN_ANNOTATION_TXDB_PACKAGE="$_BA_TXDB_PKG"; RUN_ANNOTATION_TXDB_VERSION="$_BA_TXDB_VERSION"
    RUN_ANNOTATION_TXDB_SQLITE="$txdb_sqlite"; RUN_ANNOTATION_TXDB_SQLITE_SHA256="$txdb_sha"
    RUN_ANNOTATION_ORGDB_PACKAGE="$_BA_ORGDB_PKG"; RUN_ANNOTATION_ORGDB_VERSION="$_BA_ORGDB_VERSION"
    RUN_ANNOTATION_ORGDB_SQLITE="$orgdb_sqlite"; RUN_ANNOTATION_ORGDB_SQLITE_SHA256="$orgdb_sha"
    RUN_ANNOTATION_TSS_BED="$tss_bed"; RUN_ANNOTATION_TSS_BED_SHA256="$tss_sha"
    RUN_ANNOTATION_FEATURE_BED="$feature_bed"; RUN_ANNOTATION_FEATURE_BED_SHA256="$feature_sha"
    RUN_ANNOTATION_CHROM_SIZES="$_BA_CHROM_SIZES"; RUN_ANNOTATION_CHROM_SIZES_SHA256="$_BA_CHROM_SIZES_SHA256"
    return 0
}

# ensure_builtin_annotation_assets GENOME CHROM_SIZES_FILE
# Builds (or reuses) the TSS/feature-annotation BEDs and frozen TxDb/OrgDb
# SQLite databases a built-in genome needs for --TSS-name/--anno-name,
# using the same validated builder custom profiles already use. On
# success, populates RUN_ANNOTATION_*. On failure, either dies or -- only
# with the explicit override -- records a degraded state for the caller to
# act on. Never called for custom genome profiles; those resolve
# RUN_ANNOTATION_* directly from their own CUSTOM_* arrays instead (see
# resolve_run_annotation_assets).
ensure_builtin_annotation_assets() {
    local genome="$1" chrom_sizes="$2"
    local txdb_pkg="${BUILTIN_TXDB_PKG[$genome]:-}" orgdb_pkg="${BUILTIN_ORGDB_PKG[$genome]:-}"

    if [[ -z "$txdb_pkg" || -z "$orgdb_pkg" ]]; then
        _fail_builtin_annotation "$genome" "No TxDb/OrgDb package mapping is defined in this runner for '$genome'."
        return
    fi

    header "Annotation QC Assets"
    echo -e "  Building annotation QC assets for $genome"
    echo -e "    ${DIM}TxDb: $txdb_pkg${RESET}"
    echo -e "    ${DIM}OrgDb: $orgdb_pkg${RESET}"
    echo -e "    ${DIM}This is a one-time build for this annotation version.${RESET}"

    local lock_dir="$HOME/pepatac_genomes/${genome}"
    mkdir -p "$lock_dir"
    if ! command_exists flock; then
        _fail_builtin_annotation "$genome" "flock is required but was not found. Run install.sh first."
        return
    fi
    eval "exec 201>\"$lock_dir/.annotation.lock\""
    if ! flock -w 7200 201; then
        exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Timed out after 2 hours waiting for the annotation build lock."
        return
    fi

    if ! ensure_profile_annotation_packages "$txdb_pkg" "$orgdb_pkg"; then
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Could not install/verify required annotation packages ($txdb_pkg, $orgdb_pkg)."
        return
    fi

    local versions
    versions="$(in_env_clean Rscript --vanilla - "$txdb_pkg" "$orgdb_pkg" <<'RVER'
a <- commandArgs(trailingOnly=TRUE)
cat(as.character(packageVersion(a[1])), as.character(packageVersion(a[2])), "\n")
RVER
    2>/dev/null || true)"
    local txdb_version orgdb_version
    txdb_version="$(awk '{print $1}' <<<"$versions")"
    orgdb_version="$(awk '{print $2}' <<<"$versions")"
    if [[ -z "$txdb_version" || -z "$orgdb_version" ]]; then
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Could not resolve installed package versions for $txdb_pkg / $orgdb_pkg."
        return
    fi

    local chrom_sha
    chrom_sha="$(_soft_sha256_of "$chrom_sizes")"
    if [[ -z "$chrom_sha" ]]; then
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Could not hash chrom.sizes: $chrom_sizes"
        return
    fi

    local fingerprint
    fingerprint="$(compute_annotation_fingerprint "$genome" "$chrom_sha" "$txdb_pkg" "$txdb_version" "$orgdb_pkg" "$orgdb_version")"
    if [[ -z "$fingerprint" ]]; then
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Could not compute an annotation fingerprint."
        return
    fi

    _BA_GENOME="$genome"
    _BA_TXDB_PKG="$txdb_pkg"; _BA_TXDB_VERSION="$txdb_version"
    _BA_ORGDB_PKG="$orgdb_pkg"; _BA_ORGDB_VERSION="$orgdb_version"
    _BA_CHROM_SIZES="$chrom_sizes"; _BA_CHROM_SIZES_SHA256="$chrom_sha"; _BA_FINGERPRINT="$fingerprint"

    local target_dir="$HOME/pepatac_genomes/${genome}/annotations/${fingerprint}"

    if [[ -d "$target_dir" ]] && _bind_builtin_annotation_assets "$target_dir"; then
        ok "Annotation QC assets already validated for this exact TxDb/OrgDb version -- reusing: $target_dir"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        return 0
    fi
    if [[ -d "$target_dir" ]]; then
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Existing annotation directory failed hash/validation: $target_dir"
        return
    fi

    local staging_root="$HOME/pepatac_genomes/.staging"
    mkdir -p "$staging_root"
    local stage="$staging_root/annotation_${genome}_${RUN_ID}_$$"
    rm -rf "$stage"
    # Tracked so cleanup_all_staging (registered on the script's EXIT trap)
    # removes this directory if the script is interrupted mid-build. Safe
    # to set once and leave set: every exit path below either rm -rf's
    # $stage directly (the tracked path is then simply gone, so the
    # trap-time -d check is already false and does nothing) or mv's it to
    # $target_dir (a different path, so the trap-time check on the old
    # $stage path is also already false) -- no path re-clears this.
    ACTIVE_BUILTIN_ANNOTATION_STAGE="$stage"

    if ! build_profile_annotation_assets "$genome" "$txdb_pkg" "$orgdb_pkg" "$chrom_sizes" "$stage"; then
        rm -rf "$stage"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "build_profile_annotation_assets failed for $txdb_pkg / $orgdb_pkg."
        return
    fi
    if ! validate_frozen_annotation_profile "$genome" "$stage/txdb.sqlite" "$stage/orgdb.sqlite" "$chrom_sizes"; then
        rm -rf "$stage"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Frozen TxDb/OrgDb failed validation for $genome."
        return
    fi
    local bad=false _bf _bl _be
    for _bf in tss.bed:false features.bed:false genes.bed:false exons.bed:false introns.bed:true promoters.bed:false; do
        _bl="${_bf%%:*}"; _be="${_bf##*:}"
        _soft_validate_bed "$stage/$_bl" "$chrom_sizes" "$_be" || { bad=true; break; }
    done
    if $bad; then
        rm -rf "$stage"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "One or more generated annotation BEDs failed validation."
        return
    fi

    # Permanent integrity manifest, written into staging before publish --
    # this is what every future reuse verifies bytes against (see
    # _bind_builtin_annotation_assets), rather than trusting that a
    # directory with the right name still holds the right content.
    if ! ( cd "$stage" && sha256sum txdb.sqlite orgdb.sqlite tss.bed features.bed genes.bed exons.bed introns.bed promoters.bed profile_build.conf > asset_manifest.sha256 ); then
        rm -rf "$stage"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        _fail_builtin_annotation "$genome" "Failed to write the annotation asset manifest."
        return
    fi

    mkdir -p "$(dirname "$target_dir")"
    if [[ -e "$target_dir" ]]; then
        # Published by a concurrent process while we were building in
        # staging -- compare byte-for-byte rather than trusting the name,
        # mirroring how the custom-profile publish path handles the same
        # race for its own fingerprinted annotation directories.
        local f2 mismatch=false
        for f2 in txdb.sqlite orgdb.sqlite tss.bed features.bed genes.bed exons.bed introns.bed promoters.bed profile_build.conf asset_manifest.sha256; do
            cmp -s "$stage/$f2" "$target_dir/$f2" || { mismatch=true; break; }
        done
        rm -rf "$stage"
        if $mismatch; then
            flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
            _fail_builtin_annotation "$genome" "A concurrently-published annotation directory has different bytes for the same fingerprint: $target_dir"
            return
        fi
    else
        mv "$stage" "$target_dir"
    fi

    if _bind_builtin_annotation_assets "$target_dir"; then
        ok "Published built-in annotation QC assets: $target_dir"
        flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
        return 0
    fi
    flock -u 201 2>/dev/null || true; exec 201>&- 2>/dev/null || true
    _fail_builtin_annotation "$genome" "Published annotation directory failed hash verification: $target_dir"
}

# resolve_run_annotation_assets GENOME CHROM_SIZES_FILE
# Populates the generic RUN_ANNOTATION_* scalars exactly once per run, from
# whichever source applies -- the existing CUSTOM_* arrays for a custom
# profile (already validated by verify_custom_profile_assets earlier in
# ensure_genome_assets), or a freshly built/reused fingerprinted directory
# for a built-in genome. Called once from ensure_genome_assets(); the
# pepatac.py command builder and write_reference_snapshot only ever read
# the result afterward, never re-resolve it. NOTE: on a resume where
# annotation was already resolved and recorded in resume state, the call
# site (near RUN_WIDE_FINGERPRINT below) skips calling this entirely in
# favor of verify_resumed_annotation_assets() -- see that function.
resolve_run_annotation_assets() {
    local genome="$1" chrom_sizes="$2"
    if $GENOME_IS_CUSTOM; then
        RUN_ANNOTATION_STATUS="validated"
        RUN_ANNOTATION_SOURCE="custom_profile"
        RUN_ANNOTATION_FINGERPRINT="${CUSTOM_ANNOTATION_FINGERPRINT[$genome]:-}"
        RUN_ANNOTATION_TXDB_PACKAGE="${CUSTOM_TXDB_PKG[$genome]:-}"; RUN_ANNOTATION_TXDB_VERSION="${CUSTOM_TXDB_VERSION[$genome]:-}"
        RUN_ANNOTATION_TXDB_SQLITE="${CUSTOM_TXDB_SQLITE[$genome]:-}"; RUN_ANNOTATION_TXDB_SQLITE_SHA256="${CUSTOM_TXDB_SQLITE_SHA256[$genome]:-}"
        RUN_ANNOTATION_ORGDB_PACKAGE="${CUSTOM_ORGDB_PKG[$genome]:-}"; RUN_ANNOTATION_ORGDB_VERSION="${CUSTOM_ORGDB_VERSION[$genome]:-}"
        RUN_ANNOTATION_ORGDB_SQLITE="${CUSTOM_ORGDB_SQLITE[$genome]:-}"; RUN_ANNOTATION_ORGDB_SQLITE_SHA256="${CUSTOM_ORGDB_SQLITE_SHA256[$genome]:-}"
        RUN_ANNOTATION_TSS_BED="${CUSTOM_TSS_BED[$genome]:-}"; RUN_ANNOTATION_TSS_BED_SHA256="${CUSTOM_TSS_BED_SHA256[$genome]:-}"
        RUN_ANNOTATION_FEATURE_BED="${CUSTOM_FEATURE_BED[$genome]:-}"; RUN_ANNOTATION_FEATURE_BED_SHA256="${CUSTOM_FEATURE_BED_SHA256[$genome]:-}"
        # LOCAL_CHROM_SIZES/_SHA256 are the resolved chrom.sizes for this
        # custom genome (set during its build/verification earlier in
        # ensure_genome_assets) -- same value diff_analysis.sh already
        # gets via CUSTOM_CHROM_SIZES_PATH, just also exposed generically.
        RUN_ANNOTATION_CHROM_SIZES="${LOCAL_CHROM_SIZES:-}"; RUN_ANNOTATION_CHROM_SIZES_SHA256="${LOCAL_CHROM_SIZES_SHA256:-}"
        [[ -n "$RUN_ANNOTATION_TSS_BED" && -n "$RUN_ANNOTATION_FEATURE_BED" ]] \
            || die "Internal error: custom profile '$genome' has no TSS/feature BED recorded after asset verification."
    else
        ensure_builtin_annotation_assets "$genome" "$chrom_sizes"
    fi
}

# verify_resumed_annotation_assets
# Used only when resuming AND resume state already recorded annotation
# info (RUN_ANNOTATION_STATUS non-empty after load_resume_state). Trusts
# what was already resolved and used for this run's earlier samples,
# re-verifying the recorded files still match their recorded hashes rather
# than re-resolving from whatever TxDb/OrgDb happens to be installed now.
# Deliberately does NOT recompute a fresh fingerprint or touch the live
# registry/installed-package state -- doing so is exactly how a resume
# could silently swap a later-installed package version out from under
# already-completed samples while their reference_snapshot gets
# overwritten to match, and their QC files stay built from the original.
# Uses the existing die-on-mismatch verify_file_sha256: a hash mismatch
# here means the file changed since this run started, which is a real
# integrity problem to stop on, not something to silently re-resolve past.
verify_resumed_annotation_assets() {
    case "$RUN_ANNOTATION_STATUS" in
        validated)
            verify_file_sha256 "$RUN_ANNOTATION_TXDB_SQLITE" "$RUN_ANNOTATION_TXDB_SQLITE_SHA256" "resumed run TxDb SQLite"
            verify_file_sha256 "$RUN_ANNOTATION_ORGDB_SQLITE" "$RUN_ANNOTATION_ORGDB_SQLITE_SHA256" "resumed run OrgDb SQLite"
            verify_file_sha256 "$RUN_ANNOTATION_TSS_BED" "$RUN_ANNOTATION_TSS_BED_SHA256" "resumed run TSS BED"
            verify_file_sha256 "$RUN_ANNOTATION_FEATURE_BED" "$RUN_ANNOTATION_FEATURE_BED_SHA256" "resumed run feature BED"
            # Older resume-state files predate this field -- blank means
            # "not recorded," not "corrupted," so it's skipped rather than
            # failed.
            [[ -n "$RUN_ANNOTATION_CHROM_SIZES" ]] && \
                verify_file_sha256 "$RUN_ANNOTATION_CHROM_SIZES" "$RUN_ANNOTATION_CHROM_SIZES_SHA256" "resumed run annotation chrom.sizes"
            ok "Resumed run's annotation assets verified unchanged since this run started ($RUN_ANNOTATION_SOURCE, fingerprint ${RUN_ANNOTATION_FINGERPRINT:-<none>})."
            ;;
        disabled_by_override)
            ok "Resumed run previously recorded annotation QC as disabled by override -- continuing without --TSS-name/--anno-name, unchanged from the original run."
            ;;
        *)
            die "Resume state has an unrecognized RUN_ANNOTATION_STATUS ('${RUN_ANNOTATION_STATUS}'). Refusing to guess -- start a new run, or check for a corrupted resume state file."
            ;;
    esac
}


verify_custom_profile_assets() {
    local genome="$1"
    verify_file_sha256 "$LOCAL_CUSTOM_FASTA" "$LOCAL_CUSTOM_FASTA_SHA256" "profile FASTA"
    verify_file_sha256 "$LOCAL_CHROM_SIZES" "$LOCAL_CHROM_SIZES_SHA256" "profile chrom.sizes"
    verify_bowtie2_index_sha256 "$LOCAL_BT2_INDEX" "$LOCAL_BT2_INDEX_SHA256"
    validate_fai_matches_chrom_sizes "$LOCAL_CUSTOM_FASTA" "$LOCAL_CHROM_SIZES"
    validate_bowtie2_matches_fasta "$LOCAL_BT2_INDEX" "$LOCAL_CUSTOM_FASTA"
    local key label
    for key in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
        local path_var="CUSTOM_${key}[$genome]" sha_var="CUSTOM_${key}_SHA256[$genome]" path digest
        path="${!path_var:-}"
        digest="${!sha_var:-}"
        label="custom ${key,,}"; verify_file_sha256 "$path" "$digest" "$label"
    done
    for key in TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED; do local pv="CUSTOM_${key}[$genome]" allow_empty=false; [[ "$key" == "INTRONS_BED" ]] && allow_empty=true; validate_profile_bed "${!pv}" "$LOCAL_CHROM_SIZES" "${key,,}" "$allow_empty"; done
    validate_frozen_annotation_profile "${CUSTOM_UCSC_ASSEMBLY[$genome]:-$genome}" "${CUSTOM_TXDB_SQLITE[$genome]}" "${CUSTOM_ORGDB_SQLITE[$genome]}" "$LOCAL_CHROM_SIZES"
    ok "Validated genome profile passed sequence, annotation, and hash checks."
}

verify_registered_custom_reference() {
    local genome="$1" assembly expected_fasta expected_chrom
    [[ "$CUSTOM_REG_SCHEMA" -ge 3 ]] || die "This registration predates validated TxDb/OrgDb profiles. Re-register '$genome' with runner v1.25 or later."
    assembly="${CUSTOM_REG_UCSC_ASSEMBLY:-$genome}"
    expected_fasta="$(ucsc_profile_fasta_url "$assembly")"
    expected_chrom="$(ucsc_profile_chrom_sizes_url "$assembly")"
    [[ "$CUSTOM_REG_FASTA_SOURCE" == "$expected_fasta" && "$CUSTOM_REG_CHROM_SOURCE" == "$expected_chrom" ]] \
        || die "Registered profile '$genome' is not tied to the official UCSC sequence sources for '$assembly'. Re-register it with runner v1.25 or later."
    [[ "$CUSTOM_REG_TXDB_PKG" == *".UCSC.${assembly}."* ]] \
        || die "Registered TxDb package '$CUSTOM_REG_TXDB_PKG' does not identify UCSC assembly '$assembly'. Re-register this profile."
    verify_or_restore_registered_file "$CUSTOM_REG_FASTA_CACHE" "$CUSTOM_REG_FASTA_SOURCE" "$CUSTOM_REG_FASTA_SHA256" "cached FASTA"
    verify_or_restore_registered_file "$CUSTOM_REG_CHROM_SIZES" "$CUSTOM_REG_CHROM_SOURCE" "$CUSTOM_REG_CHROM_SIZES_SHA256" "cached chrom.sizes"
    verify_bowtie2_index_sha256 "$CUSTOM_REG_BT2_INDEX" "$CUSTOM_REG_BT2_INDEX_SHA256"
    LOCAL_CUSTOM_FASTA="$CUSTOM_REG_FASTA_CACHE"; LOCAL_CUSTOM_FASTA_SHA256="$CUSTOM_REG_FASTA_SHA256"
    LOCAL_CHROM_SIZES="$CUSTOM_REG_CHROM_SIZES"; LOCAL_CHROM_SIZES_SHA256="$CUSTOM_REG_CHROM_SIZES_SHA256"
    LOCAL_BT2_INDEX="$CUSTOM_REG_BT2_INDEX"; LOCAL_BT2_INDEX_SHA256="$CUSTOM_REG_BT2_INDEX_SHA256"
    CUSTOM_FASTA_URL["$genome"]="$CUSTOM_REG_FASTA_SOURCE"; CUSTOM_CHROMSIZES_URL["$genome"]="$CUSTOM_REG_CHROM_SOURCE"
    CUSTOM_UCSC_ASSEMBLY["$genome"]="$CUSTOM_REG_UCSC_ASSEMBLY"; CUSTOM_TXDB_PKG["$genome"]="$CUSTOM_REG_TXDB_PKG"; CUSTOM_TXDB_VERSION["$genome"]="$CUSTOM_REG_TXDB_VERSION"
    CUSTOM_ORGDB_PKG["$genome"]="$CUSTOM_REG_ORGDB_PKG"; CUSTOM_ORGDB_VERSION["$genome"]="$CUSTOM_REG_ORGDB_VERSION"; CUSTOM_KEGG_ORG["$genome"]="$CUSTOM_REG_KEGG_ORG"
    CUSTOM_PROFILE_DIR["$genome"]="$CUSTOM_REG_PROFILE_DIR"; CUSTOM_ANNOTATION_FINGERPRINT["$genome"]="$CUSTOM_REG_ANNOTATION_FINGERPRINT"
    local key
    for key in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
        local path_var="CUSTOM_${key}[$genome]" sha_var="CUSTOM_${key}_SHA256[$genome]" reg_path="CUSTOM_REG_${key}" reg_sha="CUSTOM_REG_${key}_SHA256"
        printf -v "$path_var" '%s' "${!reg_path:-}"; printf -v "$sha_var" '%s' "${!reg_sha:-}"
    done
    CUSTOM_GENOME_SIZE["$genome"]="$CUSTOM_REG_GENOME_SIZE"; CUSTOM_ASSEMBLY_DIR="${CUSTOM_REG_ASSEMBLY_DIR:-$(dirname "$LOCAL_CUSTOM_FASTA")}"; RUN_GENOME_SIZE="$CUSTOM_REG_EFFECTIVE_GENOME_SIZE"; GENOME_SIZE_METHOD="$CUSTOM_REG_GENOME_SIZE_METHOD"
    verify_custom_profile_assets "$genome"; GENOME_MODE="local_explicit"
}

bytes_for_file() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo 0
        return 0
    fi
    stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || wc -c < "$file"
}

sum_selected_fastq_bytes() {
    local total=0
    local i r1 r2 size

    for i in "${!SAMPLE_NAMES[@]}"; do
        r1="${R1_FILES[$i]}"
        r2="${R2_FILES_ARR[$i]}"

        size=$(bytes_for_file "$r1")
        total=$(( total + size ))

        if [[ -n "$r2" ]]; then
            size=$(bytes_for_file "$r2")
            total=$(( total + size ))
        fi
    done

    echo "$total"
}

human_bytes() {
    local bytes="${1:-0}"
    awk -v b="$bytes" 'BEGIN {
        split("B KB MB GB TB PB", unit)
        i = 1
        while (b >= 1024 && i < 6) {
            b = b / 1024
            i++
        }
        printf "%.1f %s", b, unit[i]
    }'
}

nearest_existing_path() {
    local target="$1"

    if [[ -z "$target" ]]; then
        echo "."
        return 0
    fi

    while [[ ! -e "$target" ]]; do
        target=$(dirname "$target")
        if [[ -z "$target" || "$target" == "." ]]; then
            echo "."
            return 0
        fi
        if [[ "$target" == "/" ]]; then
            echo "/"
            return 0
        fi
    done

    echo "$target"
}

available_bytes_for_path() {
    local target="$1"
    local check_path avail_kb

    check_path=$(nearest_existing_path "$target")
    avail_kb=$(df -Pk "$check_path" | awk 'NR==2 {print $4}')

    if [[ -z "$avail_kb" ]]; then
        echo 0
    else
        echo $(( avail_kb * 1024 ))
    fi
}

csv_escape() {
    local value="${1:-}"
    value="${value//\"/\"\"}"
    printf '"%s"' "$value"
}

csv_line() {
    local first=true
    local field

    for field in "$@"; do
        if $first; then
            first=false
        else
            printf ','
        fi
        csv_escape "$field"
    done
    printf '\n'
}

gzip_file_valid() {
    local file="$1"

    [[ -s "$file" ]] || return 1

    if [[ "$file" == *.gz ]]; then
        command_exists gzip || die "gzip is required to verify compressed files, but gzip was not found."
        gzip -t "$file" >/dev/null 2>&1
    else
        return 0
    fi
}

# append_command LABEL CMD [ARGS...]
# Writes a labelled command to the reproducibility log.
append_command() {
    local label="$1"
    shift
    local block
    block="$(printf '# %s\n' "$label"; printf '%q ' "$@"; printf '\n\n')"
    {
        flock -x 9
        printf '%s' "$block" >&9
    } 9>> "$COMMANDS_FILE"
}

# ─────────────────────────────────────────────────────────────
# Resume state
#
# Stored as one `KEY=base64(value)` pair per line and parsed by reading
# and splitting each line -- never sourced or eval'd. The previous format
# was a shell snippet loaded with `source`, which meant a hand-edited or
# corrupted resume file could execute arbitrary shell code. Base64-encoding
# every value also sidesteps any need to shell-quote paths/names that
# might contain spaces, quotes, or other special characters.
# ─────────────────────────────────────────────────────────────

RESUME_MODE=false
SRA_MODE=false
SRA_ACCESSION_USED=""
RESUME_STATE_SCHEMA=4

# _rs_b64 VALUE
# Base64-encodes a single value for storage. The loader reverses this with
# `base64 -d` -- neither side ever interprets the value as shell syntax.
_rs_b64() {
    printf '%s' "$1" | base64 -w0
}

write_resume_state() {
    local state_file="$OUTPUT_DIR/.pepatac_resume_state.kv"
    local tmp_file="${state_file}.tmp.$$"
    local i

    {
        echo "# Auto-generated by PEPATAC_run.sh -- do not hand-edit."
        echo "# key=base64(value) per line. Parsed as plain data only --"
        echo "# never sourced or eval'd as shell."
        printf 'SCHEMA_VERSION=%s\n' "$(_rs_b64 "$RESUME_STATE_SCHEMA")"
        printf 'INPUT_DIR=%s\n'          "$(_rs_b64 "$INPUT_DIR")"
        printf 'SRA_ACCESSION_USED=%s\n' "$(_rs_b64 "${SRA_ACCESSION_USED:-}")"
        printf 'OUTPUT_DIR=%s\n'         "$(_rs_b64 "$OUTPUT_DIR")"
        printf 'GENOME=%s\n'             "$(_rs_b64 "$GENOME")"
        printf 'GENOME_MODE=%s\n'        "$(_rs_b64 "${GENOME_MODE:-refgenie}")"
        printf 'GENOME_IS_CUSTOM=%s\n'   "$(_rs_b64 "${GENOME_IS_CUSTOM:-false}")"
        printf 'GENOME_BUILD_METHOD=%s\n' "$(_rs_b64 "${GENOME_BUILD_METHOD:-refgenie}")"
        printf 'LOCAL_BT2_INDEX=%s\n'    "$(_rs_b64 "${LOCAL_BT2_INDEX:-}")"
        printf 'LOCAL_CHROM_SIZES=%s\n'  "$(_rs_b64 "${LOCAL_CHROM_SIZES:-}")"
        printf 'LOCAL_CUSTOM_FASTA=%s\n' "$(_rs_b64 "${LOCAL_CUSTOM_FASTA:-}")"
        printf 'LOCAL_CUSTOM_FASTA_SHA256=%s\n' "$(_rs_b64 "${LOCAL_CUSTOM_FASTA_SHA256:-}")"
        printf 'LOCAL_CHROM_SIZES_SHA256=%s\n' "$(_rs_b64 "${LOCAL_CHROM_SIZES_SHA256:-}")"
        printf 'LOCAL_BT2_INDEX_SHA256=%s\n' "$(_rs_b64 "${LOCAL_BT2_INDEX_SHA256:-}")"
        printf 'CUSTOM_ASSEMBLY_DIR=%s\n' "$(_rs_b64 "${CUSTOM_ASSEMBLY_DIR:-}")"
        printf 'EFFECTIVE_GENOME_SIZE=%s\n' "$(_rs_b64 "${RUN_GENOME_SIZE:-}")"
        printf 'GENOME_SIZE_METHOD=%s\n' "$(_rs_b64 "${GENOME_SIZE_METHOD:-}")"
        # Generic annotation fields, written for every genome (built-in or
        # custom), not just inside the GENOME_IS_CUSTOM block below. This
        # is what lets a resumed run bind the exact annotation assets it
        # already used instead of re-resolving from whatever TxDb/OrgDb
        # happens to be installed now -- see load_resume_state's use of
        # these and the "trust resume state" check at the call site.
        printf 'ANNOTATION_STATUS=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_STATUS:-}")"
        printf 'ANNOTATION_SOURCE=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_SOURCE:-}")"
        printf 'ANNOTATION_FINGERPRINT=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_FINGERPRINT:-}")"
        printf 'ANNOTATION_QC_FAILURE_REASON=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_QC_FAILURE_REASON:-}")"
        printf 'ANNOTATION_TXDB_PACKAGE=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TXDB_PACKAGE:-}")"
        printf 'ANNOTATION_TXDB_PACKAGE_VERSION=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TXDB_VERSION:-}")"
        printf 'ANNOTATION_TXDB_SQLITE=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TXDB_SQLITE:-}")"
        printf 'ANNOTATION_TXDB_SQLITE_SHA256=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TXDB_SQLITE_SHA256:-}")"
        printf 'ANNOTATION_ORGDB_PACKAGE=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_ORGDB_PACKAGE:-}")"
        printf 'ANNOTATION_ORGDB_PACKAGE_VERSION=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_ORGDB_VERSION:-}")"
        printf 'ANNOTATION_ORGDB_SQLITE=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_ORGDB_SQLITE:-}")"
        printf 'ANNOTATION_ORGDB_SQLITE_SHA256=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_ORGDB_SQLITE_SHA256:-}")"
        printf 'ANNOTATION_TSS_BED=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TSS_BED:-}")"
        printf 'ANNOTATION_TSS_BED_SHA256=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_TSS_BED_SHA256:-}")"
        printf 'ANNOTATION_FEATURE_BED=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_FEATURE_BED:-}")"
        printf 'ANNOTATION_FEATURE_BED_SHA256=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_FEATURE_BED_SHA256:-}")"
        printf 'ANNOTATION_CHROM_SIZES=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_CHROM_SIZES:-}")"
        printf 'ANNOTATION_CHROM_SIZES_SHA256=%s\n' "$(_rs_b64 "${RUN_ANNOTATION_CHROM_SIZES_SHA256:-}")"
        if ${GENOME_IS_CUSTOM:-false}; then
            printf 'CUSTOM_UCSC_ASSEMBLY=%s\n' "$(_rs_b64 "${CUSTOM_UCSC_ASSEMBLY[$GENOME]:-}")"
            printf 'CUSTOM_TXDB_PACKAGE=%s\n' "$(_rs_b64 "${CUSTOM_TXDB_PKG[$GENOME]:-}")"
            printf 'CUSTOM_TXDB_PACKAGE_VERSION=%s\n' "$(_rs_b64 "${CUSTOM_TXDB_VERSION[$GENOME]:-}")"
            printf 'CUSTOM_ORGDB_PACKAGE=%s\n' "$(_rs_b64 "${CUSTOM_ORGDB_PKG[$GENOME]:-}")"
            printf 'CUSTOM_ORGDB_PACKAGE_VERSION=%s\n' "$(_rs_b64 "${CUSTOM_ORGDB_VERSION[$GENOME]:-}")"
            printf 'CUSTOM_KEGG_ORG=%s\n' "$(_rs_b64 "${CUSTOM_KEGG_ORG[$GENOME]:-}")"
            printf 'CUSTOM_PROFILE_DIR=%s\n' "$(_rs_b64 "${CUSTOM_PROFILE_DIR[$GENOME]:-}")"
            printf 'CUSTOM_ANNOTATION_FINGERPRINT=%s\n' "$(_rs_b64 "${CUSTOM_ANNOTATION_FINGERPRINT[$GENOME]:-}")"
            for _k in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
                _pv="CUSTOM_${_k}[$GENOME]"; _sv="CUSTOM_${_k}_SHA256[$GENOME]"
                printf 'CUSTOM_%s=%s\n' "$_k" "$(_rs_b64 "${!_pv:-}")"
                printf 'CUSTOM_%s_SHA256=%s\n' "$_k" "$(_rs_b64 "${!_sv:-}")"
            done
        fi
        printf 'THREADS=%s\n'            "$(_rs_b64 "$THREADS")"
        printf 'PARALLEL_SAMPLES=%s\n'   "$(_rs_b64 "${PARALLEL_SAMPLES:-1}")"
        printf 'PAIRED=%s\n'             "$(_rs_b64 "$PAIRED")"
        printf 'USE_BLACKLIST=%s\n'      "$(_rs_b64 "$USE_BLACKLIST")"
        printf 'BLACKLIST_PATH=%s\n'     "$(_rs_b64 "${BLACKLIST_PATH:-}")"
        printf 'BLACKLIST_SOURCE=%s\n'   "$(_rs_b64 "${BLACKLIST_SOURCE:-none}")"
        printf 'RUN_ID=%s\n'             "$(_rs_b64 "$RUN_ID")"
        printf 'RUN_STARTED=%s\n'        "$(_rs_b64 "$RUN_STARTED")"
        printf 'SAMPLE_COUNT=%s\n'       "$(_rs_b64 "${#SAMPLE_NAMES[@]}")"
        for i in "${!SAMPLE_NAMES[@]}"; do
            printf 'SAMPLE_NAME_%d=%s\n' "$i" "$(_rs_b64 "${SAMPLE_NAMES[$i]}")"
            printf 'SAMPLE_R1_%d=%s\n'   "$i" "$(_rs_b64 "${R1_FILES[$i]}")"
            printf 'SAMPLE_R2_%d=%s\n'   "$i" "$(_rs_b64 "${R2_FILES_ARR[$i]:-}")"
        done
    } > "$tmp_file"
    mv -f "$tmp_file" "$state_file"
}

load_resume_state() {
    local resume_dir="$1"
    local state_file="$resume_dir/.pepatac_resume_state.kv"
    local legacy_file="$resume_dir/.pepatac_resume_state.sh"

    if [[ ! -f "$state_file" ]]; then
        if [[ -f "$legacy_file" ]]; then
            err "Found a resume state file from an older version of this script:"
            err "  $legacy_file"
            err "That format was a sourced shell script. This version never"
            err "sources resume state, by design, so it can't read it."
            err "Re-run that sample set from scratch with this version instead."
            return 1
        fi
        err "No resume state found in: $resume_dir"
        err "(expected: $state_file)"
        err "This folder may not be a previous run of this script, or it"
        err "predates resume support being added."
        return 1
    fi

    local -A rs=()
    local line key val
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        key="${line%%=*}"
        val="${line#*=}"
        if [[ ! "$key" =~ ^[A-Za-z0-9_]+$ ]]; then
            warn "Ignoring malformed resume state line: $line"
            continue
        fi
        rs["$key"]="$(printf '%s' "$val" | base64 -d 2>/dev/null || echo '')"
    done < "$state_file"

    local loaded_schema="${rs[SCHEMA_VERSION]:-0}"
    case "$loaded_schema" in
        2|3|4) ;;
        *)
            err "Resume state schema '$loaded_schema' is not supported by this script (supported: 2, 3, and $RESUME_STATE_SCHEMA)."
            return 1
            ;;
    esac

    INPUT_DIR="${rs[INPUT_DIR]:-}"
    SRA_ACCESSION_USED="${rs[SRA_ACCESSION_USED]:-}"
    OUTPUT_DIR="${rs[OUTPUT_DIR]:-}"
    GENOME="${rs[GENOME]:-}"
    GENOME_MODE="${rs[GENOME_MODE]:-refgenie}"
    [[ "$GENOME_MODE" == "custom" || "$GENOME_MODE" == "local_ucsc" ]] && GENOME_MODE="local_explicit"
    GENOME_IS_CUSTOM="${rs[GENOME_IS_CUSTOM]:-false}"
    GENOME_BUILD_METHOD="${rs[GENOME_BUILD_METHOD]:-refgenie}"
    [[ "$GENOME_BUILD_METHOD" == "ucsc" ]] && GENOME_BUILD_METHOD="local_build"
    LOCAL_CUSTOM_FASTA="${rs[LOCAL_CUSTOM_FASTA]:-}"
    LOCAL_CUSTOM_FASTA_SHA256="${rs[LOCAL_CUSTOM_FASTA_SHA256]:-}"
    LOCAL_BT2_INDEX="${rs[LOCAL_BT2_INDEX]:-}"
    LOCAL_BT2_INDEX_SHA256="${rs[LOCAL_BT2_INDEX_SHA256]:-}"
    LOCAL_CHROM_SIZES="${rs[LOCAL_CHROM_SIZES]:-}"
    LOCAL_CHROM_SIZES_SHA256="${rs[LOCAL_CHROM_SIZES_SHA256]:-}"
    CUSTOM_ASSEMBLY_DIR="${rs[CUSTOM_ASSEMBLY_DIR]:-}"
    if [[ "$GENOME_IS_CUSTOM" == "true" ]]; then
        [[ "$loaded_schema" -ge 3 ]] || { err "This custom-genome resume state is the older GTF-based format. Start a new run and register a validated UCSC/TxDb/OrgDb profile."; return 1; }
        CUSTOM_UCSC_ASSEMBLY["$GENOME"]="${rs[CUSTOM_UCSC_ASSEMBLY]:-$GENOME}"
        CUSTOM_TXDB_PKG["$GENOME"]="${rs[CUSTOM_TXDB_PACKAGE]:-}"; CUSTOM_TXDB_VERSION["$GENOME"]="${rs[CUSTOM_TXDB_PACKAGE_VERSION]:-}"
        CUSTOM_ORGDB_PKG["$GENOME"]="${rs[CUSTOM_ORGDB_PACKAGE]:-}"; CUSTOM_ORGDB_VERSION["$GENOME"]="${rs[CUSTOM_ORGDB_PACKAGE_VERSION]:-}"; CUSTOM_KEGG_ORG["$GENOME"]="${rs[CUSTOM_KEGG_ORG]:-}"
        CUSTOM_PROFILE_DIR["$GENOME"]="${rs[CUSTOM_PROFILE_DIR]:-}"; CUSTOM_ANNOTATION_FINGERPRINT["$GENOME"]="${rs[CUSTOM_ANNOTATION_FINGERPRINT]:-}"
        for _k in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
            _pv="CUSTOM_${_k}[$GENOME]"; _sv="CUSTOM_${_k}_SHA256[$GENOME]"; printf -v "$_pv" '%s' "${rs[CUSTOM_${_k}]:-}"; printf -v "$_sv" '%s' "${rs[CUSTOM_${_k}_SHA256]:-}"
        done
    fi
    RUN_GENOME_SIZE="${rs[EFFECTIVE_GENOME_SIZE]:-}"
    GENOME_SIZE_METHOD="${rs[GENOME_SIZE_METHOD]:-}"
    # Only present from schema 4 onward -- blank on an older resume state,
    # which is exactly what makes the "trust resume state" check at the
    # call site correctly fall through to fresh resolution for those.
    RUN_ANNOTATION_STATUS="${rs[ANNOTATION_STATUS]:-}"
    RUN_ANNOTATION_SOURCE="${rs[ANNOTATION_SOURCE]:-}"
    RUN_ANNOTATION_FINGERPRINT="${rs[ANNOTATION_FINGERPRINT]:-}"
    RUN_ANNOTATION_QC_FAILURE_REASON="${rs[ANNOTATION_QC_FAILURE_REASON]:-}"
    RUN_ANNOTATION_TXDB_PACKAGE="${rs[ANNOTATION_TXDB_PACKAGE]:-}"
    RUN_ANNOTATION_TXDB_VERSION="${rs[ANNOTATION_TXDB_PACKAGE_VERSION]:-}"
    RUN_ANNOTATION_TXDB_SQLITE="${rs[ANNOTATION_TXDB_SQLITE]:-}"
    RUN_ANNOTATION_TXDB_SQLITE_SHA256="${rs[ANNOTATION_TXDB_SQLITE_SHA256]:-}"
    RUN_ANNOTATION_ORGDB_PACKAGE="${rs[ANNOTATION_ORGDB_PACKAGE]:-}"
    RUN_ANNOTATION_ORGDB_VERSION="${rs[ANNOTATION_ORGDB_PACKAGE_VERSION]:-}"
    RUN_ANNOTATION_ORGDB_SQLITE="${rs[ANNOTATION_ORGDB_SQLITE]:-}"
    RUN_ANNOTATION_ORGDB_SQLITE_SHA256="${rs[ANNOTATION_ORGDB_SQLITE_SHA256]:-}"
    RUN_ANNOTATION_TSS_BED="${rs[ANNOTATION_TSS_BED]:-}"
    RUN_ANNOTATION_TSS_BED_SHA256="${rs[ANNOTATION_TSS_BED_SHA256]:-}"
    RUN_ANNOTATION_FEATURE_BED="${rs[ANNOTATION_FEATURE_BED]:-}"
    RUN_ANNOTATION_FEATURE_BED_SHA256="${rs[ANNOTATION_FEATURE_BED_SHA256]:-}"
    RUN_ANNOTATION_CHROM_SIZES="${rs[ANNOTATION_CHROM_SIZES]:-}"
    RUN_ANNOTATION_CHROM_SIZES_SHA256="${rs[ANNOTATION_CHROM_SIZES_SHA256]:-}"
    THREADS="${rs[THREADS]:-}"
    PARALLEL_SAMPLES="${rs[PARALLEL_SAMPLES]:-1}"
    PAIRED="${rs[PAIRED]:-}"
    USE_BLACKLIST="${rs[USE_BLACKLIST]:-}"
    BLACKLIST_PATH="${rs[BLACKLIST_PATH]:-}"
    BLACKLIST_SOURCE="${rs[BLACKLIST_SOURCE]:-none}"
    RUN_ID="${rs[RUN_ID]:-}"
    RUN_STARTED="${rs[RUN_STARTED]:-}"

    local sample_count="${rs[SAMPLE_COUNT]:-0}"
    if [[ ! "$sample_count" =~ ^[0-9]+$ ]]; then
        err "Resume state has an invalid SAMPLE_COUNT ('$sample_count')."
        return 1
    fi

    SAMPLE_NAMES=()
    R1_FILES=()
    R2_FILES_ARR=()
    local i
    for (( i=0; i<sample_count; i++ )); do
        SAMPLE_NAMES+=("${rs[SAMPLE_NAME_${i}]:-}")
        R1_FILES+=("${rs[SAMPLE_R1_${i}]:-}")
        R2_FILES_ARR+=("${rs[SAMPLE_R2_${i}]:-}")
    done

    return 0
}

# read_sample_status SAMPLE
# Returns the on-disk status for a sample from a prior run.
# Echoes PASS, FAIL, or NOT_RUN.
read_sample_status() {
    local sample="$1"
    local status_file="$LOG_DIR/${sample}.status"

    if [[ -f "$status_file" ]]; then
        cat "$status_file"
    else
        echo "NOT_RUN"
    fi
}

# ─────────────────────────────────────────────────────────────
# Run signature
#
# A PASS status file alone only proves that pepatac.py exited 0 at some
# point in the past -- it says nothing about whether this invocation's
# FASTQs, genome assets, blacklist, or pipeline version are still the same
# ones that produced it. Build a fingerprint of everything that can change
# a sample's actual output, so a later run only skips a PASS when nothing
# relevant has changed since it was recorded.
#
# Cheap fingerprints (size+mtime) are used for large files that are
# expensive to hash (FASTQs, Bowtie2 index); real content hashes are used
# for small files where hashing is essentially free (chrom.sizes,
# blacklist BED). Neither proves byte-for-byte content is identical in
# every possible case, but both catch the realistic ways these files
# change between runs (replaced, regenerated, or edited).
# ─────────────────────────────────────────────────────────────

# compute_run_wide_fingerprint
# Everything here is identical for every sample in a single invocation, so
# this is computed exactly once per run (see call site near genome asset
# resolution) rather than once per sample.
compute_run_wide_fingerprint() {
    local pepatac_commit script_hash bt2_dir bt2_base bt2_fp chrom_fp bl_fp

    pepatac_commit="$(git -C "$PEPATAC_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
    script_hash="$(sha256sum "${BASH_SOURCE[0]:-$0}" 2>/dev/null | awk '{print $1}')"
    [[ -z "$script_hash" ]] && script_hash="unknown"

    bt2_dir="$(dirname -- "$RUN_GENOME_INDEX")"
    bt2_base="$(basename -- "$RUN_GENOME_INDEX")"
    bt2_fp="$(find "$bt2_dir" -maxdepth 1 -name "${bt2_base}*" -printf '%f %s %T@\n' 2>/dev/null | sort | sha256sum | awk '{print $1}')"
    [[ -z "$bt2_fp" ]] && bt2_fp="missing"

    if [[ -f "$RUN_CHROM_SIZES" ]]; then
        chrom_fp="$(sha256sum "$RUN_CHROM_SIZES" 2>/dev/null | awk '{print $1}')"
    else
        chrom_fp="missing"
    fi

    if $USE_BLACKLIST && [[ -n "${BLACKLIST_PATH:-}" && -f "$BLACKLIST_PATH" ]]; then
        bl_fp="$(sha256sum "$BLACKLIST_PATH" 2>/dev/null | awk '{print $1}')"
    else
        bl_fp="none"
    fi

    printf 'pepatac_commit=%s|script_version=%s|script_hash=%s|genome=%s|genome_size=%s|bt2_index_fp=%s|chrom_sizes_fp=%s|blacklist_fp=%s|paired=%s|annotation_status=%s|annotation_fp=%s|txdb_sha=%s|orgdb_sha=%s|tss_bed_sha=%s|feature_bed_sha=%s' \
        "$pepatac_commit" "$SCRIPT_VERSION" "$script_hash" "$GENOME" "${RUN_GENOME_SIZE:-}" \
        "$bt2_fp" "$chrom_fp" "$bl_fp" "$PAIRED" \
        "${RUN_ANNOTATION_STATUS:-}" "${RUN_ANNOTATION_FINGERPRINT:-}" \
        "${RUN_ANNOTATION_TXDB_SQLITE_SHA256:-}" "${RUN_ANNOTATION_ORGDB_SQLITE_SHA256:-}" \
        "${RUN_ANNOTATION_TSS_BED_SHA256:-}" "${RUN_ANNOTATION_FEATURE_BED_SHA256:-}"
}

# sample_signature INDEX
# Combines the run-wide fingerprint with this sample's own FASTQ
# fingerprint into one signature. Requires RUN_WIDE_FINGERPRINT to already
# be set (computed once, after genome assets are resolved).
sample_signature() {
    local i="$1"
    local r1="${R1_FILES[$i]}" r2="${R2_FILES_ARR[$i]:-}"
    local r1_fp r2_fp

    if [[ -f "$r1" ]]; then
        r1_fp="$(stat -c '%s %Y' "$r1" 2>/dev/null || echo missing)"
    else
        r1_fp="missing"
    fi
    if [[ -n "$r2" ]]; then
        if [[ -f "$r2" ]]; then
            r2_fp="$(stat -c '%s %Y' "$r2" 2>/dev/null || echo missing)"
        else
            r2_fp="missing"
        fi
    else
        r2_fp="none"
    fi

    printf '%s|r1=%s:%s|r2=%s:%s' "$RUN_WIDE_FINGERPRINT" "$r1" "$r1_fp" "$r2" "$r2_fp" \
        | sha256sum | awk '{print $1}'
}

# write_sample_status SAMPLE STATUS [INDEX]
# When STATUS is PASS and INDEX is given, also records the run signature
# that produced this PASS (see sample_signature above), so a later
# invocation can tell whether inputs, genome assets, blacklist, or
# pipeline version have changed since. A FAIL clears any old PASS
# signature so a stale one is never left behind claiming to describe a
# PASS that no longer holds. The separate .attempt.sig (see
# write_attempt_signature below) is intentionally NOT touched here on
# FAIL -- it needs to survive a failure so the next run can compare
# against it -- but is cleaned up here on PASS, since a successful run's
# own .sig supersedes it.
write_sample_status() {
    local sample="$1"
    local status="$2"
    local idx="${3:-}"
    local sig_file="$LOG_DIR/${sample}.sig"

    echo "$status" > "$LOG_DIR/${sample}.status"

    if [[ "$status" == "PASS" && -n "$idx" ]]; then
        sample_signature "$idx" > "$sig_file"
        rm -f "$LOG_DIR/${sample}.attempt.sig"
    else
        rm -f "$sig_file"
    fi
}

# write_attempt_signature INDEX
# Records the run signature THIS ATTEMPT is about to use, before
# pepatac.py even starts. Unlike the PASS signature above, this is
# written unconditionally at the start of every attempt and, critically,
# is left in place if the attempt fails (write_sample_status only ever
# removes it on PASS). That means a later run, finding a prior FAIL, can
# compare its current signature against what the failed attempt actually
# used: if they match, resuming with -R (recover) is safe; if they
# differ -- FASTQ replaced, genome changed, pipeline updated -- resuming
# with -R would let pypiper pick up from checkpoints built against inputs
# that no longer apply, exactly as unsafe as resuming a stale PASS
# without -N. See the FAIL branch in Pass 1 below.
write_attempt_signature() {
    local idx="$1"
    local sample="${SAMPLE_NAMES[$idx]}"
    sample_signature "$idx" > "$LOG_DIR/${sample}.attempt.sig"
}

# ─────────────────────────────────────────────────────────────
# Run output initialisation
# ─────────────────────────────────────────────────────────────

initialize_run_outputs() {
    header "Run Records"

    mkdir -p "$OUTPUT_DIR"
    LOG_DIR="$OUTPUT_DIR/logs"
    mkdir -p "$LOG_DIR"

    RUN_MANIFEST="$OUTPUT_DIR/run_manifest.txt"
    COMMANDS_FILE="$OUTPUT_DIR/pepatac_commands.sh"
    R_SAMPLE_SHEET="$OUTPUT_DIR/samples_for_R_template.csv"
    QC_SUMMARY="$OUTPUT_DIR/qc_summary.csv"
    R_DETECTED_SAMPLE_SHEET="$OUTPUT_DIR/samples_for_R_autodetected.csv"

    write_run_manifest
    write_resume_state
    write_r_sample_sheet_template
    initialize_commands_file

    ok "Run manifest: $RUN_MANIFEST"
    ok "R sample sheet template: $R_SAMPLE_SHEET"
    ok "Command log: $COMMANDS_FILE"
}

write_run_manifest() {
    local input_bytes library_label
    input_bytes=$(sum_selected_fastq_bytes)
    library_label=$( $PAIRED && echo "Paired-end" || echo "Single-end" )

    {
        echo "PEPATAC RUN MANIFEST"
        echo "===================="
        echo "Run ID: $RUN_ID"
        echo "Script version: $SCRIPT_VERSION"
        echo "Run started: $RUN_STARTED"
        echo "Host: $(hostname 2>/dev/null || echo unknown)"
        echo "User: ${USER:-unknown}"
        echo "Operating system: $(uname -a 2>/dev/null || echo unknown)"
        echo ""
        echo "RUN SETTINGS"
        echo "------------"
        echo "Input folder: $INPUT_DIR"
        if [[ -n "${SRA_ACCESSION_USED:-}" ]]; then
            echo "Input source: Downloaded from SRA/ENA, accession $SRA_ACCESSION_USED"
        fi
        echo "Output folder: $OUTPUT_DIR"
        echo "Genome: $GENOME"
        echo "Genome mode: ${GENOME_MODE:-refgenie}"
        echo "Genome custom: $( ${GENOME_IS_CUSTOM:-false} && echo yes || echo no )"
        if $GENOME_IS_CUSTOM; then
            echo "User-defined genome profile registry: $(custom_genome_registry_file "$GENOME")"
        fi
        if [[ "${GENOME_MODE:-refgenie}" == "local_explicit" ]]; then
            echo "Local explicit FASTA:         ${LOCAL_CUSTOM_FASTA:-}"
            echo "Local explicit FASTA SHA256:  ${LOCAL_CUSTOM_FASTA_SHA256:-}"
            echo "Local explicit Bowtie2 index: ${LOCAL_BT2_INDEX:-}"
            echo "Local explicit BT2 SHA256:    ${LOCAL_BT2_INDEX_SHA256:-}"
            echo "Local explicit chrom sizes:   ${LOCAL_CHROM_SIZES:-}"
            echo "Local chrom sizes SHA256:     ${LOCAL_CHROM_SIZES_SHA256:-}"
            if $GENOME_IS_CUSTOM; then
                echo "UCSC assembly:                ${CUSTOM_UCSC_ASSEMBLY[$GENOME]:-}"
                echo "TxDb source package:          ${CUSTOM_TXDB_PKG[$GENOME]:-} ${CUSTOM_TXDB_VERSION[$GENOME]:-}"
                echo "OrgDb source package:         ${CUSTOM_ORGDB_PKG[$GENOME]:-} ${CUSTOM_ORGDB_VERSION[$GENOME]:-}"
                echo "Frozen TxDb SQLite:           ${CUSTOM_TXDB_SQLITE[$GENOME]:-}"
                echo "Frozen TxDb SHA256:           ${CUSTOM_TXDB_SQLITE_SHA256[$GENOME]:-}"
                echo "Frozen OrgDb SQLite:          ${CUSTOM_ORGDB_SQLITE[$GENOME]:-}"
                echo "Frozen OrgDb SHA256:          ${CUSTOM_ORGDB_SQLITE_SHA256[$GENOME]:-}"
                echo "PEPATAC TSS BED:              ${CUSTOM_TSS_BED[$GENOME]:-}"
                echo "PEPATAC feature BED:          ${CUSTOM_FEATURE_BED[$GENOME]:-}"
            fi
        fi
        echo "Effective genome size: ${RUN_GENOME_SIZE:-<not resolved yet>}"
        echo "Genome size method: ${GENOME_SIZE_METHOD:-<not resolved yet>}"
        echo "Threads: $THREADS"
        echo "Concurrent samples: ${PARALLEL_SAMPLES:-1}"
        echo "Maximum requested cores: $(( THREADS * ${PARALLEL_SAMPLES:-1} ))"
        echo "Library: $library_label"
        echo "Blacklist: $( $USE_BLACKLIST && echo "$BLACKLIST_PATH" || echo No )"
        echo "Sample count: ${#SAMPLE_NAMES[@]}"
        echo "Selected FASTQ size: $(human_bytes "$input_bytes") ($input_bytes bytes)"
        echo ""
        echo "ANNOTATION PROFILE"
        echo "------------------"
        # Generic for every genome, built-in or custom -- previously this
        # information only appeared in the custom/local section above, so
        # a built-in genome's annotation status, fingerprint, frozen
        # paths, hashes, and any override reason were invisible in the
        # human-readable manifest even after run.sh started tracking them
        # in the machine-readable snapshot.
        echo "Status: ${RUN_ANNOTATION_STATUS:-<not resolved yet>}"
        echo "Source: ${RUN_ANNOTATION_SOURCE:-<not resolved yet>}"
        echo "Fingerprint: ${RUN_ANNOTATION_FINGERPRINT:-<none>}"
        echo "Builder schema: ${ANNOTATION_BUILDER_SCHEMA:-<none>}"
        echo "TxDb package/version: ${RUN_ANNOTATION_TXDB_PACKAGE:-<none>} ${RUN_ANNOTATION_TXDB_VERSION:-}"
        echo "TxDb SQLite: ${RUN_ANNOTATION_TXDB_SQLITE:-<none>}"
        echo "TxDb SHA-256: ${RUN_ANNOTATION_TXDB_SQLITE_SHA256:-<none>}"
        echo "OrgDb package/version: ${RUN_ANNOTATION_ORGDB_PACKAGE:-<none>} ${RUN_ANNOTATION_ORGDB_VERSION:-}"
        echo "OrgDb SQLite: ${RUN_ANNOTATION_ORGDB_SQLITE:-<none>}"
        echo "OrgDb SHA-256: ${RUN_ANNOTATION_ORGDB_SQLITE_SHA256:-<none>}"
        echo "TSS BED: ${RUN_ANNOTATION_TSS_BED:-<none>}"
        echo "TSS BED SHA-256: ${RUN_ANNOTATION_TSS_BED_SHA256:-<none>}"
        echo "Feature BED: ${RUN_ANNOTATION_FEATURE_BED:-<none>}"
        echo "Feature BED SHA-256: ${RUN_ANNOTATION_FEATURE_BED_SHA256:-<none>}"
        echo "Chrom sizes: ${RUN_ANNOTATION_CHROM_SIZES:-<none>}"
        echo "Chrom sizes SHA-256: ${RUN_ANNOTATION_CHROM_SIZES_SHA256:-<none>}"
        if [[ "${RUN_ANNOTATION_STATUS:-}" == "disabled_by_override" ]]; then
            echo "QC failure reason: ${RUN_ANNOTATION_QC_FAILURE_REASON:-<not recorded>}"
            echo "NOTE: --TSS-name/--anno-name were omitted for this run via PEPATAC_ALLOW_MISSING_ANNOTATION_QC=1."
        fi
        echo ""
        echo "PATHS"
        echo "-----"
        echo "PEPATAC directory: $PEPATAC_DIR"
        echo "PEPATAC pipeline: $PIPELINE"
        echo "Refgenie config: $REFGENIE_CONFIG"
        echo "Blacklist directory: $BLACKLIST_DIR"
        echo "Logs directory: $LOG_DIR"
        echo "R sample sheet template: $R_SAMPLE_SHEET"
        echo "QC summary: $QC_SUMMARY"
        echo "Autodetected R sample sheet: $R_DETECTED_SAMPLE_SHEET"
        echo "Command log: $COMMANDS_FILE"
        echo ""
        echo "SAMPLES"
        echo "-------"
        local i r2_label
        for i in "${!SAMPLE_NAMES[@]}"; do
            r2_label="${R2_FILES_ARR[$i]}"
            [[ -z "$r2_label" ]] && r2_label="NA"
            echo "$((i+1)). ${SAMPLE_NAMES[$i]}"
            echo "   R1: ${R1_FILES[$i]}"
            echo "   R2: $r2_label"
        done
    } > "$RUN_MANIFEST"
}

initialize_commands_file() {
    {
        echo "#!/usr/bin/env bash"
        echo "# Commands generated by PEPATAC_run.sh $SCRIPT_VERSION"
        echo "# Run ID: $RUN_ID"
        echo "# Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "# Reproducibility log. Can be rerun if Conda and the same input files are still available."
        echo "set -euo pipefail"
        printf 'ENV_NAME=%q\n' "$ENV_NAME"
        echo 'CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"'
        echo '[[ -f "$CONDA_SH" ]] && source "$CONDA_SH"'
        echo 'in_env() { conda run --no-capture-output -n "$ENV_NAME" "$@"; }'
        echo ""
    } > "$COMMANDS_FILE"
    chmod +x "$COMMANDS_FILE"
}

write_r_sample_sheet_template() {
    local i sample r1 r2 sample_out log_file library_label
    library_label=$( $PAIRED && echo "paired" || echo "single" )

    {
        csv_line "SampleID" "Condition" "Replicate" "R1" "R2" "Library" "Genome" "OutputDir" "LogFile" "bamReads" "Peaks" "PeakCaller"
        for i in "${!SAMPLE_NAMES[@]}"; do
            sample="${SAMPLE_NAMES[$i]}"
            r1="${R1_FILES[$i]}"
            r2="${R2_FILES_ARR[$i]}"
            sample_out="$OUTPUT_DIR/$sample"
            log_file="$LOG_DIR/${sample}.log"
            csv_line "$sample" "" "" "$r1" "$r2" "$library_label" "$GENOME" "$sample_out" "$log_file" "" "" "narrow"
        done
    } > "$R_SAMPLE_SHEET"
}

# ─────────────────────────────────────────────────────────────
# Disk space
# ─────────────────────────────────────────────────────────────

check_disk_space() {
    header "Disk Space Check"

    local input_bytes available_bytes recommended_bytes home_available_bytes
    input_bytes=$(sum_selected_fastq_bytes)
    available_bytes=$(available_bytes_for_path "$OUTPUT_DIR")
    home_available_bytes=$(available_bytes_for_path "$HOME")

    # 6× FASTQ size accounts for BAMs, filtered BAMs, bigWigs,
    # peak files, trimmed FASTQs, and intermediate outputs.
    recommended_bytes=$(( input_bytes * 6 ))

    echo -e "  ${BOLD}Selected FASTQ size:${RESET}   $(human_bytes "$input_bytes")"
    echo -e "  ${BOLD}Free space for output:${RESET}  $(human_bytes "$available_bytes")"
    echo -e "  ${BOLD}Free space under HOME:${RESET}  $(human_bytes "$home_available_bytes")"
    echo -e "  ${DIM}Rule of thumb used here: output free space should be at least ~6× selected FASTQ size.${RESET}"
    echo -e "  ${DIM}Genome downloads/builds may need extra space under $HOME/refgenie.${RESET}"
    blank

    if [[ "$input_bytes" -gt 0 && "$available_bytes" -lt "$recommended_bytes" ]]; then
        warn "Output disk space may be low for this run."
        echo -e "  ${DIM}Recommended free output space: $(human_bytes "$recommended_bytes")${RESET}"
        echo "  Continue anyway? [y/N]"
        read -r -p "  > " DISK_CONFIRM
        [[ "${DISK_CONFIRM,,}" == "y" ]] || die "Aborted due to low disk space."
    else
        ok "Output disk space looks reasonable."
    fi
}

# ─────────────────────────────────────────────────────────────
# QC output helpers
# ─────────────────────────────────────────────────────────────

# Ambiguity side-channel for detect_best_bam()/detect_best_peak()
#
# Both functions are called as `result="$(detect_best_bam ...)"` at every
# existing call site. Command substitution always runs its command in a
# subshell, so a plain variable assignment made inside the function would
# vanish the instant that subshell exits -- it can never reach the caller.
# A file survives that boundary. Both functions are only ever called
# sequentially in the main process (never inside a background/parallel
# worker -- see Pass 2/Pass 3 below), so one fixed path, overwritten on
# each call and read back immediately after, is safe with no locking.
#
# read_detect_ambiguous
# Prints one candidate path per line from the most recent detect_best_bam()
# or detect_best_peak() call. Empty output means that call was NOT
# ambiguous (it either resolved to exactly one file or found nothing).
# Must be read before the next detect_* call, which overwrites it.
read_detect_ambiguous() {
    local ambig_file="${LOG_DIR:-/tmp}/.detect_ambiguous.tmp"
    [[ -f "$ambig_file" ]] && cat "$ambig_file"
    return 0
}

detect_best_bam() {
    local sample_out="$1"
    local pattern
    local -a matches
    local ambig_file="${LOG_DIR:-/tmp}/.detect_ambiguous.tmp"

    : > "$ambig_file"

    for pattern in "*final*.bam" "*filtered*.bam" "*dedup*.bam" "*sort*.bam" "*.bam"; do
        mapfile -t matches < <(find "$sample_out" -type f -name "$pattern" 2>/dev/null | sort)
        if [[ ${#matches[@]} -eq 1 ]]; then
            echo "${matches[0]}"
            return 0
        elif [[ ${#matches[@]} -gt 1 ]]; then
            printf '%s\n' "${matches[@]}" > "$ambig_file"
            echo ""
            return 0
        fi
        # Zero matches at this tier -- fall through to the next, broader pattern.
    done

    echo ""
}

detect_best_peak() {
    local sample_out="$1"
    local pattern
    local -a matches
    local ambig_file="${LOG_DIR:-/tmp}/.detect_ambiguous.tmp"

    : > "$ambig_file"

    # PEPATAC writes *_rmBlacklist.narrowPeak alongside the raw *_peaks.narrowPeak
    # whenever blacklist filtering is enabled. The raw file's name still ends in
    # "peaks.narrowPeak" and would otherwise match first, silently discarding the
    # cleaned output and feeding blacklisted/repeat-artifact regions downstream.
    for pattern in "*_rmBlacklist.narrowPeak" "*peaks.narrowPeak" "*.narrowPeak" "*summits.bed" "*.bed"; do
        mapfile -t matches < <(find "$sample_out" -type f -name "$pattern" 2>/dev/null | sort)
        if [[ ${#matches[@]} -eq 1 ]]; then
            echo "${matches[0]}"
            return 0
        elif [[ ${#matches[@]} -gt 1 ]]; then
            printf '%s\n' "${matches[@]}" > "$ambig_file"
            echo ""
            return 0
        fi
    done

    echo ""
}

# peak_caller_for PEAK_FILE
# Reports the actual format of the detected peak file instead of assuming
# "narrow" -- detect_best_peak() can fall back to a summits/generic BED file
# when no proper narrowPeak exists, and DiffBind cannot parse that as narrow
# format (wrong column expectations).
peak_caller_for() {
    local peak_file="$1"
    case "$(basename "${peak_file:-}")" in
        *.narrowPeak) echo "narrow" ;;
        "") echo "narrow" ;;   # nothing detected -- keep the prior default
        *)  echo "bed" ;;
    esac
}

count_peak_lines() {
    local peak_file="$1"

    if [[ -z "$peak_file" || ! -f "$peak_file" ]]; then
        echo "NA"
        return 0
    fi

    if [[ "$peak_file" == *.gz ]]; then
        gzip -cd "$peak_file" 2>/dev/null | wc -l | awk '{print $1}'
    else
        wc -l < "$peak_file" | awk '{print $1}'
    fi
}

write_qc_summary() {
    local i sample status exit_code log_file sample_out bam_file peak_file peak_count bam_size bam_size_label

    header "End-of-Run QC Summary"

    {
        csv_line "SampleID" "Status" "ExitCode" "PeakCount" "BamSizeBytes" "BamSizeHuman" "BamFile" "PeakFile" "OutputDir" "LogFile" "R1" "R2"
        for i in "${!SAMPLE_NAMES[@]}"; do
            sample="${SAMPLE_NAMES[$i]}"
            status="${SAMPLE_STATUS[$i]:-NOT_RUN}"
            exit_code="${SAMPLE_EXIT_CODES[$i]:-NA}"
            log_file="${SAMPLE_LOGS[$i]:-$LOG_DIR/${sample}.log}"
            sample_out="${SAMPLE_OUT_DIRS[$i]:-$OUTPUT_DIR/$sample}"
            bam_file="${SAMPLE_BAMS[$i]:-}"
            peak_file="${SAMPLE_PEAKS[$i]:-}"
            peak_count="${SAMPLE_PEAK_COUNTS[$i]:-NA}"
            bam_size="${SAMPLE_BAM_SIZES[$i]:-0}"
            bam_size_label=$(human_bytes "$bam_size")
            csv_line "$sample" "$status" "$exit_code" "$peak_count" "$bam_size" "$bam_size_label" "$bam_file" "$peak_file" "$sample_out" "$log_file" "${R1_FILES[$i]}" "${R2_FILES_ARR[$i]}"
        done
    } > "$QC_SUMMARY"

    printf "  %-32s  %-8s  %-10s  %-10s  %s\n" "Sample" "Status" "Peaks" "BAM size" "Log"
    printf "  %-32s  %-8s  %-10s  %-10s  %s\n" "------" "------" "-----" "--------" "---"

    for i in "${!SAMPLE_NAMES[@]}"; do
        sample="${SAMPLE_NAMES[$i]}"
        status="${SAMPLE_STATUS[$i]:-NOT_RUN}"
        peak_count="${SAMPLE_PEAK_COUNTS[$i]:-NA}"
        bam_size="${SAMPLE_BAM_SIZES[$i]:-0}"
        bam_size_label=$(human_bytes "$bam_size")
        log_file="${SAMPLE_LOGS[$i]:-$LOG_DIR/${sample}.log}"
        printf "  %-32s  %-8s  %-10s  %-10s  %s\n" "$sample" "$status" "$peak_count" "$bam_size_label" "$log_file"
    done

    blank
    ok "QC summary written: $QC_SUMMARY"
}

write_detected_r_sample_sheet() {
    local i sample sample_out log_file bam_file peak_file

    {
        csv_line "SampleID" "Condition" "Replicate" "bamReads" "Peaks" "PeakCaller" "Status" "PeakCount" "LogFile"
        for i in "${!SAMPLE_NAMES[@]}"; do
            sample="${SAMPLE_NAMES[$i]}"
            sample_out="${SAMPLE_OUT_DIRS[$i]:-$OUTPUT_DIR/$sample}"
            log_file="${SAMPLE_LOGS[$i]:-$LOG_DIR/${sample}.log}"
            bam_file="${SAMPLE_BAMS[$i]:-}"
            peak_file="${SAMPLE_PEAKS[$i]:-}"
            csv_line "$sample" "" "" "$bam_file" "$peak_file" "$(peak_caller_for "$peak_file")" "${SAMPLE_STATUS[$i]:-NOT_RUN}" "${SAMPLE_PEAK_COUNTS[$i]:-NA}" "$log_file"
        done
    } > "$R_DETECTED_SAMPLE_SHEET"

    ok "Autodetected R sample sheet written: $R_DETECTED_SAMPLE_SHEET"
}

# ─────────────────────────────────────────────────────────────
# Conda / launcher
# ─────────────────────────────────────────────────────────────

source_conda_or_die() {
    local conda_sh="$HOME/miniconda3/etc/profile.d/conda.sh"

    if [[ -f "$conda_sh" ]]; then
        # shellcheck source=/dev/null
        source "$conda_sh"
    else
        die "Conda not found. Run install.sh first. Expected: $conda_sh"
    fi
}

check_launcher_requirements() {
    header "Launcher Check"

    source_conda_or_die
    ok "Conda startup script found."

    if command_exists setsid; then
        ok "setsid found for strong Ctrl+C cleanup."
    else
        warn "setsid not found; Ctrl+C cleanup will use a less robust fallback."
    fi

    if ! conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
        die "Conda environment '$ENV_NAME' not found. Run install.sh first."
    fi
    ok "Conda environment '$ENV_NAME' found."

    if [[ ! -f "$PIPELINE" ]]; then
        die "Pipeline not found at $PIPELINE. Run install.sh first."
    fi
    ok "Pipeline found: $PIPELINE"

    if [[ ! -f "$REFGENIE_CONFIG" ]]; then
        die "Refgenie config not found at $REFGENIE_CONFIG. Run install.sh first."
    fi
    ok "Refgenie config found: $REFGENIE_CONFIG"

    export REFGENIE="$REFGENIE_CONFIG"
}

# ─────────────────────────────────────────────────────────────
# Genome assets (Refgenie-managed path)
# ─────────────────────────────────────────────────────────────

pull_refgenie_asset_once() {
    local asset="$1"
    local seek="$2"

    if asset_ok "$seek"; then
        ok "$seek already present."
    else
        label "Pulling $asset..."
        # shellcheck disable=SC2086
        # NOTE: refgenie pull exits 0 even when the asset isn't available on any
        # subscribed server — it only prints a warning line. So the command's own
        # exit code can't be trusted; verify the asset actually landed on disk.
        in_env refgenie pull -c "$REFGENIE_CONFIG" "$asset" || true
        if asset_ok "$seek"; then
            ok "$asset installed."
        else
            warn "$asset was not obtained (likely unavailable on this server)."
            return 1
        fi
    fi
}

ensure_genome_assets() {
    header "Reference Sequence Assets"

    if [[ "${GENOME_MODE:-refgenie}" == "local_explicit" ]]; then
        if $GENOME_IS_CUSTOM; then
            if $RESUME_MODE && load_reference_snapshot_into_runner; then
                :
            else
                _acquire_custom_genome_build_lock "$GENOME"
                assert_custom_registry_unchanged_since_selection "$GENOME"
                load_custom_genome_registry "$GENOME" \
                    || die "User-defined genome profile registry disappeared before reuse: $GENOME"
                verify_registered_custom_reference "$GENOME"
                _release_custom_genome_build_lock
            fi
            if [[ -z "${RUN_GENOME_SIZE:-}" || -z "${GENOME_SIZE_METHOD:-}" ]]; then
                resolve_effective_genome_size "$GENOME" "$LOCAL_CHROM_SIZES"
            fi
            write_reference_snapshot "$GENOME"
        else
            validate_genome_reference_integrity "$LOCAL_BT2_INDEX" "$LOCAL_CHROM_SIZES"
        fi
        ok "Explicit local genome sequence assets are ready."
        return 0
    fi

    if [[ "${GENOME_BUILD_METHOD:-refgenie}" == "local_build" ]]; then
        if $GENOME_IS_CUSTOM; then
            _acquire_custom_genome_build_lock "$GENOME"
            assert_custom_registry_unchanged_since_selection "$GENOME"
        fi
        _build_local_genome "$GENOME"
        validate_genome_reference_integrity "$LOCAL_BT2_INDEX" "$LOCAL_CHROM_SIZES"
        if $GENOME_IS_CUSTOM; then
            resolve_effective_genome_size "$GENOME" "$LOCAL_CHROM_SIZES"
            write_custom_genome_registry "$GENOME"
            write_reference_snapshot "$GENOME"
            _release_custom_genome_build_lock
        fi
        ok "Explicit local genome sequence assets are ready for $GENOME."
        return 0
    fi

    local required_assets=(
        "${GENOME}/fasta.fasta"
        "${GENOME}/bowtie2_index.bowtie2_index"
        "${GENOME}/fasta.chrom_sizes"
    )

    label "Checking required Refgenie sequence assets for $GENOME..."
    local missing_required=false
    local asset
    for asset in "${required_assets[@]}"; do
        if asset_ok "$asset"; then
            ok "$asset"
        else
            warn "$asset missing."
            missing_required=true
        fi
    done

    if $missing_required; then
        blank
        warn "One or more required sequence assets are missing for $GENOME."
        label "Installing missing Refgenie assets now because the run has been confirmed."
        blank
        _pull_genome_from_refgenie "$GENOME"
    else
        blank
        ok "All required Refgenie sequence assets already exist for $GENOME."
    fi

    blank
    label "Verifying required Refgenie sequence assets..."
    for asset in "${required_assets[@]}"; do
        asset_ok "$asset" || die "$asset is still missing after installation attempt."
        ok "$asset"
    done

    ok "Refgenie sequence assets ready for $GENOME."
}

# ─────────────────────────────────────────────────────────────
# _pull_genome_from_refgenie GENOME
# Downloads pre-built assets from the standard Refgenie server
# After pulling, calls _repair_refgenie_registration to ensure all
# assets are properly registered in the config yaml — refgenie pull
# sometimes downloads and extracts files without updating the config.
# ─────────────────────────────────────────────────────────────
_pull_genome_from_refgenie() {
    local genome="$1"

    # Refgenie is deliberately limited to sequence/alignment assets here.
    # Gene models and peak annotation are handled downstream by an
    # assembly-matched Bioconductor TxDb plus a species-matched OrgDb.
    pull_refgenie_asset_once "${genome}/fasta"         "${genome}/fasta.fasta"
    pull_refgenie_asset_once "${genome}/fasta"         "${genome}/fasta.chrom_sizes"
    pull_refgenie_asset_once "${genome}/bowtie2_index" "${genome}/bowtie2_index.bowtie2_index"

    # Repair config registration if files downloaded but seek keys were not written.
    _repair_refgenie_registration "$genome"
}

# ─────────────────────────────────────────────────────────────
# _repair_refgenie_registration GENOME
# Checks whether each required asset seek path resolves correctly.
# If not, uses the alias directory (which always has genome-named symlinks)
# to find the actual data path and re-registers assets.
# This handles a known refgenie pull bug where assets are extracted
# but not written to the config yaml.
# ─────────────────────────────────────────────────────────────
_repair_refgenie_registration() {
    local genome="$1"
    local alias_dir="$HOME/refgenie/alias/${genome}"
    local data_dir
    data_dir="$(dirname "$REFGENIE_CONFIG")/data"

    label "Verifying Refgenie config registration for $genome..."

    # ── fasta + chrom_sizes ───────────────────────────────────────────────────
    if ! asset_ok "${genome}/fasta.fasta" || ! asset_ok "${genome}/fasta.chrom_sizes"; then
        # The alias directory always contains genome-named symlinks regardless
        # of how the underlying data files are named (digest vs genome name).
        local fasta_alias_dir="${alias_dir}/fasta/default"
        if [[ -d "$fasta_alias_dir" ]]; then
            # Find the real data directory by resolving the alias directory's parent
            local real_fasta_dir
            real_fasta_dir="$(realpath "$fasta_alias_dir" 2>/dev/null || readlink -f "$fasta_alias_dir" 2>/dev/null || echo "")"

            # If realpath fails (alias is not a symlink but a real dir), use it directly
            if [[ -z "$real_fasta_dir" ]] || [[ "$real_fasta_dir" == "$fasta_alias_dir" ]]; then
                real_fasta_dir="$fasta_alias_dir"
            fi

            # Find the genome-named or digest-named fasta file
            local fa_file=""
            local chrom_file=""
            local fai_file=""
            fa_file="$(ls "$fasta_alias_dir/"*.fa 2>/dev/null | head -1 || true)"
            [[ -n "$fa_file" ]] && fa_file="$(basename "$fa_file")"
            chrom_file="$(ls "$fasta_alias_dir/"*.chrom.sizes 2>/dev/null | head -1 || true)"
            [[ -n "$chrom_file" ]] && chrom_file="$(basename "$chrom_file")"
            fai_file="$(ls "$fasta_alias_dir/"*.fa.fai 2>/dev/null | head -1 || true)"
            [[ -n "$fai_file" ]] && fai_file="$(basename "$fai_file")"

            if [[ -n "$fa_file" && -n "$chrom_file" ]]; then
                # Get path relative to refgenie genome_folder
                local rel_path
                rel_path="data/$(realpath --relative-to="$data_dir" "$fasta_alias_dir" 2>/dev/null || \
                    python3 -c "import os; print(os.path.relpath('$fasta_alias_dir', '$data_dir'))" 2>/dev/null || \
                    echo "$(basename "$(dirname "$(dirname "$fasta_alias_dir")")")/fasta/default")"

                local seek_json="{\"fasta\": \"${fa_file}\", \"chrom_sizes\": \"${chrom_file}\""
                [[ -n "$fai_file" ]] && seek_json+=", \"fai\": \"${fai_file}\""
                seek_json+="}"

                in_env refgenie add -c "$REFGENIE_CONFIG" \
                    -p "$rel_path" \
                    -s "$seek_json" \
                    --force \
                    "${genome}/fasta:default" \
                    && ok "Re-registered: ${genome}/fasta.fasta and ${genome}/fasta.chrom_sizes" \
                    || warn "Could not re-register ${genome}/fasta"
            else
                warn "Could not find fasta files in $fasta_alias_dir — skipping fasta repair."
            fi
        else
            warn "Alias directory not found: $fasta_alias_dir — skipping fasta repair."
        fi
    fi

    # ── bowtie2_index ─────────────────────────────────────────────────────────
    if ! asset_ok "${genome}/bowtie2_index.bowtie2_index"; then
        local bt2_alias_dir="${alias_dir}/bowtie2_index/default"
        if [[ -d "$bt2_alias_dir" ]]; then
            # Find the index prefix — use ls to avoid glob expansion issues
            local bt2_prefix=""
            local bt2_file
            bt2_file="$(ls "$bt2_alias_dir/"*.1.bt2 2>/dev/null | head -1 || ls "$bt2_alias_dir/"*.1.bt2l 2>/dev/null | head -1 || true)"
            if [[ -n "$bt2_file" ]]; then
                bt2_prefix="$(basename "$bt2_file" | sed 's/\.1\.bt2l\?$//')"
            fi

            if [[ -n "$bt2_prefix" ]]; then
                local rel_path
                rel_path="data/$(realpath --relative-to="$data_dir" "$bt2_alias_dir" 2>/dev/null || \
                    python3 -c "import os; print(os.path.relpath('$bt2_alias_dir', '$data_dir'))" 2>/dev/null || \
                    echo "$(basename "$(dirname "$(dirname "$bt2_alias_dir")")")/bowtie2_index/default")"

                in_env refgenie add -c "$REFGENIE_CONFIG" \
                    -p "$rel_path" \
                    -s "{\"bowtie2_index\": \"${bt2_prefix}\"}" \
                    --force \
                    "${genome}/bowtie2_index:default" \
                    && ok "Re-registered: ${genome}/bowtie2_index.bowtie2_index" \
                    || warn "Could not re-register ${genome}/bowtie2_index"
            else
                warn "Could not find bowtie2 index files in $bt2_alias_dir — skipping bowtie2_index repair."
            fi
        fi
    fi

    ok "Registration repair complete for $genome."
}

# _fetch_sequence_asset SOURCE DEST_FILE LABEL
# Resolves SOURCE as either an http(s) URL or an existing local file and
# lands the (decompressed, if gzipped) result at DEST_FILE. Shared by both
# the FASTA and chrom-sizes steps below so a user-defined genome profile's user-supplied
# path/URL is handled identically to a built-in UCSC one.
#
# Compression is detected from the downloaded/copied bytes themselves (the
# gzip magic number, 1f 8b) rather than from a ".gz" suffix on the URL or
# filename. A suffix check breaks on perfectly ordinary URLs like
# "assembly.fa.gz?download=1" or a redirect/CDN link with no .gz in it at
# all -- either would have been treated as already-decompressed and fed
# straight into bowtie2-build/samtools as gzip bytes.
_fetch_sequence_asset() {
    local source="$1"
    local dest="$2"
    local label_text="$3"
    local raw_tmp="${dest}.download"

    rm -f "$raw_tmp"
    if [[ "$source" == http://* || "$source" == https://* ]]; then
        label "Downloading $label_text from $source ..."
        wget -q --show-progress -O "$raw_tmp" "$source" \
            || die "Failed to download $label_text from $source"
    elif [[ -f "$source" ]]; then
        label "Using local $label_text: $source"
        cp "$source" "$raw_tmp" || die "Failed to copy $source"
    else
        die "$label_text source is neither a URL nor an existing local file: $source"
    fi

    local magic
    magic="$(head -c2 "$raw_tmp" 2>/dev/null | od -An -tx1 | tr -d ' \n')"
    if [[ "$magic" == "1f8b" ]]; then
        gunzip -c "$raw_tmp" > "$dest" || die "Failed to decompress $label_text"
        rm -f "$raw_tmp"
    else
        mv -f "$raw_tmp" "$dest"
    fi
}

# ─────────────────────────────────────────────────────────────
# _build_local_genome GENOME
# Builds either a built-in UCSC-backed local reference or a custom
# hash-addressed immutable reference.
# ─────────────────────────────────────────────────────────────

samtools_for_reference_build() {
    local bin="$HOME/miniconda3/envs/${ENV_NAME}/bin/samtools"
    [[ -x "$bin" ]] || bin="$(command -v samtools 2>/dev/null || true)"
    [[ -n "$bin" && -x "$bin" ]] || die "samtools not found. Run install.sh first."
    printf '%s\n' "$bin"
}

_build_builtin_local_genome() {
    local genome="$1"
    local work_dir="$HOME/pepatac_genomes/${genome}"
    mkdir -p "$work_dir"
    local fasta_source="${UCSC_FASTA_URL[$genome]:-}"
    local chromsizes_source="${UCSC_CHROMSIZES_URL[$genome]:-}"
    [[ -n "$fasta_source" ]] || die "No FASTA source defined for $genome."

    local fasta_fa="$work_dir/${genome}.fa"
    local chromsizes="$work_dir/${genome}.chrom.sizes"
    local bt2_prefix="$work_dir/bowtie2/${genome}"
    [[ -f "$fasta_fa" ]] || _fetch_sequence_asset "$fasta_source" "$fasta_fa" "$genome FASTA"

    local samtools_bin
    samtools_bin="$(samtools_for_reference_build)"
    [[ -s "${fasta_fa}.fai" ]] || "$samtools_bin" faidx "$fasta_fa"

    if [[ ! -f "$chromsizes" ]]; then
        if [[ -n "$chromsizes_source" ]]; then
            _fetch_sequence_asset "$chromsizes_source" "$chromsizes" "$genome chrom sizes"
        else
            cut -f1,2 "${fasta_fa}.fai" > "$chromsizes"
        fi
    fi

    mkdir -p "$(dirname "$bt2_prefix")"
    if ! (bowtie2_index_files "$bt2_prefix" >/dev/null 2>&1); then
        rm -f "${bt2_prefix}".*.bt2 "${bt2_prefix}".*.bt2l
        in_env bowtie2-build --threads "$THREADS" "$fasta_fa" "$bt2_prefix" \
            || die "bowtie2-build failed for $genome"
    fi

    LOCAL_CUSTOM_FASTA="$fasta_fa"
    LOCAL_BT2_INDEX="$bt2_prefix"
    LOCAL_CHROM_SIZES="$chromsizes"
    GENOME_MODE="local_explicit"
    validate_fai_matches_chrom_sizes "$fasta_fa" "$chromsizes"
    validate_bowtie2_matches_fasta "$bt2_prefix" "$fasta_fa"
}

_build_custom_reference_staged() {
    local genome="$1"
    local assembly="${CUSTOM_UCSC_ASSEMBLY[$genome]:-$genome}"
    local fasta_source="${CUSTOM_FASTA_URL[$genome]:-}" chrom_source="${CUSTOM_CHROMSIZES_URL[$genome]:-}"
    local expected_fasta expected_chrom
    expected_fasta="$(ucsc_profile_fasta_url "$assembly")"
    expected_chrom="$(ucsc_profile_chrom_sizes_url "$assembly")"
    local txdb_pkg="${CUSTOM_TXDB_PKG[$genome]:-}" orgdb_pkg="${CUSTOM_ORGDB_PKG[$genome]:-}"
    [[ "$fasta_source" == "$expected_fasta" && "$chrom_source" == "$expected_chrom" ]] \
        || die "Validated profiles must use the official UCSC FASTA and chrom.sizes paths for '$assembly'."$'\n'"       FASTA: $expected_fasta"$'\n'"       Chrom sizes: $expected_chrom"
    [[ -n "$txdb_pkg" && -n "$orgdb_pkg" ]] || die "Both TxDb and OrgDb package names are required for '$genome'."
    [[ "$txdb_pkg" == *".UCSC.${assembly}."* ]] \
        || die "The TxDb package name must identify the selected UCSC assembly '$assembly' (expected '.UCSC.${assembly}.' in the package name)."
    ensure_profile_annotation_packages "$txdb_pkg" "$orgdb_pkg"
    mkdir -p "$CUSTOM_GENOME_STAGING_ROOT/$RUN_ID"; ACTIVE_CUSTOM_STAGE="$CUSTOM_GENOME_STAGING_ROOT/$RUN_ID/${genome}_$$"; rm -rf "$ACTIVE_CUSTOM_STAGE"
    mkdir -p "$ACTIVE_CUSTOM_STAGE/bowtie2" "$ACTIVE_CUSTOM_STAGE/annotation_build"
    local stage_fasta="$ACTIVE_CUSTOM_STAGE/genome.fa" stage_chrom="$ACTIVE_CUSTOM_STAGE/genome.chrom.sizes" stage_bt2="$ACTIVE_CUSTOM_STAGE/bowtie2/${genome}" samtools_bin
    samtools_bin="$(samtools_for_reference_build)"
    _fetch_sequence_asset "$fasta_source" "$stage_fasta" "$assembly UCSC FASTA"
    "$samtools_bin" faidx "$stage_fasta" || die "samtools faidx failed on staged FASTA"
    _fetch_sequence_asset "$chrom_source" "$stage_chrom" "$assembly UCSC chrom.sizes"
    validate_fai_matches_chrom_sizes "$stage_fasta" "$stage_chrom"
    in_env bowtie2-build --threads "$THREADS" "$stage_fasta" "$stage_bt2" || die "bowtie2-build failed for genome profile '$genome'"
    validate_bowtie2_matches_fasta "$stage_bt2" "$stage_fasta"
    build_profile_annotation_assets "$assembly" "$txdb_pkg" "$orgdb_pkg" "$stage_chrom" "$ACTIVE_CUSTOM_STAGE/annotation_build"
    validate_frozen_annotation_profile "$assembly" \
        "$ACTIVE_CUSTOM_STAGE/annotation_build/txdb.sqlite" \
        "$ACTIVE_CUSTOM_STAGE/annotation_build/orgdb.sqlite" \
        "$stage_chrom"
    local staged_bed staged_label staged_allow_empty
    for staged_bed in tss.bed features.bed genes.bed exons.bed introns.bed promoters.bed; do
        staged_label="${staged_bed%.bed}"; staged_allow_empty=false
        [[ "$staged_bed" == "introns.bed" ]] && staged_allow_empty=true
        validate_profile_bed "$ACTIVE_CUSTOM_STAGE/annotation_build/$staged_bed" "$stage_chrom" "$staged_label" "$staged_allow_empty"
    done
    local build_conf="$ACTIVE_CUSTOM_STAGE/annotation_build/profile_build.conf" line key val
    [[ -s "$build_conf" ]] || die "Annotation profile builder did not produce profile_build.conf"
    local -A pb=(); while IFS= read -r line || [[ -n "$line" ]]; do [[ -z "$line" || "$line" == \#* ]] && continue; key="${line%%=*}"; val="${line#*=}"; pb["$key"]="$val"; done < "$build_conf"
    local fasta_sha annotation_manifest annotation_fingerprint
    fasta_sha="$(require_sha256 "$stage_fasta" "staged FASTA")"; annotation_manifest="$(mktemp)"
    local file
    for file in txdb.sqlite orgdb.sqlite tss.bed features.bed genes.bed exons.bed introns.bed promoters.bed profile_build.conf; do printf '%s  %s\n' "$(require_sha256 "$ACTIVE_CUSTOM_STAGE/annotation_build/$file" "$file")" "$file"; done | sort -k2,2 > "$annotation_manifest"
    annotation_fingerprint="$(require_sha256 "$annotation_manifest" "annotation fingerprint manifest")"; rm -f "$annotation_manifest"
    mkdir -p "$ACTIVE_CUSTOM_STAGE/annotations"; mv "$ACTIVE_CUSTOM_STAGE/annotation_build" "$ACTIVE_CUSTOM_STAGE/annotations/$annotation_fingerprint"
    local target_dir target_profile
    target_dir="$(custom_genome_assembly_dir "$genome" "$fasta_sha")"
    target_profile="$target_dir/annotations/$annotation_fingerprint"
    mkdir -p "$(dirname "$target_dir")"
    if [[ ! -e "$target_dir" ]]; then mv "$ACTIVE_CUSTOM_STAGE" "$target_dir"; ACTIVE_CUSTOM_STAGE=""
    else
        verify_file_sha256 "$target_dir/genome.fa" "$fasta_sha" "existing immutable FASTA"
        validate_fai_matches_chrom_sizes "$target_dir/genome.fa" "$target_dir/genome.chrom.sizes"
        validate_bowtie2_matches_fasta "$target_dir/bowtie2/${genome}" "$target_dir/genome.fa"
        mkdir -p "$target_dir/annotations"
        if [[ -e "$target_profile" ]]; then
            for file in txdb.sqlite orgdb.sqlite tss.bed features.bed genes.bed exons.bed introns.bed promoters.bed profile_build.conf; do cmp -s "$ACTIVE_CUSTOM_STAGE/annotations/$annotation_fingerprint/$file" "$target_profile/$file" || die "Existing immutable annotation fingerprint has different bytes: $target_profile/$file"; done
        else mv "$ACTIVE_CUSTOM_STAGE/annotations/$annotation_fingerprint" "$target_profile"; fi
        cleanup_custom_reference_staging
    fi
    LOCAL_CUSTOM_FASTA="$target_dir/genome.fa"; LOCAL_CUSTOM_FASTA_SHA256="$fasta_sha"
    LOCAL_CHROM_SIZES="$target_dir/genome.chrom.sizes"; LOCAL_CHROM_SIZES_SHA256="$(require_sha256 "$LOCAL_CHROM_SIZES" "published chrom.sizes")"
    LOCAL_BT2_INDEX="$target_dir/bowtie2/${genome}"; LOCAL_BT2_INDEX_SHA256="$(bowtie2_index_sha256 "$LOCAL_BT2_INDEX")"
    CUSTOM_ASSEMBLY_DIR="$target_dir"; CUSTOM_PROFILE_DIR["$genome"]="$target_profile"; CUSTOM_ANNOTATION_FINGERPRINT["$genome"]="$annotation_fingerprint"
    CUSTOM_TXDB_PKG["$genome"]="$txdb_pkg"; CUSTOM_TXDB_VERSION["$genome"]="${pb[CUSTOM_TXDB_PACKAGE_VERSION]:-}"
    CUSTOM_ORGDB_PKG["$genome"]="$orgdb_pkg"; CUSTOM_ORGDB_VERSION["$genome"]="${pb[CUSTOM_ORGDB_PACKAGE_VERSION]:-}"
    local key2 file2
    for key2 in TXDB_SQLITE ORGDB_SQLITE TSS_BED FEATURE_BED GENES_BED EXONS_BED INTRONS_BED PROMOTERS_BED PROFILE_BUILD_CONF; do
        case "$key2" in FEATURE_BED) file2="features.bed";; TXDB_SQLITE) file2="txdb.sqlite";; ORGDB_SQLITE) file2="orgdb.sqlite";; PROFILE_BUILD_CONF) file2="profile_build.conf";; *) file2="${key2,,}"; file2="${file2%_bed}.bed";; esac
        local path_var="CUSTOM_${key2}[$genome]" sha_var="CUSTOM_${key2}_SHA256[$genome]" full="$target_profile/$file2"
        printf -v "$path_var" '%s' "$full"; printf -v "$sha_var" '%s' "$(require_sha256 "$full" "$file2")"
    done
    GENOME_MODE="local_explicit"; verify_custom_profile_assets "$genome"; ok "Published immutable validated genome profile: $target_dir"
}

_build_local_genome() {
    local genome="$1"
    if $GENOME_IS_CUSTOM; then
        _build_custom_reference_staged "$genome"
    else
        _build_builtin_local_genome "$genome"
    fi
}

# ─────────────────────────────────────────────────────────────
# collect_custom_genome
# Interactive Step-3 branch for an assembly outside the five built-ins.
# A user-defined profile must identify a UCSC assembly plus matching,
# installable TxDb and OrgDb packages. KEGG remains optional. No arbitrary
# GTF/GFF is accepted and no runtime TxDb is fabricated from user annotation.
#
# On success, sets the profile globals used by the staged immutable builder.
# The FASTA/chrom.sizes download, package installation, compatibility checks,
# Bowtie2 build, frozen database creation, and publication all occur only
# after the final run confirmation.
# ─────────────────────────────────────────────────────────────
collect_custom_genome() {
    blank
    echo -e "  ${BOLD}Validated user-defined genome profile${RESET}"
    echo -e "  ${DIM}This mode requires a UCSC assembly plus matching TxDb and OrgDb packages.${RESET}"
    echo -e "  ${DIM}Arbitrary FASTA/GTF combinations are not accepted.${RESET}"
    echo -e "  ${DIM}No files or packages are changed until after the final run confirmation.${RESET}"
    blank
    local CG_NAME
    while true; do
        echo "  UCSC assembly name (for example galGal6):"
        read -r -p "  > " CG_NAME
        CG_NAME="${CG_NAME%$'\r'}"
        if ! is_valid_custom_genome_name "$CG_NAME"; then err "Use only letters, numbers, and underscores (no spaces)."; continue; fi
        if is_supported_genome "$CG_NAME"; then err "'$CG_NAME' is built in. Choose it from the main list."; continue; fi
        break
    done
    local found_registry=false
    CUSTOM_SELECTED_REGISTRY_SHA256="__ABSENT__"
    _acquire_custom_genome_build_lock "$CG_NAME"
    if load_custom_genome_registry "$CG_NAME"; then
        found_registry=true
        CUSTOM_SELECTED_REGISTRY_SHA256="$CUSTOM_REGISTRY_SHA256"
    fi
    _release_custom_genome_build_lock
    if $found_registry; then
        blank; ok "Found an existing validated profile for '$CG_NAME':"
        echo -e "    ${DIM}UCSC assembly:${RESET}  ${CUSTOM_REG_UCSC_ASSEMBLY:-<missing>}"
        echo -e "    ${DIM}FASTA SHA-256:${RESET}  ${CUSTOM_REG_FASTA_SHA256:-<missing>}"
        echo -e "    ${DIM}TxDb:${RESET}           ${CUSTOM_REG_TXDB_PKG:-<missing>} ${CUSTOM_REG_TXDB_VERSION:-}"
        echo -e "    ${DIM}OrgDb:${RESET}          ${CUSTOM_REG_ORGDB_PKG:-<missing>} ${CUSTOM_REG_ORGDB_VERSION:-}"
        echo -e "    ${DIM}KEGG:${RESET}           ${CUSTOM_REG_KEGG_ORG:-<not provided>}"
        blank
        echo "  Reuse this registered profile? [Y/n]"
        read -r -p "  > " CG_REUSE
        if [[ "${CG_REUSE,,}" != "n" ]]; then
            [[ "$CUSTOM_REG_SCHEMA" -ge 3 ]] || die "This is an older GTF-based registration. Re-register it once as a validated UCSC/TxDb/OrgDb profile."
            GENOME="$CG_NAME"; GENOME_IS_CUSTOM=true; GENOME_MODE="local_explicit"; GENOME_BUILD_METHOD="local_build"
            CUSTOM_FASTA_URL["$CG_NAME"]="$CUSTOM_REG_FASTA_SOURCE"; CUSTOM_CHROMSIZES_URL["$CG_NAME"]="$CUSTOM_REG_CHROM_SOURCE"; CUSTOM_UCSC_ASSEMBLY["$CG_NAME"]="$CUSTOM_REG_UCSC_ASSEMBLY"
            CUSTOM_TXDB_PKG["$CG_NAME"]="$CUSTOM_REG_TXDB_PKG"; CUSTOM_TXDB_VERSION["$CG_NAME"]="$CUSTOM_REG_TXDB_VERSION"; CUSTOM_ORGDB_PKG["$CG_NAME"]="$CUSTOM_REG_ORGDB_PKG"; CUSTOM_ORGDB_VERSION["$CG_NAME"]="$CUSTOM_REG_ORGDB_VERSION"; CUSTOM_KEGG_ORG["$CG_NAME"]="$CUSTOM_REG_KEGG_ORG"; CUSTOM_GENOME_SIZE["$CG_NAME"]="$CUSTOM_REG_GENOME_SIZE"
            ok "Selected registered validated genome profile: $CG_NAME"; return 0
        fi
        warn "A new version will be staged after final confirmation. Existing immutable versions remain untouched."
    fi
    local CG_FASTA_SOURCE CG_CHROM_SOURCE
    CG_FASTA_SOURCE="$(ucsc_profile_fasta_url "$CG_NAME")"
    CG_CHROM_SOURCE="$(ucsc_profile_chrom_sizes_url "$CG_NAME")"
    echo -e "  ${DIM}Official UCSC FASTA:       $CG_FASTA_SOURCE${RESET}"
    echo -e "  ${DIM}Official UCSC chrom.sizes: $CG_CHROM_SOURCE${RESET}"
    echo -e "  ${DIM}These official UCSC paths are required; arbitrary FASTA/GTF inputs are not accepted.${RESET}"
    local CG_TXDB CG_ORGDB CG_KEGG="" CG_GENOME_SIZE=""
    while true; do
        echo "  Matching TxDb package (required), e.g. TxDb.Hsapiens.UCSC.${CG_NAME}.knownGene:"
        read -r -p "  > " CG_TXDB
        CG_TXDB="${CG_TXDB%$'\r'}"
        if [[ "$CG_TXDB" =~ ^TxDb\.[A-Za-z0-9_.]+$ && "$CG_TXDB" == *".UCSC.${CG_NAME}."* ]]; then break; fi
        err "Enter a UCSC TxDb package for this exact assembly, containing '.UCSC.${CG_NAME}.', such as TxDb.Ggallus.UCSC.galGal6.refGene."
    done
    while true; do
        echo "  Matching OrgDb package (required), e.g. org.Hs.eg.db:"
        read -r -p "  > " CG_ORGDB
        CG_ORGDB="${CG_ORGDB%$'\r'}"
        [[ "$CG_ORGDB" =~ ^org\.[A-Za-z0-9_.]+\.db$ ]] && break
        err "Enter an OrgDb package name such as org.Gg.eg.db."
    done
    while true; do
        echo "  KEGG organism code [optional -- leave blank to skip], e.g. hsa:"
        read -r -p "  > " CG_KEGG
        CG_KEGG="${CG_KEGG%$'\r'}"
        [[ -z "$CG_KEGG" || "$CG_KEGG" =~ ^[a-z]{3,5}$ ]] && break
        err "Use a 3-5 letter lowercase KEGG organism code."
    done
    while true; do
        echo "  Effective genome size in bp [optional -- leave blank for a 90% estimate]:"
        read -r -p "  > " CG_GENOME_SIZE
        CG_GENOME_SIZE="${CG_GENOME_SIZE%$'\r'}"
        [[ -z "$CG_GENOME_SIZE" || ( "$CG_GENOME_SIZE" =~ ^[0-9]+$ && "$CG_GENOME_SIZE" -gt 0 ) ]] && break
        err "Enter a positive whole number or leave blank."
    done
    GENOME="$CG_NAME"; GENOME_IS_CUSTOM=true; GENOME_MODE="refgenie"; GENOME_BUILD_METHOD="local_build"
    CUSTOM_UCSC_ASSEMBLY["$CG_NAME"]="$CG_NAME"; CUSTOM_FASTA_URL["$CG_NAME"]="$CG_FASTA_SOURCE"; CUSTOM_CHROMSIZES_URL["$CG_NAME"]="$CG_CHROM_SOURCE"
    CUSTOM_TXDB_PKG["$CG_NAME"]="$CG_TXDB"; CUSTOM_ORGDB_PKG["$CG_NAME"]="$CG_ORGDB"; CUSTOM_KEGG_ORG["$CG_NAME"]="$CG_KEGG"; CUSTOM_GENOME_SIZE["$CG_NAME"]="$CG_GENOME_SIZE"
    blank; ok "Genome-profile inputs collected. Package installation, validation, and reference publication begin only after final confirmation."
}

# ─────────────────────────────────────────────────────────────
# Blacklist
# ─────────────────────────────────────────────────────────────

ensure_blacklist_if_requested() {
    if ! $USE_BLACKLIST; then
        return 0
    fi

    header "Blacklist Download"

    if [[ "$BLACKLIST_SOURCE" == "standard" ]]; then
        mkdir -p "$BLACKLIST_DIR"

        local blacklist_tmp
        blacklist_tmp="${BLACKLIST_PATH}.part"

        local expected_sha256="${BLACKLIST_SHA256[$GENOME]:-}"

        if [[ -f "$BLACKLIST_PATH" ]]; then
            if gzip_file_valid "$BLACKLIST_PATH" && \
               [[ -z "$expected_sha256" || "$(sha256sum "$BLACKLIST_PATH" | awk '{print $1}')" == "$expected_sha256" ]]; then
                ok "Blacklist already downloaded and verified: $BLACKLIST_PATH"
                return 0
            else
                warn "Existing blacklist appears incomplete, corrupted, or doesn't match the pinned version. Re-downloading."
                rm -f "$BLACKLIST_PATH"
            fi
        fi

        command_exists wget || die "wget is required to download the blacklist, but wget was not found."

        rm -f "$blacklist_tmp"
        label "Downloading $GENOME blacklist from Boyle Lab (pinned commit $BLACKLIST_COMMIT)..."

        if wget -q -O "$blacklist_tmp" "${BLACKLIST_URLS[$GENOME]}"; then
            if ! gzip_file_valid "$blacklist_tmp"; then
                rm -f "$blacklist_tmp"
                die "Blacklist download finished, but the file failed gzip verification. Please rerun the script."
            fi
            if [[ -n "$expected_sha256" ]]; then
                local actual_sha256
                actual_sha256="$(sha256sum "$blacklist_tmp" | awk '{print $1}')"
                if [[ "$actual_sha256" != "$expected_sha256" ]]; then
                    rm -f "$blacklist_tmp"
                    die "Blacklist checksum mismatch for $GENOME. Expected $expected_sha256, got $actual_sha256. Refusing to use an unverified blacklist -- the pinned commit's content may have changed, or the download was corrupted/intercepted. Re-run, or if the file has legitimately changed, update BLACKLIST_COMMIT/BLACKLIST_SHA256 above after manually verifying the new content."
                fi
            fi
            mv "$blacklist_tmp" "$BLACKLIST_PATH"
            ok "Downloaded and verified (gzip + SHA256): $BLACKLIST_PATH"
        else
            rm -f "$blacklist_tmp"
            die "Blacklist download failed. Please check internet connection and rerun the script."
        fi
    elif [[ "$BLACKLIST_SOURCE" == "custom" ]]; then
        if [[ -f "$BLACKLIST_PATH" ]]; then
            ok "Custom blacklist found: $BLACKLIST_PATH"
        else
            die "Custom blacklist file not found: $BLACKLIST_PATH"
        fi
    else
        die "Internal error: unknown BLACKLIST_SOURCE='$BLACKLIST_SOURCE'"
    fi
}

# validate_custom_blacklist BLACKLIST CHROM_SIZES
# BED coordinates are 0-based half-open. Every record must be numeric,
# positive-width, on a known sequence, and within chromosome bounds.
validate_custom_blacklist() {
    local blacklist="$1" chrom_sizes="$2"
    [[ -s "$blacklist" ]] || die "Custom blacklist is missing or empty: $blacklist"
    [[ -s "$chrom_sizes" ]] || die "Cannot validate custom blacklist without chrom.sizes: $chrom_sizes"

    local report reader
    report="$(mktemp)"
    if [[ "$(head -c2 "$blacklist" 2>/dev/null | od -An -tx1 | tr -d ' \n')" == "1f8b" ]]; then
        reader=(gzip -cd "$blacklist")
    else
        reader=(cat "$blacklist")
    fi

    if ! "${reader[@]}" | awk -F'\t' -v sizes="$chrom_sizes" -v report="$report" '
        BEGIN {
            while ((getline line < sizes) > 0) {
                split(line, a, "\t")
                if (a[1] != "" && a[2] ~ /^[0-9]+$/) len[a[1]]=a[2]
            }
            close(sizes)
        }
        /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ || NF==0 {next}
        {
            if (NF < 3) {print "line " NR ": fewer than 3 BED columns" > report; bad=1; exit}
            if (!($1 in len)) {print "line " NR ": unknown chromosome " $1 > report; bad=1; exit}
            if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {print "line " NR ": non-numeric start/end" > report; bad=1; exit}
            if ($2 < 0) {print "line " NR ": start < 0" > report; bad=1; exit}
            if ($3 <= $2) {print "line " NR ": end <= start" > report; bad=1; exit}
            if ($3 > len[$1]) {print "line " NR ": end exceeds chromosome length" > report; bad=1; exit}
            matched++
        }
        END {
            if (!bad && matched==0) {print "no valid blacklist records matched the assembly" > report; exit 1}
            if (bad) exit 1
        }
    '; then
        local why
        why="$(cat "$report" 2>/dev/null || true)"
        rm -f "$report"
        die "Custom blacklist validation failed: ${why:-unknown BED error}."
    fi
    rm -f "$report"
    ok "Custom blacklist coordinates match the selected assembly."
}

# ─────────────────────────────────────────────────────────────
# Preflight
# ─────────────────────────────────────────────────────────────

run_preflight_checks() {
    header "Preflight Check"

    if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
        ok "Conda environment '$ENV_NAME' found."
    else
        die "Conda environment '$ENV_NAME' not found. Run install.sh first."
    fi

    [[ -f "$PIPELINE" ]] || die "Pipeline not found at $PIPELINE. Run install.sh first."
    ok "Pipeline found."

    [[ -f "$REFGENIE_CONFIG" ]] || die "Refgenie config not found at $REFGENIE_CONFIG. Run install.sh first."
    ok "Refgenie config found."

    export REFGENIE="$REFGENIE_CONFIG"

    local genome_index chrom_sizes

    if [[ "${GENOME_MODE:-refgenie}" == "local_explicit" ]]; then
        genome_index="$LOCAL_BT2_INDEX"
        chrom_sizes="$LOCAL_CHROM_SIZES"
        [[ -n "$genome_index" ]] || die "Local explicit Bowtie2 index path is empty."
        [[ -n "$chrom_sizes"  ]] || die "Local explicit chrom sizes path is empty."
    else
        genome_index=$(in_env refgenie seek -c "$REFGENIE_CONFIG" "${GENOME}/bowtie2_index.bowtie2_index") || die "Cannot resolve Bowtie2 index for $GENOME."
        chrom_sizes=$(in_env refgenie seek -c "$REFGENIE_CONFIG" "${GENOME}/fasta.chrom_sizes") || die "Cannot resolve chrom sizes for $GENOME."
        [[ -n "$genome_index" ]] || die "Refgenie returned an empty Bowtie2 index path for $GENOME."
        [[ -n "$chrom_sizes"  ]] || die "Refgenie returned an empty chrom sizes path for $GENOME."
        ok "Bowtie2 index resolved."
        ok "Chrom sizes resolved."
    fi

    validate_genome_reference_integrity "$genome_index" "$chrom_sizes"
    if $GENOME_IS_CUSTOM; then
        verify_file_sha256 "$LOCAL_CUSTOM_FASTA" "$LOCAL_CUSTOM_FASTA_SHA256" "custom FASTA"
        verify_file_sha256 "$LOCAL_CHROM_SIZES" "$LOCAL_CHROM_SIZES_SHA256" "custom chrom.sizes"
        verify_custom_profile_assets "$GENOME"
    fi

    local i r1 r2
    for i in "${!SAMPLE_NAMES[@]}"; do
        r1="${R1_FILES[$i]}"
        r2="${R2_FILES_ARR[$i]}"

        [[ -f "$r1" ]] || die "Missing R1 FASTQ for sample ${SAMPLE_NAMES[$i]}: $r1"

        if $PAIRED && [[ -n "$r2" ]]; then
            [[ -f "$r2" ]] || die "Missing R2 FASTQ for sample ${SAMPLE_NAMES[$i]}: $r2"
        fi
    done
    ok "FASTQ files verified."

    mkdir -p "$OUTPUT_DIR"
    [[ -d "$OUTPUT_DIR" ]] || die "Could not create output directory: $OUTPUT_DIR"
    ok "Output directory ready: $OUTPUT_DIR"

    if $USE_BLACKLIST; then
        [[ -f "$BLACKLIST_PATH" ]] || die "Blacklist requested but file not found: $BLACKLIST_PATH"
        if [[ "$BLACKLIST_SOURCE" == "custom" ]]; then
            validate_custom_blacklist "$BLACKLIST_PATH" "$chrom_sizes"
        fi
        ok "Blacklist ready: $BLACKLIST_PATH"
    else
        ok "Blacklist filtering not requested."
    fi


    blank
    echo -e "  ${GREEN}${BOLD}Preflight passed. Ready for takeoff.${RESET}"
}

# ─────────────────────────────────────────────────────────────
# SRA / ENA accession-based FASTQ acquisition
#
# sra_lookup()   — Step 1b: prompt for accession, query ENA metadata
#                  only (no FASTQ bytes), let user pick runs, detect
#                  paired/single from metadata.
# sra_download() — After final confirmation: download FASTQs, verify
#                  checksums, rename to _R1/_R2 convention, set INPUT_DIR.
# ─────────────────────────────────────────────────────────────

normalize_sra_accession() {
    local raw="$1"
    raw="$(echo "$raw" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')"
    echo "$raw"
}

fetch_ena_filereport() {
    local accession="$1"
    local out_tsv="$2"
    local url="${ENA_FILEREPORT_BASE}?accession=${accession}&result=read_run&fields=${ENA_FIELDS}&format=tsv"

    command_exists curl || die "curl is required to query ENA but was not found. Run install.sh first."

    label "Querying ENA for: $accession"
    if ! curl -fsSL "$url" -o "$out_tsv" 2>/dev/null; then
        die "Could not reach ENA's API. Check your internet connection and try again."
    fi

    if [[ ! -s "$out_tsv" ]] || [[ "$(wc -l < "$out_tsv")" -le 1 ]]; then
        die "No sequencing runs found at ENA for accession: $accession
        Double-check the accession number. Note this only works for public,
        already-released data — embargoed/private SRA records won't resolve here."
    fi
}

sra_sample_label() {
    local sample_title="$1" experiment_title="$2" run_accession="$3"
    local label="$sample_title"
    [[ -z "$label" || "$label" == "NA" ]] && label="$experiment_title"
    [[ -z "$label" || "$label" == "NA" ]] && label="$run_accession"
    label=$(echo "$label" | tr -c '[:alnum:]._-' '_' | sed 's/__*/_/g; s/^_//; s/_$//')
    echo "${label}_${run_accession}"
}

sra_lookup() {
    local accession_raw report_tsv

    header "Accession Lookup"
    echo -e "  Enter a public SRA study accession (SRP) or BioProject accession (PRJNA)."
    echo -e "  ${DIM}Example: SRP240350 or PRJNA599930${RESET}"
    blank

    while true; do
        echo "  Accession:"
        read -r -p "  > " accession_raw
        if [[ -z "$accession_raw" ]]; then
            err "Please enter an accession."
            continue
        fi
        SRA_ACCESSION=$(normalize_sra_accession "$accession_raw")
        break
    done

    blank
    local default_dl_base="${DEFAULT_SRA_DOWNLOAD_DIR:-$DEFAULT_OUTPUT_BASE/sra_downloads}"
    local default_dl_preview="$default_dl_base/${SRA_ACCESSION}_${RUN_ID}"
    echo -e "  Where should the downloaded FASTQ files be saved?"
    echo -e "  ${DIM}Default: $default_dl_preview${RESET}"
    echo -e "  ${DIM}Tip: add -newdefault after a path to save it as the default download folder.${RESET}"
    blank

    echo "  Download folder [press Enter for default]:"
    IFS= read -r -e -p "  > " SRA_DOWNLOAD_DIR_INPUT
    if strip_newdefault_marker SRA_DOWNLOAD_DIR_INPUT; then
        normalize_path_var SRA_DOWNLOAD_DIR_INPUT
        [[ -z "$SRA_DOWNLOAD_DIR_INPUT" ]] && SRA_DOWNLOAD_DIR_INPUT="$default_dl_base"
        save_user_default DEFAULT_SRA_DOWNLOAD_DIR "$SRA_DOWNLOAD_DIR_INPUT"
        DEFAULT_SRA_DOWNLOAD_DIR="$SRA_DOWNLOAD_DIR_INPUT"
        SRA_DOWNLOAD_DIR="$SRA_DOWNLOAD_DIR_INPUT/${SRA_ACCESSION}_${RUN_ID}"
    else
        normalize_path_var SRA_DOWNLOAD_DIR_INPUT
        if [[ -z "$SRA_DOWNLOAD_DIR_INPUT" ]]; then
            SRA_DOWNLOAD_DIR="$default_dl_preview"
        else
            SRA_DOWNLOAD_DIR="$SRA_DOWNLOAD_DIR_INPUT"
        fi
    fi

    ok "Download folder: $SRA_DOWNLOAD_DIR"
    mkdir -p "$SRA_DOWNLOAD_DIR"
    report_tsv="$SRA_DOWNLOAD_DIR/.ena_filereport.tsv"

    fetch_ena_filereport "$SRA_ACCESSION" "$report_tsv"

    # Parse TSV into parallel arrays. Column order from ENA_FIELDS:
    # run_accession, fastq_ftp, fastq_md5, fastq_bytes, library_layout,
    # library_strategy, instrument_platform, sample_title, experiment_title,
    # read_count, base_count
    local line run ftp md5 bytes layout strategy platform stitle etitle rcount bcount
    SRA_RUN_IDS=()
    SRA_RUN_LAYOUT=()
    SRA_RUN_FTP=()
    SRA_RUN_MD5=()
    SRA_RUN_BYTES=()
    SRA_RUN_LABEL=()
    SRA_RUN_READ_LEN=()
    while IFS=$'\t' read -r run ftp md5 bytes layout strategy platform stitle etitle rcount bcount; do
        [[ "$run" == "run_accession" ]] && continue
        [[ -z "$run" ]] && continue
        SRA_RUN_IDS+=("$run")
        SRA_RUN_LAYOUT+=("$layout")
        SRA_RUN_FTP+=("$ftp")
        SRA_RUN_MD5+=("$md5")
        SRA_RUN_BYTES+=("$bytes")
        SRA_RUN_LABEL+=("$(sra_sample_label "$stitle" "$etitle" "$run")")
        local this_len=""
        if [[ "$rcount" =~ ^[0-9]+$ && "$bcount" =~ ^[0-9]+$ && "$rcount" -gt 0 ]]; then
            if [[ "$layout" == "PAIRED" ]]; then
                this_len=$(( bcount / rcount / 2 ))
            else
                this_len=$(( bcount / rcount ))
            fi
        fi
        SRA_RUN_READ_LEN+=("$this_len")
    done < "$report_tsv"

    if [[ ${#SRA_RUN_IDS[@]} -eq 0 ]]; then
        die "ENA returned a report but no usable runs could be parsed from it: $report_tsv"
    fi

    blank
    echo -e "  ${BOLD}Found ${#SRA_RUN_IDS[@]} run(s) for $SRA_ACCESSION:${RESET}"
    blank
    local i total_bytes=0 n_paired=0 n_single=0
    SRA_RUN_TOTAL_BYTES=()
    for i in "${!SRA_RUN_IDS[@]}"; do
        local n_files run_bytes
        n_files=$(awk -F';' '{print NF}' <<< "${SRA_RUN_FTP[$i]}")
        run_bytes=$(awk -F';' '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}' <<< "${SRA_RUN_BYTES[$i]}")
        SRA_RUN_TOTAL_BYTES+=("$run_bytes")
        printf "    ${CYAN}%2d. %-15s${RESET}  %-10s  %2s file(s)  %10s  %s\n" \
            "$((i+1))" "${SRA_RUN_IDS[$i]}" "${SRA_RUN_LAYOUT[$i]:-NA}" "$n_files" \
            "$(human_bytes "$run_bytes")" "${SRA_RUN_LABEL[$i]}"
        total_bytes=$((total_bytes + run_bytes))
        if [[ "${SRA_RUN_LAYOUT[$i]}" == "PAIRED" ]]; then
            n_paired=$((n_paired+1))
        else
            n_single=$((n_single+1))
        fi
    done
    blank
    echo -e "  ${BOLD}Total size (all runs):${RESET}  $(human_bytes "$total_bytes")"
    echo -e "  ${BOLD}Library layout:${RESET}         $n_paired paired-end run(s), $n_single single-end run(s)"
    blank

    echo -e "  ${BOLD}Which runs do you want to download?${RESET}"
    echo -e "  ${DIM}Examples: 'all' · '3' · '1,3,5' · '1-4' · '1-3,6'${RESET}"
    blank

    SRA_SELECTED_IDX=()
    while true; do
        echo "  Runs to download [default: all]:"
        read -r -p "  > " RUN_SELECTION
        RUN_SELECTION="${RUN_SELECTION:-all}"
        RUN_SELECTION="$(echo "$RUN_SELECTION" | tr -d '[:space:]')"

        if [[ "${RUN_SELECTION,,}" == "all" ]]; then
            SRA_SELECTED_IDX=()
            for i in "${!SRA_RUN_IDS[@]}"; do
                SRA_SELECTED_IDX+=("$i")
            done
            break
        fi

        SRA_SELECTED_IDX=()
        local token range_start range_end n bad_token=""
        IFS=',' read -r -a tokens <<< "$RUN_SELECTION"
        for token in "${tokens[@]}"; do
            if [[ "$token" =~ ^[0-9]+$ ]]; then
                n="$token"
                if [[ "$n" -lt 1 || "$n" -gt "${#SRA_RUN_IDS[@]}" ]]; then
                    bad_token="$token"; break
                fi
                SRA_SELECTED_IDX+=("$((n-1))")
            elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
                range_start="${BASH_REMATCH[1]}"
                range_end="${BASH_REMATCH[2]}"
                if [[ "$range_start" -lt 1 || "$range_end" -gt "${#SRA_RUN_IDS[@]}" || "$range_start" -gt "$range_end" ]]; then
                    bad_token="$token"; break
                fi
                for ((n = range_start; n <= range_end; n++)); do
                    SRA_SELECTED_IDX+=("$((n-1))")
                done
            else
                bad_token="$token"; break
            fi
        done

        if [[ -n "$bad_token" ]]; then
            err "Couldn't understand '$bad_token' — use a number 1-${#SRA_RUN_IDS[@]}, a range like 1-4, or 'all'."
            continue
        fi

        # Dedupe while preserving order.
        declare -A seen_idx=()
        declare -a deduped_idx=()
        for i in "${SRA_SELECTED_IDX[@]}"; do
            if [[ -z "${seen_idx[$i]:-}" ]]; then
                deduped_idx+=("$i")
                seen_idx[$i]=1
            fi
        done
        SRA_SELECTED_IDX=("${deduped_idx[@]}")
        unset seen_idx deduped_idx

        if [[ ${#SRA_SELECTED_IDX[@]} -eq 0 ]]; then
            err "No runs selected — enter at least one number, a range, or 'all'."
            continue
        fi

        # Reject mixed paired/single-end selections here, before they're
        # accepted. Sample discovery only runs in one library mode per
        # invocation, and a paired-end run's *_1.fastq.gz is indistinguishable
        # from a single-end sample's own naming at that point -- there is no
        # reliable way to "skip the wrong-layout runs" after the fact, only a
        # way to silently misprocess some of them. Make the user pick one
        # layout per run of this script instead.
        local sel_n_paired=0 sel_n_single=0
        for i in "${SRA_SELECTED_IDX[@]}"; do
            if [[ "${SRA_RUN_LAYOUT[$i]}" == "PAIRED" ]]; then
                sel_n_paired=$((sel_n_paired+1))
            else
                sel_n_single=$((sel_n_single+1))
            fi
        done
        if [[ "$sel_n_paired" -gt 0 && "$sel_n_single" -gt 0 ]]; then
            blank
            err "This selection mixes $sel_n_paired paired-end run(s) and $sel_n_single single-end run(s)."
            err "This script processes one library type per run. Select only the"
            err "paired-end runs or only the single-end runs, then run this script"
            err "again for the other layout."
            continue
        fi

        break
    done

    # Recompute totals for selected subset.
    local selected_bytes=0
    n_paired=0; n_single=0
    for i in "${SRA_SELECTED_IDX[@]}"; do
        selected_bytes=$((selected_bytes + SRA_RUN_TOTAL_BYTES[i]))
        if [[ "${SRA_RUN_LAYOUT[$i]}" == "PAIRED" ]]; then
            n_paired=$((n_paired+1))
        else
            n_single=$((n_single+1))
        fi
    done

    blank
    echo -e "  ${BOLD}Selected ${#SRA_SELECTED_IDX[@]} of ${#SRA_RUN_IDS[@]} run(s):${RESET}"
    for i in "${SRA_SELECTED_IDX[@]}"; do
        echo -e "    ${CYAN}${SRA_RUN_IDS[$i]}${RESET}  ${SRA_RUN_LABEL[$i]}"
    done
    blank
    echo -e "  ${BOLD}Selected download size:${RESET}  $(human_bytes "$selected_bytes")"
    SRA_TOTAL_BYTES="$selected_bytes"

    # Defensive backstop only: the selection loop above already rejects
    # mixed layouts and makes the user re-select, so this should be
    # unreachable. If it's ever hit, fail loudly rather than silently
    # continuing with a run-wide library mode that doesn't match every run.
    if [[ "$n_paired" -gt 0 && "$n_single" -gt 0 ]]; then
        die "Internal error: mixed-layout selection reached past the selection guard ($n_paired paired, $n_single single)."
    fi

    local avail_bytes
    avail_bytes=$(available_bytes_for_path "$SRA_DOWNLOAD_DIR")
    if [[ "$avail_bytes" -gt 0 && "$selected_bytes" -gt "$avail_bytes" ]]; then
        warn "Estimated download ($(human_bytes "$selected_bytes")) exceeds free space at destination ($(human_bytes "$avail_bytes"))."
    fi

    if [[ "$n_paired" -gt 0 && "$n_single" -eq 0 ]]; then
        SRA_AUTO_PAIRED=true
    elif [[ "$n_single" -gt 0 && "$n_paired" -eq 0 ]]; then
        SRA_AUTO_PAIRED=false
    else
        SRA_AUTO_PAIRED=""
    fi

    blank
    ok "Run selection ready: ${#SRA_SELECTED_IDX[@]} run(s), $(human_bytes "$selected_bytes")."
    echo -e "  ${DIM}These will be downloaded after the rest of the run settings are confirmed.${RESET}"
    blank
}

sra_download() {
    local run

    header "Downloading FASTQ Files from ENA"
    echo -e "  Accession: ${BOLD}$SRA_ACCESSION${RESET}"
    echo -e "  Download destination: ${BOLD}$SRA_DOWNLOAD_DIR${RESET}"
    echo -e "  Runs: ${#SRA_SELECTED_IDX[@]}   Size: $(human_bytes "$SRA_TOTAL_BYTES")"
    blank
    echo "  Download these ${#SRA_SELECTED_IDX[@]} run(s) now? [Y/n]"
    read -r -p "  > " CONFIRM_DL
    [[ "${CONFIRM_DL,,}" == "n" ]] && die "Aborted by user."

    command_exists wget || die "wget is required to download FASTQs but was not found. Run install.sh first."

    local i url ok_count=0 fail_count=0
    local -a SRA_SUCCESSFUL_IDX=()
    for i in "${SRA_SELECTED_IDX[@]}"; do
        run="${SRA_RUN_IDS[$i]}"
        IFS=';' read -r -a urls <<< "${SRA_RUN_FTP[$i]}"
        IFS=';' read -r -a md5s <<< "${SRA_RUN_MD5[$i]}"

        if [[ ${#urls[@]} -eq 0 || -z "${urls[0]}" ]]; then
            warn "No FASTQ URL listed for $run (often a controlled-access or SRA-only record) — skipping."
            fail_count=$((fail_count+1))
            continue
        fi

        local j dest run_ok=true
        for j in "${!urls[@]}"; do
            url="${urls[$j]}"
            [[ "$url" =~ ^https?:// || "$url" =~ ^ftp:// ]] || url="https://${url}"
            dest="$SRA_DOWNLOAD_DIR/$(basename "$url")"
            local md5sum_expected="${md5s[$j]:-}"

            if [[ -s "$dest" ]]; then
                ok "Already present: $(basename "$dest")"
            else
                label "Downloading: $(basename "$dest")  ($run, file $((j+1))/${#urls[@]})"
                if ! wget -q --show-progress -O "${dest}.part" "$url"; then
                    err "Download failed: $url"
                    rm -f "${dest}.part"
                    run_ok=false
                    continue
                fi
                mv "${dest}.part" "$dest"
            fi

            if [[ -n "$md5sum_expected" ]] && command_exists md5sum; then
                local actual_md5
                actual_md5=$(md5sum "$dest" | awk '{print $1}')
                if [[ "$actual_md5" != "$md5sum_expected" ]]; then
                    err "Checksum mismatch for $(basename "$dest") — expected $md5sum_expected, got $actual_md5"
                    rm -f "$dest"
                    run_ok=false
                else
                    ok "Checksum verified: $(basename "$dest")"
                fi
            fi
        done

        if $run_ok; then
            ok_count=$((ok_count+1))
            SRA_SUCCESSFUL_IDX+=("$i")
        else
            fail_count=$((fail_count+1))
        fi
    done

    blank
    ok "Downloaded $ok_count run(s) successfully."
    if [[ "$fail_count" -gt 0 ]]; then
        warn "$fail_count run(s) had a download or checksum problem — see above."
        echo "  Continue anyway with the runs that succeeded? [y/N]"
        read -r -p "  > " CONTINUE_PARTIAL
        [[ "${CONTINUE_PARTIAL,,}" != "y" ]] && die "Aborted by user after partial download failure."
    fi

    # Rename ENA's run-accession filenames to _R1/_R2 convention so
    # sample discovery picks them up exactly like locally-named FASTQs.
    # Only successful downloads -- a failed run's corrupt/partial file (now
    # deleted above on checksum mismatch, or never completed on download
    # failure) must never be promoted into the naming convention that
    # discover_samples() treats as ready-to-process input.
    label "Renaming downloaded files to sample-labeled R1/R2 convention..."
    for i in "${SRA_SUCCESSFUL_IDX[@]}"; do
        run="${SRA_RUN_IDS[$i]}"
        local clean_label="${SRA_RUN_LABEL[$i]}"
        if [[ -f "$SRA_DOWNLOAD_DIR/${run}_1.fastq.gz" ]]; then
            mv -n "$SRA_DOWNLOAD_DIR/${run}_1.fastq.gz" "$SRA_DOWNLOAD_DIR/${clean_label}_R1.fastq.gz" 2>/dev/null || true
        fi
        if [[ -f "$SRA_DOWNLOAD_DIR/${run}_2.fastq.gz" ]]; then
            mv -n "$SRA_DOWNLOAD_DIR/${run}_2.fastq.gz" "$SRA_DOWNLOAD_DIR/${clean_label}_R2.fastq.gz" 2>/dev/null || true
        fi
        if [[ -f "$SRA_DOWNLOAD_DIR/${run}.fastq.gz" ]]; then
            mv -n "$SRA_DOWNLOAD_DIR/${run}.fastq.gz" "$SRA_DOWNLOAD_DIR/${clean_label}.fastq.gz" 2>/dev/null || true
        fi
    done
    ok "Files renamed."

    INPUT_DIR="$SRA_DOWNLOAD_DIR"
    SRA_ACCESSION_USED="$SRA_ACCESSION"

    blank
    ok "Data source ready: $INPUT_DIR"
    blank
}

# ─────────────────────────────────────────────────────────────
# Banner
# ─────────────────────────────────────────────────────────────

# Test/library mode: define all helpers without starting the interactive runner.
if [[ "${PEPATAC_LIBRARY_ONLY:-0}" == "1" ]]; then
    return 0 2>/dev/null || exit 0
fi

clear 2>/dev/null || true
echo ""
echo -e "${BOLD}${MAGENTA}"
cat << 'BANNER'
  ██████╗ ███████╗██████╗  █████╗ ████████╗ █████╗  ██████╗
  ██╔══██╗██╔════╝██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗██╔════╝
  ██████╔╝█████╗  ██████╔╝███████║   ██║   ███████║██║
  ██╔═══╝ ██╔══╝  ██╔═══╝ ██╔══██║   ██║   ██╔══██║██║
  ██║     ███████╗██║     ██║  ██║   ██║   ██║  ██║╚██████╗
  ╚═╝     ╚══════╝╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝
BANNER
echo -e "${RESET}"
echo -e "  ${DIM}ATAC-seq Pipeline Runner  ·  Local Sequential Mode${RESET}"
echo -e "  ${DIM}Configure first → download/check → take off${RESET}"
echo ""

# ─────────────────────────────────────────────────────────────
# STEP 0 — New run or resume a previous one
# ─────────────────────────────────────────────────────────────

header "New Run or Resume"

echo -e "  1.  New run  ${DIM}(default)${RESET}"
echo -e "  2.  Resume a previous run"
echo -e "      ${DIM}Skips samples that already PASSED. Retries FAILed and NOT_RUN samples.${RESET}"
blank

echo "  Choice [1/2, default: 1]:"
read -r -p "  > " NEW_OR_RESUME_CHOICE
if [[ "${NEW_OR_RESUME_CHOICE:-1}" == "2" ]]; then
    RESUME_MODE=true

    blank
    echo -e "  Enter the OUTPUT folder of the run you want to resume."
    echo -e "  ${DIM}It should contain run_manifest.txt and a logs/ subfolder.${RESET}"
    blank

    while true; do
        echo "  Previous run's output folder:"
        IFS= read -r -e -p "  > " RESUME_DIR_INPUT
        normalize_path_var RESUME_DIR_INPUT

        if [[ -z "$RESUME_DIR_INPUT" ]]; then
            err "Please enter a path."
            continue
        fi
        if [[ ! -d "$RESUME_DIR_INPUT" ]]; then
            err "Directory not found: $RESUME_DIR_INPUT"
            continue
        fi
        if load_resume_state "$RESUME_DIR_INPUT"; then
            ok "Loaded settings from previous run: $OUTPUT_DIR"
            if $GENOME_IS_CUSTOM; then
                if ! load_reference_snapshot_into_runner false; then
                    warn "This run has no schema-3 validated genome-profile snapshot."
                    warn "Falling back to its resume-state paths; re-run from scratch before publication if this is a legacy run."
                fi
            fi
            break
        else
            err "Could not load a previous run from that folder (see above)."
        fi
    done

    blank
    echo -e "  ${BOLD}Restored settings:${RESET}"
    echo -e "    Input folder:  $INPUT_DIR"
    echo -e "    Output folder: $OUTPUT_DIR"
    echo -e "    Genome:        $GENOME$( $GENOME_IS_CUSTOM && echo " (custom)" || echo "" )"
    echo -e "    Genome mode:   ${GENOME_MODE:-refgenie}"
    echo -e "    Samples:       ${#SAMPLE_NAMES[@]}"
    blank

    # Show per-sample status from previous run.
    LOG_DIR="$OUTPUT_DIR/logs"
    echo -e "  ${BOLD}Sample status from previous run:${RESET}"
    skipped=0
    retry=0
    for sname in "${SAMPLE_NAMES[@]}"; do
        prev_status=""
        prev_status=$(read_sample_status "$sname")
        if [[ "$prev_status" == "PASS" ]]; then
            echo -e "    ${GREEN}✔${RESET}  $sname  ${DIM}(PASS — will skip)${RESET}"
            skipped=$((skipped+1))
        else
            echo -e "    ${YELLOW}↺${RESET}  $sname  ${DIM}($prev_status — will retry)${RESET}"
            retry=$((retry+1))
        fi
    done
    blank
    ok "$skipped sample(s) will be skipped. $retry sample(s) will be retried."
    blank

    echo "  Proceed with this resume? [Y/n]"
    read -r -p "  > " RESUME_CONFIRM
    [[ "${RESUME_CONFIRM,,}" == "n" ]] && die "Aborted by user."

    # Allow overriding the thread/concurrency settings from the previous run.
    blank
    AVAIL_CORES=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

    echo -e "  ${BOLD}Concurrent samples from previous run:${RESET} ${PARALLEL_SAMPLES:-1}  ${DIM}(system has $AVAIL_CORES logical cores)${RESET}"
    echo "  Change concurrent-samples count? Leave blank to keep [${PARALLEL_SAMPLES:-1}]:"
    read -r -p "  > " PARALLEL_OVERRIDE
    if [[ -n "$PARALLEL_OVERRIDE" ]]; then
        if [[ "$PARALLEL_OVERRIDE" =~ ^[0-9]+$ ]] && [[ "$PARALLEL_OVERRIDE" -ge 1 ]]; then
            PARALLEL_SAMPLES="$PARALLEL_OVERRIDE"
            ok "Concurrent samples updated to $PARALLEL_SAMPLES."
        else
            warn "Invalid input — keeping previous concurrent-samples count: ${PARALLEL_SAMPLES:-1}."
        fi
    else
        ok "Keeping concurrent samples: ${PARALLEL_SAMPLES:-1}."
    fi

    echo -e "  ${BOLD}Threads from previous run:${RESET} $THREADS"
    echo "  Change thread count? Leave blank to keep [$THREADS]:"
    read -r -p "  > " THREADS_OVERRIDE
    if [[ -n "$THREADS_OVERRIDE" ]]; then
        if [[ "$THREADS_OVERRIDE" =~ ^[0-9]+$ ]] && [[ "$THREADS_OVERRIDE" -ge 1 ]]; then
            THREADS="$THREADS_OVERRIDE"
            ok "Thread count updated to $THREADS."
        else
            warn "Invalid input — keeping previous thread count: $THREADS."
        fi
    else
        ok "Keeping thread count: $THREADS."
    fi

    REQUESTED_CORES=$(( THREADS * ${PARALLEL_SAMPLES:-1} ))
    if [[ "$REQUESTED_CORES" -gt "$AVAIL_CORES" ]]; then
        warn "$THREADS threads x ${PARALLEL_SAMPLES:-1} concurrent samples = $REQUESTED_CORES cores requested, but only $AVAIL_CORES are available."
    fi
    blank

    # In resume mode: skip all the interactive steps below and jump straight
    # to the download/preflight/processing sequence.
fi

if ! $RESUME_MODE; then

# ─────────────────────────────────────────────────────────────
# STEP 1 — Data source
# ─────────────────────────────────────────────────────────────

header "Step 1 · Data Source"

echo -e "  Where should FASTQ files come from?"
blank
echo -e "    ${CYAN}1${RESET}.  My own data  ${DIM}(point at a local folder)${RESET}"
echo -e "    ${CYAN}2${RESET}.  Download from SRA  ${DIM}(give an SRP or PRJNA accession)${RESET}"
blank

echo "  Choice [1/2, default: 1]:"
read -r -p "  > " DATA_SOURCE_CHOICE
DATA_SOURCE_CHOICE="${DATA_SOURCE_CHOICE:-1}"

SRA_AUTO_PAIRED=""
SRA_ACCESSION_USED=""
SRA_MODE=false

if [[ "$DATA_SOURCE_CHOICE" == "2" ]]; then
    SRA_MODE=true
    sra_lookup
    # sra_lookup sets SRA_AUTO_PAIRED for Step 5 to use as a default.
    # The actual download happens right before Step 9 (sample discovery),
    # after all other settings are confirmed.
else
    echo -e "  Enter the path to your folder containing raw FASTQ files."
    echo -e "  ${DIM}Supports: sample_R1.fastq.gz / sample_R2.fastq.gz${RESET}"
echo -e "  ${DIM}Paste Linux paths, C:\\... / F:\\..., or \\\\wsl.localhost\\... paths directly.${RESET}"
    echo -e "  ${DIM}Tip: add -newdefault after a path to save it as the default output base folder.${RESET}"
    blank

    while true; do
        echo "  Input folder:"
        IFS= read -r -e -p "  > " INPUT_DIR
        normalize_path_var INPUT_DIR

        if [[ -z "$INPUT_DIR" ]]; then
            err "Please enter a path."
        elif [[ ! -d "$INPUT_DIR" ]]; then
            err "Directory not found: $INPUT_DIR"
        else
            FASTQ_COUNT=$(find "$INPUT_DIR" -maxdepth 1 \( \
                -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \
                -o -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" \
                -o -name "*_1.fastq.gz"  -o -name "*_1.fq.gz" \
                -o -name "*_R1_*.fastq"  -o -name "*_R1.fastq" \
                -o -name "*_1.fastq" \
            \) 2>/dev/null | wc -l)

            if [[ "$FASTQ_COUNT" -eq 0 ]]; then
                FASTQ_COUNT=$(find "$INPUT_DIR" -maxdepth 1 \( \
                    -name "*.fastq.gz" -o -name "*.fq.gz" \
                    -o -name "*.fastq" -o -name "*.fq" \
                \) 2>/dev/null | wc -l)
            fi

            if [[ "$FASTQ_COUNT" -eq 0 ]]; then
                warn "No FASTQ files found in $INPUT_DIR"
                echo "  Use this folder anyway? [y/N]"
                read -r -p "  > " CONFIRM
                [[ "${CONFIRM,,}" == "y" ]] && break
            else
                ok "Found $FASTQ_COUNT FASTQ file(s) in $INPUT_DIR"
                break
            fi
        fi
    done
fi

# ─────────────────────────────────────────────────────────────
# STEP 2 — Output folder
# ─────────────────────────────────────────────────────────────

header "Step 2 · Output Folder"

DEFAULT_OUTPUT="$DEFAULT_OUTPUT_BASE/$RUN_ID"
echo -e "  Where should results be written?"
echo -e "  ${DIM}Default: $DEFAULT_OUTPUT${RESET}"
echo -e "  ${DIM}Tip: add -newdefault after a path to save it as the default output base for future runs.${RESET}"
blank

echo "  Output folder [press Enter for default]:"
IFS= read -r -e -p "  > " OUTPUT_DIR
if strip_newdefault_marker OUTPUT_DIR; then
    normalize_path_var OUTPUT_DIR
    [[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$DEFAULT_OUTPUT_BASE"
    save_user_default DEFAULT_OUTPUT_BASE "$OUTPUT_DIR"
    DEFAULT_OUTPUT_BASE="$OUTPUT_DIR"
    OUTPUT_DIR="$DEFAULT_OUTPUT_BASE/$RUN_ID"
else
    normalize_path_var OUTPUT_DIR
    [[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$DEFAULT_OUTPUT"
fi

ok "Output directory selected: $OUTPUT_DIR"
echo -e "  ${DIM}It will be created during preflight after confirmation.${RESET}"

# ─────────────────────────────────────────────────────────────
# STEP 3 — Genome selection
# ─────────────────────────────────────────────────────────────

header "Step 3 · Reference Genome"

# Five assemblies are supported out of the box. Refgenie supplies
# sequence/alignment assets when available. For rn7 and danRer11, the same
# assembly is downloaded from UCSC and indexed locally. A user-defined profile option extends this to any UCSC assembly with matching TxDb and OrgDb packages.
GENOME_MODE="refgenie"
GENOME_IS_CUSTOM=false
LOCAL_BT2_INDEX=""
LOCAL_BT2_INDEX_SHA256=""
LOCAL_CHROM_SIZES=""
LOCAL_CHROM_SIZES_SHA256=""
LOCAL_CUSTOM_FASTA=""
LOCAL_CUSTOM_FASTA_SHA256=""
CUSTOM_ASSEMBLY_DIR=""
RUN_GENOME_SIZE=""
GENOME_SIZE_METHOD=""

echo -e "  Supported genomes:"
echo -e "  ${DIM}(R) = available on Refgenie server — fast automatic download${RESET}"
echo -e "  ${DIM}(U) = UCSC sequence download + local Bowtie2 build${RESET}"
blank
for i in "${!SUPPORTED_GENOMES[@]}"; do
    g="${SUPPORTED_GENOMES[$i]}"
    if [[ "${GENOME_ON_REFGENIE[$g]:-no}" == "yes" ]]; then
        tag="${GREEN}(R)${RESET}"
    else
        tag="${YELLOW}(U)${RESET}"
    fi
    printf "    ${CYAN}%2d${RESET}.  %-12s %b\n" "$((i+1))" "$g" "$tag"
done
printf "    ${CYAN}%2s${RESET}.  %-12s %b\n" "C" "User profile" "${DIM}(UCSC + TxDb + OrgDb)${RESET}"
blank
echo -e "  ${DIM}User-defined profiles require a UCSC assembly, matching TxDb, and matching OrgDb.${RESET}"
echo -e "  ${DIM}KEGG remains optional; required combinations are validated after confirmation.${RESET}"
echo -e "  ${DIM}Genome assets are checked or created after the final run summary.${RESET}"
blank

while true; do
    echo "  Enter genome name, number, or C for custom [default: 1/mm10]:"
    read -r -p "  > " GENOME_INPUT
    GENOME_INPUT="${GENOME_INPUT:-1}"

    if [[ "${GENOME_INPUT,,}" == "c" || "${GENOME_INPUT,,}" == "custom" ]]; then
        collect_custom_genome
        break
    fi

    if [[ "$GENOME_INPUT" =~ ^[0-9]+$ ]]; then
        IDX=$((GENOME_INPUT - 1))
        if [[ $IDX -ge 0 && $IDX -lt ${#SUPPORTED_GENOMES[@]} ]]; then
            GENOME="${SUPPORTED_GENOMES[$IDX]}"
        else
            err "Invalid number. Enter 1-${#SUPPORTED_GENOMES[@]}, or C for custom."
            continue
        fi
    else
        GENOME="$GENOME_INPUT"
    fi

    if is_supported_genome "$GENOME"; then
        ok "Genome selected: $GENOME"
        break
    fi

    err "Unsupported genome: $GENOME"
    echo -e "  ${DIM}Choose one of: ${SUPPORTED_GENOMES[*]}, or C for a user-defined genome profile.${RESET}"
done

# ─────────────────────────────────────────────────────────────
# STEP 3b — Download method
# ─────────────────────────────────────────────────────────────
# Only shown for the five built-in genomes. A user-defined genome profile always builds
# locally (collect_custom_genome already set GENOME_BUILD_METHOD="local_build"),
# so this whole step is skipped for it -- there's no Refgenie server entry
# to choose between.
# If genome is not on the Refgenie server, UCSC build is required
# regardless of preference — user must confirm before proceeding.
# ─────────────────────────────────────────────────────────────

if ! $GENOME_IS_CUSTOM; then

GENOME_BUILD_METHOD=""   # "refgenie" | "ucsc"

    blank
    header "Step 3b · Genome Download Method"

    _on_refgenie="${GENOME_ON_REFGENIE[$GENOME]:-no}"

    if [[ "$_on_refgenie" == "yes" ]]; then
        echo -e "  How would you like to obtain the ${BOLD}$GENOME${RESET} genome assets?"
        blank
        echo -e "    ${CYAN}1${RESET}.  ${GREEN}Refgenie download${RESET}  ${DIM}(recommended — fast, automatic, ~minutes)${RESET}"
        echo -e "        Downloads a pre-built index from the Refgenie server."
        blank
        echo -e "    ${CYAN}2${RESET}.  ${YELLOW}UCSC build${RESET}          ${DIM}(slower — downloads FASTA and builds index locally, ~1-2 hrs)${RESET}"
        echo -e "        Downloads the raw FASTA from UCSC and builds the Bowtie2 index"
        echo -e "        on your machine. Gene annotation is handled later by Bioconductor."
        blank

        while true; do
            echo "  Choice [1/2, default: 1]:"
            read -r -p "  > " _method_choice
            _method_choice="${_method_choice:-1}"
            case "$_method_choice" in
                1) GENOME_BUILD_METHOD="refgenie"; ok "Using Refgenie download."; break ;;
                2) GENOME_BUILD_METHOD="local_build";     ok "Using local build from UCSC.";        break ;;
                *) err "Enter 1 or 2." ;;
            esac
        done
    else
        # Genome not on Refgenie server — UCSC build is the only option.
        echo -e "  ${YELLOW}⚠  $GENOME is not available on the standard Refgenie server.${RESET}"
        echo -e "  A ${BOLD}local explicit build${RESET} is required regardless of preference."
        blank
        echo -e "  The runner will automatically:"
        echo -e "    1.  Download the $GENOME FASTA from UCSC (~several GB)"
        echo -e "    2.  Build a Bowtie2 index                (~1-2 hours)"
        echo -e "    3.  Reuse the local files on future runs"
        blank
        echo -e "  ${DIM}This only needs to happen once. Subsequent runs reuse the local cache.${RESET}"
        blank

        while true; do
            echo "  Understood — proceed with UCSC build? [y/N]"
            read -r -p "  > " _ucsc_confirm
            case "${_ucsc_confirm,,}" in
                y|yes) GENOME_BUILD_METHOD="local_build"; ok "Local build from UCSC confirmed."; break ;;
                n|no|"") die "Genome build cancelled. Choose a genome available on Refgenie (marked with (R))." ;;
                *) err "Enter y or n." ;;
            esac
        done
    fi

fi   # ! $GENOME_IS_CUSTOM

# ─────────────────────────────────────────────────────────────
# STEP 4 — Threads
# ─────────────────────────────────────────────────────────────

header "Step 4 · Concurrent Samples"

AVAIL_CORES=$(nproc 2>/dev/null || echo 4)
CPU_BUDGET=$(( AVAIL_CORES > 8 ? AVAIL_CORES - 8 : AVAIL_CORES ))

echo -e "  Many pipeline steps (peak calling in particular) do not scale with"
echo -e "  ${BOLD}--threads${RESET} -- they use roughly one core no matter how many threads"
echo -e "  a single sample is given. Processing multiple samples ${BOLD}at once${RESET} uses"
echo -e "  the rest of the machine instead of leaving it idle between those steps."
blank

DEFAULT_PARALLEL=$(( CPU_BUDGET < 4 ? CPU_BUDGET : 4 ))
[[ "$DEFAULT_PARALLEL" -lt 1 ]] && DEFAULT_PARALLEL=1
echo -e "  Available CPU cores: ${BOLD}$AVAIL_CORES${RESET}"
echo -e "  ${DIM}Suggested: $DEFAULT_PARALLEL concurrent sample(s)${RESET}"
blank

while true; do
    echo "  Samples to process concurrently [default: $DEFAULT_PARALLEL]:"
    read -r -p "  > " PARALLEL_INPUT
    if [[ -z "$PARALLEL_INPUT" ]]; then
        PARALLEL_SAMPLES="$DEFAULT_PARALLEL"
        break
    elif [[ "$PARALLEL_INPUT" =~ ^[0-9]+$ ]] && [[ "$PARALLEL_INPUT" -ge 1 ]]; then
        PARALLEL_SAMPLES="$PARALLEL_INPUT"
        break
    else
        err "Enter a positive integer (1 = sequential)."
    fi
done

if [[ "$PARALLEL_SAMPLES" -gt 1 ]]; then
    ok "Processing up to $PARALLEL_SAMPLES samples concurrently."
    echo -e "  ${DIM}Per-sample console output is suppressed at this concurrency level to avoid${RESET}"
    echo -e "  ${DIM}interleaved output -- each sample still logs live to its own file in logs/.${RESET}"
else
    ok "Processing samples sequentially (concurrency disabled)."
fi

# ─────────────────────────────────────────────────────────────
# STEP 4b — Threads per sample
# ─────────────────────────────────────────────────────────────

header "Step 4b · CPU Threads"

DEFAULT_THREADS=$(( CPU_BUDGET / PARALLEL_SAMPLES ))
[[ "$DEFAULT_THREADS" -lt 1 ]] && DEFAULT_THREADS=1

echo -e "  Available CPU cores: ${BOLD}$AVAIL_CORES${RESET}   Concurrent samples: ${BOLD}$PARALLEL_SAMPLES${RESET}"
echo -e "  ${DIM}Suggested: $DEFAULT_THREADS thread(s)/sample ($DEFAULT_THREADS x $PARALLEL_SAMPLES = $(( DEFAULT_THREADS * PARALLEL_SAMPLES )) cores)${RESET}"
blank

while true; do
    echo "  Threads to use per sample [default: $DEFAULT_THREADS]:"
    read -r -p "  > " THREADS_INPUT
    if [[ -z "$THREADS_INPUT" ]]; then
        THREADS="$DEFAULT_THREADS"
        break
    elif [[ "$THREADS_INPUT" =~ ^[0-9]+$ ]] && [[ "$THREADS_INPUT" -ge 1 ]]; then
        THREADS="$THREADS_INPUT"
        break
    else
        err "Enter a positive integer."
    fi
done

ok "Using $THREADS threads per sample."

REQUESTED_CORES=$(( THREADS * PARALLEL_SAMPLES ))
if [[ "$REQUESTED_CORES" -gt "$AVAIL_CORES" ]]; then
    warn "$THREADS threads x $PARALLEL_SAMPLES concurrent samples = $REQUESTED_CORES cores requested, but only $AVAIL_CORES are available."
    echo -e "  ${DIM}This will oversubscribe the CPU -- samples will contend for cores rather${RESET}"
    echo -e "  ${DIM}than speed things up. Consider lowering threads or concurrency.${RESET}"
fi

# ─────────────────────────────────────────────────────────────
# STEP 5 — Paired-end or single-end
# ─────────────────────────────────────────────────────────────

header "Step 5 · Library Type"

echo -e "  Is this paired-end or single-end data?"
blank
echo -e "    ${CYAN}1${RESET}.  Paired-end  (R1 + R2)"
echo -e "    ${CYAN}2${RESET}.  Single-end  (R1 only)"
blank

PE_DEFAULT=1
if [[ -n "${SRA_AUTO_PAIRED:-}" ]]; then
    if $SRA_AUTO_PAIRED; then
        PE_DEFAULT=1
        ok "Detected from SRA/ENA metadata ($SRA_ACCESSION): all selected runs are paired-end."
    else
        PE_DEFAULT=2
        ok "Detected from SRA/ENA metadata ($SRA_ACCESSION): all selected runs are single-end."
    fi
    echo -e "  ${DIM}Press Enter to accept, or choose the other option if this looks wrong.${RESET}"
    blank
fi

while true; do
    echo "  Choice [1/2, default: $PE_DEFAULT]:"
    read -r -p "  > " PE_INPUT
    PE_INPUT="${PE_INPUT:-$PE_DEFAULT}"
    case "$PE_INPUT" in
        1) PAIRED=true;  ok "Paired-end mode.";  break ;;
        2) PAIRED=false; ok "Single-end mode."; break ;;
        *) err "Enter 1 or 2." ;;
    esac
done

# ─────────────────────────────────────────────────────────────
# STEP 6 — Blacklist filtering
# ─────────────────────────────────────────────────────────────

header "Step 7 · Blacklist Filtering"

echo -e "  Exclude known problematic genomic regions from peak calls?"
echo -e "  ${DIM}Removes artifactual signal from repetitive/telomeric regions.${RESET}"
echo -e "  ${DIM}Recommended for publication-quality results.${RESET}"
echo -e "  ${DIM}Blacklist files are downloaded after the final run summary.${RESET}"
blank

USE_BLACKLIST=false
BLACKLIST_PATH=""
BLACKLIST_SOURCE="none"

if [[ -n "${BLACKLIST_URLS[$GENOME]+x}" ]]; then
    BLACKLIST_FILE="$BLACKLIST_DIR/${GENOME}-blacklist.v2.bed.gz"

    echo "  Enable blacklist filtering for $GENOME? [Y/n]"
    read -r -p "  > " BL_INPUT
    if [[ "${BL_INPUT,,}" != "n" ]]; then
        USE_BLACKLIST=true
        BLACKLIST_PATH="$BLACKLIST_FILE"
        BLACKLIST_SOURCE="standard"

        if [[ -f "$BLACKLIST_FILE" ]]; then
            ok "Blacklist already present: $BLACKLIST_FILE"
        else
            ok "Blacklist queued for download after final confirmation."
            echo -e "  ${DIM}Will download to: $BLACKLIST_FILE${RESET}"
        fi
    else
        ok "Skipping blacklist filtering."
    fi
else
    warn "No standard blacklist available for $GENOME."
    echo -e "  ${DIM}You can provide a custom blacklist BED file if you have one.${RESET}"
    blank
    echo "  Path to custom blacklist BED file (or press Enter to skip):"
    IFS= read -r -e -p "  > " CUSTOM_BL
    normalize_path_var CUSTOM_BL

    if [[ -n "$CUSTOM_BL" ]]; then
        USE_BLACKLIST=true
        BLACKLIST_PATH="$CUSTOM_BL"
        BLACKLIST_SOURCE="custom"

        if [[ -f "$CUSTOM_BL" ]]; then
            ok "Custom blacklist selected: $BLACKLIST_PATH"
        else
            warn "Custom blacklist path was entered but does not currently exist."
            echo -e "  ${DIM}Preflight will stop before launch if this file is still missing.${RESET}"
        fi
    else
        ok "Skipping blacklist filtering."
    fi
fi

# ─────────────────────────────────────────────────────────────
# STEP 8 — Sample discovery
# ─────────────────────────────────────────────────────────────

# In SRA mode files don't exist yet — discovery runs again after the
# actual download (right after final confirmation). Show a metadata
# preview here so the run summary is still useful.

discover_samples() {
    mapfile -t R1_FILES < <(
        find "$INPUT_DIR" -maxdepth 1 \( \
            -name "*_R1_*.fastq.gz" -o \
            -name "*_R1_*.fq.gz"    -o \
            -name "*_R1.fastq.gz"   -o \
            -name "*_R1.fq.gz"      -o \
            -name "*_1.fastq.gz"    -o \
            -name "*_1.fq.gz"       -o \
            -name "*_R1_*.fastq"    -o \
            -name "*_R1.fastq"      -o \
            -name "*_1.fastq" \
        \) | sort
    )

    if [[ ${#R1_FILES[@]} -eq 0 ]] && ! $PAIRED; then
        mapfile -t R1_FILES < <(
            find "$INPUT_DIR" -maxdepth 1 \( \
                -name "*.fastq.gz" -o -name "*.fq.gz" \
                -o -name "*.fastq" -o -name "*.fq" \
            \) | sort
        )
    fi

    [[ ${#R1_FILES[@]} -gt 0 ]] || die "No FASTQ files found in $INPUT_DIR"

    SAMPLE_NAMES=()
    R2_FILES_ARR=()

    for R1 in "${R1_FILES[@]}"; do
        BASENAME=$(basename "$R1")
        SAMPLE=$(echo "$BASENAME" \
            | sed -E 's/_R1_[^.]+//' \
            | sed -E 's/_R1\./\./'   \
            | sed -E 's/_1\./\./'    \
            | sed -E 's/\.(fastq|fq)(\.gz)?$//' \
        )
        SAMPLE_NAMES+=("$SAMPLE")

        # "_1"/"_2" is SRA/ENA's own mate-pair suffix, so it's matched here
        # on purpose. But the same suffix is also a very common way to label
        # independent replicates (e.g. Control_1.fastq.gz, Control_2.fastq.gz
        # as two separate single-end samples, not mates of one sample). If
        # the match came from the unambiguous _R1 pattern there's no risk;
        # if it only came from the bare _1 pattern, flag it here so it's
        # visible at the "Proceed with these samples?" gate below instead of
        # silently collapsing two real samples into one paired-end sample.
        ambiguous_pairing=false
        if [[ "$BASENAME" != *_R1_* && "$BASENAME" != *_R1.* ]] && [[ "$BASENAME" == *_1.* ]]; then
            ambiguous_pairing=true
        fi

        if $PAIRED; then
            R2=$(echo "$R1" \
                | sed 's/_R1_/_R2_/' \
                | sed 's/_R1\./_R2./' \
                | sed 's/_1\./_2./' \
            )
            if [[ -f "$R2" ]]; then
                R2_FILES_ARR+=("$R2")
                if $ambiguous_pairing; then
                    printf "    ${YELLOW}⚠${RESET}  %-40s  +  %s  ${YELLOW}(ambiguous _1/_2 suffix — confirm these are true mates, not two independent replicates)${RESET}\n" "$SAMPLE" "$(basename "$R2")"
                else
                    printf "    ${GREEN}✔${RESET}  %-40s  +  %s\n" "$SAMPLE" "$(basename "$R2")"
                fi
            else
                R2_FILES_ARR+=("")
                printf "    ${YELLOW}⚠${RESET}  %-40s  (R2 not found — will run single-end)\n" "$SAMPLE"
            fi
        else
            R2_FILES_ARR+=("")
            printf "    ${GREEN}✔${RESET}  %s\n" "$SAMPLE"
        fi
    done
}

# assert_unique_sample_names
# The _R1_/_1 stripping in discover_samples() derives a sample name from each
# FASTQ's basename. Two different files can reduce to the identical derived
# name -- most commonly separate lane files from the same instrument run
# (Sample_R1_L001.fastq.gz and Sample_R1_L002.fastq.gz both become "Sample").
# There is no lane-merging step in this script: if that happens, two workers
# would later write into the same sample output directory, the same log
# file, and the same status/exit-code records, and the result would be an
# undefined mixture of both jobs. Fail loudly here, before any output paths
# are created, instead of letting that collision happen silently downstream.
assert_unique_sample_names() {
    local -A seen_at=()
    local collided=false
    local i name

    for i in "${!SAMPLE_NAMES[@]}"; do
        name="${SAMPLE_NAMES[$i]}"
        if [[ -n "${seen_at[$name]:-}" ]]; then
            collided=true
            err "Duplicate derived sample name: '$name'"
            err "  ${seen_at[$name]}"
            err "  ${R1_FILES[$i]}"
        else
            seen_at[$name]="${R1_FILES[$i]}"
        fi
    done

    if $collided; then
        blank
        err "Two or more FASTQ files reduce to the same sample name above."
        err "This is almost always separate lane files (e.g. _L001/_L002) from"
        err "the same sequencing run. This script has no lane-merging step, so"
        err "running them as-is would let concurrent jobs collide on the same"
        err "output directory, log file, and status records."
        err ""
        err "Fix by either renaming the files so each biological sample has a"
        err "unique name, or concatenating each sample's lane FASTQs into one"
        err "file per sample/mate before running this script."
        die "Aborted: duplicate derived sample names (see above)."
    fi
}

# confirm_missing_r2
# In paired mode, silently treating a sample with no R2 mate as single-end
# changes alignment, duplicate handling, fragment interpretation, and peak
# calling for that one sample relative to the rest of the batch -- and mixes
# single- and paired-end results into what's supposed to be one comparison.
# Require an explicit choice instead of quietly proceeding.
confirm_missing_r2() {
    $PAIRED || return 0

    local missing_idx=()
    local i
    for i in "${!SAMPLE_NAMES[@]}"; do
        if [[ -z "${R2_FILES_ARR[$i]}" ]]; then
            missing_idx+=("$i")
        fi
    done

    [[ ${#missing_idx[@]} -eq 0 ]] && return 0

    blank
    warn "This run is paired-end, but ${#missing_idx[@]} sample(s) have no R2 mate:"
    for i in "${missing_idx[@]}"; do
        echo -e "    ${YELLOW}✘${RESET}  ${SAMPLE_NAMES[$i]}"
    done
    echo -e "  ${DIM}Running these as single-end changes alignment, duplicate handling,${RESET}"
    echo -e "  ${DIM}fragment interpretation, and peak calling for just these samples.${RESET}"
    blank
    echo "  Exclude these samples and continue with the rest? [Y/n]"
    read -r -p "  > " R2_CHOICE
    if [[ "${R2_CHOICE,,}" == "n" ]]; then
        die "Stopped: resolve the missing R2 files (or rerun in single-end mode) and try again."
    fi

    local new_names=() new_r1=() new_r2=()
    for i in "${!SAMPLE_NAMES[@]}"; do
        local skip=false m
        for m in "${missing_idx[@]}"; do
            [[ "$i" == "$m" ]] && skip=true && break
        done
        if ! $skip; then
            new_names+=("${SAMPLE_NAMES[$i]}")
            new_r1+=("${R1_FILES[$i]}")
            new_r2+=("${R2_FILES_ARR[$i]}")
        fi
    done
    SAMPLE_NAMES=("${new_names[@]}")
    R1_FILES=("${new_r1[@]}")
    R2_FILES_ARR=("${new_r2[@]}")

    [[ ${#SAMPLE_NAMES[@]} -eq 0 ]] && die "No samples remain after excluding those missing R2."
    ok "Excluded ${#missing_idx[@]} sample(s). Continuing with ${#SAMPLE_NAMES[@]} sample(s)."
}

header "Step 9 · Sample Discovery"

if $SRA_MODE; then
    echo -e "  ${DIM}Files have not been downloaded yet — showing the planned sample list${RESET}"
    echo -e "  ${DIM}from SRA/ENA metadata. Discovery re-runs for real after downloading.${RESET}"
    blank
    for i in "${SRA_SELECTED_IDX[@]}"; do
        if [[ "${SRA_RUN_LAYOUT[$i]}" == "PAIRED" ]]; then
            printf "    ${GREEN}✔${RESET}  %-40s  +  %s\n" "${SRA_RUN_LABEL[$i]}" "(R2 expected)"
        else
            printf "    ${GREEN}✔${RESET}  %s\n" "${SRA_RUN_LABEL[$i]}"
        fi
    done
    SRA_PLANNED_SAMPLE_COUNT="${#SRA_SELECTED_IDX[@]}"
    SAMPLE_NAMES=()
    R1_FILES=()
    R2_FILES_ARR=()
else
    declare -a SAMPLE_NAMES=()
    declare -a R1_FILES=()
    declare -a R2_FILES_ARR=()
    discover_samples
    assert_unique_sample_names
    confirm_missing_r2
fi

blank
echo "  Proceed with these samples? [Y/n]"
read -r -p "  > " CONFIRM_SAMPLES
[[ "${CONFIRM_SAMPLES,,}" == "n" ]] && die "Aborted by user."

# ─────────────────────────────────────────────────────────────
# STEP 10 — Run summary and final confirmation
# ─────────────────────────────────────────────────────────────

header "Step 10 · Run Summary"

if $SRA_MODE; then
    input_bytes="$SRA_TOTAL_BYTES"
    sample_count_display="$SRA_PLANNED_SAMPLE_COUNT (planned)"
else
    input_bytes=$(sum_selected_fastq_bytes)
    sample_count_display="${#SAMPLE_NAMES[@]}"
fi

echo -e "  ${BOLD}Input folder:${RESET}    $( $SRA_MODE && echo "(downloads pending — $SRA_DOWNLOAD_DIR)" || echo "$INPUT_DIR" )"
echo -e "  ${BOLD}Output folder:${RESET}   $OUTPUT_DIR"
echo -e "  ${BOLD}Genome:${RESET}          $GENOME$( $GENOME_IS_CUSTOM && echo " (custom)" || echo "" )"
echo -e "  ${BOLD}Genome mode:${RESET}     ${GENOME_MODE:-refgenie}"
if $GENOME_IS_CUSTOM; then
    echo -e "  ${BOLD}Asset source:${RESET}    validated UCSC genome profile"
    echo -e "  ${BOLD}UCSC assembly:${RESET}   ${CUSTOM_UCSC_ASSEMBLY[$GENOME]:-}"
    echo -e "  ${BOLD}TxDb package:${RESET}    ${CUSTOM_TXDB_PKG[$GENOME]:-}"
    echo -e "  ${BOLD}OrgDb package:${RESET}   ${CUSTOM_ORGDB_PKG[$GENOME]:-}"
    echo -e "  ${BOLD}KEGG code:${RESET}       ${CUSTOM_KEGG_ORG[$GENOME]:-<none -- KEGG enrichment skipped>}"
else
    echo -e "  ${BOLD}Asset source:${RESET}    ${GENOME_BUILD_METHOD:-refgenie}$( [[ "${GENOME_BUILD_METHOD:-refgenie}" == "local_build" ]] && echo " (local explicit build)" || echo " (Refgenie server)" )"
fi
echo -e "  ${BOLD}Threads:${RESET}         $THREADS"
echo -e "  ${BOLD}Library:${RESET}         $( $PAIRED && echo "Paired-end" || echo "Single-end" )"
echo -e "  ${BOLD}Blacklist:${RESET}       $( $USE_BLACKLIST && echo "$BLACKLIST_PATH" || echo "No" )"
echo -e "  ${BOLD}Samples:${RESET}         $sample_count_display"
echo -e "  ${BOLD}FASTQ size:${RESET}      $(human_bytes "$input_bytes")"
blank

STEP_NUM=1
echo -e "  ${DIM}After this confirmation, the runner will:${RESET}"
if $SRA_MODE; then
    echo -e "  ${DIM}  $STEP_NUM. Download FASTQ files from SRA/ENA (ENA FTP, with checksum verification)${RESET}"
    STEP_NUM=$((STEP_NUM+1))
fi
echo -e "  ${DIM}  $STEP_NUM. Write a run manifest, resume state, and R sample sheet template${RESET}"
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM}  $STEP_NUM. Check available disk space${RESET}"
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM}  $STEP_NUM. Check the launcher installation${RESET}"
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM}  $STEP_NUM. Download/build missing sequence assets once (Refgenie or explicit local reference)${RESET}"
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM}  $STEP_NUM. Safely download and verify the blacklist once, if requested${RESET}"
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM}  $STEP_NUM. Run preflight checks${RESET}"
STEP_NUM=$((STEP_NUM+1))
if [[ "${PARALLEL_SAMPLES:-1}" -gt 1 ]]; then
    echo -e "  ${DIM}  $STEP_NUM. Process up to $PARALLEL_SAMPLES samples concurrently${RESET}"
else
    echo -e "  ${DIM}  $STEP_NUM. Process samples sequentially${RESET}"
fi
STEP_NUM=$((STEP_NUM+1))
echo -e "  ${DIM} $STEP_NUM. Write end-of-run QC tables${RESET}"
blank

echo "  Start downloads, preflight, and processing? [Y/n]"
read -r -p "  > " FINAL_CONFIRM
[[ "${FINAL_CONFIRM,,}" == "n" ]] && die "Aborted by user."

fi   # end of "if ! $RESUME_MODE"

# ─────────────────────────────────────────────────────────────
# STEP 11 — Downloads, index, preflight, process
# ─────────────────────────────────────────────────────────────

# In SRA mode: all prompts answered, confirmation given — download now.
if ! $RESUME_MODE && $SRA_MODE; then
    sra_download
    header "Step 9 (continued) · Sample Discovery"
    declare -a SAMPLE_NAMES=()
    declare -a R1_FILES=()
    declare -a R2_FILES_ARR=()
    discover_samples
    assert_unique_sample_names
    confirm_missing_r2
fi

initialize_run_outputs
check_disk_space
check_launcher_requirements
ensure_genome_assets
ensure_blacklist_if_requested

# Resolve and freeze every run-wide reference value before preflight and
# before the final manifest/resume-state rewrite.
if [[ "${GENOME_MODE:-refgenie}" == "local_explicit" ]]; then
    RUN_GENOME_INDEX="$LOCAL_BT2_INDEX"
    RUN_CHROM_SIZES="$LOCAL_CHROM_SIZES"
else
    RUN_GENOME_INDEX=$(in_env refgenie seek -c "$REFGENIE_CONFIG" "${GENOME}/bowtie2_index.bowtie2_index")
    RUN_CHROM_SIZES=$(in_env refgenie seek -c "$REFGENIE_CONFIG" "${GENOME}/fasta.chrom_sizes")
fi
resolve_effective_genome_size "$GENOME" "$RUN_CHROM_SIZES"
ok "Resolved genome assets once for this run (effective genome size: $RUN_GENOME_SIZE; method: $GENOME_SIZE_METHOD)."

# Same "resolve once here" treatment for --TSS-name/--anno-name: custom
# profiles already have validated frozen TxDb/OrgDb assets by this point
# (verify_custom_profile_assets ran inside ensure_genome_assets); built-in
# genomes get theirs built or reused here for the first time. Either way
# this runs exactly once, uniformly, regardless of which of the three
# ensure_genome_assets() branches produced RUN_CHROM_SIZES above.
#
# On a resume where annotation was already resolved for this run (resume
# state schema 4+, RUN_ANNOTATION_STATUS non-empty from load_resume_state),
# trust and re-verify that instead of resolving fresh -- otherwise a
# resumed run could silently pick up a newer TxDb/OrgDb installed since
# the run started, overwrite the snapshot to match it, while
# already-completed samples' QC output still reflects the original.
if $RESUME_MODE && [[ -n "${RUN_ANNOTATION_STATUS:-}" ]]; then
    verify_resumed_annotation_assets
else
    resolve_run_annotation_assets "$GENOME" "$RUN_CHROM_SIZES"
fi

# write_reference_snapshot now runs for every genome, not just custom
# profiles -- schema 4 records the generic ANNOTATION_* fields (built from
# RUN_ANNOTATION_* above) for built-ins too, which is what lets
# diff_analysis.sh load a frozen, fingerprinted TxDb/OrgDb for a built-in
# genome instead of whatever happens to be installed when it runs.
write_reference_snapshot "$GENOME"
if $GENOME_IS_CUSTOM; then
    if $USE_BLACKLIST && [[ "$BLACKLIST_SOURCE" == "custom" ]]; then
        validate_custom_blacklist "$BLACKLIST_PATH" "$RUN_CHROM_SIZES"
    fi
fi

write_run_manifest
write_resume_state
run_preflight_checks

# ─────────────────────────────────────────────────────────────
# STEP 12 — Process samples
# ─────────────────────────────────────────────────────────────

header "Processing Samples"

TOTAL=${#SAMPLE_NAMES[@]}
PASS=0
FAIL=0
SKIPPED=0
FAIL_LIST=()

declare -a SAMPLE_STATUS=()
declare -a SAMPLE_EXIT_CODES=()
declare -a SAMPLE_LOGS=()
declare -a SAMPLE_OUT_DIRS=()
declare -a SAMPLE_BAMS=()
declare -a SAMPLE_PEAKS=()
declare -a SAMPLE_PEAK_COUNTS=()
declare -a SAMPLE_BAM_SIZES=()

LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# Computed once here (immediately after the assets it fingerprints are
# resolved) and reused by every sample_signature() call below -- see the
# "Run signature" section for why this exists.
RUN_WIDE_FINGERPRINT="$(compute_run_wide_fingerprint)"

if [[ "${PARALLEL_SAMPLES:-1}" -gt 1 ]]; then
    EST_MEM_PER_SAMPLE_MB=$(estimate_memory_per_sample_mb "$RUN_GENOME_INDEX")
    EST_TOTAL_MEM_MB=$(( EST_MEM_PER_SAMPLE_MB * PARALLEL_SAMPLES ))
    AVAIL_MEM_MB=$(free -m 2>/dev/null | awk '/^Mem:/{print $7}')
    if [[ -n "$AVAIL_MEM_MB" && "$EST_TOTAL_MEM_MB" -gt "$AVAIL_MEM_MB" ]]; then
        warn "Estimated memory need (~${EST_MEM_PER_SAMPLE_MB} MB/sample x $PARALLEL_SAMPLES = ~${EST_TOTAL_MEM_MB} MB) exceeds available memory (~${AVAIL_MEM_MB} MB)."
        echo -e "  ${DIM}This is a rough estimate based on genome index size, not an exact${RESET}"
        echo -e "  ${DIM}measurement -- but running this concurrently risks swapping or${RESET}"
        echo -e "  ${DIM}out-of-memory failures. Consider lowering concurrent samples.${RESET}"
    fi
fi

# process_one_sample INDEX
# Self-contained: safe to run in a backgrounded subshell for concurrent
# processing. Does NOT write to the SAMPLE_* summary arrays -- array writes
# made inside a background subshell never propagate back to the parent
# shell, so results are collected afterward by collect_sample_result()
# instead, reading back from disk.
process_one_sample() {
    local i="$1"
    local SAMPLE="${SAMPLE_NAMES[$i]}"
    local R1="${R1_FILES[$i]}"
    local R2="${R2_FILES_ARR[$i]}"
    local SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"
    local LOG_FILE="$LOG_DIR/${SAMPLE}.log"

    rm -f "${LOG_FILE%.log}.exitcode" "${LOG_FILE%.log}.pid"
    mkdir -p "$SAMPLE_OUT"
    write_attempt_signature "$i"

    local local_genome_index="$RUN_GENOME_INDEX"
    local local_chrom_sizes="$RUN_CHROM_SIZES"
    local GENOME_SIZE="$RUN_GENOME_SIZE"

    local CMD=(
        "${CONDA_EXE:-conda}"
        run
        --no-capture-output
        -n "$ENV_NAME"
        env REFGENIE="$REFGENIE_CONFIG"
        python "$PIPELINE"
        -S "$SAMPLE"
        -I "$R1"
        -G "$GENOME"
        --genome-index "$local_genome_index"
        --chrom-sizes  "$local_chrom_sizes"
        -gs "$GENOME_SIZE"
        -O "$OUTPUT_DIR"
        -P "$THREADS"
    )

    if $PAIRED && [[ -n "$R2" ]]; then
        CMD+=(-Q paired -I2 "$R2")
    else
        CMD+=(-Q single)
    fi

    if $USE_BLACKLIST && [[ -n "$BLACKLIST_PATH" ]]; then
        CMD+=(--blacklist "$BLACKLIST_PATH")
    fi

    # Unconditional on genome type -- conditional only on the deliberately
    # recorded RUN_ANNOTATION_STATUS from resolve_run_annotation_assets().
    # Hashes are re-verified here, at the point of use, on every single
    # sample launch -- not just once at build time -- so a file edited or
    # replaced on disk between samples (or between runs, for a reused
    # fingerprinted built-in annotation) is caught before it reaches
    # PEPATAC rather than trusted because the path still exists.
    case "$RUN_ANNOTATION_STATUS" in
        validated)
            verify_file_sha256 "$RUN_ANNOTATION_TSS_BED" "$RUN_ANNOTATION_TSS_BED_SHA256" "run TSS BED"
            verify_file_sha256 "$RUN_ANNOTATION_FEATURE_BED" "$RUN_ANNOTATION_FEATURE_BED_SHA256" "run feature BED"
            CMD+=(--TSS-name "$RUN_ANNOTATION_TSS_BED" --anno-name "$RUN_ANNOTATION_FEATURE_BED")
            ;;
        disabled_by_override)
            : # PEPATAC_ALLOW_MISSING_ANNOTATION_QC=1 -- omit both flags, already warned and recorded.
            ;;
        *)
            die "Annotation QC assets were not resolved for $GENOME (RUN_ANNOTATION_STATUS='${RUN_ANNOTATION_STATUS:-<unset>}')."$'\n'"       This should not be reachable -- resolve_run_annotation_assets() should have set a status or died already."
            ;;
    esac

    # See Pass 1 (where SAMPLE_RERUN_MODE is populated) for why each mode
    # gets its specific flag: a stale prior PASS cannot trust pypiper's own
    # checkpoints (-N forces it to disregard all of them and genuinely
    # start over); an ordinary prior failure can resume past its leftover
    # lock file cheaply (-R); a sample that was never attempted needs
    # neither.
    local rerun_mode="${SAMPLE_RERUN_MODE[$i]:-}"
    local rerun_note=""
    case "$rerun_mode" in
        restart) CMD+=(-N); rerun_note=" [full restart: -N]" ;;
        recover) CMD+=(-R); rerun_note=" [recovering: -R]" ;;
    esac

    append_command "$SAMPLE" "${CMD[@]}"

    local EXIT_CODE
    if [[ "${PARALLEL_SAMPLES:-1}" -gt 1 ]]; then
        echo -e "  ${CYAN}▶${RESET}  Started: $SAMPLE${rerun_note}  (gs=$GENOME_SIZE, log: $LOG_FILE)"
        if run_sample_command_with_cleanup "$LOG_FILE" --quiet "${CMD[@]}"; then
            EXIT_CODE=0
        else
            EXIT_CODE="$?"
        fi
        if [[ "$EXIT_CODE" -eq 0 ]]; then
            echo -e "  ${GREEN}✔${RESET}  Finished: $SAMPLE  (exit 0)"
        else
            echo -e "  ${RED}✘${RESET}  Finished: $SAMPLE  (exit $EXIT_CODE) -- see $LOG_FILE"
        fi
    else
        blank
        echo -e "  ${DIM}Started: $(date '+%Y-%m-%d %H:%M:%S')${rerun_note}${RESET}"
        echo -e "  ${DIM}Effective genome size (-gs): $GENOME_SIZE${RESET}"
        echo -e "  ${DIM}Log: $LOG_FILE${RESET}"
        echo ""
        CURRENT_SAMPLE="$SAMPLE"
        if run_sample_command_with_cleanup "$LOG_FILE" "${CMD[@]}"; then
            EXIT_CODE=0
        else
            EXIT_CODE="$?"
        fi
        CURRENT_SAMPLE=""
        if [[ "$EXIT_CODE" -eq 0 ]]; then
            ok "pepatac.py exited 0 for $SAMPLE  (finished: $(date '+%H:%M:%S'))"
        else
            err "pepatac.py exited $EXIT_CODE for $SAMPLE"
            err "    See log: $LOG_FILE"
        fi
    fi

    # The exit code is the only thing that needs to survive this subshell --
    # collect_sample_result() reads it back to do full output validation.
    echo "$EXIT_CODE" > "${LOG_FILE%.log}.exitcode"
}

# collect_sample_result INDEX
# Runs serially after process_one_sample has finished for this index (either
# just now, or in a previous run being skipped). Reads results back from
# disk and does the actual PASS/FAIL validation exactly once, in exactly one
# place, regardless of whether the sample was skipped, run sequentially, or
# run concurrently.
collect_sample_result() {
    local i="$1"
    local SAMPLE="${SAMPLE_NAMES[$i]}"
    local SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"
    local LOG_FILE="$LOG_DIR/${SAMPLE}.log"
    local exitcode_file="${LOG_FILE%.log}.exitcode"

    local EXIT_CODE
    if [[ -f "$exitcode_file" ]]; then
        EXIT_CODE="$(cat "$exitcode_file" 2>/dev/null || echo 125)"
        rm -f "$exitcode_file"
    elif [[ "${WAS_SKIPPED[$i]:-0}" == "1" ]]; then
        EXIT_CODE=0
    else
        EXIT_CODE=125
    fi

    SAMPLE_EXIT_CODES[$i]="$EXIT_CODE"
    SAMPLE_LOGS[$i]="$LOG_FILE"
    SAMPLE_OUT_DIRS[$i]="$SAMPLE_OUT"
    SAMPLE_BAMS[$i]="$(detect_best_bam "$SAMPLE_OUT")"
    local -a bam_ambiguous
    mapfile -t bam_ambiguous < <(read_detect_ambiguous)
    SAMPLE_PEAKS[$i]="$(detect_best_peak "$SAMPLE_OUT")"
    local -a peak_ambiguous
    mapfile -t peak_ambiguous < <(read_detect_ambiguous)
    SAMPLE_PEAK_COUNTS[$i]="$(count_peak_lines "${SAMPLE_PEAKS[$i]}")"

    if [[ -n "${SAMPLE_BAMS[$i]}" ]]; then
        SAMPLE_BAM_SIZES[$i]="$(bytes_for_file "${SAMPLE_BAMS[$i]}")"
    else
        SAMPLE_BAM_SIZES[$i]="0"
    fi

    local fail_reason=""
    if [[ "$EXIT_CODE" == "125" ]]; then
        fail_reason="worker ended without recording a result (crashed, killed, or interrupted before finishing)"
    elif [[ "$EXIT_CODE" -ne 0 ]]; then
        fail_reason="pepatac.py exited $EXIT_CODE"
    elif [[ "${#bam_ambiguous[@]}" -gt 1 ]]; then
        fail_reason="multiple candidate BAM files matched at the same priority tier -- refusing to guess which is correct: ${bam_ambiguous[*]}"
    elif [[ -z "${SAMPLE_BAMS[$i]}" ]]; then
        fail_reason="pepatac.py exited 0 but no BAM was found"
    elif [[ ! -s "${SAMPLE_BAMS[$i]}" ]]; then
        fail_reason="pepatac.py exited 0 but the BAM is empty"
    elif ! in_env samtools quickcheck "${SAMPLE_BAMS[$i]}" 2>/dev/null; then
        fail_reason="pepatac.py exited 0 but the BAM failed samtools quickcheck"
    elif [[ "${#peak_ambiguous[@]}" -gt 1 ]]; then
        fail_reason="multiple candidate peak files matched at the same priority tier -- refusing to guess which is correct: ${peak_ambiguous[*]}"
    elif [[ -z "${SAMPLE_PEAKS[$i]}" ]]; then
        fail_reason="pepatac.py exited 0 but no peak file was found"
    elif [[ ! -s "${SAMPLE_PEAKS[$i]}" ]]; then
        fail_reason="pepatac.py exited 0 but the peak file is empty"
    fi

    if [[ -z "$fail_reason" ]]; then
        ok "Sample PASSED: $SAMPLE  (validated: BAM intact, peak file present)"
        SAMPLE_STATUS[$i]="PASS"
        write_sample_status "$SAMPLE" "PASS" "$i"
        (( PASS++ )) || true
    else
        err "Sample FAILED validation: $SAMPLE  ($fail_reason)"
        err "    See log: $LOG_FILE"
        SAMPLE_STATUS[$i]="FAIL"
        write_sample_status "$SAMPLE" "FAIL"
        FAIL_LIST+=("$SAMPLE")
        (( FAIL++ )) || true
    fi
}

# ── Pass 1: fast skip-determination (cheap filesystem checks only) ──────
NEEDS_PROCESSING=()
declare -A WAS_SKIPPED=()
# index -> "restart" (stale PASS, needs -N) | "recover" (prior FAIL, needs
# -R) | unset (never attempted -- no special flag needed). Populated below,
# read by process_one_sample() in Pass 2.
declare -A SAMPLE_RERUN_MODE=()
for i in "${!SAMPLE_NAMES[@]}"; do
    SAMPLE="${SAMPLE_NAMES[$i]}"
    NUM=$(( i + 1 ))

    blank
    echo -e "  ${BOLD}${CYAN}[$NUM/$TOTAL]${RESET}  ${BOLD}$SAMPLE${RESET}"

    PRIOR_STATUS=$(read_sample_status "$SAMPLE")
    if [[ "$PRIOR_STATUS" == "PASS" ]]; then
        SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"
        candidate_bam="$(detect_best_bam "$SAMPLE_OUT")"
        mapfile -t bam_ambiguous < <(read_detect_ambiguous)
        candidate_peak="$(detect_best_peak "$SAMPLE_OUT")"
        mapfile -t peak_ambiguous < <(read_detect_ambiguous)

        # A bare "PASS" string doesn't confirm the actual output still exists
        # or is intact -- FASTQs can be replaced, the genome can change, or a
        # prior output directory can be left with stale/partial files without
        # ever updating this status file. Re-validate before trusting it.
        stale_reason=""
        if [[ "${#bam_ambiguous[@]}" -gt 1 ]]; then
            stale_reason="multiple candidate BAM files matched at the same priority tier -- reprocessing instead of guessing: ${bam_ambiguous[*]}"
        elif [[ -z "$candidate_bam" ]]; then
            stale_reason="no BAM found"
        elif [[ ! -s "$candidate_bam" ]]; then
            stale_reason="BAM is empty"
        elif ! in_env samtools quickcheck "$candidate_bam" 2>/dev/null; then
            stale_reason="BAM failed samtools quickcheck"
        elif [[ "${#peak_ambiguous[@]}" -gt 1 ]]; then
            stale_reason="multiple candidate peak files matched at the same priority tier -- reprocessing instead of guessing: ${peak_ambiguous[*]}"
        elif [[ -z "$candidate_peak" ]]; then
            stale_reason="no peak file found"
        elif [[ ! -s "$candidate_peak" ]]; then
            stale_reason="peak file is empty"
        else
            # Outputs look intact, but intact isn't the same as current --
            # they could have been produced from FASTQs, a genome build, a
            # blacklist, or a pipeline version that this invocation is no
            # longer using. Only trust the skip when the recorded signature
            # from PASS time matches what this invocation would use now.
            sample_sig_file="$LOG_DIR/${SAMPLE}.sig"
            if [[ ! -f "$sample_sig_file" ]]; then
                stale_reason="no run signature on record for this PASS (from before signature tracking was added) -- cannot confirm it matches current inputs/genome/settings"
            else
                recorded_sig="$(cat "$sample_sig_file" 2>/dev/null || echo '')"
                current_sig="$(sample_signature "$i")"
                if [[ "$recorded_sig" != "$current_sig" ]]; then
                    stale_reason="run signature changed -- FASTQ, genome assets, blacklist, or pipeline version differ from the run that produced this PASS"
                fi
            fi
        fi

        if [[ -z "$stale_reason" ]]; then
            ok "Skipping — already PASSED in a previous run."
            WAS_SKIPPED["$i"]=1
            (( SKIPPED++ )) || true
            continue
        else
            # A stale PASS means pepatac.py's own checkpoint files on disk
            # reflect the OLD inputs, not the current ones. Reprocessing
            # this sample with no special flag would leave pypiper's
            # default automatic checkpointing free to silently skip any
            # stage whose checkpoint file still exists from that stale run
            # -- exactly the failure mode this staleness check exists to
            # catch in the first place. -N ("start over and run every
            # command even if output exists," per pypiper/PEPATAC's own
            # docs) forces every checkpoint to be disregarded, so
            # "reprocessing" here can't quietly end up mixing old-stage
            # outputs with new ones.
            warn "Prior PASS for $SAMPLE looks stale ($stale_reason)."
            warn "Forcing a full fresh restart (-N) instead of a plain re-run, so no stage"
            warn "from the stale run can be silently skipped via its old checkpoint file."
            SAMPLE_RERUN_MODE["$i"]="restart"
        fi
    elif [[ "$PRIOR_STATUS" == "FAIL" ]]; then
        # A genuine failure (crash, timeout, interruption -- not staleness)
        # typically leaves a pypiper lock file on whatever stage was in
        # progress when it died. -R ("recover") tells pypiper to overwrite
        # that lock and continue from there, per PEPATAC's own documented
        # recovery workflow -- rather than either risking an error against
        # the leftover lock (no flag) or wastefully discarding legitimately
        # completed upstream stages by forcing a full restart (-N, which
        # is reserved for cases where prior checkpoints can't be trusted).
        #
        # But -R is only safe if the failed attempt's own inputs still
        # match what this invocation would use now. A FAIL status alone
        # doesn't confirm that -- the FASTQ could have been replaced, or
        # the genome/blacklist/pipeline version could have changed, since
        # that attempt started. Resuming with -R against different inputs
        # than the failed attempt used would let pypiper pick up from
        # checkpoints built against inputs that no longer apply -- exactly
        # the same risk a stale PASS carries without -N. Compare against
        # the attempt signature written at the START of that failed
        # attempt (see write_attempt_signature), not the PASS signature,
        # which doesn't exist for a FAIL.
        attempt_sig_file="$LOG_DIR/${SAMPLE}.attempt.sig"
        if [[ ! -f "$attempt_sig_file" ]]; then
            warn "Prior attempt for $SAMPLE failed, but no attempt signature is on record"
            warn "(from before this tracking was added) -- forcing a full restart (-N) to be safe."
            SAMPLE_RERUN_MODE["$i"]="restart"
        else
            recorded_attempt_sig="$(cat "$attempt_sig_file" 2>/dev/null || echo '')"
            current_sig="$(sample_signature "$i")"
            if [[ "$recorded_attempt_sig" == "$current_sig" ]]; then
                warn "Prior attempt for $SAMPLE failed -- resuming with recovery mode (-R)."
                SAMPLE_RERUN_MODE["$i"]="recover"
            else
                warn "Prior attempt for $SAMPLE failed, and inputs changed since then --"
                warn "forcing a full restart (-N) instead of recovery mode, since resuming"
                warn "against different inputs than the failed attempt used is unsafe."
                SAMPLE_RERUN_MODE["$i"]="restart"
            fi
        fi
    else
        # NOT_RUN (no status file at all) doesn't guarantee nothing was
        # ever attempted here. If the runner itself was killed -- WSL or
        # the machine stopping, for instance -- after pepatac.py started
        # but before collect_sample_result ever ran to write a status
        # file, the sample's output directory can already hold partial
        # output or lock files with no status record at all. Running such
        # a sample with no special flag would leave pypiper free to pick
        # up from that unrecorded partial state. Detect non-empty leftover
        # output and force a clean restart instead of trusting it.
        SAMPLE_OUT_CHECK="$OUTPUT_DIR/$SAMPLE"
        if [[ -d "$SAMPLE_OUT_CHECK" ]] && \
           [[ -n "$(find "$SAMPLE_OUT_CHECK" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
            warn "$SAMPLE has no run record, but $SAMPLE_OUT_CHECK already has content --"
            warn "a prior attempt may have been interrupted before it could record a result."
            warn "Forcing a full restart (-N) rather than trusting pypiper's checkpoints in"
            warn "an unknown state."
            SAMPLE_RERUN_MODE["$i"]="restart"
        fi
    fi

    NEEDS_PROCESSING+=("$i")
done

# ── Pass 2: run whatever needs processing, bounded by PARALLEL_SAMPLES ──
if [[ ${#NEEDS_PROCESSING[@]} -gt 0 ]]; then
    blank
    if [[ "${PARALLEL_SAMPLES:-1}" -gt 1 ]]; then
        header "Processing ${#NEEDS_PROCESSING[@]} sample(s), up to $PARALLEL_SAMPLES at a time"
    fi

    running=0
    for i in "${NEEDS_PROCESSING[@]}"; do
        process_one_sample "$i" &
        running=$(( running + 1 ))
        if [[ "$running" -ge "${PARALLEL_SAMPLES:-1}" ]]; then
            wait -n || true
            running=$(( running - 1 ))
        fi
    done
    wait || true
fi

# ── Pass 3: collect + validate results for every sample ─────────────────
blank
for i in "${!SAMPLE_NAMES[@]}"; do
    collect_sample_result "$i"
done

write_qc_summary
write_detected_r_sample_sheet

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────

header "Run Complete"

echo -e "  ${BOLD}Samples total:${RESET}            $TOTAL"
echo -e "  ${GREEN}${BOLD}Passed (this run):${RESET}        $((PASS - SKIPPED))"
if [[ "$SKIPPED" -gt 0 ]]; then
    echo -e "  ${GREEN}${BOLD}Skipped (prior PASS):${RESET}     $SKIPPED"
fi

if [[ $FAIL -gt 0 ]]; then
    echo -e "  ${RED}${BOLD}Failed:${RESET}                   $FAIL"
    blank
    echo -e "  ${RED}Failed samples:${RESET}"
    for s in "${FAIL_LIST[@]}"; do
        echo -e "    ${RED}✘${RESET}  $s"
    done
fi

blank
echo -e "  ${BOLD}Results written to:${RESET}       $OUTPUT_DIR"
echo -e "  ${BOLD}Logs:${RESET}                     $LOG_DIR"
echo -e "  ${BOLD}Run manifest:${RESET}             $RUN_MANIFEST"
echo -e "  ${BOLD}QC summary:${RESET}               $QC_SUMMARY"
echo -e "  ${BOLD}R template sheet:${RESET}         $R_SAMPLE_SHEET"
echo -e "  ${BOLD}R autodetected sheet:${RESET}     $R_DETECTED_SAMPLE_SHEET"
echo -e "  ${BOLD}Command log:${RESET}              $COMMANDS_FILE"
blank

if [[ $FAIL -eq 0 ]]; then
    echo -e "  ${GREEN}${BOLD}All samples completed successfully.${RESET}"
else
    echo -e "  ${YELLOW}Some samples failed. Re-run with resume mode to retry them.${RESET}"
fi
blank

# A nonzero FAIL count must produce a nonzero exit status. Without this, a
# partially failed run looks identical to a fully successful one to any
# caller checking $? -- CI, workflow managers, or a wrapper script -- which
# can let a broken run silently pass as complete.
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
