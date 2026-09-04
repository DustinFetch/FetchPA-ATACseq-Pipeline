#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# PEPATAC EXPLORER
# Interactive, menu-driven exploration of a completed
# PEPATAC differential analysis run.
#
# Reads explorer_bundle.rds (written by PEPATAC_diff_analysis.sh)
# and lets you generate, on demand, in any order, repeatedly:
#   1. PCA plot           — any subset of detected samples
#   2. Re-plot contrast   — volcano + MA at new FDR/LFC thresholds
#   3. Peak boxplot       — accessibility at specific peaks (coords or index)
#   4. Sample correlation heatmap
#   5. Motif summary      — top known TF motifs for a contrast
#   6. GO/KEGG summary    — re-draw enrichment dot plot
#   7. Tornado plot       — deepTools heatmap centered on diff peaks
#                           (peak center/start/end via deepTools TSS/TES mode ·
#                           window · sort · Up/Down split) for a contrast
#   8. Composite tornado plot — same, but averaged into user-defined
#                           replicate groups with +/- SD shading on the trace
#
# Flow:
#   1. Point at a diff-analysis output folder (loads explorer_bundle.rds)
#   2. Choose the explorer output root folder
#   3. Load bundle, print samples / groups / contrasts
#   4. Menu loop — each analysis run gets its own organized subfolder
#
# Requires: PEPATAC_diff_analysis.sh v1.1-txdb-chipseeker-annotation or later
#           (earlier versions did not write explorer_bundle.rds)
# Best with: PEPATAC_diff_analysis.sh v1.3-modern-contrast-design or later
#           (earlier versions normalized once globally rather than per
#           contrast; this explorer still works against their bundles via
#           its whole-experiment/legacy-filename fallback tiers, just
#           without a contrast-specific normalization factor to prefer)
# ============================================================

SCRIPT_VERSION="0.6-annotation-integrity-audit-compatible"
RUN_ID="$(date '+%Y%m%d_%H%M%S')"
ENV_NAME="FetchPA"

# Tornado plots (analyses 8/9): deepTools computeMatrix writes its
# intermediate matrix.gz using a fixed "%.6f" text format regardless of
# value magnitude. PEPATAC's smoothShift bigWigs carry normalized ATAC
# signal in the ~1e-7..1e-9 range, which is below that 6-decimal floor --
# real signal gets silently serialized as "0.000000" and is unrecoverable
# once written. Confirmed directly against real data: --scale multiplies
# values BEFORE computeMatrix's own text serialization (not just a header
# annotation applied later), so every tornado matrix is computed pre-scaled
# by this constant and PLOTTED IN THOSE SCALED UNITS PERMANENTLY -- values
# are never divided back down. Every plot title/axis/output that shows
# signal must say so (see run_tornado()). This is simpler and less
# failure-prone than scaling up, writing, then unscaling and rewriting a
# second matrix file; the only cost is that plotted values read as
# "smoothShift signal x 1,000,000" instead of raw units, which is called
# out everywhere it appears.
TORNADO_SIGNAL_SCALE=1000000

# explorer_bundle.rds's schema_version this explorer was written against and
# verified to understand. Checked (not just read) right after the bundle
# loads -- see the SCHEMA_VERSION compatibility check in Step 3.
EXPLORER_SUPPORTED_SCHEMA=3

# ─────────────────────────────────────────────────────────────
# Terminal colors / formatting  (identical to diff_analysis.sh)
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

label()  { echo -e "  ${BOLD}${BLUE}▸${RESET} $*"; }
ok()     { echo -e "  ${GREEN}✔${RESET}  $*"; }
warn()   { echo -e "  ${YELLOW}⚠${RESET}  $*"; }
err()    { echo -e "  ${RED}✘${RESET}  $*"; }
die()    { err "$*"; exit 1; }
blank()  { echo ""; }

# ─────────────────────────────────────────────────────────────
# Safe Bash -> R string embedding
# (ported verbatim from PEPATAC_diff_analysis.sh, which fixed the same
# class of bug: raw "$VAR" interpolation into generated R source breaks
# or corrupts the script if the value contains a `"` or `\`. Every value
# embedded into an R heredoc below goes through one of these instead.)
# Defined here, before Step 1, because BUNDLE_PATH_R is computed at
# top level as soon as the bundle is found -- a bash function must be
# defined before any top-level call reaches it, not just before the
# call syntactically appears later in the file.
# ─────────────────────────────────────────────────────────────

# url_encode STRING
# Percent-encodes STRING to [A-Za-z0-9._~-] plus %XX escapes, so the
# encoded output can never contain a `"`, `\`, or anything else with
# meaning in R syntax. Decoded on the R side by .pepatac_url_decode().
url_encode() {
    local LC_ALL=C
    local string="${1:-}"
    local strlen=${#string} pos c hex
    local encoded=""
    for (( pos = 0; pos < strlen; pos++ )); do
        c="${string:pos:1}"
        case "$c" in
            [a-zA-Z0-9.~_-]) encoded+="$c" ;;
            *) printf -v hex '%%%02X' "'$c"
               encoded+="$hex" ;;
        esac
    done
    printf '%s' "$encoded"
}

# r_string_literal VALUE
# Emits an R expression that safely evaluates to VALUE as a character
# scalar: .pepatac_url_decode("<percent-encoded VALUE>").
r_string_literal() {
    printf '.pepatac_url_decode("%s")' "$(url_encode "$1")"
}

# r_vector_literal VALUE...
# Emits an R expression that safely evaluates to a character vector of
# the given values, in order: .pepatac_decode_vec(c("enc1", "enc2", ...)).
# Empty input emits character(0) rather than an empty c().
r_vector_literal() {
    if [[ $# -eq 0 ]]; then
        printf '.pepatac_decode_vec(character(0))'
        return
    fi
    local out="" first=true enc
    for val in "$@"; do
        enc="$(url_encode "$val")"
        if $first; then out="\"$enc\""; first=false
        else out="$out, \"$enc\""; fi
    done
    printf '.pepatac_decode_vec(c(%s))' "$out"
}

# ─────────────────────────────────────────────────────────────
# Path normalisation  (shared with run.sh / diff_analysis.sh / install.sh)
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

in_env() {
    conda run --no-capture-output -n "$ENV_NAME" "$@"
}

in_env_clean() {
    conda run --no-capture-output -n "$ENV_NAME" \
        env -u R_ARCH -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$@"
}

source_conda_or_die() {
    local conda_sh="$HOME/miniconda3/etc/profile.d/conda.sh"
    if [[ -f "$conda_sh" ]]; then
        # shellcheck source=/dev/null
        source "$conda_sh"
    else
        die "Conda not found at $conda_sh. Run PEPATAC_install.sh first."
    fi
}

# ─────────────────────────────────────────────────────────────
# Abort cleanup
# ─────────────────────────────────────────────────────────────

ABORTING=false

abort_cleanup() {
    local signal="${1:-INT}"
    if $ABORTING; then exit 130; fi
    ABORTING=true
    trap - INT TERM HUP
    blank
    err "Abort requested ($signal). Stopping child processes..."
    pkill -TERM -P $$ >/dev/null 2>&1 || true
    sleep 2
    pkill -KILL -P $$ >/dev/null 2>&1 || true
    err "Explorer aborted."
    exit 130
}

trap 'abort_cleanup INT'  INT
trap 'abort_cleanup TERM' TERM
trap 'abort_cleanup HUP'  HUP

# ─────────────────────────────────────────────────────────────
# Banner
# ─────────────────────────────────────────────────────────────

clear 2>/dev/null || true
echo ""
echo -e "${BOLD}${MAGENTA}"
cat << 'BANNER'
  ██████╗ ██╗███████╗███████╗
  ██╔══██╗██║██╔════╝██╔════╝
  ██║  ██║██║█████╗  █████╗
  ██║  ██║██║██╔══╝  ██╔══╝
  ██████╔╝██║██║     ██║
  ╚═════╝ ╚═╝╚═╝     ╚═╝
  Interactive Explorer
BANNER
echo -e "${RESET}"
echo -e "  ${DIM}PEPATAC Interactive Explorer  ·  PCA · Re-plot · Boxplot · Motifs · GO · Tornado${RESET}"
echo ""

# ─────────────────────────────────────────────────────────────
# STEP 1 — Locate diff-analysis output folder
# ─────────────────────────────────────────────────────────────

header "Step 1 · Diff-Analysis Output Folder"

echo -e "  Point this at the output folder from a completed PEPATAC_diff_analysis.sh run."
echo -e "  ${DIM}Looking for: explorer_bundle.rds${RESET}"
blank

while true; do
    IFS= read -r -e -p "  Diff-analysis output folder: " ANALYSIS_DIR
    normalize_path_var ANALYSIS_DIR
    ANALYSIS_DIR="${ANALYSIS_DIR%/}"

    if [[ -z "$ANALYSIS_DIR" ]]; then
        err "Please enter a path."
        continue
    fi
    if [[ ! -d "$ANALYSIS_DIR" ]]; then
        err "Folder not found: $ANALYSIS_DIR"
        continue
    fi

    BUNDLE_PATH="$ANALYSIS_DIR/explorer_bundle.rds"
    if [[ ! -f "$BUNDLE_PATH" ]]; then
        err "No explorer_bundle.rds found in that folder."
        warn "This file is written by PEPATAC_diff_analysis.sh (v0.7-motif-go-analysis or later)."
        warn "If this folder came from an older version, re-run the diff analysis to generate it."
        continue
    fi

    ok "Found bundle: $BUNDLE_PATH"
    break
done

# Computed once, used by every analysis below wherever the bundle path is
# embedded into a generated R script (see url_encode/r_string_literal
# above for why this is needed instead of raw "$BUNDLE_PATH").
BUNDLE_PATH_R="$(r_string_literal "$BUNDLE_PATH")"

# ─────────────────────────────────────────────────────────────
# STEP 2 — Output folder for explorer plots
# ─────────────────────────────────────────────────────────────

header "Step 2 · Output Folder"

echo -e "  Where should explorer plots be saved?"
echo -e "  ${DIM}A timestamped subfolder of the diff output is the default.${RESET}"
blank

DEFAULT_EXPLORE_OUT="$ANALYSIS_DIR/explorer_$RUN_ID"
echo -e "  ${DIM}Default: $DEFAULT_EXPLORE_OUT${RESET}"
blank

IFS= read -r -e -p "  Output folder [press Enter for default]: " EXPLORE_OUT
normalize_path_var EXPLORE_OUT
[[ -z "$EXPLORE_OUT" ]] && EXPLORE_OUT="$DEFAULT_EXPLORE_OUT"
mkdir -p "$EXPLORE_OUT"
ok "Explorer output root: $EXPLORE_OUT"

# Keep session-level helper files out of the visible analysis folders.
EXPLORE_TMP="$EXPLORE_OUT/.tmp"
mkdir -p "$EXPLORE_TMP"

# ─────────────────────────────────────────────────────────────
# STEP 3 — Load bundle metadata (via a fast R snippet)
# ─────────────────────────────────────────────────────────────

header "Step 3 · Loading Bundle"

source_conda_or_die

INSPECT_R="$EXPLORE_TMP/.inspect_bundle.R"
INSPECT_TSV="$EXPLORE_TMP/.bundle_inventory.tsv"

cat > "$INSPECT_R" << 'RCHECK'
b <- readRDS(commandArgs(trailingOnly=TRUE)[1])
tsv_path <- commandArgs(trailingOnly=TRUE)[2]

# Write sample/condition table.
con <- file(tsv_path, "w")
writeLines("SAMPLE\tCONDITION", con)
if (!is.null(b$col_data) && nrow(b$col_data) > 0) {
    for (i in seq_len(nrow(b$col_data))) {
        writeLines(paste(b$col_data$SampleID[i],
                         b$col_data$Condition[i], sep="\t"), con)
    }
}
close(con)

cat("OK\n")
cat("SCHEMA_VERSION=",    if (!is.null(b$schema_version)) b$schema_version else 0, "\n", sep="")
cat("N_SAMPLES=",         if (!is.null(b$col_data)) nrow(b$col_data) else 0, "\n", sep="")
cat("N_CONTRASTS=",       length(b$contrast_labels), "\n", sep="")
cat("HAS_VST=",           !is.null(b$vst_matrix), "\n", sep="")
cat("HAS_ORGDB=",         !is.null(b$orgdb), "\n", sep="")
cat("HAS_GO=",            if (!is.null(b$has_go)) isTRUE(b$has_go) else !is.null(b$orgdb), "\n", sep="")
cat("HAS_KEGG=",          if (!is.null(b$has_kegg)) isTRUE(b$has_kegg) else !is.null(b$orgdb), "\n", sep="")
cat("HAS_TXDB=",          !is.null(b$txdb), "\n", sep="")
cat("GENOME=",            if (!is.null(b$genome)) b$genome else "unknown", "\n", sep="")
cat("FDR_CUTOFF=",        if (!is.null(b$fdr_cutoff)) b$fdr_cutoff else 0.05, "\n", sep="")
cat("FC_CUTOFF=",         if (!is.null(b$fc_cutoff)) b$fc_cutoff else 0.585, "\n", sep="")
cat("ORGDB=",             if (!is.null(b$orgdb)) b$orgdb else "none", "\n", sep="")
cat("TXDB=",              if (!is.null(b$txdb)) b$txdb else "none", "\n", sep="")
cat("TXDB_KEYTYPE=",      if (!is.null(b$txdb_gene_keytype)) b$txdb_gene_keytype else "unknown", "\n", sep="")
cat("HAS_HOMER=",         local({
    # Whether HOMER produced usable, already-saved results -- not whether
    # the HOMER binary itself is still installed at its recorded path.
    # run_motif_summary() only ever reads this pre-computed CSV; HOMER
    # never needs to re-run for the explorer to show it. Checking the
    # binary path instead would hide a fully usable feature after moving
    # the explorer to a machine without HOMER, or after an env rebuild.
    if (is.null(b$diff_out) || is.null(b$contrast_labels) || length(b$contrast_labels) == 0) {
        FALSE
    } else {
        any(vapply(b$contrast_labels, function(lbl) {
            file.exists(file.path(b$diff_out, lbl, "motifs",
                                  paste0(lbl, "_known_motifs_top50.csv")))
        }, logical(1)))
    }
}), "\n", sep="")
cat("CONTRASTS=",         paste(b$contrast_labels, collapse="|"), "\n", sep="")
cat("GENERATED=",         if (!is.null(b$generated)) b$generated else "unknown", "\n", sep="")
cat("SCRIPT_VERSION=",    if (!is.null(b$script_version)) b$script_version else "unknown", "\n", sep="")
RCHECK

label "Reading metadata from bundle..."
INSPECT_OUT=$(in_env_clean Rscript --vanilla "$INSPECT_R" \
    "$BUNDLE_PATH" "$INSPECT_TSV" 2>&1) || {
    err "Failed to read explorer_bundle.rds. Output:"
    echo "$INSPECT_OUT"
    die "Cannot continue without a readable bundle."
}

# Parse key-value output lines.
_get() { echo "$INSPECT_OUT" | grep "^${1}=" | head -1 | cut -d= -f2-; }

BUNDLE_SCHEMA=$(_get SCHEMA_VERSION)

# Compatibility check -- previously this value was read into a variable and
# never looked at again. Asymmetric on purpose: older or missing schema
# degrades gracefully (this script already guards most individual fields
# defensively -- HAS_VST, HAS_ORGDB, etc. -- so a missing field reports
# unavailable per-feature rather than failing here). A NEWER schema is a
# different situation: a field could have changed shape or meaning in a
# way this explorer has no way to detect, so silently reading it and
# hoping for the best is the riskier direction. That gets a hard stop
# instead of a warning that's easy to miss.
if [[ -z "$BUNDLE_SCHEMA" || "$BUNDLE_SCHEMA" == "0" ]]; then
    warn "Bundle has no schema_version (written by an older diff_analysis.sh)."
    warn "Some fields this explorer expects may be missing; features that need"
    warn "them will report unavailable individually rather than failing here."
elif [[ "$BUNDLE_SCHEMA" -gt "$EXPLORER_SUPPORTED_SCHEMA" ]]; then
    err "Bundle schema_version ($BUNDLE_SCHEMA) is newer than this explorer"
    err "understands (supports up to $EXPLORER_SUPPORTED_SCHEMA)."
    err "Fields may have changed shape or meaning in a way this explorer"
    err "cannot detect -- continuing could silently misinterpret them."
    die "Update PEPATAC_explore_organized_outputs.sh to a version that supports schema $BUNDLE_SCHEMA, or re-run against a bundle written by a matching diff_analysis.sh version."
elif [[ "$BUNDLE_SCHEMA" -lt "$EXPLORER_SUPPORTED_SCHEMA" ]]; then
    label "Bundle schema_version ($BUNDLE_SCHEMA) is older than current ($EXPLORER_SUPPORTED_SCHEMA) -- should still work; older bundles may just be missing newer fields."
fi

N_SAMPLES=$(_get N_SAMPLES)
N_CONTRASTS=$(_get N_CONTRASTS)
HAS_VST=$(_get HAS_VST)
HAS_ORGDB=$(_get HAS_ORGDB)
HAS_GO=$(_get HAS_GO)
HAS_KEGG=$(_get HAS_KEGG)
HAS_TXDB=$(_get HAS_TXDB)
HAS_HOMER=$(_get HAS_HOMER)
BUNDLE_GENOME=$(_get GENOME)
BUNDLE_FDR=$(_get FDR_CUTOFF)
BUNDLE_FC=$(_get FC_CUTOFF)
BUNDLE_ORGDB=$(_get ORGDB)
BUNDLE_TXDB=$(_get TXDB)
BUNDLE_TXDB_KEYTYPE=$(_get TXDB_KEYTYPE)
BUNDLE_GENERATED=$(_get GENERATED)
BUNDLE_SCRIPT_VER=$(_get SCRIPT_VERSION)
CONTRAST_LIST_RAW=$(_get CONTRASTS)   # pipe-separated

# Parse samples from TSV.
declare -a SAMPLE_NAMES=()
declare -a SAMPLE_CONDITIONS=()
declare -A SAMPLE_TO_CONDITION

if [[ -s "$INSPECT_TSV" ]]; then
    {
        read -r _header
        while IFS=$'\t' read -r s c; do
            [[ -z "$s" ]] && continue
            SAMPLE_NAMES+=("$s")
            SAMPLE_CONDITIONS+=("$c")
            SAMPLE_TO_CONDITION["$s"]="$c"
        done
    } < "$INSPECT_TSV"
fi

N_SAMPLES_LOADED=${#SAMPLE_NAMES[@]}

# Build ordered unique group list.
declare -a GROUP_LIST=()
declare -A GROUP_SEEN
for c in "${SAMPLE_CONDITIONS[@]:-}"; do
    if [[ -n "$c" && -z "${GROUP_SEEN[$c]:-}" ]]; then
        GROUP_LIST+=("$c")
        GROUP_SEEN["$c"]=1
    fi
done

# Parse contrast list.
declare -a CONTRAST_NAMES=()
if [[ -n "$CONTRAST_LIST_RAW" ]]; then
    IFS='|' read -ra CONTRAST_NAMES <<< "$CONTRAST_LIST_RAW"
fi

# Print bundle summary.
blank
ok "Bundle loaded from: $ANALYSIS_DIR"
ok "Generated:  $BUNDLE_GENERATED  (diff script $BUNDLE_SCRIPT_VER)"
ok "Genome:     $BUNDLE_GENOME"
if [[ "$HAS_TXDB" == "TRUE" ]]; then
    ok "Coordinates:  $BUNDLE_TXDB"
    ok "Gene mapping: $BUNDLE_ORGDB ($BUNDLE_TXDB_KEYTYPE)"
fi
ok "Original thresholds:  FDR < $BUNDLE_FDR  |  |log2FC| ≥ $BUNDLE_FC"
blank

echo -e "  ${BOLD}Samples (${N_SAMPLES_LOADED}):${RESET}"
for (( i=0; i<${#SAMPLE_NAMES[@]}; i++ )); do
    printf "    ${CYAN}%3d${RESET}  %-34s  ${DIM}(%s)${RESET}\n" \
        "$((i+1))" "${SAMPLE_NAMES[$i]}" "${SAMPLE_CONDITIONS[$i]}"
done
blank

echo -e "  ${BOLD}Contrasts (${#CONTRAST_NAMES[@]}):${RESET}"
for (( i=0; i<${#CONTRAST_NAMES[@]}; i++ )); do
    printf "    ${CYAN}%3d${RESET}  %s\n" "$((i+1))" "${CONTRAST_NAMES[$i]}"
done
blank

if [[ "$HAS_VST" != "TRUE" ]]; then
    warn "Bundle has no VST matrix — PCA and correlation heatmap will be unavailable."
    warn "(Original PCA step was skipped, likely < 3 samples.)"
fi
if [[ "$HAS_HOMER" != "TRUE" ]]; then
    warn "No saved HOMER motif results found in any contrast folder."
    warn "Motif summary will be unavailable. Menu item 5 will be skipped."
fi
if [[ "$HAS_GO" != "TRUE" ]]; then
    warn "No saved GO enrichment results were found in this bundle."
fi
if [[ "$HAS_KEGG" != "TRUE" ]]; then
    warn "No saved KEGG enrichment results were found; KEGG choices will be hidden."
fi
if [[ "$HAS_TXDB" != "TRUE" ]]; then
    warn "No TxDb provenance recorded; this appears to be an older analysis bundle."
fi

# ─────────────────────────────────────────────────────────────
# Shared bash helpers
# ─────────────────────────────────────────────────────────────

# resolve_contrast INPUT -> sets RESOLVED_CONTRAST or returns 1
resolve_contrast() {
    local inp="$1"
    RESOLVED_CONTRAST=""
    if [[ "$inp" =~ ^[0-9]+$ ]]; then
        local idx=$(( inp - 1 ))
        if [[ "$idx" -lt 0 || "$idx" -ge "${#CONTRAST_NAMES[@]}" ]]; then
            err "Number out of range (1–${#CONTRAST_NAMES[@]})."
            return 1
        fi
        RESOLVED_CONTRAST="${CONTRAST_NAMES[$idx]}"
    else
        for c in "${CONTRAST_NAMES[@]}"; do
            if [[ "$c" == "$inp" ]]; then
                RESOLVED_CONTRAST="$c"
                break
            fi
        done
        if [[ -z "$RESOLVED_CONTRAST" ]]; then
            err "Contrast not found: $inp"
            return 1
        fi
    fi
    return 0
}

show_contrast_menu() {
    echo -e "  Available contrasts:"
    for (( i=0; i<${#CONTRAST_NAMES[@]}; i++ )); do
        printf "    ${CYAN}%3d${RESET}  %s\n" "$((i+1))" "${CONTRAST_NAMES[$i]}"
    done
}

show_sample_menu() {
    echo -e "  Samples:"
    for (( i=0; i<${#SAMPLE_NAMES[@]}; i++ )); do
        printf "    ${CYAN}%3d${RESET}  %-34s  ${DIM}(%s)${RESET}\n" \
            "$((i+1))" "${SAMPLE_NAMES[$i]}" "${SAMPLE_CONDITIONS[$i]}"
    done
    blank
    echo -e "  Groups:"
    for g in "${GROUP_LIST[@]}"; do
        local cnt=0
        for c in "${SAMPLE_CONDITIONS[@]}"; do
            [[ "$c" == "$g" ]] && (( cnt++ )) || true
        done
        printf "    ${CYAN}%-14s${RESET}  (%d samples)\n" "$g" "$cnt"
    done
}

# prompt_sample_selection PROMPT MIN_COUNT -> sets RESOLVED_SAMPLES array
prompt_sample_selection() {
    local prompt_text="$1"
    local min_count="${2:-1}"
    RESOLVED_SAMPLES=()
    while true; do
        blank
        show_sample_menu
        blank
        echo -e "  ${DIM}Enter numbers, sample names, or group names (space-separated).${RESET}"
        echo -e "  ${DIM}Type 'all' for every sample.${RESET}"
        blank
        read -p "  $prompt_text " SEL_INPUT

        if [[ -z "$SEL_INPUT" ]]; then
            err "No input given."
            continue
        fi

        if [[ "${SEL_INPUT,,}" == "all" ]]; then
            RESOLVED_SAMPLES=("${SAMPLE_NAMES[@]}")
        else
            RESOLVED_SAMPLES=()
            local ok_flag=true
            IFS=' ' read -ra _tokens <<< "$SEL_INPUT"
            for tok in "${_tokens[@]}"; do
                if [[ "$tok" =~ ^[0-9]+$ ]]; then
                    local idx=$(( tok - 1 ))
                    if [[ "$idx" -lt 0 || "$idx" -ge "${#SAMPLE_NAMES[@]}" ]]; then
                        err "Sample number $tok out of range."
                        ok_flag=false; break
                    fi
                    RESOLVED_SAMPLES+=("${SAMPLE_NAMES[$idx]}")
                else
                    # Try exact sample name.
                    local found_sample=""
                    for s in "${SAMPLE_NAMES[@]}"; do
                        [[ "$s" == "$tok" ]] && found_sample="$s" && break
                    done
                    if [[ -n "$found_sample" ]]; then
                        RESOLVED_SAMPLES+=("$found_sample")
                        continue
                    fi
                    # Try group name — expand to all members.
                    local found_group=false
                    for g in "${GROUP_LIST[@]}"; do
                        if [[ "$g" == "$tok" ]]; then
                            found_group=true
                            for (( ii=0; ii<${#SAMPLE_NAMES[@]}; ii++ )); do
                                [[ "${SAMPLE_CONDITIONS[$ii]}" == "$g" ]] && \
                                    RESOLVED_SAMPLES+=("${SAMPLE_NAMES[$ii]}")
                            done
                            break
                        fi
                    done
                    if ! $found_group && [[ -z "$found_sample" ]]; then
                        err "Not recognized as a sample number, name, or group: $tok"
                        ok_flag=false; break
                    fi
                fi
            done
            $ok_flag || continue
            # De-duplicate, preserve order.
            local -a deduped=()
            local -A _seen_s
            for s in "${RESOLVED_SAMPLES[@]}"; do
                if [[ -z "${_seen_s[$s]:-}" ]]; then
                    deduped+=("$s")
                    _seen_s["$s"]=1
                fi
            done
            RESOLVED_SAMPLES=("${deduped[@]}")
        fi

        if [[ "${#RESOLVED_SAMPLES[@]}" -lt "$min_count" ]]; then
            err "Need at least $min_count sample(s); got ${#RESOLVED_SAMPLES[@]}."
            continue
        fi

        ok "Selected ${#RESOLVED_SAMPLES[@]} sample(s): ${RESOLVED_SAMPLES[*]}"
        break
    done
}

run_r_script() {
    local r_script="$1"
    local log_file="$2"
    set +e
    in_env_clean Rscript --vanilla "$r_script" 2>&1 | tee "$log_file"
    local r_exit=${PIPESTATUS[0]}
    set -e
    return "$r_exit"
}

# ─────────────────────────────────────────────────────────────
# Analysis 1 — PCA
# ─────────────────────────────────────────────────────────────

run_pca_analysis() {
    header "PCA Plot"

    if [[ "$HAS_VST" != "TRUE" ]]; then
        err "Bundle has no VST matrix — PCA is unavailable."
        err "(Original diff run had < 3 samples or the VST step failed.)"
        return
    fi

    echo -e "  Choose which samples to include."
    echo -e "  ${DIM}You can use any subset — e.g. exclude an outlier to see the rest more clearly.${RESET}"
    prompt_sample_selection "Samples for PCA:" 3

    blank
    read -p "  Top N most variable peaks to use [default: 500, 'all' for every peak]: " PCA_NTOP
    PCA_NTOP="${PCA_NTOP:-500}"
    if [[ "${PCA_NTOP,,}" == "all" ]]; then
        PCA_NTOP_R="Inf"
    elif [[ "$PCA_NTOP" =~ ^[0-9]+$ ]]; then
        PCA_NTOP_R="$PCA_NTOP"
    else
        warn "Invalid — using default of 500."
        PCA_NTOP_R="500"
        PCA_NTOP="500"
    fi

    read -p "  Label points with sample names? [Y/n]: " PCA_LABELS
    PCA_LABEL_R=$([[ "${PCA_LABELS,,}" == "n" ]] && echo "FALSE" || echo "TRUE")

    local TAG="pca_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/01_PCA/PCA_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_PDF="$ANALYSIS_OUT_DIR/PCA_${TS}.pdf"
    local OUT_PNG="${OUT_PDF%.pdf}.png"
    local OUT_CSV="${OUT_PDF%.pdf}_coordinates.csv"

    local R_SAMPLES
    R_SAMPLES=$(r_vector_literal "${RESOLVED_SAMPLES[@]}")
    local OUT_CSV_R OUT_PDF_R OUT_PNG_R
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"
    OUT_PDF_R="$(r_string_literal "$OUT_PDF")"
    OUT_PNG_R="$(r_string_literal "$OUT_PNG")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(matrixStats)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b           <- readRDS(${BUNDLE_PATH_R})
vst_mat     <- b\$vst_matrix
sel_samples <- $R_SAMPLES
ntop        <- $PCA_NTOP_R
label_pts   <- $PCA_LABEL_R

avail <- intersect(sel_samples, colnames(vst_mat))
if (length(avail) < length(sel_samples)) {
    missing <- setdiff(sel_samples, colnames(vst_mat))
    warning("Samples not in VST matrix, skipping: ", paste(missing, collapse=", "))
}
if (length(avail) < 3) stop("Need at least 3 samples in the VST matrix to run PCA.")

mat <- vst_mat[, avail, drop=FALSE]
rv  <- matrixStats::rowVars(mat)
ntop_use <- if (is.infinite(ntop)) nrow(mat) else min(ntop, nrow(mat))
selected <- order(rv, decreasing=TRUE)[seq_len(ntop_use)]
mat      <- mat[selected, , drop=FALSE]

pca      <- prcomp(t(mat), scale.=FALSE)
pct_var  <- round(100 * pca\$sdev^2 / sum(pca\$sdev^2), 1)

conditions <- if (!is.null(b\$col_data)) {
    setNames(as.character(b\$col_data\$Condition), b\$col_data\$SampleID)
} else {
    setNames(rep("unknown", length(avail)), avail)
}

pca_df <- data.frame(
    PC1       = pca\$x[, 1],
    PC2       = pca\$x[, 2],
    sample    = avail,
    condition = conditions[avail],
    stringsAsFactors = FALSE
)

write.csv(
    cbind(pca_df, pca\$x[, seq_len(min(5, ncol(pca\$x))), drop=FALSE]),
    ${OUT_CSV_R}, row.names=FALSE
)

p <- ggplot(pca_df, aes(PC1, PC2, color=condition, label=sample)) +
    geom_point(size=4, alpha=0.85) +
    theme_bw(base_size=14) +
    labs(
        x      = paste0("PC1 (", pct_var[1], "% variance)"),
        y      = paste0("PC2 (", pct_var[2], "% variance)"),
        color  = "Group",
        title  = "PCA - ATAC-seq Consensus Peaks",
        subtitle = paste0(length(avail), " samples, top ", ntop_use, " variable peaks (DESeq2 VST)")
    )

if (label_pts) {
    p <- p + ggrepel::geom_text_repel(size=3.5, show.legend=FALSE, max.overlaps=20)
}

ggsave(${OUT_PDF_R}, plot=p, width=9, height=7)
ggsave(${OUT_PNG_R}, plot=p, width=9, height=7, dpi=150)
cat("PCA plot saved:", ${OUT_PDF_R}, "\n")
cat("PCA coordinates saved:", ${OUT_CSV_R}, "\n")
RSCRIPT_EOF

    label "Generating PCA plot..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "PCA plot:        $OUT_PDF"
        ok "PCA PNG:         $OUT_PNG"
        ok "PCA coordinates: $OUT_CSV"
    else
        err "PCA failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 2 — Re-plot contrast (volcano + MA at new thresholds)
# ─────────────────────────────────────────────────────────────

run_replot_analysis() {
    header "Re-plot Contrast · Volcano + MA"

    if [[ "${#CONTRAST_NAMES[@]}" -eq 0 ]]; then
        err "No contrasts found in this bundle."
        return
    fi

    blank
    show_contrast_menu
    blank

    while true; do
        read -p "  Choose contrast (name or number): " CINPUT
        resolve_contrast "$CINPUT" && break
    done

    ok "Contrast: $RESOLVED_CONTRAST"
    blank

    echo -e "  Effect-size estimate to use for this re-plot:"
    echo -e "    ${CYAN}1${RESET}.  DiffBind shrunken Fold ${DIM}(default)${RESET}"
    echo -e "        ${DIM}Stabilized estimate -- good for ranking peaks or general${RESET}"
    echo -e "        ${DIM}exploration. Shrinkage strength depends on the same uncertainty${RESET}"
    echo -e "        ${DIM}that drives significance, which can make effect size and${RESET}"
    echo -e "        ${DIM}significance look more tightly related than they really are.${RESET}"
    echo -e "    ${CYAN}2${RESET}.  DESeq2 MLE ${DIM}(unshrunken)${RESET}"
    echo -e "        ${DIM}The actual coefficient the significance test evaluates. Noisier,${RESET}"
    echo -e "        ${DIM}but the conventional volcano-plot choice, and avoids conflating a${RESET}"
    echo -e "        ${DIM}regularized effect size with an unregularized test statistic.${RESET}"
    echo -e "        ${DIM}Recommended for publication figures. Requires the diff_analysis.sh${RESET}"
    echo -e "        ${DIM}run to have retrieved it (falls back to option 1 if unavailable).${RESET}"
    blank
    read -p "  Choice [1/2, default 1]: " LFC_CHOICE_INPUT
    LFC_CHOICE_INPUT="${LFC_CHOICE_INPUT:-1}"
    local USE_MLE_LFC=false
    case "$LFC_CHOICE_INPUT" in
        2) USE_MLE_LFC=true;  ok "Using DESeq2 MLE (unshrunken) log2FoldChange." ;;
        *) USE_MLE_LFC=false; ok "Using DiffBind's shrunken Fold." ;;
    esac
    local USE_MLE_LFC_R
    USE_MLE_LFC_R=$($USE_MLE_LFC && echo "TRUE" || echo "FALSE")

    echo -e "  Original run thresholds:  FDR < ${BUNDLE_FDR}  |  |log2FC| ≥ ${BUNDLE_FC}"
    echo -e "  ${DIM}Enter new values to override, or press Enter to keep the originals.${RESET}"
    blank

    read -p "  FDR cutoff [default: $BUNDLE_FDR]: " NEW_FDR
    NEW_FDR="${NEW_FDR:-$BUNDLE_FDR}"
    [[ "$NEW_FDR" =~ ^[0-9]*\.?[0-9]+$ ]] || { warn "Invalid — using $BUNDLE_FDR."; NEW_FDR="$BUNDLE_FDR"; }

    read -p "  |log2FC| cutoff [default: $BUNDLE_FC]: " NEW_FC
    NEW_FC="${NEW_FC:-$BUNDLE_FC}"
    [[ "$NEW_FC" =~ ^[0-9]*\.?[0-9]+$ ]] || { warn "Invalid — using $BUNDLE_FC."; NEW_FC="$BUNDLE_FC"; }

    ok "Using FDR < $NEW_FDR  |  |log2FC| ≥ $NEW_FC"

    read -p "  Volcano Y-axis max (-log10 padj) [default: automatic]: " NEW_Y_MAX
    if [[ -n "$NEW_Y_MAX" ]]; then
        if [[ "$NEW_Y_MAX" =~ ^[0-9]*\.?[0-9]+$ ]] && awk "BEGIN{exit !($NEW_Y_MAX > 0)}"; then
            ok "Y-axis capped at $NEW_Y_MAX"
        else
            warn "Invalid — using automatic Y-axis scaling."
            NEW_Y_MAX=""
        fi
    fi
    local Y_MAX_R
    Y_MAX_R="${NEW_Y_MAX:-NA}"

    read -p "  Volcano X-axis max (|log2FC|) [default: automatic, symmetric]: " NEW_X_MAX
    if [[ -n "$NEW_X_MAX" ]]; then
        if [[ "$NEW_X_MAX" =~ ^[0-9]*\.?[0-9]+$ ]] && awk "BEGIN{exit !($NEW_X_MAX > 0)}"; then
            ok "X-axis capped at ±$NEW_X_MAX"
        else
            warn "Invalid — using automatic (symmetric) X-axis scaling."
            NEW_X_MAX=""
        fi
    fi
    local X_MAX_R
    X_MAX_R="${NEW_X_MAX:-NA}"

    blank
    echo -e "  Optionally label the most significant peaks on the volcano plot by"
    echo -e "  coordinate (leave blank to skip)."
    echo -e "  ${DIM}Ranked by significance (adjusted p-value), separately for up- and${RESET}"
    echo -e "  ${DIM}down-regulated peaks.${RESET}"
    blank
    read -p "  Label top N most significant up/down peaks [default: 0]: " TOPN_INPUT
    TOPN_INPUT="${TOPN_INPUT:-0}"
    local TOPN_LABEL
    if [[ "$TOPN_INPUT" =~ ^[0-9]+$ ]]; then
        TOPN_LABEL="$TOPN_INPUT"
    else
        warn "Not a valid number, using 0 (off)."
        TOPN_LABEL="0"
    fi
    [[ "$TOPN_LABEL" -gt 0 ]] && ok "Will label the top $TOPN_LABEL most significant up peaks and top $TOPN_LABEL down peaks."

    local TAG="replot_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/02_Contrast_Replots/${SAFE_LABEL}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    # The mle/shrunken filename tag is NOT decided here in bash -- only R,
    # after checking whether this bundle actually has MLE_log2FoldChange for
    # this contrast, knows whether the requested estimate is available or
    # falls back to shrunken. Building the tag here (as an earlier version
    # of this feature did) meant a requested-but-unavailable MLE plot got
    # saved with "_mle_" in the filename despite silently containing the
    # shrunken data -- exactly the kind of mislabeling this whole feature
    # exists to avoid. R constructs out_vol_pdf/out_vol_png/out_ma_pdf/
    # out_ma_png itself below, after the fallback check, so the name on
    # disk always matches what was actually plotted.
    local OUT_CSV="$ANALYSIS_OUT_DIR/${SAFE_LABEL}_replot_${TS}.csv"

    local RESOLVED_CONTRAST_R OUT_CSV_R OUT_DIR_R SAFE_LABEL_R TS_TAG_R FDR_TAG_R FC_TAG_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"
    OUT_DIR_R="$(r_string_literal "$ANALYSIS_OUT_DIR")"
    SAFE_LABEL_R="$(r_string_literal "$SAFE_LABEL")"
    TS_TAG_R="$(r_string_literal "$TS")"
    FDR_TAG_R="$(r_string_literal "$NEW_FDR")"
    FC_TAG_R="$(r_string_literal "$NEW_FC")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b           <- readRDS(${BUNDLE_PATH_R})
label       <- ${RESOLVED_CONTRAST_R}
fdr_cutoff  <- $NEW_FDR
fc_cutoff   <- $NEW_FC
y_max_override <- $Y_MAX_R   # NA = automatic; user-set value caps the volcano's Y axis
x_max_override <- $X_MAX_R   # NA = automatic (symmetric to data); user-set value caps |log2FC|
top_n_label    <- $TOPN_LABEL   # 0 = no labels; else label top N most significant peaks per direction

res <- b\$all_results[[label]]
if (is.null(res)) stop("Contrast '", label, "' not found in bundle.")

# Which effect-size estimate to plot/threshold on -- chosen interactively
# above. DiffBind's own Fold column may be apeglm/ashr-shrunk; MLE_* columns
# (if diff_analysis.sh successfully retrieved them for this contrast -- see
# that script's per-contrast analysis loop) hold the actual unshrunken
# DESeq2 coefficient the significance test evaluates. Substituting into
# res\$log2FoldChange in place, rather than branching the rest of this
# script, means every downstream step (thresholding, plotting, the CSV
# export) automatically uses whichever estimate was chosen with no special-
# casing needed past this point.
use_mle_lfc <- $USE_MLE_LFC_R
lfc_source_label <- "DiffBind shrunken Fold"
if (use_mle_lfc) {
    if ("MLE_log2FoldChange" %in% names(res)) {
        res\$log2FoldChange <- res\$MLE_log2FoldChange
        lfc_source_label <- "DESeq2 MLE (unshrunken)"
        cat("Using unshrunken MLE log2FoldChange for this re-plot.\n")
    } else {
        cat("[WARN] MLE_log2FoldChange not found in this bundle (older diff_analysis.sh run,\n")
        cat("       or extraction failed for this contrast) -- falling back to DiffBind's\n")
        cat("       shrunken Fold.\n")
        use_mle_lfc <- FALSE
    }
}

# Filenames are built HERE, not in bash, specifically because use_mle_lfc
# has just been finalized (post-fallback-check) -- building them in bash
# before this point could only ever encode what was requested, not what
# actually happened.
lfc_tag <- if (use_mle_lfc) "mle" else "shrunken"
out_dir    <- ${OUT_DIR_R}
safe_label <- ${SAFE_LABEL_R}
ts_tag     <- ${TS_TAG_R}
fdr_tag    <- ${FDR_TAG_R}
fc_tag     <- ${FC_TAG_R}
out_vol_pdf <- file.path(out_dir, sprintf("%s_volcano_%s_fdr%s_lfc%s_%s.pdf",
                                           safe_label, lfc_tag, fdr_tag, fc_tag, ts_tag))
out_vol_png <- sub("\\\\.pdf$", ".png", out_vol_pdf)
out_ma_pdf  <- sub("_volcano_", "_MA_", out_vol_pdf, fixed=TRUE)
out_ma_png  <- sub("\\\\.pdf$", ".png", out_ma_pdf)

# Re-apply thresholds.
if (fc_cutoff > 0) {
    res\$Sig <- !is.na(res\$padj) & res\$padj < fdr_cutoff &
                abs(res\$log2FoldChange) >= fc_cutoff
} else {
    res\$Sig <- !is.na(res\$padj) & res\$padj < fdr_cutoff
}
res\$Direction <- ifelse(res\$Sig & res\$log2FoldChange > 0, "Up",
                 ifelse(res\$Sig & res\$log2FoldChange < 0, "Down", "NS"))

n_up   <- sum(res\$Direction == "Up",   na.rm=TRUE)
n_down <- sum(res\$Direction == "Down", na.rm=TRUE)
cat(sprintf("Contrast: %s\nUp: %d  Down: %d\n", label, n_up, n_down))

write.csv(res, ${OUT_CSV_R}, row.names=FALSE)

# ── Volcano ──────────────────────────────────────────────────
res_plot <- res[!is.na(res\$padj), ]
res_plot\$neg_log10_padj <- -log10(res_plot\$padj + 1e-300)

# Top-N most significant peaks per direction (0 = no labels), ranked by
# adjusted p-value -- separately for up and down so one direction dominating
# on significance doesn't crowd the other off the plot.
if (top_n_label > 0) {
    up_labeled   <- head(res_plot[res_plot\$Direction == "Up",   ][order(res_plot[res_plot\$Direction == "Up",   ]\$padj), ],   top_n_label)
    down_labeled <- head(res_plot[res_plot\$Direction == "Down", ][order(res_plot[res_plot\$Direction == "Down", ]\$padj), ], top_n_label)
    label_peaks <- rbind(up_labeled, down_labeled)
    label_peaks\$peak_id <- paste0(label_peaks\$Chr, ":", label_peaks\$Start, "-", label_peaks\$End)
    cat(sprintf("Labeling top %d up peak(s) and top %d down peak(s) (%d total, by significance).\n",
                nrow(up_labeled), nrow(down_labeled), nrow(label_peaks)))
} else {
    label_peaks <- res_plot[0, ]
    label_peaks\$peak_id <- character(0)
}

subtitle_vol <- if (fc_cutoff > 0) {
    sprintf("%d up, %d down  (FDR < %.4g  |  |log2FC| >= %.3g)  --  x-axis: %s",
            n_up, n_down, fdr_cutoff, fc_cutoff, lfc_source_label)
} else {
    sprintf("%d up, %d down  (FDR < %.4g)  --  x-axis: %s", n_up, n_down, fdr_cutoff, lfc_source_label)
}

# Symmetric x-axis: cosmetic framing only, centered on 0 so the plot reads
# left/right at a glance -- does not filter, alter, or hide any point, unless
# the user explicitly asked for a smaller cap via x_max_override (in which
# case points outside it are clipped, same as any manual axis limit).
# Automatic case: max() over the actual plotted data means every point still
# fits inside the frame; this only changes how far the frame extends past
# zero on whichever side happens to have less spread.
volcano_x_lim <- if (!is.na(x_max_override)) {
    x_max_override
} else {
    max(abs(res_plot\$log2FoldChange), na.rm=TRUE)
}

p_vol <- ggplot(res_plot, aes(x=log2FoldChange, y=neg_log10_padj, color=Direction)) +
    geom_point(alpha=0.5, size=1.2) +
    { if (top_n_label > 0)
        geom_text_repel(data=label_peaks, aes(label=peak_id),
                        size=2.2, max.overlaps=15, show.legend=FALSE)
      else list() } +
    geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed",
               color="grey50", linewidth=0.5) +
    { if (fc_cutoff > 0)
        geom_vline(xintercept=c(-fc_cutoff, fc_cutoff), linetype="dashed",
                   color="grey50", linewidth=0.5)
      else list() } +
    scale_color_manual(values=c(Up="firebrick3", Down="steelblue3", NS="grey70")) +
    # oob=scales::squish (instead of plain xlim()/ylim(), whose default
    # oob=censor turns out-of-range points into NA and silently drops them):
    # any point outside the frame gets clamped to sit right at the edge
    # instead of vanishing, so a manual cap narrower than the data still
    # shows where those points piled up. Inert in the automatic case (NA
    # upper bound, or a limit computed from the data's own max) since
    # nothing ever falls outside a limit derived from the data itself.
    scale_x_continuous(limits=c(-volcano_x_lim, volcano_x_lim), oob=scales::squish) +
    scale_y_continuous(limits=c(0, y_max_override), oob=scales::squish) +
    labs(title=paste0("Volcano: ", label),
         subtitle=subtitle_vol,
         x=paste0("log2 Fold Change (", lfc_source_label, ")"),
         y="-log10(adjusted p-value)") +
    theme_bw(base_size=13)

ggsave(out_vol_pdf, p_vol, width=8, height=6)
ggsave(out_vol_png, p_vol, width=8, height=6, dpi=150)
cat("Volcano saved:", out_vol_pdf, "\n")

# ── MA plot ───────────────────────────────────────────────────
conc_col <- grep("^Conc", names(res), value=TRUE)[1]
if (!is.na(conc_col)) {
    subtitle_ma <- if (fc_cutoff > 0) {
        sprintf("FDR < %.4g  |  |log2FC| >= %.3g", fdr_cutoff, fc_cutoff)
    } else {
        sprintf("FDR < %.4g", fdr_cutoff)
    }
    p_ma <- ggplot(res, aes(x=.data[[conc_col]], y=log2FoldChange, color=Direction)) +
        geom_point(alpha=0.4, size=1.2) +
        geom_hline(yintercept=0, linetype="solid", color="grey40") +
        scale_color_manual(values=c(Up="firebrick3", Down="steelblue3", NS="grey70")) +
        labs(title=paste0("MA: ", label),
             subtitle=subtitle_ma,
             x=paste0("Mean Accessibility (", conc_col, ", log2 concentration)"),
             y=paste0("log2 Fold Change (", lfc_source_label, ")")) +
        theme_bw(base_size=13)
    ggsave(out_ma_pdf, p_ma, width=8, height=6)
    ggsave(out_ma_png, p_ma, width=8, height=6, dpi=150)
    cat("MA plot saved:", out_ma_pdf, "\n")
} else {
    cat("  [NOTE] MA plot skipped: no 'Conc' column in results.\n")
}
RSCRIPT_EOF

    label "Running re-plot..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        # Discovered by globbing the output dir, not assumed from a bash-side
        # guess -- R decides the mle/shrunken tag at runtime (see above), so
        # bash has no reliable way to know the exact filename in advance.
        local found_vol_pdf found_vol_png found_ma_pdf found_ma_png
        found_vol_pdf=$(find "$ANALYSIS_OUT_DIR" -maxdepth 1 -name "*_volcano_*.pdf" | head -1)
        found_vol_png=$(find "$ANALYSIS_OUT_DIR" -maxdepth 1 -name "*_volcano_*.png" | head -1)
        found_ma_pdf=$(find "$ANALYSIS_OUT_DIR" -maxdepth 1 -name "*_MA_*.pdf" | head -1)
        found_ma_png=$(find "$ANALYSIS_OUT_DIR" -maxdepth 1 -name "*_MA_*.png" | head -1)
        [[ -n "$found_vol_pdf" ]] && ok "Volcano PDF: $found_vol_pdf"
        [[ -n "$found_vol_png" ]] && ok "Volcano PNG: $found_vol_png"
        [[ -n "$found_ma_pdf" ]] && ok "MA PDF: $found_ma_pdf"
        [[ -n "$found_ma_png" ]] && ok "MA PNG: $found_ma_png"
        ok "Results CSV: $OUT_CSV"
    else
        err "Re-plot failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 3 — Peak accessibility boxplot
# ─────────────────────────────────────────────────────────────

run_peak_boxplot_analysis() {
    header "Peak Accessibility Boxplot"

    blank
    echo -e "  Choose how to specify peaks:"
    echo -e "    ${CYAN}1${RESET}.  Genomic coordinates  (e.g. chr1:1000000-1001000)"
    echo -e "    ${CYAN}2${RESET}.  Peak index            (row number from a contrast result CSV)"
    blank

    local PEAK_MODE
    while true; do
        read -p "  Choice [1/2]: " PEAK_MODE_INPUT
        case "$PEAK_MODE_INPUT" in
            1) PEAK_MODE="coords"; break ;;
            2) PEAK_MODE="index";  break ;;
            *) err "Enter 1 or 2." ;;
        esac
    done

    blank
    local PEAK_INPUT_STR=""
    if [[ "$PEAK_MODE" == "coords" ]]; then
        echo -e "  Enter one or more genomic coordinates, space-separated."
        echo -e "  ${DIM}Format: chr:start-end   e.g.  chr1:1000000-1001000 chr5:55000-56000${RESET}"
        blank
        read -p "  Coordinates: " PEAK_INPUT_STR
    else
        blank
        show_contrast_menu
        blank
        while true; do
            read -p "  Which contrast's result CSV to use for row index lookup? " CINPUT
            resolve_contrast "$CINPUT" && break
        done
        blank
        echo -e "  Enter one or more peak row indices (1-based), space-separated."
        echo -e "  ${DIM}These are row numbers from ${RESOLVED_CONTRAST}_all_peaks.csv${RESET}"
        blank
        read -p "  Row indices: " PEAK_INPUT_STR
    fi

    if [[ -z "$PEAK_INPUT_STR" ]]; then
        err "No peaks entered."
        return
    fi

    blank
    echo -e "  Choose which samples to include."
    prompt_sample_selection "Samples to plot:" 2

    blank
    read -p "  Overlay individual sample points on boxplots? [Y/n]: " BOX_PTS
    BOX_PTS_R=$([[ "${BOX_PTS,,}" == "n" ]] && echo "FALSE" || echo "TRUE")

    local TAG="boxplot_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/03_Peak_Boxplots/Peak_Boxplot_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_PDF="$ANALYSIS_OUT_DIR/peak_boxplot_${TS}.pdf"
    local OUT_PNG="${OUT_PDF%.pdf}.png"
    local OUT_CSV="${OUT_PDF%.pdf}_values.csv"

    local R_SAMPLES
    R_SAMPLES=$(r_vector_literal "${RESOLVED_SAMPLES[@]}")

    # Peak input is raw, free-form keyboard input -- the highest-risk value
    # in this whole script if embedded raw, so it goes through the same
    # safe encoding as everything else here.
    local PEAK_INPUT_R
    PEAK_INPUT_R="$(r_string_literal "$PEAK_INPUT_STR")"
    local INDEX_CONTRAST_R
    INDEX_CONTRAST_R="$(r_string_literal "${RESOLVED_CONTRAST:-}")"
    local OUT_CSV_R OUT_PDF_R OUT_PNG_R
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"
    OUT_PDF_R="$(r_string_literal "$OUT_PDF")"
    OUT_PNG_R="$(r_string_literal "$OUT_PNG")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b            <- readRDS(${BUNDLE_PATH_R})
peak_mode    <- "$PEAK_MODE"
peak_input   <- $PEAK_INPUT_R
sel_samples  <- $R_SAMPLES
show_pts     <- $BOX_PTS_R
diff_out     <- b\$diff_out
idx_contrast <- $INDEX_CONTRAST_R

# ── Locate the raw count matrix saved by the diff run ────────
counts_rds <- file.path(diff_out, "diagnostics", "consensus_peak_raw_counts.rds")
if (!file.exists(counts_rds)) {
    stop("Raw count matrix not found at: ", counts_rds,
         "\nExpected diagnostics/consensus_peak_raw_counts.rds from the diff run.")
}
count_mat <- readRDS(counts_rds)
cat("Loaded count matrix:", nrow(count_mat), "peaks x", ncol(count_mat), "columns\n")

# ── Identify coordinate columns ──────────────────────────────
coord_cols <- intersect(c("Chr","chr","CHR"), names(count_mat))
start_cols <- intersect(c("Start","start","START"), names(count_mat))
end_cols   <- intersect(c("End","end","END"), names(count_mat))

if (length(coord_cols) == 0 || length(start_cols) == 0 || length(end_cols) == 0) {
    stop("Cannot find Chr/Start/End columns in count matrix. Available: ",
         paste(names(count_mat), collapse=", "))
}
chr_col   <- coord_cols[1]
start_col <- start_cols[1]
end_col   <- end_cols[1]

peak_ids_all <- paste0(count_mat[[chr_col]], ":",
                       count_mat[[start_col]], "-",
                       count_mat[[end_col]])
rownames(count_mat) <- make.unique(peak_ids_all)

# ── Identify sample count columns ────────────────────────────
meta_cols  <- c(chr_col, start_col, end_col,
                "Chr","Start","End","chr","start","end",
                "CHR","START","END","Conc","width","strand","Score","score","Name","name")
count_cols <- setdiff(names(count_mat), meta_cols)

# Further restrict to numeric columns only.
count_cols <- count_cols[vapply(count_mat[, count_cols, drop=FALSE], is.numeric, logical(1))]

cat("Sample count columns found:", length(count_cols), "\n")
if (length(count_cols) == 0) {
    stop("No numeric sample columns identified in the count matrix.")
}

# ── Select requested samples ─────────────────────────────────
# Try exact match, then make.names(), then punctuation-insensitive.
resolve_col <- function(sample_id, cols) {
    if (sample_id %in% cols) return(sample_id)
    mn <- make.names(sample_id)
    if (mn %in% make.names(cols)) return(cols[make.names(cols) == mn][1])
    norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
    hits <- cols[norm(cols) == norm(sample_id)]
    if (length(hits) == 1) return(hits)
    NA_character_
}

resolved_cols <- vapply(sel_samples, resolve_col, character(1), cols=count_cols)
missing_s     <- sel_samples[is.na(resolved_cols)]
if (length(missing_s) > 0) {
    warning("Could not find count columns for samples: ", paste(missing_s, collapse=", "))
}
resolved_cols <- na.omit(resolved_cols)
if (length(resolved_cols) < 2) {
    stop("Fewer than 2 sample columns resolved — cannot draw boxplots.")
}
sample_names_matched <- sel_samples[!is.na(vapply(sel_samples, resolve_col, character(1), cols=count_cols))]

raw_sub <- as.matrix(count_mat[, resolved_cols, drop=FALSE])
colnames(raw_sub) <- sample_names_matched

# ── Resolve requested peaks ───────────────────────────────────
tokens <- trimws(strsplit(peak_input, "\\\\s+")[[1]])
tokens <- tokens[nzchar(tokens)]
peak_rows <- integer(0)

if (peak_mode == "coords") {
    for (tok in tokens) {
        # Parse chr:start-end, tolerating spaces around : and -
        m <- regmatches(tok, regexec("^([^:]+):([0-9]+)-([0-9]+)$", tok))[[1]]
        if (length(m) != 4) {
            warning("Could not parse coordinate: ", tok, " — skipping.")
            next
        }
        q_chr   <- m[2]
        q_start <- as.integer(m[3])
        q_end   <- as.integer(m[4])
        # Find peaks that overlap the query interval.
        hits <- which(
            count_mat[[chr_col]]   == q_chr &
            as.integer(count_mat[[end_col]])   >= q_start &
            as.integer(count_mat[[start_col]]) <= q_end
        )
        if (length(hits) == 0) {
            warning("No peaks overlap ", tok, " — skipping.")
        } else {
            cat(sprintf("  %s -> %d matching peak(s)\n", tok, length(hits)))
            peak_rows <- c(peak_rows, hits)
        }
    }
} else {
    # Index mode: use the row numbers from the specified contrast CSV.
    c_dir    <- file.path(diff_out, idx_contrast)
    csv_path <- file.path(c_dir, paste0(idx_contrast, "_all_peaks.csv"))
    if (!file.exists(csv_path)) {
        stop("Result CSV not found: ", csv_path)
    }
    ref_df <- read.csv(csv_path)
    for (tok in tokens) {
        idx <- suppressWarnings(as.integer(tok))
        if (is.na(idx) || idx < 1 || idx > nrow(ref_df)) {
            warning("Row index out of range or invalid: ", tok, " — skipping.")
            next
        }
        # Match this row back into the count matrix by coordinates.
        q_chr   <- as.character(ref_df\$Chr[idx])
        q_start <- as.integer(ref_df\$Start[idx])
        q_end   <- as.integer(ref_df\$End[idx])
        hits <- which(
            count_mat[[chr_col]]   == q_chr &
            as.integer(count_mat[[start_col]]) == q_start &
            as.integer(count_mat[[end_col]])   == q_end
        )
        if (length(hits) == 0) {
            warning("Row ", idx, " (", q_chr, ":", q_start, "-", q_end,
                    ") not found in count matrix — skipping.")
        } else {
            cat(sprintf("  Row %d -> %s:%d-%d (%d match(es))\n",
                        idx, q_chr, q_start, q_end, length(hits)))
            peak_rows <- c(peak_rows, hits[1])
        }
    }
}

peak_rows <- unique(peak_rows)
if (length(peak_rows) == 0) {
    stop("No peaks were resolved from the input.")
}
cat(sprintf("Plotting %d peak(s).\n", length(peak_rows)))

raw_counts <- raw_sub[peak_rows, , drop=FALSE]

# ── Normalized counts ──────────────────────────────────────────
# Prefer the REAL per-sample normalization factors DiffBind/DESeq2 used.
# diff_analysis.sh normalizes each contrast independently, on a subset
# containing only its own two groups (DiffBind's modern design-based
# contrast mode) -- so there is no longer one single "real" factor per
# sample across the whole experiment; a sample used in two different
# contrasts can legitimately have two different real factors. Preference
# order:
#   1. This contrast's own <label>/diagnostics/diffbind_normalization_factors.tsv
#      -- the actual factors that specific comparison's DESeq2 model used.
#      Only tried when a contrast is known (index mode).
#   2. diagnostics/diffbind_normalization_factors_whole_experiment.tsv --
#      real RLE/background= factors (same method as any contrast), just
#      computed once across all samples rather than scoped to one
#      comparison. Not what any specific contrast used, but still a real
#      DiffBind normalization, not an approximation.
#   3. diagnostics/diffbind_normalization_factors.tsv (no "_whole_experiment"
#      suffix) -- the filename used before this per-contrast/whole-experiment
#      split existed. An older schema-2 bundle can have real factors sitting
#      right here under the old name; without this fallback such a bundle
#      would skip straight past two genuinely real sources to the ad hoc
#      estimate for no reason other than a filename rename.
#   4. An ad hoc poscounts refit on just the selected samples, clearly
#      labeled as such -- only reached if none of the above resolve.
norm_label  <- "Normalized counts"
norm_counts <- NULL
norm_source <- NULL
norm_source_is_contrast_specific <- FALSE

if (nzchar(idx_contrast)) {
    contrast_nf_path <- file.path(diff_out, idx_contrast, "diagnostics",
                                  "diffbind_normalization_factors.tsv")
    if (file.exists(contrast_nf_path)) {
        norm_source <- contrast_nf_path
        norm_source_is_contrast_specific <- TRUE
    }
}
if (is.null(norm_source)) {
    whole_exp_nf_path <- file.path(diff_out, "diagnostics",
                                   "diffbind_normalization_factors_whole_experiment.tsv")
    if (file.exists(whole_exp_nf_path)) {
        norm_source <- whole_exp_nf_path
    }
}
if (is.null(norm_source)) {
    legacy_nf_path <- file.path(diff_out, "diagnostics",
                                "diffbind_normalization_factors.tsv")
    if (file.exists(legacy_nf_path)) {
        norm_source <- legacy_nf_path
    }
}

if (!is.null(norm_source)) {
    nf_df <- tryCatch(
        read.delim(norm_source, stringsAsFactors=FALSE, check.names=FALSE),
        error=function(e) NULL
    )
    if (!is.null(nf_df) && all(c("Sample","NormFactor") %in% colnames(nf_df))) {
        # Resolve each selected sample to a row in the norm-factors table,
        # tolerating the same punctuation/casing drift as resolve_col()
        # above (exact match, then make.names(), then alnum-only).
        resolve_nf <- function(sample_id) {
            if (sample_id %in% nf_df\$Sample) {
                return(nf_df\$NormFactor[nf_df\$Sample == sample_id][1])
            }
            mn  <- make.names(sample_id)
            hit <- nf_df\$Sample[make.names(nf_df\$Sample) == mn]
            if (length(hit) == 1) return(nf_df\$NormFactor[nf_df\$Sample == hit][1])
            norm_str <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
            hit2 <- nf_df\$Sample[norm_str(nf_df\$Sample) == norm_str(sample_id)]
            if (length(hit2) == 1) return(nf_df\$NormFactor[nf_df\$Sample == hit2][1])
            NA_real_
        }
        factors <- vapply(sample_names_matched, resolve_nf, numeric(1))
        if (!anyNA(factors) && all(factors > 0)) {
            norm_counts <- sweep(raw_counts, 2, factors[colnames(raw_counts)], "/")
            if (norm_source_is_contrast_specific) {
                cat(sprintf("Normalized counts: using contrast '%s's actual DiffBind normalization factors from\n",
                            idx_contrast))
                cat(sprintf("  %s\n", norm_source))
            } else {
                norm_label <- "Normalized counts (whole-experiment factors, not contrast-specific)"
                cat("Normalized counts: using real whole-experiment DiffBind normalization factors\n")
                cat("(not scoped to one contrast -- pick a contrast in index mode for that contrast's\n")
                cat("own factors instead):\n")
                cat(sprintf("  %s\n", norm_source))
            }
        } else {
            cat("Could not match all selected samples to", basename(norm_source),
                "-- falling back to an ad hoc estimate for display only.\n")
        }
    }
}

col_data_r <- if (!is.null(b\$col_data)) {
    df <- b\$col_data[b\$col_data\$SampleID %in% sample_names_matched, , drop=FALSE]
    rownames(df) <- df\$SampleID
    df[sample_names_matched, , drop=FALSE]
} else {
    data.frame(SampleID=sample_names_matched,
               Condition="unknown",
               row.names=sample_names_matched)
}

if (is.null(norm_counts)) {
    norm_label <- "Normalized counts (ad hoc)"
    norm_counts <- tryCatch({
        full_raw <- as.matrix(count_mat[, resolved_cols, drop=FALSE])
        colnames(full_raw) <- sample_names_matched
        full_raw[!is.finite(full_raw) | full_raw < 0] <- 0
        full_raw <- round(full_raw)
        storage.mode(full_raw) <- "integer"
        dds_tmp <- DESeqDataSetFromMatrix(
            countData = full_raw,
            colData   = col_data_r,
            design    = ~1
        )
        dds_tmp <- estimateSizeFactors(dds_tmp, type="poscounts")
        cat("Normalized counts: real normalization factors were not available (or didn't\n")
        cat("cover the selected samples) -- estimating ad hoc (poscounts, this sample\n")
        cat("subset only) for display purposes. These numbers do NOT match the pipeline's\n")
        cat("actual statistics; use them to eyeball relative accessibility only, not as\n")
        cat("the analysis result.\n")
        counts(dds_tmp, normalized=TRUE)[peak_rows, , drop=FALSE]
    }, error = function(e) {
        warning("Normalization failed (", conditionMessage(e), ") — showing raw counts only.")
        NULL
    })
}

# ── Conditions for x-axis ────────────────────────────────────
conditions_map <- if (!is.null(b\$col_data)) {
    setNames(as.character(b\$col_data\$Condition), b\$col_data\$SampleID)
} else {
    setNames(rep("unknown", length(sample_names_matched)), sample_names_matched)
}

peak_labels <- rownames(raw_counts)

make_long <- function(mat, count_type) {
    df <- as.data.frame(t(mat))
    df\$sample <- rownames(df)
    tidyr::pivot_longer(df, -sample, names_to="peak", values_to="count") %>%
        dplyr::mutate(
            condition  = conditions_map[sample],
            count_type = count_type
        )
}

long_raw  <- make_long(raw_counts, "Raw counts")
long_norm <- if (!is.null(norm_counts)) make_long(norm_counts, norm_label) else NULL

plot_df <- dplyr::bind_rows(long_raw, long_norm)
plot_df\$count_type <- factor(plot_df\$count_type,
                              levels=c("Raw counts", norm_label))

write.csv(plot_df, ${OUT_CSV_R}, row.names=FALSE)

n_peaks  <- length(unique(plot_df\$peak))
n_panels <- if (!is.null(norm_counts)) 2 else 1
n_rows   <- ceiling(n_peaks / min(4, n_peaks))

subtitle_txt <- paste0(n_peaks, " peak(s) - raw and normalized counts")
if (grepl("ad hoc", norm_label, fixed=TRUE)) {
    subtitle_txt <- paste0(subtitle_txt,
        " (normalized panel is an ad hoc estimate, not the pipeline's real normalization)")
} else if (grepl("whole-experiment", norm_label, fixed=TRUE)) {
    subtitle_txt <- paste0(subtitle_txt,
        " (normalized panel uses whole-experiment factors, not this specific contrast's own)")
}

p <- ggplot(plot_df, aes(condition, count, fill=condition)) +
    geom_boxplot(outlier.shape=NA, alpha=0.7) +
    facet_grid(count_type ~ peak, scales="free_y") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="none",
          strip.text.x=element_text(size=7)) +
    labs(x=NULL, y="Count", title="Peak Accessibility by Group",
         subtitle=subtitle_txt)

if (show_pts) {
    p <- p + geom_jitter(width=0.15, size=1.4, alpha=0.7, color="black")
}

plot_w <- max(8, n_peaks * 2.5 + 2)
plot_h <- max(5, n_panels * 3.5)
ggsave(${OUT_PDF_R}, plot=p, width=plot_w, height=plot_h)
ggsave(${OUT_PNG_R}, plot=p, width=plot_w, height=plot_h, dpi=150)
cat("Boxplot saved:", ${OUT_PDF_R}, "\n")
cat("Values CSV saved:", ${OUT_CSV_R}, "\n")
RSCRIPT_EOF

    label "Generating peak boxplot..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "Boxplot PDF: $OUT_PDF"
        ok "Boxplot PNG: $OUT_PNG"
        ok "Values CSV: $OUT_CSV"
    else
        err "Peak boxplot failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 4 — Sample correlation heatmap
# ─────────────────────────────────────────────────────────────

run_correlation_heatmap() {
    header "Sample Correlation Heatmap"

    if [[ "$HAS_VST" != "TRUE" ]]; then
        err "Bundle has no VST matrix — correlation heatmap is unavailable."
        return
    fi

    echo -e "  Choose which samples to include."
    prompt_sample_selection "Samples for correlation heatmap:" 2

    blank
    echo -e "  Correlation method:"
    echo -e "    ${CYAN}1${RESET}.  Pearson"
    echo -e "    ${CYAN}2${RESET}.  Spearman"
    blank
    local CORR_METHOD="pearson"
    while true; do
        read -p "  Choice [1/2]: " CORR_IN
        case "$CORR_IN" in
            1) CORR_METHOD="pearson";  ok "Pearson selected.";  break ;;
            2) CORR_METHOD="spearman"; ok "Spearman selected."; break ;;
            *) err "Enter 1 or 2." ;;
        esac
    done

    local TAG="corrheat_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/04_Sample_Correlation/${CORR_METHOD}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_PDF="$ANALYSIS_OUT_DIR/sample_correlation_${CORR_METHOD}_${TS}.pdf"
    local OUT_PNG="${OUT_PDF%.pdf}.png"
    local OUT_CSV="${OUT_PDF%.pdf}_matrix.csv"

    local R_SAMPLES
    R_SAMPLES=$(r_vector_literal "${RESOLVED_SAMPLES[@]}")
    local OUT_CSV_R OUT_PDF_R OUT_PNG_R
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"
    OUT_PDF_R="$(r_string_literal "$OUT_PDF")"
    OUT_PNG_R="$(r_string_literal "$OUT_PNG")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(pheatmap)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b           <- readRDS(${BUNDLE_PATH_R})
vst_mat     <- b\$vst_matrix
sel_samples <- $R_SAMPLES
method      <- "$CORR_METHOD"

avail <- intersect(sel_samples, colnames(vst_mat))
if (length(avail) < 2) stop("Need at least 2 samples in the VST matrix.")

mat        <- vst_mat[, avail, drop=FALSE]
cor_matrix <- cor(mat, method=method)
write.csv(cor_matrix, ${OUT_CSV_R})

conditions <- if (!is.null(b\$col_data)) {
    setNames(as.character(b\$col_data\$Condition), b\$col_data\$SampleID)
} else {
    setNames(rep("unknown", length(avail)), avail)
}

col_ann <- data.frame(Group=conditions[avail], row.names=avail)
grp_colors <- setNames(
    colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628"))(
        length(unique(col_ann\$Group))),
    unique(as.character(col_ann\$Group))
)
ann_colors <- list(Group=grp_colors)

cor_breaks <- c(seq(-1, 0, length.out=100), seq(0, 1, length.out=100)[-1])
cor_colors <- c(
    colorRampPalette(c("darkblue","blue","cornflowerblue","white"))(99),
    colorRampPalette(c("white","lightyellow","orange","red","darkred"))(99)
)

w <- max(6, ncol(cor_matrix) * 0.9 + 2)
h <- max(5, ncol(cor_matrix) * 0.8 + 2)

pdf(${OUT_PDF_R}, width=w, height=h)
pheatmap(cor_matrix,
    cluster_rows=TRUE, cluster_cols=TRUE,
    display_numbers=TRUE, number_format="%.2f", fontsize_number=8,
    color=cor_colors, breaks=cor_breaks,
    annotation_col=col_ann, annotation_colors=ann_colors,
    main=paste0("Sample Correlation (", tools::toTitleCase(method), ")"))
dev.off()

# PNG via a temporary PDF-to-image approach is not available in this env;
# re-draw with ggsave wrapper instead.
png(${OUT_PNG_R}, width=w*150, height=h*150, res=150)
pheatmap(cor_matrix,
    cluster_rows=TRUE, cluster_cols=TRUE,
    display_numbers=TRUE, number_format="%.2f", fontsize_number=8,
    color=cor_colors, breaks=cor_breaks,
    annotation_col=col_ann, annotation_colors=ann_colors,
    main=paste0("Sample Correlation (", tools::toTitleCase(method), ")"))
dev.off()

cat("Correlation heatmap saved:", ${OUT_PDF_R}, "\n")
cat("Correlation matrix saved:", ${OUT_CSV_R}, "\n")
RSCRIPT_EOF

    label "Generating correlation heatmap..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "Heatmap PDF: $OUT_PDF"
        ok "Heatmap PNG: $OUT_PNG"
        ok "Matrix CSV:  $OUT_CSV"
    else
        err "Correlation heatmap failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 5 — Motif summary
# ─────────────────────────────────────────────────────────────

run_motif_summary() {
    header "Motif Summary"

    if [[ "$HAS_HOMER" != "TRUE" ]]; then
        err "No saved HOMER motif result CSV found in any contrast folder."
        err "Expected: <diff_out>/<contrast>/motifs/<contrast>_known_motifs_top50.csv"
        err "Either motif analysis wasn't run/enabled, or every contrast had too few"
        err "significant peaks for it to produce output."
        return
    fi

    blank
    show_contrast_menu
    blank

    while true; do
        read -p "  Choose contrast (name or number): " CINPUT
        resolve_contrast "$CINPUT" && break
    done
    ok "Contrast: $RESOLVED_CONTRAST"

    blank
    read -p "  How many top motifs to show? [default: 10]: " N_MOTIFS
    N_MOTIFS="${N_MOTIFS:-10}"
    [[ "$N_MOTIFS" =~ ^[0-9]+$ ]] || { warn "Invalid — using 10."; N_MOTIFS=10; }

    local TAG="motifs_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/05_Motif_Summaries/${SAFE_LABEL}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_CSV="$ANALYSIS_OUT_DIR/${SAFE_LABEL}_motif_summary_top${N_MOTIFS}_${TS}.csv"

    local RESOLVED_CONTRAST_R OUT_CSV_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b        <- readRDS(${BUNDLE_PATH_R})
label    <- ${RESOLVED_CONTRAST_R}
n_top    <- $N_MOTIFS
diff_out <- b\$diff_out

motif_csv <- file.path(diff_out, label, "motifs",
                       paste0(label, "_known_motifs_top50.csv"))

if (!file.exists(motif_csv)) {
    cat(sprintf(
        "  No significant peaks for contrast '%s', or motif analysis was not run.\n",
        label))
    cat(sprintf("  Expected: %s\n", motif_csv))
    quit(status=0)
}

df <- read.csv(motif_csv, check.names=FALSE)
if (nrow(df) == 0) {
    cat("  Motif CSV found but is empty.\n")
    quit(status=0)
}

n_show <- min(n_top, nrow(df))

# Identify key columns by partial name match (HOMER column names vary).
name_col <- grep("Motif.Name|motif.name|Name", colnames(df), ignore.case=TRUE, value=TRUE)[1]
pval_col <- grep("P-value|p.value|pvalue", colnames(df), ignore.case=TRUE, value=TRUE)[1]
pct_col  <- grep("% of Target|Target%|target.percent", colnames(df), ignore.case=TRUE, value=TRUE)[1]
qval_col <- grep("q.value|Q-value|qvalue|FDR", colnames(df), ignore.case=TRUE, value=TRUE)[1]

cat(sprintf("\n  Top %d known motifs for: %s\n", n_show, label))
cat(sprintf("  (from %s)\n\n", basename(motif_csv)))
cat(sprintf("  %-4s  %-42s  %-12s  %-12s  %s\n",
            "#", "Motif Name", "P-value", "q-value", "% Target"))
cat(paste0(rep("-", 85), collapse=""), "\n")

for (i in seq_len(n_show)) {
    mname <- if (!is.na(name_col)) as.character(df[i, name_col]) else "?"
    mpval <- if (!is.na(pval_col)) as.character(df[i, pval_col]) else "?"
    mqval <- if (!is.na(qval_col)) as.character(df[i, qval_col]) else "?"
    mpct  <- if (!is.na(pct_col))  as.character(df[i, pct_col])  else "?"
    cat(sprintf("  %-4d  %-42s  %-12s  %-12s  %s\n",
                i, substr(mname, 1, 42), mpval, mqval, mpct))
}

out_df <- df[seq_len(n_show), , drop=FALSE]
out_csv_path <- ${OUT_CSV_R}
write.csv(out_df, out_csv_path, row.names=FALSE)
cat(sprintf("\n  Saved top %d motifs to: %s\n", n_show, basename(out_csv_path)))
RSCRIPT_EOF

    label "Loading motif summary..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "Motif summary CSV: $OUT_CSV"
    else
        err "Motif summary failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 6 — Annotation distribution plots
# ─────────────────────────────────────────────────────────────
# Reads <label>_annotated_peaks.tsv (written by diff_analysis.sh's
# ChIPseeker step) and redraws the genomic-feature / TSS / TES /
# gene-body / width / chromosome distributions on demand. These plots
# used to be generated once during the diff run; diff_analysis.sh now
# only computes and saves the annotation table, so this is the only
# place they get drawn.

run_annotation_plots() {
    header "Annotation Distribution Plots"

    blank
    show_contrast_menu
    blank

    while true; do
        read -p "  Choose contrast (name or number): " CINPUT
        resolve_contrast "$CINPUT" && break
    done
    ok "Contrast: $RESOLVED_CONTRAST"

    blank
    echo -e "  Available plots (from ChIPseeker annotation of significant peaks):"
    echo -e "    ${CYAN}1${RESET}.  Genomic feature composition (100% stacked bar)"
    echo -e "    ${CYAN}2${RESET}.  Distance to nearest TSS (histogram)"
    echo -e "    ${CYAN}3${RESET}.  Distance to nearest-gene TES (histogram)"
    echo -e "    ${CYAN}4${RESET}.  Position within gene body (histogram)"
    echo -e "    ${CYAN}5${RESET}.  Peak width distribution (histogram)"
    echo -e "    ${CYAN}6${RESET}.  Chromosome distribution (bar)"
    echo -e "    ${CYAN}7${RESET}.  All of the above"
    blank

    local -a PLOT_KEYS
    while true; do
        read -p "  Choice (space-separated, e.g. '1 3 6') [default: 7]: " ANNPLOT_IN
        ANNPLOT_IN="${ANNPLOT_IN:-7}"
        PLOT_KEYS=()
        local bad=false all_selected=false
        IFS=' ' read -ra _sel <<< "$ANNPLOT_IN"
        for tok in "${_sel[@]}"; do
            case "$tok" in
                1) PLOT_KEYS+=("feature") ;;
                2) PLOT_KEYS+=("tss") ;;
                3) PLOT_KEYS+=("tes") ;;
                4) PLOT_KEYS+=("body") ;;
                5) PLOT_KEYS+=("width") ;;
                6) PLOT_KEYS+=("chr") ;;
                7) PLOT_KEYS=(feature tss tes body width chr); all_selected=true ;;
                *) err "Unrecognized option: $tok"; bad=true ;;
            esac
            $all_selected && break
        done
        $bad && continue
        if [[ ${#PLOT_KEYS[@]} -eq 0 ]]; then
            err "Choose at least one option."
            continue
        fi
        break
    done

    local TAG="annplots_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/06_Annotation_Plots/${SAFE_LABEL}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"

    local PLOTS_WANTED_R
    PLOTS_WANTED_R="$(r_vector_literal "${PLOT_KEYS[@]}")"
    local RESOLVED_CONTRAST_R ANALYSIS_OUT_DIR_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    ANALYSIS_OUT_DIR_R="$(r_string_literal "$ANALYSIS_OUT_DIR")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(ggplot2)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b            <- readRDS(${BUNDLE_PATH_R})
label        <- ${RESOLVED_CONTRAST_R}
diff_out     <- b\$diff_out
txdb_pkg     <- if (!is.null(b\$txdb)) b\$txdb else "TxDb"
out_dir      <- ${ANALYSIS_OUT_DIR_R}
plots_wanted <- $PLOTS_WANTED_R

ann_path <- file.path(diff_out, label, paste0(label, "_annotated_peaks.tsv"))

if (!file.exists(ann_path)) {
    cat(sprintf("  No annotated-peaks table found for contrast '%s'.\n", label))
    cat(sprintf("  Expected: %s\n", ann_path))
    quit(status=0)
}

ann <- tryCatch(
    read.delim(ann_path, header=TRUE, stringsAsFactors=FALSE,
               quote="", comment.char="", check.names=FALSE),
    error=function(e) {
        cat(sprintf("  Could not read annotation table: %s\n", conditionMessage(e)))
        NULL
    }
)
if (is.null(ann) || nrow(ann) == 0) {
    cat("  Annotation table found but has no rows (no significant peaks were annotated for this contrast).\n")
    quit(status=0)
}
cat(sprintf("  Loaded %d annotated peak(s) for '%s'.\n", nrow(ann), label))

save_both <- function(p, path_base, width, height) {
    ggsave(paste0(path_base, ".pdf"), p, width=width, height=height)
    ggsave(paste0(path_base, ".png"), p, width=width, height=height, dpi=150)
}

feature_order <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron",
                   "Downstream", "Distal Intergenic", "Other", "Unannotated")

if ("feature" %in% plots_wanted) {
    feat <- factor(ann\$Feature, levels=feature_order)
    feature_df <- as.data.frame(table(feat), stringsAsFactors=FALSE)
    colnames(feature_df) <- c("Feature", "Count")
    feature_df <- feature_df[feature_df\$Count > 0, , drop=FALSE]
    if (nrow(feature_df) > 0) {
        feature_df\$Percent <- 100 * feature_df\$Count / sum(feature_df\$Count)
        write.csv(feature_df, file.path(out_dir, "genomic_features_counts.csv"), row.names=FALSE)
        # Single horizontal 100% stacked bar: one bar spanning 0-100%, split
        # into segments by Feature share, white borders act as the dividing
        # lines between segments. Keep feature_order's ordering in the stack
        # (not table()'s alphabetical default) and only label segments large
        # enough to hold readable text -- small slivers still show via the
        # legend and the CSV.
        feature_df\$Feature <- factor(feature_df\$Feature, levels=feature_order)
        feature_df\$Label <- ifelse(feature_df\$Percent >= 3, sprintf("%.1f%%", feature_df\$Percent), "")
        # y holds the dummy single-category value and x holds Percent
        # directly (rather than mapping x=category/y=Percent and calling
        # coord_flip()) -- geom_col auto-detects the horizontal orientation
        # from which aesthetic is discrete, so axis.text.y/axis.ticks.y/
        # panel.grid.major.y below unambiguously target the dummy category
        # axis with no coord_flip theme-element ambiguity to get backwards.
        # Legend doubles as the "table on the side": each key's label is
        # "Feature (xx.x%)", one per row, in a single right-hand column --
        # this is what makes slivers too thin to hold their own on-bar
        # text (e.g. 3' UTR at ~1%) still readable, without a second
        # plot/table object (which would need gridExtra/patchwork/cowplot --
        # none of which this script currently depends on). The labeller is
        # a function keyed by name (match(breaks, feature_df\$Feature)), not
        # a plain positional vector, so it stays correctly aligned to
        # whichever levels ggplot actually keeps after dropping unused
        # feature_order categories with zero peaks.
        pct_by_feature <- setNames(feature_df\$Percent, as.character(feature_df\$Feature))
        p_feat <- ggplot(feature_df, aes(y="", x=Percent, fill=Feature)) +
            geom_col(width=0.6, color="white", linewidth=0.6) +
            geom_text(aes(label=Label), position=position_stack(vjust=0.5),
                      color="white", fontface="bold", size=3.5) +
            scale_x_continuous(limits=c(0, 100), expand=c(0, 0),
                               breaks=seq(0, 100, 25), labels=function(x) paste0(x, "%")) +
            scale_fill_discrete(labels=function(breaks) {
                sprintf("%s (%.1f%%)", breaks, pct_by_feature[breaks])
            }) +
            labs(title=paste0(label, " - Genomic Feature Annotation"),
                 subtitle=paste0("ChIPseeker / ", txdb_pkg),
                 y=NULL, x="Share of significant peaks", fill="Feature") +
            theme_minimal(base_size=11) +
            theme(plot.title=element_text(face="bold"),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  panel.grid.major.y=element_blank(),
                  panel.grid.minor=element_blank(),
                  legend.position="right") +
            guides(fill=guide_legend(ncol=1))
        save_both(p_feat, file.path(out_dir, "genomic_features_bar"), 10, 4.2)
        cat("  Saved: genomic_features_bar (pdf/png) + genomic_features_counts.csv\n")
    } else {
        cat("  Feature composition: no usable Feature values found.\n")
    }
}

if ("tss" %in% plots_wanted) {
    tss_kb <- suppressWarnings(as.numeric(ann\$DistanceToTSS) / 1000)
    tss_kb <- tss_kb[is.finite(tss_kb) & abs(tss_kb) <= 100]
    if (length(tss_kb) > 0) {
        p_tss <- ggplot(data.frame(DistanceKb=tss_kb), aes(x=DistanceKb)) +
            geom_histogram(bins=80, boundary=0, color="white", linewidth=0.15) +
            geom_vline(xintercept=0, linetype="dashed") +
            labs(title=paste0(label, " - Distance to Nearest TSS"),
                 subtitle=sprintf("%d annotated peaks shown within ±100 kb", length(tss_kb)),
                 x="Signed distance to TSS (kb)", y="Peak count") +
            theme_bw(base_size=12) +
            theme(plot.title=element_text(face="bold"))
        save_both(p_tss, file.path(out_dir, "tss_distance_histogram"), 9, 5.5)
        cat("  Saved: tss_distance_histogram (pdf/png)\n")
    } else {
        cat("  TSS distance: no usable DistanceToTSS values within ±100 kb.\n")
    }
}

if ("tes" %in% plots_wanted) {
    tes_kb <- suppressWarnings(as.numeric(ann\$DistanceToTES) / 1000)
    tes_kb <- tes_kb[is.finite(tes_kb) & abs(tes_kb) <= 100]
    if (length(tes_kb) > 0) {
        p_tes <- ggplot(data.frame(DistanceKb=tes_kb), aes(x=DistanceKb)) +
            geom_histogram(bins=80, boundary=0, color="white", linewidth=0.15) +
            geom_vline(xintercept=0, linetype="dashed") +
            labs(title=paste0(label, " - Distance to Nearest-Gene TES"),
                 subtitle=sprintf("%d annotated peaks shown within ±100 kb", length(tes_kb)),
                 x="Signed distance to TES (kb)", y="Peak count") +
            theme_bw(base_size=12) +
            theme(plot.title=element_text(face="bold"))
        save_both(p_tes, file.path(out_dir, "tes_distance_histogram"), 9, 5.5)
        cat("  Saved: tes_distance_histogram (pdf/png)\n")
    } else {
        cat("  TES distance: no usable DistanceToTES values within ±100 kb.\n")
    }
}

if ("body" %in% plots_wanted) {
    body_pos <- suppressWarnings(as.numeric(ann\$GeneBodyPosition))
    body_pos <- body_pos[is.finite(body_pos) & body_pos >= 0 & body_pos <= 1]
    if (length(body_pos) > 0) {
        p_body <- ggplot(data.frame(Position=body_pos), aes(x=Position)) +
            geom_histogram(bins=50, boundary=0, color="white", linewidth=0.15) +
            scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25),
                               labels=c("TSS", "25%", "50%", "75%", "TES")) +
            labs(title=paste0(label, " - Position Within Nearest Gene Body"),
                 subtitle=sprintf("%d peaks physically inside the assigned TxDb gene", length(body_pos)),
                 x="Relative transcriptional position", y="Peak count") +
            theme_bw(base_size=12) +
            theme(plot.title=element_text(face="bold"))
        save_both(p_body, file.path(out_dir, "gene_body_position"), 9, 5.5)
        cat("  Saved: gene_body_position (pdf/png)\n")
    } else {
        cat("  Gene body position: no peaks fall inside their assigned gene body.\n")
    }
}

if ("width" %in% plots_wanted) {
    widths <- suppressWarnings(as.numeric(ann\$Width))
    widths <- widths[is.finite(widths) & widths > 0]
    if (length(widths) > 0) {
        upper <- as.numeric(stats::quantile(widths, 0.99, na.rm=TRUE))
        width_plot <- widths[widths <= upper]
        p_width <- ggplot(data.frame(Width=width_plot), aes(x=Width)) +
            geom_histogram(bins=60, color="white", linewidth=0.15) +
            labs(title=paste0(label, " - Significant Peak Widths"),
                 subtitle=sprintf("99th percentile display limit: %.0f bp", upper),
                 x="Peak width (bp)", y="Peak count") +
            theme_bw(base_size=12) +
            theme(plot.title=element_text(face="bold"))
        save_both(p_width, file.path(out_dir, "peak_width_distribution"), 9, 5.5)
        cat("  Saved: peak_width_distribution (pdf/png)\n")
    } else {
        cat("  Peak width: no usable Width values found.\n")
    }
}

if ("chr" %in% plots_wanted) {
    chr_df <- as.data.frame(table(ann\$Chr), stringsAsFactors=FALSE)
    colnames(chr_df) <- c("Chr", "Count")
    chr_df <- chr_df[chr_df\$Count > 0, , drop=FALSE]
    if (nrow(chr_df) > 0) {
        chr_df\$Chr <- factor(chr_df\$Chr, levels=unique(chr_df\$Chr))
        p_chr <- ggplot(chr_df, aes(x=Chr, y=Count)) +
            geom_col(color="white", linewidth=0.2) +
            labs(title=paste0(label, " - Chromosome Distribution"),
                 subtitle=sprintf("%d significant peaks across %d sequences",
                                  nrow(ann), nrow(chr_df)),
                 x="Chromosome", y="Peak count") +
            theme_bw(base_size=12) +
            theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
                  plot.title=element_text(face="bold"))
        save_both(p_chr, file.path(out_dir, "chromosome_distribution"), 10, 5)
        cat("  Saved: chromosome_distribution (pdf/png)\n")
    } else {
        cat("  Chromosome distribution: no usable Chr values found.\n")
    }
}

cat(sprintf("\n  Outputs written to: %s\n", out_dir))
RSCRIPT_EOF

    label "Generating annotation distribution plot(s) for $RESOLVED_CONTRAST..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "Outputs saved under: $ANALYSIS_OUT_DIR"
    else
        err "Annotation plot generation failed. See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 7 — GO/KEGG re-plot
# ─────────────────────────────────────────────────────────────

run_go_kegg_replot() {
    header "GO/KEGG Enrichment Re-plot"

    if [[ "$HAS_GO" != "TRUE" && "$HAS_KEGG" != "TRUE" ]]; then
        err "This bundle contains no saved GO or KEGG enrichment results to re-plot."
        return
    fi

    blank
    show_contrast_menu
    blank

    while true; do
        read -p "  Choose contrast (name or number): " CINPUT
        resolve_contrast "$CINPUT" && break
    done
    ok "Contrast: $RESOLVED_CONTRAST"

    blank
    echo -e "  Available result sets:"
    if [[ "$HAS_GO" == "TRUE" ]]; then
        echo -e "    ${CYAN}1${RESET}.  GO Biological Process (BP)"
        echo -e "    ${CYAN}2${RESET}.  GO Molecular Function (MF)"
        echo -e "    ${CYAN}3${RESET}.  GO Cellular Component (CC)"
    fi
    if [[ "$HAS_KEGG" == "TRUE" ]]; then
        echo -e "    ${CYAN}4${RESET}.  KEGG Pathways"
    fi
    echo -e "    ${CYAN}5${RESET}.  Run all available"
    blank

    local -a ONT_CHOICES ONT_SUFFIXES ONT_LABELS
    while true; do
        read -p "  Choice [1/2/3/4/5]: " ONT_IN
        case "$ONT_IN" in
            1) [[ "$HAS_GO" == "TRUE" ]] || { err "GO results are unavailable."; continue; }
               ONT_CHOICES=(BP); ONT_SUFFIXES=(Biological_Process); ONT_LABELS=("GO Biological Process"); break ;;
            2) [[ "$HAS_GO" == "TRUE" ]] || { err "GO results are unavailable."; continue; }
               ONT_CHOICES=(MF); ONT_SUFFIXES=(Molecular_Function); ONT_LABELS=("GO Molecular Function"); break ;;
            3) [[ "$HAS_GO" == "TRUE" ]] || { err "GO results are unavailable."; continue; }
               ONT_CHOICES=(CC); ONT_SUFFIXES=(Cellular_Component); ONT_LABELS=("GO Cellular Component"); break ;;
            4) [[ "$HAS_KEGG" == "TRUE" ]] || { err "KEGG results are unavailable."; continue; }
               ONT_CHOICES=(KEGG); ONT_SUFFIXES=(KEGG_pathways); ONT_LABELS=("KEGG Pathways"); break ;;
            5) ONT_CHOICES=(); ONT_SUFFIXES=(); ONT_LABELS=()
               if [[ "$HAS_GO" == "TRUE" ]]; then
                   ONT_CHOICES+=(BP MF CC)
                   ONT_SUFFIXES+=(Biological_Process Molecular_Function Cellular_Component)
                   ONT_LABELS+=("GO Biological Process" "GO Molecular Function" "GO Cellular Component")
               fi
               if [[ "$HAS_KEGG" == "TRUE" ]]; then
                   ONT_CHOICES+=(KEGG)
                   ONT_SUFFIXES+=(KEGG_pathways)
                   ONT_LABELS+=("KEGG Pathways")
               fi
               break ;;
            *) err "Enter 1, 2, 3, 4, or 5." ;;
        esac
    done

    blank
    read -p "  Max terms to show in dot plot [default: 20]: " SHOW_N
    SHOW_N="${SHOW_N:-20}"
    [[ "$SHOW_N" =~ ^[0-9]+$ ]] || { warn "Invalid — using 20."; SHOW_N=20; }

    local n_ont="${#ONT_CHOICES[@]}"
    local oi
    for oi in "${!ONT_CHOICES[@]}"; do
        if [[ "$n_ont" -gt 1 ]]; then
            blank
            label "[$((oi + 1))/$n_ont] ${ONT_LABELS[$oi]}"
        fi
        _run_go_kegg_one "${ONT_CHOICES[$oi]}" "${ONT_SUFFIXES[$oi]}" "${ONT_LABELS[$oi]}" "$SHOW_N"
    done
}

# Runs a single ontology's GO/KEGG re-plot for the already-resolved contrast.
# Split out from run_go_kegg_replot so "Run all" (option 5) can call this
# once per ontology without duplicating the R heredoc four times.
_run_go_kegg_one() {
    local ONT_CHOICE="$1"
    local ONT_FILE_SUFFIX="$2"
    local ONT_LABEL="$3"
    local SHOW_N="$4"

    local TAG="go_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/07_GO_KEGG/${SAFE_LABEL}_${ONT_FILE_SUFFIX}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_PDF="$ANALYSIS_OUT_DIR/${SAFE_LABEL}_${ONT_FILE_SUFFIX}_top${SHOW_N}_${TS}.pdf"
    local OUT_PNG="${OUT_PDF%.pdf}.png"
    local ALL_TERMS_CSV="$ANALYSIS_OUT_DIR/${SAFE_LABEL}_${ONT_FILE_SUFFIX}_all_terms_${TS}.csv"

    local RESOLVED_CONTRAST_R ALL_TERMS_CSV_R OUT_PDF_R OUT_PNG_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    ALL_TERMS_CSV_R="$(r_string_literal "$ALL_TERMS_CSV")"
    OUT_PDF_R="$(r_string_literal "$OUT_PDF")"
    OUT_PNG_R="$(r_string_literal "$OUT_PNG")"

    cat > "$R_SCRIPT" << RSCRIPT_EOF
suppressPackageStartupMessages({
    library(ggplot2)
    library(clusterProfiler)
    library(enrichplot)
})

.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b            <- readRDS(${BUNDLE_PATH_R})
label        <- ${RESOLVED_CONTRAST_R}
ont_choice   <- "$ONT_CHOICE"
ont_suffix   <- "$ONT_FILE_SUFFIX"
ont_label    <- "$ONT_LABEL"
show_n       <- $SHOW_N
diff_out     <- b\$diff_out

go_dir  <- file.path(diff_out, label, "GO")
csv_path <- file.path(go_dir, paste0(label, "_", ont_suffix, ".csv"))
all_terms_path <- ${ALL_TERMS_CSV_R}

if (!file.exists(csv_path)) {
    cat(sprintf(
        "  No %s results file found for contrast '%s' (GO/KEGG may not have been run, or no genes were annotated to differential peaks).\n",
        ont_label, label))
    cat(sprintf("  Expected: %s\n", csv_path))
    quit(status=0)
}

df <- read.csv(csv_path, check.names=FALSE)

# Identify the significance column once — used both to sort the full table
# and to decide which terms count as "significant" below.
sig_col <- if ("p.adjust" %in% names(df)) "p.adjust" else names(df)[grep("p.adjust|padj|FDR", names(df), ignore.case=TRUE)[1]]

if (nrow(df) == 0) {
    cat(sprintf("  %s CSV found but contains zero tested terms.\n", ont_label))
    quit(status=0)
}

# Sort by significance (best terms first) so the saved table is useful to
# skim even when nothing clears the FDR threshold.
if (!is.na(sig_col)) {
    df <- df[order(df[[sig_col]]), , drop=FALSE]
}

# Always save the FULL results table — every term that was tested, with
# its stats (p-value, p.adjust, gene ratio/count, etc.) — regardless of
# whether anything is significant. This is what to open and browse by
# hand when the dot plot below turns out empty.
write.csv(df, all_terms_path, row.names=FALSE)
cat(sprintf("  Full %s results table saved (%d terms tested): %s\n",
            ont_label, nrow(df), all_terms_path))

# Filter to significant terms (p.adjust < 0.05) for the dot plot only.
if (!is.na(sig_col)) {
    df_sig <- df[!is.na(df[[sig_col]]) & df[[sig_col]] < 0.05, , drop=FALSE]
} else {
    df_sig <- df
}

if (nrow(df_sig) == 0) {
    cat(sprintf("  No significant %s terms (FDR < 0.05) — skipping dot plot.\n", ont_label))
    if (!is.na(sig_col)) {
        best <- df[1, ]
        best_desc <- if ("Description" %in% names(best)) best[["Description"]] else "(no Description column)"
        cat(sprintf("  Closest term by %s: \"%s\" (%s = %.4g). Full ranked table has %d terms — see the CSV above.\n",
                    sig_col, best_desc, sig_col, best[[sig_col]], nrow(df)))
    }
    quit(status=0)
}

n_show <- min(show_n, nrow(df_sig))
cat(sprintf("  Plotting top %d / %d significant %s terms.\n",
            n_show, nrow(df_sig), ont_label))

# Re-build an enrichResult-like object for dotplot, or fall back to manual ggplot.
p <- tryCatch({
    # Reconstruct a minimal enrichResult so enrichplot::dotplot works.
    er <- new("enrichResult", result=df_sig, readable=TRUE,
              organism="UNKNOWN", ontology=ont_choice,
              keytype="ENTREZID", gene=character(0), universe=character(0),
              geneSets=list(), pvalueCutoff=0.05, pAdjustMethod="BH",
              qvalueCutoff=0.2, minGSSize=1L, maxGSSize=500L,
              gene2Symbol=character(0), params=list())
    dotplot(er, showCategory=n_show) +
        labs(title=paste0(label, " - ", ont_label),
             subtitle=sprintf("Top %d significant terms (FDR < 0.05)", n_show)) +
        theme_bw(base_size=11) +
        theme(axis.text.y=element_text(size=8))
}, error = function(e) {
    # Fallback: manual dot plot from the CSV columns.
    cat("  [NOTE] enrichResult reconstruction failed — using direct ggplot fallback.\n")
    desc_col  <- intersect(c("Description","description"), names(df_sig))[1]
    count_col <- intersect(c("Count","count","GeneNum"), names(df_sig))[1]
    ratio_col <- intersect(c("GeneRatio","generatio","gene_ratio"), names(df_sig))[1]
    if (is.na(desc_col) || is.na(sig_col)) stop("Cannot identify required columns.")
    plot_df <- df_sig[seq_len(n_show), ]
    plot_df[[desc_col]] <- factor(plot_df[[desc_col]],
                                   levels=rev(plot_df[[desc_col]]))
    p2 <- ggplot(plot_df, aes(x=if (!is.na(count_col)) .data[[count_col]] else seq_len(nrow(plot_df)),
                              y=.data[[desc_col]],
                              color=-log10(.data[[sig_col]]),
                              size=if (!is.na(count_col)) .data[[count_col]] else 1)) +
        geom_point() +
        scale_color_continuous(name="-log10(FDR)", low="grey80", high="darkred") +
        labs(title=paste0(label, " - ", ont_label),
             subtitle=sprintf("Top %d significant terms (FDR < 0.05)", n_show),
             x=if (!is.na(count_col)) "Gene count" else "Rank",
             y=NULL) +
        theme_bw(base_size=11) +
        theme(axis.text.y=element_text(size=8))
    p2
})

plot_h <- max(6, n_show * 0.4 + 2)
ggsave(${OUT_PDF_R}, p, width=10, height=plot_h)
ggsave(${OUT_PNG_R}, p, width=10, height=plot_h, dpi=150)
cat("GO/KEGG dot plot saved:", ${OUT_PDF_R}, "\n")
RSCRIPT_EOF

    label "Generating GO/KEGG dot plot ($ONT_LABEL)..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        if [[ -s "$OUT_PDF" ]]; then
            ok "Dot plot PDF: $OUT_PDF"
            ok "Dot plot PNG: $OUT_PNG"
        fi
        if [[ -s "$ALL_TERMS_CSV" ]]; then
            ok "Full results table (all terms tested, sorted by significance): $ALL_TERMS_CSV"
        fi
        if [[ ! -s "$OUT_PDF" && ! -s "$ALL_TERMS_CSV" ]]; then
            warn "No $ONT_LABEL results available for this contrast. See log: $LOG_FILE"
        elif [[ ! -s "$OUT_PDF" ]]; then
            warn "No significant $ONT_LABEL terms — dot plot skipped, but the full table above is there to browse."
        fi
    else
        err "GO/KEGG re-plot failed ($ONT_LABEL). See log: $LOG_FILE"
    fi
}

# ─────────────────────────────────────────────────────────────
# Analysis 8 — Tornado plots (deepTools)
# ─────────────────────────────────────────────────────────────

run_tornado() {
    # "individual" (one heatmap row/trace line per replicate) or
    # "composite" (mean per user-defined group, +/- SD shading on the
    # trace) -- chosen by which menu option the user picked, not by an
    # in-flow prompt. Kept as two separate menu entries (rather than one
    # option with a mode switch buried mid-prompt-chain) because the two
    # modes diverge in which downstream questions are even relevant
    # (sort-reference is a sample picker vs. a group picker; sample
    # ordering doesn't apply to composite at all) -- interleaving both
    # made the prompt sequence confusing to follow.
    local T_DISPLAY_MODE="${1:-individual}"
    if [[ "$T_DISPLAY_MODE" != "individual" && "$T_DISPLAY_MODE" != "composite" ]]; then
        err "Internal error: run_tornado called with unknown mode '$T_DISPLAY_MODE'."
        return
    fi

    if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
        header "Composite Tornado Plot"
    else
        header "Tornado Plot"
    fi

    # ── deepTools check ──────────────────────────────────────────
    # computeMatrix consumes the PEPATAC smoothShift bigWig tracks directly.
    # bamCoverage is not needed — PEPATAC already produced normalised bigWigs.
    if ! conda run --no-capture-output -n "$ENV_NAME" bash -c \
            "command -v computeMatrix && command -v plotHeatmap" >/dev/null 2>&1; then
        err "deepTools (computeMatrix / plotHeatmap) not found in the FetchPA environment."
        err "Run PEPATAC_install.sh to install deeptools."
        return
    fi

    # ── Contrast selection ───────────────────────────────────────
    blank
    show_contrast_menu
    blank

    while true; do
        read -p "  Choose contrast (name or number): " CINPUT
        resolve_contrast "$CINPUT" && break
    done
    ok "Contrast: $RESOLVED_CONTRAST"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"

    # Output dir/log set up now (not later) -- the composite-mode
    # auto-grouping step below needs somewhere to write its own small R
    # query before the rest of the prompt chain runs.
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local TAG="tornado_${RUN_ID}_$$"
    # Composite gets its own numbered top-level folder, matching every
    # other analysis's one-folder-per-menu-option convention (01_PCA ..
    # 07_GO_KEGG) -- #8 and #9 are separate menu entries now, so they
    # shouldn't share 08_Tornado_Plots just because they used to be one
    # option with a mode switch.
    local TORNADO_OUT
    if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
        TORNADO_OUT="$EXPLORE_OUT/09_Composite_Tornado_Plots/${SAFE_LABEL}_${TS}"
    else
        TORNADO_OUT="$EXPLORE_OUT/08_Tornado_Plots/${SAFE_LABEL}_${TS}"
    fi
    mkdir -p "$TORNADO_OUT"
    # Every intermediate/scratch file (prep scripts, logs, matrices, sorted
    # BEDs, group TSVs) lives one level deeper here, so TORNADO_OUT itself
    # only ever shows the final heatmap/trace PNGs and PDFs -- Windows
    # Explorer (unlike a Unix `ls`) does not hide dot-prefixed files, so
    # these were all showing up as clutter when browsed from the Windows
    # side even though most are already dot-prefixed.
    local TORNADO_COMPUTE_DIR="$TORNADO_OUT/files_for_computation"
    mkdir -p "$TORNADO_COMPUTE_DIR"
    local PREP_R="$TORNADO_COMPUTE_DIR/.${TAG}_prep.R"
    local PREP_OUT="$TORNADO_COMPUTE_DIR/.${TAG}_prep.tsv"
    local LOG_FILE="$TORNADO_COMPUTE_DIR/${TAG}.log"

    # ── Provenance: record the scaling workaround up front ─────────
    # This run compensates for a known deepTools computeMatrix limitation
    # (see TORNADO_SIGNAL_SCALE above) -- recorded here, not just applied
    # silently, since it materially affects the plotted values' units.
    local DEEPTOOLS_VERSION
    DEEPTOOLS_VERSION="$(conda run --no-capture-output -n "$ENV_NAME" computeMatrix --version 2>&1 | head -1)"
    {
        echo "[tornado-precision] deepTools version: ${DEEPTOOLS_VERSION:-unknown}"
        echo "[tornado-precision] scale factor applied (computeMatrix --scale, kept in plotted units, never unscaled): ${TORNADO_SIGNAL_SCALE}"
    } >> "$LOG_FILE"

    # Appended to every plot title below so scaled units are never mistaken
    # for raw smoothShift signal -- every tornado/profile plot from this
    # point on is plotted directly from the scaled matrix, permanently.
    local TORNADO_UNITS_NOTE=" (smoothShift signal x ${TORNADO_SIGNAL_SCALE})"

    # matplotlib's "Reds" colormap does NOT start at true white -- its
    # lightest step is #fff5f0 (RGB 255,245,240), a pale pink, since it's
    # built from ColorBrewer's sequential 9-class Reds palette. This is
    # that exact same palette with only the first (near-white) stop
    # replaced by true white, so the 0 floor renders as actual white while
    # every other step -- including the darkest red at the top end -- is
    # unchanged from what --colorMap Reds already produced. --colorList
    # overrides --colorMap entirely (deepTools does not support using both).
    local TORNADO_HEATMAP_COLORLIST="white,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d"

    # ── Centering ────────────────────────────────────────────────
    # IMPORTANT: "Gene TSS"/"Gene TES" below are NOT the same operation as
    # running deepTools --referencePoint TSS/TES directly on the peak BED.
    # A differential-peak interval is not a gene: its BED start/end are just
    # the two edges of that (typically sub-kb) peak, and this pipeline's own
    # peak BEDs carry no strand ("." in column 6, so deepTools would treat
    # every peak as "+"). Feeding peaks straight into TSS/TES mode produces
    # two panels that are both just "distance from a peak edge" -- for a
    # promoter mark like H3K4me3 those look nearly identical to each other,
    # which is NOT the same thing as seeing true, strand-correct TSS/TES
    # enrichment. Real gene anchoring requires real gene coordinates + real
    # strand, which is why the two gene-anchored options below instead read
    # the nearest-gene assignment (GeneChr/GeneStart/GeneEnd/GeneStrand)
    # that diff_analysis.sh's ChIPseeker/TxDb step already wrote per
    # significant peak into <contrast>_annotated_peaks.tsv, and build a
    # proper stranded gene BED from that.
    blank
    echo -e "  ${BOLD}Center on:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  Peak centers  ${DIM}(differential peak interval itself)${RESET}"
    echo -e "    ${CYAN}2${RESET}.  Gene TSS      ${DIM}(nearest gene per peak; true strand-aware TSS)${RESET}"
    echo -e "    ${CYAN}3${RESET}.  Gene TES      ${DIM}(nearest gene per peak; true strand-aware TES)${RESET}"
    blank
    local T_CENTER T_CENTER_LABEL T_GENE_ANCHOR
    while true; do
        read -p "  Choice [1/2/3, default 1]: " T_C
        T_C="${T_C:-1}"
        case "$T_C" in
            # deepTools accepts exactly TSS, TES, or center.  "midpoint" is
            # not a valid --referencePoint value.
            1) T_CENTER="center"; T_CENTER_LABEL="peak center"; T_GENE_ANCHOR=false; break ;;
            2) T_CENTER="TSS";    T_CENTER_LABEL="TSS";         T_GENE_ANCHOR=true;  break ;;
            3) T_CENTER="TES";    T_CENTER_LABEL="TES";         T_GENE_ANCHOR=true;  break ;;
            *) err "Enter 1, 2, or 3." ;;
        esac
    done
    ok "Centering: $T_CENTER_LABEL$($T_GENE_ANCHOR && echo " (gene-anchored)")"
    local T_GENE_ANCHOR_R
    T_GENE_ANCHOR_R=$($T_GENE_ANCHOR && echo "TRUE" || echo "FALSE")

    # ── Window size ──────────────────────────────────────────────
    blank
    echo -e "  ${BOLD}Window size${RESET} around center (kb each side)."
    echo -e "  ${DIM}e.g. 2 = ±2 kb (total 4 kb window).${RESET}"
    local T_KB T_BP
    while true; do
        read -p "  Window in kb [default: 2]: " T_KB
        T_KB="${T_KB:-2}"
        if [[ "$T_KB" =~ ^[0-9]+([.][0-9]+)?$ ]] && \
           awk "BEGIN{exit !($T_KB > 0 && $T_KB <= 50)}"; then
            T_BP=$(awk "BEGIN{printf \"%d\", $T_KB * 1000}")
            ok "Window: ±${T_KB} kb (${T_BP} bp each side)"
            break
        else
            err "Enter a number between 0.1 and 50."
        fi
    done

    # ── Row sorting ──────────────────────────────────────────────
    blank
    echo -e "  ${BOLD}Row sorting:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  Mean signal"
    echo -e "    ${CYAN}2${RESET}.  BED peak score"
    echo -e "    ${CYAN}3${RESET}.  Genomic position"
    echo -e "    ${CYAN}4${RESET}.  Preserve BED order"
    blank
    local T_SORT_MODE T_SORT_LABEL T_SORT_REGIONS T_SORT_USING
    while true; do
        read -p "  Choice [1/2/3/4, default 1]: " T_S
        T_S="${T_S:-1}"
        case "$T_S" in
            1) T_SORT_MODE="mean";       T_SORT_LABEL="mean signal";       T_SORT_REGIONS="descend"; T_SORT_USING="mean"; break ;;
            2) T_SORT_MODE="peak_score"; T_SORT_LABEL="BED peak score";    T_SORT_REGIONS="keep";    T_SORT_USING="";     break ;;
            3) T_SORT_MODE="genomic";    T_SORT_LABEL="genomic position";  T_SORT_REGIONS="keep";    T_SORT_USING="";     break ;;
            4) T_SORT_MODE="keep";       T_SORT_LABEL="input BED order";   T_SORT_REGIONS="keep";    T_SORT_USING="";     break ;;
            *) err "Enter 1, 2, 3, or 4." ;;
        esac
    done
    ok "Sort: $T_SORT_LABEL"

    # ── Sort reference (which samples define the shared row order) ──
    # Only meaningful for mean-based sorting (T_SORT_USING="mean") --
    # peak_score/genomic/keep modes never consult sample values at all.
    # deepTools sorts the matrix rows ONCE, producing one shared order
    # displayed across every panel (verified against its actual source:
    # one computeMatrix call + one plotHeatmap call means there is no way
    # for panels to be sorted independently of each other here). Left at
    # the default, that shared order is based on the mean across every
    # sample combined. Restricting --sortUsingSamples to a reference group
    # (e.g. WT replicates) instead orders every panel -- including the
    # mutant/treatment ones -- by what the reference looked like, which is
    # usually the more scientifically meaningful question for a two-
    # condition comparison. Reuses prompt_sample_selection() (same
    # numbers/names/group-name picker already used everywhere else in this
    # script) rather than a bespoke menu.
    local -a T_SORT_REF_SAMPLE_NAMES=()
    local -a T_COMPOSITE_GROUP_NAMES=()
    local -a T_COMPOSITE_GROUP_MEMBERS=()   # space-joined bundle sample names, one entry per group
    local T_SORT_REF_GROUP_NAME=""          # composite mode only: which group's mean sorts the rows

    if [[ "$T_DISPLAY_MODE" == "individual" ]]; then
        if [[ "$T_SORT_MODE" == "mean" ]]; then
            blank
            echo -e "  ${BOLD}Sort genes by mean signal in:${RESET}"
            echo -e "  ${DIM}Restricting this to a reference group (e.g. your WT replicates) orders${RESET}"
            echo -e "  ${DIM}every panel's rows the same way -- by what the reference looked like --${RESET}"
            echo -e "  ${DIM}instead of by all samples averaged together. Every panel still shows${RESET}"
            echo -e "  ${DIM}every sample's own signal; only the shared row ORDER changes.${RESET}"
            echo -e "  ${DIM}Type 'all' for the default (sort by all samples combined).${RESET}"
            prompt_sample_selection "Reference sample(s)/group to sort by:" 1
            if [[ "${#RESOLVED_SAMPLES[@]}" -lt "${#SAMPLE_NAMES[@]}" ]]; then
                T_SORT_REF_SAMPLE_NAMES=("${RESOLVED_SAMPLES[@]}")
                ok "Sorting by mean signal in: ${T_SORT_REF_SAMPLE_NAMES[*]}  (order applied to every panel)."
            else
                ok "Sorting by mean signal across all samples (default)."
            fi
        fi
    else
        # ── Composite groups ─────────────────────────────────────
        # The contrast picks the PEAK SET; which biological group(s) get
        # displayed over that peak set is a separate, independent choice --
        # a WT-only composite over DKO_vs_WT's differential peaks is a
        # completely valid figure on its own, with no requirement that DKO
        # also appear. So: look up this contrast's saved biological groups
        # (b$sample_ids / b$sample_groups / b$contrast_cases /
        # b$contrast_ctrls -- the same fields the tornado prep step already
        # uses) and let the user pick any ONE OR MORE of them by name or
        # number -- no need to "define" a group and re-pick its replicates
        # when the bundle already knows exactly who's in it. 'custom' opts
        # into the free-form manual builder for ad hoc groupings the
        # bundle doesn't already know about.
        local AUTO_GROUPS_R="$TORNADO_COMPUTE_DIR/.${TAG}_autogroups.R"
        local AUTO_GROUPS_TSV="$TORNADO_COMPUTE_DIR/.${TAG}_autogroups.tsv"
        local AUTO_GROUPS_LOG="$TORNADO_COMPUTE_DIR/.${TAG}_autogroups.log"
        local _RESOLVED_CONTRAST_R _AUTO_GROUPS_TSV_R
        _RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
        _AUTO_GROUPS_TSV_R="$(r_string_literal "$AUTO_GROUPS_TSV")"
        cat > "$AUTO_GROUPS_R" << RSCRIPT_EOF
.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b     <- readRDS(${BUNDLE_PATH_R})
label <- ${_RESOLVED_CONTRAST_R}
out   <- file(${_AUTO_GROUPS_TSV_R}, "w")
ci <- which(b\$contrast_labels == label)
if (length(ci) == 1 && !is.null(b\$contrast_cases) && !is.null(b\$contrast_ctrls) &&
    !is.null(b\$sample_ids) && !is.null(b\$sample_groups)) {
    case_grp <- b\$contrast_cases[ci]
    ctrl_grp <- b\$contrast_ctrls[ci]
    case_samples <- b\$sample_ids[b\$sample_groups == case_grp]
    ctrl_samples <- b\$sample_ids[b\$sample_groups == ctrl_grp]
    writeLines(paste("case", case_grp, paste(case_samples, collapse=" "), sep="\t"), out)
    writeLines(paste("ctrl", ctrl_grp, paste(ctrl_samples, collapse=" "), sep="\t"), out)
}
close(out)
cat("OK\n")
RSCRIPT_EOF
        run_r_script "$AUTO_GROUPS_R" "$AUTO_GROUPS_LOG" > /dev/null 2>&1 || true

        local T_AUTO_ROLE T_AUTO_NAME T_AUTO_MEMBERS
        local -a T_AUTO_GROUP_NAMES=() T_AUTO_GROUP_MEMBERS=()
        if [[ -s "$AUTO_GROUPS_TSV" ]]; then
            while IFS=$'\t' read -r T_AUTO_ROLE T_AUTO_NAME T_AUTO_MEMBERS; do
                [[ -z "$T_AUTO_NAME" || -z "$T_AUTO_MEMBERS" ]] && continue
                T_AUTO_GROUP_NAMES+=("$T_AUTO_NAME")
                T_AUTO_GROUP_MEMBERS+=("$T_AUTO_MEMBERS")
            done < "$AUTO_GROUPS_TSV"
        fi

        # No manual/custom grouping path: composite mode always uses
        # exactly what the bundle already knows about this contrast, one
        # independent run per group. If that's not available, this mode
        # simply isn't usable for this contrast/bundle -- fail clearly
        # rather than fall back to an interactive group-builder.
        local -a T_RUN_LABELS=()             # one entry per independent output run
        local -a T_RUN_GROUP_INDEX_SETS=()   # space-joined 0-based indices into T_COMPOSITE_GROUP_NAMES/MEMBERS, one entry per run
        if [[ "${#T_AUTO_GROUP_NAMES[@]}" -lt 1 ]]; then
            err "Could not auto-derive this contrast's biological groups from the bundle."
            err "Composite mode requires b\$contrast_cases / b\$contrast_ctrls / b\$sample_ids / b\$sample_groups"
            err "to be present and to resolve this contrast -- re-run PEPATAC_diff_analysis.sh to regenerate"
            err "the bundle if this contrast predates those fields."
            return
        fi
        T_COMPOSITE_GROUP_NAMES=("${T_AUTO_GROUP_NAMES[@]}")
        T_COMPOSITE_GROUP_MEMBERS=("${T_AUTO_GROUP_MEMBERS[@]}")
        blank
        echo -e "  ${BOLD}Composite groups for this contrast:${RESET}"
        local _gi3
        for _gi3 in "${!T_AUTO_GROUP_NAMES[@]}"; do
            local -a _gm=(${T_AUTO_GROUP_MEMBERS[$_gi3]})
            echo -e "    ${T_AUTO_GROUP_NAMES[$_gi3]}  ${DIM}(${#_gm[@]} replicate(s))${RESET}"
        done
        echo -e "  ${DIM}Each group gets its own independent output -- its own heatmap, its own${RESET}"
        echo -e "  ${DIM}trace -- never combined into one image. This run saves output for all${RESET}"
        echo -e "  ${DIM}${#T_AUTO_GROUP_NAMES[@]} of them, exactly as if you ran this contrast once per group.${RESET}"
        ok "Will generate independent composite output for: ${T_AUTO_GROUP_NAMES[*]}"
        for _gi3 in "${!T_COMPOSITE_GROUP_NAMES[@]}"; do
            T_RUN_LABELS+=("${T_COMPOSITE_GROUP_NAMES[$_gi3]}")
            T_RUN_GROUP_INDEX_SETS+=("$_gi3")
        done
        # Every auto-mode run is single-group by construction, so which
        # group sorts the rows is unambiguous (itself) -- no prompt needed;
        # T_SORT_REF_GROUP_NAME simply stays unset.
    fi

    # ── Sample grouping ──────────────────────────────────────────
    # This only reorders which column each individual replicate lands in
    # within computeMatrix's raw matrix -- irrelevant for composite mode,
    # which regroups/relabels by the user-picked groups afterward anyway,
    # so skip asking and just leave the underlying order untouched.
    local T_SPLIT=false
    if [[ "$T_DISPLAY_MODE" == "individual" ]]; then
        blank
        echo -e "  ${BOLD}Sample ordering:${RESET}"
        echo -e "    ${CYAN}1${RESET}.  Original sample order"
        echo -e "    ${CYAN}2${RESET}.  Order by condition group (case then control)"
        blank
        while true; do
            read -p "  Choice [1/2, default 2]: " T_G
            T_G="${T_G:-2}"
            case "$T_G" in
                1) T_SPLIT=false; break ;;
                2) T_SPLIT=true;  break ;;
                *) err "Enter 1 or 2." ;;
            esac
        done
    fi

    # ── Peak subset ──────────────────────────────────────────────
    blank
    echo -e "  ${BOLD}Peak subset:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  All significant peaks"
    echo -e "    ${CYAN}2${RESET}.  Opening peaks only (Up)"
    echo -e "    ${CYAN}3${RESET}.  Closing peaks only (Down)"
    echo -e "    ${CYAN}4${RESET}.  Opening + closing as separate row groups"
    blank
    local T_SUBSET
    while true; do
        read -p "  Choice [1/2/3/4, default 1]: " T_P
        T_P="${T_P:-1}"
        case "$T_P" in
            1) T_SUBSET="all";       break ;;
            2) T_SUBSET="up";        break ;;
            3) T_SUBSET="down";      break ;;
            4) T_SUBSET="split_dir"; break ;;
            *) err "Enter 1, 2, 3, or 4." ;;
        esac
    done

    # ── Threads ──────────────────────────────────────────────────
    blank
    read -p "  deepTools threads [default: 4]: " T_THREADS
    T_THREADS="${T_THREADS:-4}"
    [[ "$T_THREADS" =~ ^[0-9]+$ ]] && [[ "$T_THREADS" -ge 1 ]] || T_THREADS=4
    ok "Threads: $T_THREADS"

    # ── Locate BAMs and peaks from explorer bundle ───────────────
    local RESOLVED_CONTRAST_R PREP_OUT_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    PREP_OUT_R="$(r_string_literal "$PREP_OUT")"

    # Use R to extract BAM paths + groups + peak BEDs from the bundle.
    cat > "$PREP_R" << RSCRIPT_EOF
.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

b         <- readRDS(${BUNDLE_PATH_R})
label     <- ${RESOLVED_CONTRAST_R}
subset    <- "$T_SUBSET"
split_grp <- as.logical("$T_SPLIT")
diff_out  <- b\$diff_out
out_tsv   <- ${PREP_OUT_R}

# Locate significant peaks BED.
c_dir    <- file.path(diff_out, label)
sig_bed  <- file.path(c_dir, paste0(label, "_sig_peaks.bed"))
sig_csv  <- file.path(c_dir, paste0(label, "_significant_peaks.csv"))

# Prefer the directional BEDs diff_analysis.sh writes natively (already in
# correct 0-based BED coordinates) over reconstructing them here. Only
# fall back to reconstruction for older runs from before diff_analysis.sh
# wrote these directly.
native_up_bed <- file.path(c_dir, paste0(label, "_up_peaks.bed"))
native_dn_bed <- file.path(c_dir, paste0(label, "_down_peaks.bed"))
recon_up_bed  <- file.path(c_dir, "tornado", paste0(label, "_up_peaks.bed"))
recon_dn_bed  <- file.path(c_dir, "tornado", paste0(label, "_down_peaks.bed"))

if (file.exists(native_up_bed) || file.exists(native_dn_bed)) {
    up_bed <- native_up_bed
    dn_bed <- native_dn_bed
} else {
    up_bed <- recon_up_bed
    dn_bed <- recon_dn_bed

    # Reconstruct direction BEDs from the significant-peaks CSV -- only
    # reached for runs from before diff_analysis.sh wrote these natively.
    # DiffBind/GRanges coordinates are 1-based; BED starts are 0-based.
    # A prior version of this reconstruction wrote the raw Start value
    # straight through with no conversion, shifting every peak here 1bp
    # from where diff_analysis.sh's own BED writer (which does apply
    # pmax(0L, Start - 1L)) puts the same peak.
    if (subset %in% c("up","down","split_dir") && file.exists(sig_csv)) {
        dir.create(dirname(up_bed), showWarnings=FALSE, recursive=TRUE)
        df <- tryCatch(read.csv(sig_csv, stringsAsFactors=FALSE), error=function(e) NULL)
        if (!is.null(df)) {
            dir_col <- grep("^Direction\$", colnames(df), value=TRUE, ignore.case=TRUE)[1]
            chr_col <- grep("^Chr\$",       colnames(df), value=TRUE, ignore.case=TRUE)[1]
            sta_col <- grep("^Start\$",     colnames(df), value=TRUE, ignore.case=TRUE)[1]
            end_col <- grep("^End\$",       colnames(df), value=TRUE, ignore.case=TRUE)[1]
            lfc_col <- grep("^log2FoldChange\$|^Fold\$", colnames(df), value=TRUE, ignore.case=TRUE)[1]
            if (!is.na(dir_col) && !is.na(chr_col) && !is.na(sta_col) && !is.na(end_col)) {
                up_rows <- df[!is.na(df[[dir_col]]) & df[[dir_col]] == "Up",   , drop=FALSE]
                dn_rows <- df[!is.na(df[[dir_col]]) & df[[dir_col]] == "Down", , drop=FALSE]
                write_bed <- function(rows, path) {
                    if (nrow(rows) > 0) {
                        score <- if (!is.na(lfc_col)) abs(as.numeric(rows[[lfc_col]])) else rep(0, nrow(rows))
                        score[!is.finite(score)] <- 0
                        starts_0based <- pmax(0L, as.integer(rows[[sta_col]]) - 1L)
                        write.table(data.frame(rows[[chr_col]], starts_0based, rows[[end_col]],
                                               ".", score, "."),
                                    path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
                    }
                }
                write_bed(up_rows, up_bed)
                write_bed(dn_rows, dn_bed)
            }
        }
    }
}

# Determine BEDs to use.
gene_anchor <- ${T_GENE_ANCHOR_R}
beds  <- character(0)
blabs <- character(0)

if (gene_anchor) {
    # Gene-anchored TSS/TES: a differential peak is NOT a gene, and this
    # pipeline's peak BEDs carry no strand ("." in column 6), so running
    # deepTools --referencePoint TSS/TES directly on a peak BED just means
    # "distance from a peak edge" -- not a real, strand-correct gene TSS/TES.
    # Instead, build a proper gene BED from the nearest-gene assignment
    # (GeneChr/GeneStart/GeneEnd/GeneStrand) that diff_analysis.sh's
    # ChIPseeker/TxDb step already computed and wrote per significant peak.
    ann_path <- file.path(c_dir, paste0(label, "_annotated_peaks.tsv"))
    gene_bed_dir <- file.path(c_dir, "tornado")
    dir.create(gene_bed_dir, showWarnings=FALSE, recursive=TRUE)

    if (!file.exists(ann_path)) {
        cat("ERROR: Gene-anchored tornado needs", basename(ann_path), "which was not found.\n")
        cat("Re-run PEPATAC_diff_analysis.sh (ChIPseeker annotation step) for this contrast.\n")
        quit(status=1)
    }
    ann <- tryCatch(
        read.delim(ann_path, header=TRUE, stringsAsFactors=FALSE,
                   quote="", comment.char="", check.names=FALSE),
        error=function(e) NULL
    )
    if (is.null(ann) || nrow(ann) == 0) {
        cat("ERROR: Could not read", ann_path, "or it is empty.\n")
        quit(status=1)
    }

    # One row per unique gene: a gene can be the nearest gene to more than
    # one significant peak, and a duplicate BED row would just repeat an
    # identical TSS/TES profile in the heatmap. Keep, as each gene's
    # representative row, whichever of its peaks has the strongest effect
    # size, so "sort by BED score" still means something for this mode.
    # Genes with no coordinate/strand ("*", ambiguous-locus genes) or no
    # gene assigned at all are dropped -- they cannot be strand-anchored.
    build_gene_bed <- function(ann_df, out_path) {
        keep <- !is.na(ann_df\$GeneID) & nzchar(ann_df\$GeneID) &
                !is.na(ann_df\$GeneChr) & !is.na(ann_df\$GeneStart) & !is.na(ann_df\$GeneEnd) &
                !is.na(ann_df\$GeneStrand) & ann_df\$GeneStrand %in% c("+", "-")
        ann_df <- ann_df[keep, , drop=FALSE]
        if (nrow(ann_df) == 0) return(0L)

        lfc_col <- grep("^log2FoldChange\$|^Fold\$", colnames(ann_df), value=TRUE, ignore.case=TRUE)[1]
        gene_score <- if (!is.na(lfc_col)) abs(as.numeric(ann_df[[lfc_col]])) else rep(0, nrow(ann_df))
        gene_score[!is.finite(gene_score)] <- 0

        ord <- order(ann_df\$GeneID, -gene_score)
        ann_df <- ann_df[ord, , drop=FALSE]
        gene_score <- gene_score[ord]
        first_idx <- !duplicated(ann_df\$GeneID)
        genes <- ann_df[first_idx, , drop=FALSE]
        genes\$.score <- gene_score[first_idx]

        starts_0based <- pmax(0L, as.integer(genes\$GeneStart) - 1L)
        write.table(
            data.frame(genes\$GeneChr, starts_0based, as.integer(genes\$GeneEnd),
                       genes\$GeneID, genes\$.score, genes\$GeneStrand),
            out_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
        )
        nrow(genes)
    }

    dir_col <- grep("^Direction\$", colnames(ann), value=TRUE, ignore.case=TRUE)[1]
    subset_ann <- function(want) {
        if (is.na(dir_col)) return(ann[0, , drop=FALSE])
        ann[!is.na(ann[[dir_col]]) & ann[[dir_col]] == want, , drop=FALSE]
    }

    all_gbed <- file.path(gene_bed_dir, paste0(label, "_genes_all.bed"))
    up_gbed  <- file.path(gene_bed_dir, paste0(label, "_genes_up.bed"))
    dn_gbed  <- file.path(gene_bed_dir, paste0(label, "_genes_down.bed"))

    gene_counts <- integer(0)
    if (subset == "all") {
        n <- build_gene_bed(ann, all_gbed)
        beds <- all_gbed; blabs <- "Genes_near_All_significant_peaks"; gene_counts <- n
    } else if (subset == "up") {
        n <- build_gene_bed(subset_ann("Up"), up_gbed)
        beds <- up_gbed; blabs <- "Genes_near_Up_peaks"; gene_counts <- n
    } else if (subset == "down") {
        n <- build_gene_bed(subset_ann("Down"), dn_gbed)
        beds <- dn_gbed; blabs <- "Genes_near_Down_peaks"; gene_counts <- n
    } else if (subset == "split_dir") {
        n_up <- build_gene_bed(subset_ann("Up"), up_gbed)
        n_dn <- build_gene_bed(subset_ann("Down"), dn_gbed)
        beds <- c(up_gbed, dn_gbed)
        blabs <- c("Genes_near_Up_peaks", "Genes_near_Down_peaks")
        gene_counts <- c(n_up, n_dn)
    }

    keep_idx <- gene_counts > 0
    if (!any(keep_idx)) {
        cat("ERROR: No significant peaks had a usable stranded gene assignment for this subset -- nothing to plot.\n")
        quit(status=1)
    }
    if (any(!keep_idx)) {
        cat("NOTE: skipping empty gene set(s):", paste(blabs[!keep_idx], collapse=", "), "\n")
    }
    beds        <- beds[keep_idx]
    blabs       <- blabs[keep_idx]
    gene_counts <- gene_counts[keep_idx]
    for (gi in seq_along(beds)) {
        cat(sprintf("  %s: %d unique gene(s) with usable coordinates/strand\n", blabs[gi], gene_counts[gi]))
    }
} else if (subset == "all") {
    beds  <- sig_bed
    blabs <- "All_significant_peaks"
} else if (subset == "up") {
    beds  <- if (file.exists(up_bed) && file.size(up_bed) > 0) up_bed else sig_bed
    blabs <- if (file.exists(up_bed) && file.size(up_bed) > 0) "Opening_peaks_Up" else "All_significant_peaks"
} else if (subset == "down") {
    beds  <- if (file.exists(dn_bed) && file.size(dn_bed) > 0) dn_bed else sig_bed
    blabs <- if (file.exists(dn_bed) && file.size(dn_bed) > 0) "Closing_peaks_Down" else "All_significant_peaks"
} else if (subset == "split_dir") {
    beds  <- c(up_bed, dn_bed)
    blabs <- c("Opening_peaks_Up", "Closing_peaks_Down")
}

# Retrieve sample metadata from the bundle.
sample_ids    <- b\$sample_ids
sample_bams   <- b\$sample_bams
sample_groups <- b\$sample_groups

if (is.null(sample_ids) || is.null(sample_bams)) {
    cat("ERROR: Bundle does not contain sample BAM paths.\n")
    cat("Re-run PEPATAC_diff_analysis.sh (v0.9 or later) to regenerate the bundle.\n")
    quit(status=1)
}

# Find case/control group for this contrast.
ci       <- which(b\$contrast_labels == label)
case_grp <- b\$contrast_cases[ci]
ctrl_grp <- b\$contrast_ctrls[ci]

# Order samples: case first then control (if requested), else original order.
if (split_grp) {
    case_idx  <- which(sample_groups == case_grp)
    ctrl_idx  <- which(sample_groups == ctrl_grp)
    idx_order <- c(case_idx, ctrl_idx)
} else {
    idx_order <- seq_along(sample_ids)
}

out_rows <- data.frame(
    SampleID = sample_ids[idx_order],
    BAM      = sample_bams[idx_order],
    Group    = sample_groups[idx_order],
    stringsAsFactors = FALSE
)

bed_meta <- data.frame(
    SampleID = paste0("__BED__", seq_along(beds)),
    BAM      = beds,
    Group    = blabs,
    stringsAsFactors = FALSE
)

write.table(rbind(bed_meta, out_rows), out_tsv,
            sep="\t", quote=FALSE, row.names=FALSE)
cat("OK\n")
RSCRIPT_EOF

    label "Preparing sample/peak metadata from bundle..."
    if ! run_r_script "$PREP_R" "$LOG_FILE"; then
        err "Failed to extract metadata from bundle. See log: $LOG_FILE"
        return
    fi

    if [[ ! -f "$PREP_OUT" ]]; then
        err "Metadata file not created. See log: $LOG_FILE"
        return
    fi

    # Parse the TSV into bash arrays.
    declare -a T_BEDS=() T_BED_LABELS=() T_BAMS=() T_SAMPLE_IDS=() T_GROUPS=()
    while IFS=$'\t' read -r sid bam grp; do
        [[ "$sid" == "SampleID" ]] && continue
        if [[ "$sid" == __BED__* ]]; then
            T_BEDS+=("$bam")
            T_BED_LABELS+=("$grp")
        else
            T_SAMPLE_IDS+=("$sid")
            T_BAMS+=("$bam")
            T_GROUPS+=("$grp")
        fi
    done < "$PREP_OUT"

    if [[ ${#T_BAMS[@]} -eq 0 ]]; then
        err "No BAM paths found in bundle for contrast '$RESOLVED_CONTRAST'."
        return
    fi

    # Verify BAMs exist and are indexed.
    declare -a VALID_BAMS=() VALID_IDS=() VALID_GROUPS=()
    local i bam
    for i in "${!T_BAMS[@]}"; do
        bam="${T_BAMS[$i]}"
        if [[ ! -f "$bam" ]]; then
            warn "BAM not found (skipping): $bam"
            continue
        fi
        if [[ ! -f "${bam}.bai" ]] && [[ ! -f "${bam%.bam}.bai" ]] && \
           [[ ! -f "${bam}.csi" ]]; then
            label "Indexing: $(basename "$bam")..."
            conda run --no-capture-output -n "$ENV_NAME" samtools index "$bam" >> "$LOG_FILE" 2>&1 || {
                warn "Could not index $(basename "$bam") — skipping."
                continue
            }
        fi
        VALID_BAMS+=("$bam")
        VALID_IDS+=("${T_SAMPLE_IDS[$i]}")
        VALID_GROUPS+=("${T_GROUPS[$i]}")
    done

    if [[ ${#VALID_BAMS[@]} -eq 0 ]]; then
        err "No valid BAM files found — aborting tornado."
        return
    fi

    # computeMatrix requires bigWig score files.  Use the smoothShift bigWigs
    # produced by PEPATAC during sample processing — these are ATAC-seq optimised
    # (Tn5 shift-corrected, read-count normalised) and live in the same aligned_*
    # folder as each BAM, named <SampleID>_smooth_shift.bw.
    declare -a VALID_BIGWIGS=() SCORE_IDS=() SCORE_GROUPS=()

    blank
    label "Locating PEPATAC smoothShift bigWig tracks..."
    for i in "${!VALID_BAMS[@]}"; do
        bam="${VALID_BAMS[$i]}"
        local bw_dir bw
        bw_dir="$(dirname "$bam")"
        bw="${bw_dir}/${VALID_IDS[$i]}_smooth_shift.bw"

        if [[ -s "$bw" ]]; then
            ok "Found: $(basename "$bw")"
        else
            err "smoothShift bigWig not found for ${VALID_IDS[$i]}: $bw"
            err "Expected PEPATAC to have written this file during sample processing."
            return
        fi

        VALID_BIGWIGS+=("$bw")
        SCORE_IDS+=("${VALID_IDS[$i]}")
        SCORE_GROUPS+=("${VALID_GROUPS[$i]}")
    done

    if [[ ${#VALID_BIGWIGS[@]} -eq 0 ]]; then
        err "No smoothShift bigWig tracks found — aborting tornado."
        return
    fi

    # Short, human-readable figure labels -- "<Group> <n>" per replicate
    # within its group (e.g. "WT 1", "WT 2", "D3a 1") instead of the full
    # raw SampleID (e.g. "WT_B1_H3K4me3_S13_SRR11785444"), which collides
    # badly across 6+ narrow deepTools panels. Full SampleIDs are still in
    # SCORE_IDS/the log for traceability -- this is purely a figure-display
    # label, built once here and reused wherever --samplesLabel is needed.
    declare -a SCORE_SHORT_LABELS=()
    declare -A GROUP_REP_COUNTER=()
    for i in "${!SCORE_IDS[@]}"; do
        local _grp="${SCORE_GROUPS[$i]:-sample}"
        GROUP_REP_COUNTER["$_grp"]=$(( ${GROUP_REP_COUNTER["$_grp"]:-0} + 1 ))
        SCORE_SHORT_LABELS+=("${_grp} ${GROUP_REP_COUNTER[$_grp]}")
    done

    # Bigwig paths in the exact order passed to --scoreFileName, one per
    # line -- the validation script below needs this to open the right
    # bigwig for each matrix sample column when independently re-deriving
    # values to check the scaled matrix against.
    local TORNADO_BIGWIGS_TXT="$TORNADO_COMPUTE_DIR/.${TAG}_bigwigs.txt"
    printf '%s\n' "${VALID_BIGWIGS[@]}" > "$TORNADO_BIGWIGS_TXT"

    # ── Composite group resolution (composite mode only) ──────────
    # T_COMPOSITE_GROUP_NAMES/MEMBERS above were collected against the
    # bundle-wide SAMPLE_NAMES before BAM/bigwig validation could drop any.
    # Resolve each group's members to their 0-based position within
    # SCORE_IDS (the order bigwigs were actually handed to computeMatrix),
    # same pattern as the individual-mode sort-reference resolution above.
    # The rule here protects USER INTENT, not a fixed group count: every
    # group the user explicitly asked to see must survive validation with
    # at least one replicate, or the whole run aborts. A single requested
    # group (e.g. "WT" alone) surviving is a complete success, not a
    # partial one -- there is no minimum group count. What must NOT happen
    # is silently dropping a group the user asked for and plotting
    # whatever's left, which is the actual failure mode to guard against.
    declare -a T_COMPOSITE_GROUP_INDICES=()  # space-joined 0-based SCORE_IDS indices, one entry per surviving group
    declare -a T_COMPOSITE_GROUP_LABELS=()   # "<name> (n=<k>)" display label, one entry per surviving group
    local T_SORT_REF_GROUP_INDEX=""          # 1-based index into the surviving groups list, or "" for none
    if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
        local gi
        local -a _dropped_groups=()
        for gi in "${!T_COMPOSITE_GROUP_NAMES[@]}"; do
            local _gname="${T_COMPOSITE_GROUP_NAMES[$gi]}"
            local -a _members=(${T_COMPOSITE_GROUP_MEMBERS[$gi]})
            local -a _idx=()
            local _m _j
            for _m in "${_members[@]}"; do
                for _j in "${!SCORE_IDS[@]}"; do
                    if [[ "${SCORE_IDS[$_j]}" == "$_m" ]]; then
                        _idx+=("$_j")
                        break
                    fi
                done
            done
            if [[ "${#_idx[@]}" -eq 0 ]]; then
                _dropped_groups+=("$_gname")
                continue
            fi
            T_COMPOSITE_GROUP_INDICES+=("${_idx[*]}")
            T_COMPOSITE_GROUP_LABELS+=("${_gname} (n=${#_idx[@]})")
            if [[ -n "$T_SORT_REF_GROUP_NAME" && "$T_SORT_REF_GROUP_NAME" == "$_gname" ]]; then
                T_SORT_REF_GROUP_INDEX="${#T_COMPOSITE_GROUP_INDICES[@]}"
            fi
        done

        if [[ "${#_dropped_groups[@]}" -gt 0 ]]; then
            err "Requested group(s) with NO replicates surviving BAM/bigwig validation: ${_dropped_groups[*]}"
            err "Aborting -- every group you asked to display must survive; none are silently dropped."
            return
        fi
        if [[ -n "$T_SORT_REF_GROUP_NAME" && -z "$T_SORT_REF_GROUP_INDEX" ]]; then
            warn "Sort-reference group '${T_SORT_REF_GROUP_NAME}' had no valid replicates -- sorting by all groups instead."
        fi
    fi

    # ── Matrix precision validation script (both modes) ─────────────
    # The matrix computeMatrix just writes (once --scale is wired in below)
    # stays in scaled units PERMANENTLY -- it is never divided back down or
    # rewritten. This script only independently re-derives a sample of
    # matrix positions straight from the source bigWigs (replicating
    # computeMatrix's own reference-point + bin-averaging +
    # --missingDataAsZero logic) and compares matrix_value / scale against
    # that, to confirm the scaled matrix actually represents the source
    # signal before anything is plotted from it. A failure here is
    # treated as fatal by the caller (see the per-BED loop below) --
    # this is a correctness gate, not an FYI. Generated once, reused per
    # BED (only --matrix changes).
    local VALIDATE_PY="$TORNADO_COMPUTE_DIR/.${TAG}_validate.py"
    cat > "$VALIDATE_PY" << 'VALIDATEPYEOF'
#!/usr/bin/env python3
import argparse, gzip, json, random, sys
import numpy as np

def read_matrix(path):
    with gzip.open(path, "rt") as fh:
        header_line = fh.readline()
        if not header_line.startswith("@"):
            sys.exit(f"ERROR: {path} does not look like a deepTools computeMatrix output.")
        header = json.loads(header_line[1:])
        rows, meta = [], []
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            meta.append(fields[:6])
            rows.append(np.array(
                [np.nan if v == "nan" else float(v) for v in fields[6:]], dtype=float
            ))
    return header, meta, (np.vstack(rows) if rows else np.empty((0, 0)))

def ref_point(mode, start, end, strand):
    if mode == "center":
        return (start + end) // 2
    if mode == "TSS":
        return start if strand == "+" else end
    return end if strand == "+" else start  # TES

def fail(msg):
    print(f"[tornado-precision] {msg}")
    print("[tornado-precision] Precision validation: FAIL")
    print("[tornado-precision] Tornado generation aborted -- matrix does not faithfully represent source bigWig signal.")
    sys.exit(1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--scale", required=True, type=float)
    ap.add_argument("--bigwigs-file", required=True)
    ap.add_argument("--center-mode", required=True, choices=["center", "TSS", "TES"])
    ap.add_argument("--n-large", type=int, default=10, help="largest-magnitude nonzero matrix bins to check")
    ap.add_argument("--n-random-nonzero", type=int, default=10, help="random nonzero matrix bins to check")
    ap.add_argument("--n-zero", type=int, default=15,
                     help="random ZERO-valued matrix bins to check -- guards against the exact "
                          "failure mode this validator exists for: computeMatrix silently "
                          "serializing real signal as 0.000000")
    ap.add_argument("--tol-rel", type=float, default=0.01)
    ap.add_argument("--tol-abs", type=float, default=2e-12,
                     help="with scale=1e6 and computeMatrix's 6-decimal serialization, the "
                          "theoretical quantization floor in real units is ~1e-12 -- this leaves "
                          "a little headroom above that, not a loose absolute allowance")
    args = ap.parse_args()

    header, meta, scaled_body = read_matrix(args.matrix)

    with open(args.bigwigs_file) as fh:
        bigwig_paths = [line.strip() for line in fh if line.strip()]

    sb = header["sample_boundaries"]
    n_samples = len(sb) - 1
    if n_samples != len(bigwig_paths):
        fail(f"sample count mismatch ({n_samples} matrix samples vs {len(bigwig_paths)} bigwig paths)")
    upstream = header.get("upstream", [0])[0]
    bin_size = header.get("bin size", [10])[0]

    import pyBigWig
    bw_cache = {}
    def get_bw(i):
        if i not in bw_cache:
            bw_cache[i] = pyBigWig.open(bigwig_paths[i])
        return bw_cache[i]

    # Stratified candidate selection -- this is the core of what this
    # validator exists to catch. Checking only nonzero matrix values would
    # prove "the numbers that survived serialization are accurate" but say
    # nothing about whether OTHER positions were real signal that got
    # silently zeroed out -- which is exactly the bug this whole precision
    # workaround was built for. So zero-valued matrix bins are sampled and
    # checked too: if the source bigWig shows real signal at a position the
    # matrix claims is zero, that is a hard failure.
    nonzero_candidates, zero_candidates = [], []
    for s in range(n_samples):
        block = scaled_body[:, sb[s]:sb[s + 1]]
        nz_rows, nz_bins = np.nonzero(block)
        for r, b in zip(nz_rows.tolist(), nz_bins.tolist()):
            nonzero_candidates.append((r, s, b, float(block[r, b])))
        z_rows, z_bins = np.where(block == 0.0)
        for r, b in zip(z_rows.tolist(), z_bins.tolist()):
            zero_candidates.append((r, s, b, 0.0))

    if not nonzero_candidates:
        fail("no nonzero values anywhere in the matrix to validate against -- cannot confirm signal integrity")

    nonzero_candidates.sort(key=lambda c: -abs(c[3]))
    n_large = min(len(nonzero_candidates), args.n_large)
    chosen = nonzero_candidates[:n_large]
    remaining_nonzero = nonzero_candidates[n_large:]
    n_rand_nz = min(len(remaining_nonzero), args.n_random_nonzero)
    if n_rand_nz > 0:
        chosen += random.sample(remaining_nonzero, n_rand_nz)

    n_zero = min(len(zero_candidates), args.n_zero)
    if n_zero > 0:
        chosen += random.sample(zero_candidates, n_zero)

    max_abs_err, max_rel_err, tested, all_pass = 0.0, 0.0, 0, True
    for row_idx, sample_idx, bin_idx, scaled_actual in chosen:
        chrom, start, end, name, score, strand = meta[row_idx]
        start, end = int(start), int(end)
        ref = ref_point(args.center_mode, start, end, strand)
        wstart, wend = ref - upstream + bin_idx * bin_size, ref - upstream + (bin_idx + 1) * bin_size
        bw = get_bw(sample_idx)
        chrom_len = bw.chroms().get(chrom)
        if chrom_len is None or wstart < 0 or wend > chrom_len:
            continue  # out-of-bounds/unknown-chrom window -- can't independently check this one
        raw = np.nan_to_num(np.array(bw.values(chrom, wstart, wend), dtype=float))  # mirrors --missingDataAsZero
        expected = float(raw.mean()) if len(raw) else 0.0  # real, unscaled units

        actual = scaled_actual / args.scale  # undo the scale for comparison only -- the matrix itself is untouched
        abs_err = abs(actual - expected)
        rel_err = (abs_err / abs(expected)) if expected != 0 else (0.0 if actual == 0 else float("inf"))
        max_abs_err = max(max_abs_err, abs_err)
        if rel_err != float("inf"):
            max_rel_err = max(max_rel_err, rel_err)
        ok_here = abs_err <= (args.tol_abs + args.tol_rel * abs(expected))
        all_pass = all_pass and ok_here
        tested += 1
        print(f"[tornado-precision] validate row={row_idx} sample={sample_idx} bin={bin_idx}: "
              f"matrix/scale={actual:.6e} pyBigWig={expected:.6e} rel_err={rel_err:.3e} {'OK' if ok_here else 'MISMATCH'}")

    if tested == 0:
        fail("every sampled position fell outside a valid chromosome window -- cannot confirm signal integrity")

    print(f"[tornado-precision] {tested} randomly selected sample/region/bin values checked")
    print(f"[tornado-precision] Scale factor: {args.scale:g}")
    print(f"[tornado-precision] Maximum relative difference: {max_rel_err:.3e}")
    print(f"[tornado-precision] Maximum absolute difference: {max_abs_err:.3e}")
    print(f"[tornado-precision] Precision validation: {'PASS' if all_pass else 'FAIL'}")
    if not all_pass:
        print("[tornado-precision] Tornado generation aborted -- matrix does not faithfully represent source bigWig signal.")
        sys.exit(1)

if __name__ == "__main__":
    main()
VALIDATEPYEOF

    # ── Composite averaging/plotting script (composite mode only) ──
    # Generated once and reused for every BED below (only --matrix/--out
    # paths change per BED). Reads the deepTools matrix.gz that computeMatrix
    # already produced (JSON header line + tab-delimited body -- the exact
    # same file individual mode feeds straight to plotHeatmap/plotProfile),
    # averages the requested sample columns per group into a new matrix.gz
    # with one sample-block per group (which then drops straight into an
    # unmodified plotHeatmap call), and separately computes each group's
    # per-bin mean +/- SD across its replicates' own region-averaged
    # profiles for the shaded composite trace plot. numpy/matplotlib are
    # deepTools' own dependencies, so nothing new is required of the
    # FetchPA environment.
    local COMPOSITE_PY=""
    if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
        COMPOSITE_PY="$TORNADO_COMPUTE_DIR/.${TAG}_composite.py"
        cat > "$COMPOSITE_PY" << 'PYEOF'
#!/usr/bin/env python3
import argparse, gzip, json, sys
import numpy as np

def read_matrix(path):
    with gzip.open(path, "rt") as fh:
        header_line = fh.readline()
        if not header_line.startswith("@"):
            sys.exit(f"ERROR: {path} does not look like a deepTools computeMatrix output.")
        header = json.loads(header_line[1:])
        rows = []
        meta = []
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            meta.append(fields[:6])
            vals = np.array(
                [np.nan if v == "nan" else float(v) for v in fields[6:]],
                dtype=float,
            )
            rows.append(vals)
    return header, meta, np.vstack(rows) if rows else np.empty((0, 0))

def write_matrix(path, header, meta, body):
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        for m, row in zip(meta, body):
            vals = ["nan" if np.isnan(v) else repr(float(v)) for v in row]
            fh.write("\t".join(list(m) + vals) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--groups-file", required=True,
                     help="TSV: label<TAB>space-joined 0-based sample indices, one group per line")
    ap.add_argument("--composite-matrix-out", required=True)
    ap.add_argument("--profile-tsv-out", required=True)
    ap.add_argument("--profile-plot-png", required=True)
    ap.add_argument("--profile-plot-pdf", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--ref-label", default="center")
    ap.add_argument("--units-note", default="",
                     help="appended to the y-axis label and TSV header comment, e.g. ' (smoothShift signal x 1000000)'")
    args = ap.parse_args()

    groups = []
    with open(args.groups_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            label, idx_str = line.rstrip("\n").split("\t")
            groups.append({"label": label, "indices": [int(x) for x in idx_str.split()]})
    header, meta, body = read_matrix(args.matrix)

    sb = header["sample_boundaries"]
    n_samples = len(sb) - 1
    n_bins = sb[1] - sb[0]
    if any((sb[k + 1] - sb[k]) != n_bins for k in range(n_samples)):
        sys.exit("ERROR: sample blocks have unequal bin counts -- cannot average across samples.")

    def sample_block(sample_idx):
        start = sb[sample_idx]
        end = sb[sample_idx + 1]
        return body[:, start:end]

    # ── Composite matrix: one averaged sample-block per group, feeds
    # straight into plotHeatmap unmodified. ──
    comp_blocks = []
    for g in groups:
        blocks = np.stack([sample_block(i) for i in g["indices"]], axis=0)  # (reps, regions, bins)
        comp_blocks.append(np.nanmean(blocks, axis=0))
    comp_body = np.hstack(comp_blocks)

    comp_header = dict(header)
    comp_header["sample_labels"] = [g["label"] for g in groups]
    new_bounds = [0]
    for _ in groups:
        new_bounds.append(new_bounds[-1] + n_bins)
    comp_header["sample_boundaries"] = new_bounds
    write_matrix(args.composite_matrix_out, comp_header, meta, comp_body)

    # ── Composite profile: per group, per bin, mean +/- SD across
    # replicates' own region-averaged profiles. ──
    upstream = header.get("upstream", [0])[0]
    bin_size = header.get("bin size", [10])[0]
    bin_centers = -upstream + (np.arange(n_bins) + 0.5) * bin_size

    # Column names carry the units note directly (rather than a leading
    # "#" comment line) so the file stays a plain, single-header TSV that
    # naive readers (pandas/R read.delim with no comment handling) parse
    # correctly without silently misreading a comment row as data.
    mean_col = "mean" + args.units_note.replace(" ", "_").replace("(", "").replace(")", "")
    sd_col = "sd" + args.units_note.replace(" ", "_").replace("(", "").replace(")", "")
    with open(args.profile_tsv_out, "w") as fh:
        fh.write(f"group\tbin_center\t{mean_col}\t{sd_col}\tn\n")
        for g in groups:
            reps = [np.nanmean(sample_block(i), axis=0) for i in g["indices"]]  # each: (bins,)
            reps = np.vstack(reps)  # (reps, bins)
            mean = np.nanmean(reps, axis=0)
            sd = np.nanstd(reps, axis=0, ddof=1) if reps.shape[0] > 1 else np.zeros(n_bins)
            for b in range(n_bins):
                # ATAC-seq bigWig signal (RPGC/CPM-normalized) is commonly
                # ~1e-9 in magnitude; ".6f" silently rounds every such value
                # to 0.000000. ".6g" keeps 6 significant digits and falls
                # back to scientific notation for small/large magnitudes.
                fh.write(f"{g['label']}\t{bin_centers[b]:.1f}\t{mean[b]:.6g}\t{sd[b]:.6g}\t{reps.shape[0]}\n")

    # ── Shaded composite trace plot ──
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.get_cmap("tab10").colors
    for gi, g in enumerate(groups):
        reps = [np.nanmean(sample_block(i), axis=0) for i in g["indices"]]
        reps = np.vstack(reps)
        mean = np.nanmean(reps, axis=0)
        sd = np.nanstd(reps, axis=0, ddof=1) if reps.shape[0] > 1 else np.zeros(n_bins)
        color = colors[gi % len(colors)]
        ax.plot(bin_centers, mean, label=g["label"], color=color, linewidth=1.8)
        ax.fill_between(bin_centers, mean - sd, mean + sd, color=color, alpha=0.2, linewidth=0)

    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel(f"Distance from {args.ref_label} (bp)")
    ax.set_ylabel(f"Mean signal{args.units_note}")
    ax.set_title(args.title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.profile_plot_png, dpi=200)
    fig.savefig(args.profile_plot_pdf)

if __name__ == "__main__":
    main()
PYEOF
    fi

    # One groups-file per RUN (not one global file) -- each run only feeds
    # COMPOSITE_PY the subset of groups it's supposed to show, per
    # T_RUN_GROUP_INDEX_SETS built above (one index set per independent
    # output; auto mode has one single-group set per run, custom mode has
    # exactly one set containing everything defined).
    declare -a T_RUN_GROUPS_TSV=()
    if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
        local _ri
        for _ri in "${!T_RUN_LABELS[@]}"; do
            local _run_tsv="$TORNADO_COMPUTE_DIR/.${TAG}_composite_groups_${_ri}.tsv"
            : > "$_run_tsv"
            local -a _run_idx_set=(${T_RUN_GROUP_INDEX_SETS[$_ri]})
            local _gidx
            for _gidx in "${_run_idx_set[@]}"; do
                printf '%s\t%s\n' "${T_COMPOSITE_GROUP_LABELS[$_gidx]}" "${T_COMPOSITE_GROUP_INDICES[$_gidx]}" \
                    >> "$_run_tsv"
            done
            T_RUN_GROUPS_TSV+=("$_run_tsv")
        done
    fi

    # ── Per-BED tornado generation ───────────────────────────────
    local bi
    for bi in "${!T_BEDS[@]}"; do
        local BED_FILE="${T_BEDS[$bi]}"
        local BED_LABEL="${T_BED_LABELS[$bi]}"

        if [[ ! -f "$BED_FILE" ]] || [[ $(wc -l < "$BED_FILE") -lt 5 ]]; then
            warn "BED for '$BED_LABEL' is missing or too small — skipping."
            continue
        fi

        # Prepare a correctly ordered BED when the requested order is based on
        # the BED itself.  Mean-signal ordering is handled by plotHeatmap.
        local BED_FOR_MATRIX="$BED_FILE"
        local ORDERED_BED=""
        case "$T_SORT_MODE" in
            peak_score)
                ORDERED_BED="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_ordered_by_peak_score_${TS}.bed"
                awk 'BEGIN{OFS="\t"} NF>=3 {s=(NF>=5 && $5 ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$/) ? $5 : 0; print s,$0}' \
                    "$BED_FILE" | sort -k1,1gr | cut -f2- > "$ORDERED_BED"
                BED_FOR_MATRIX="$ORDERED_BED"
                local N_SCORES
                N_SCORES=$(awk '{print (NF>=5 ? $5 : 0)}' "$BED_FILE" | sort -u | wc -l)
                if [[ "$N_SCORES" -le 1 ]]; then
                    warn "BED scores are all identical for '$BED_LABEL'; peak-score order will preserve ties arbitrarily."
                fi
                ;;
            genomic)
                ORDERED_BED="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_ordered_genomically_${TS}.bed"
                sort -k1,1V -k2,2n -k3,3n "$BED_FILE" > "$ORDERED_BED"
                BED_FOR_MATRIX="$ORDERED_BED"
                ;;
            keep|mean)
                ;;
        esac

        local N_PEAKS
        N_PEAKS=$(wc -l < "$BED_FOR_MATRIX")
        local MATRIX_OUT="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_matrix.gz"
        local HEATMAP_PNG="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_tornado_${TS}.png"
        local HEATMAP_PDF="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_tornado_${TS}.pdf"
        local SORTED_BED="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_sorted_regions_${TS}.bed"

        blank
        label "Computing matrix: $BED_LABEL ($N_PEAKS regions, ${#VALID_BIGWIGS[@]} bigWigs)..."
        echo -e "  ${DIM}Log: $LOG_FILE${RESET}"

        set +e
        conda run --no-capture-output -n "$ENV_NAME" \
            computeMatrix reference-point \
                --referencePoint "$T_CENTER" \
                --beforeRegionStartLength "$T_BP" \
                --afterRegionStartLength  "$T_BP" \
                --regionsFileName "$BED_FOR_MATRIX" \
                --scoreFileName   "${VALID_BIGWIGS[@]}" \
                --outFileName     "$MATRIX_OUT" \
                --numberOfProcessors "$T_THREADS" \
                --skipZeros \
                --missingDataAsZero \
                --scale "$TORNADO_SIGNAL_SCALE" \
            >> "$LOG_FILE" 2>&1
        local MATRIX_EXIT=$?
        set -e

        if [[ "$MATRIX_EXIT" -ne 0 ]]; then
            err "computeMatrix failed for '$BED_LABEL'. Log: $LOG_FILE"
            continue
        fi
        ok "Matrix written: $(basename "$MATRIX_OUT")  ${DIM}(scaled x${TORNADO_SIGNAL_SCALE} -- see below)${RESET}"

        # ── Precision validation (hard gate, not advisory) ─────────
        # The matrix above is pre-scaled by TORNADO_SIGNAL_SCALE so real
        # signal survives computeMatrix's own 6-decimal text serialization
        # (see that constant's definition up top). It is NOT divided back
        # down -- both #8 and #9 plot directly from these scaled values,
        # which is why every title/label below says so. Before anything
        # reads this file, independently re-derive a STRATIFIED sample of
        # matrix positions straight from the source bigWigs -- large
        # nonzero values, random nonzero values, AND random ZERO-valued
        # matrix positions (checking only nonzero values would prove the
        # numbers that survived are accurate but say nothing about whether
        # real signal got silently zeroed out, which is the actual bug
        # this whole workaround exists for) -- and require every one to
        # match (matrix value / scale) within tolerance. A failure here
        # aborts this BED's plot rather than silently plotting unverified
        # numbers.
        label "Validating scaled matrix against source bigWigs..."
        set +e
        local VALIDATE_OUTPUT
        VALIDATE_OUTPUT=$(conda run --no-capture-output -n "$ENV_NAME" python3 "$VALIDATE_PY" \
            --matrix "$MATRIX_OUT" \
            --scale "$TORNADO_SIGNAL_SCALE" \
            --bigwigs-file "$TORNADO_BIGWIGS_TXT" \
            --center-mode "$T_CENTER" \
            --n-large 10 \
            --n-random-nonzero 10 \
            --n-zero 15 2>&1)
        local VALIDATE_EXIT=$?
        set -e
        echo "$VALIDATE_OUTPUT" >> "$LOG_FILE"

        if [[ "$VALIDATE_EXIT" -ne 0 ]]; then
            err "Tornado generation aborted for '$BED_LABEL' — matrix does not faithfully represent source bigWig signal."
            err "See validation details in log: $LOG_FILE"
            continue
        fi
        ok "Precision validation PASS — scaled matrix confirmed against source bigWigs. See log for error bounds."

        if [[ "$T_DISPLAY_MODE" == "composite" ]]; then
            # One independent output set per run (auto mode: one run per
            # group, e.g. WT and DKO separately -- exactly as if this
            # contrast were run twice, once per group; custom mode: one
            # run containing everything the user combined). Never
            # overlaid together regardless of how many runs there are.
            local _ri
            for _ri in "${!T_RUN_LABELS[@]}"; do
                local -a _run_idx_set=(${T_RUN_GROUP_INDEX_SETS[$_ri]})
                local -a _run_labels=()
                local _gidx
                for _gidx in "${_run_idx_set[@]}"; do
                    _run_labels+=("${T_COMPOSITE_GROUP_LABELS[$_gidx]}")
                done
                local RUN_SAFE_LABEL="${T_RUN_LABELS[$_ri]//[^A-Za-z0-9_]/_}"

                # Titles name peak-set source and displayed signal separately
                # (e.g. "DKO_vs_WT significant peaks -- WT composite (n=3)")
                # so a single-group composite can never be misread as a
                # missing/failed two-group comparison -- it's a deliberate,
                # complete figure on its own.
                local _group_join=""
                local _gi4
                for _gi4 in "${!_run_labels[@]}"; do
                    if [[ "$_gi4" -eq 0 ]]; then
                        _group_join="${_run_labels[$_gi4]}"
                    else
                        _group_join="${_group_join} + ${_run_labels[$_gi4]}"
                    fi
                done
                local COMPOSITE_TITLE="${RESOLVED_CONTRAST} ${BED_LABEL//_/ } — ${_group_join} composite${TORNADO_UNITS_NOTE}"

                local COMPOSITE_MATRIX_OUT="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_matrix.gz"
                local COMPOSITE_HEATMAP_PNG="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_tornado_${TS}.png"
                local COMPOSITE_HEATMAP_PDF="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_tornado_${TS}.pdf"
                local COMPOSITE_SORTED_BED="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_sorted_regions_${TS}.bed"
                local COMPOSITE_PROFILE_TSV="$TORNADO_COMPUTE_DIR/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_profile_${TS}.tsv"
                local COMPOSITE_PROFILE_PNG="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_profile_${TS}.png"
                local COMPOSITE_PROFILE_PDF="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_${RUN_SAFE_LABEL}_composite_profile_${TS}.pdf"

                label "Averaging replicates into composite group(s): ${_group_join}..."
                set +e
                conda run --no-capture-output -n "$ENV_NAME" python3 "$COMPOSITE_PY" \
                    --matrix "$MATRIX_OUT" \
                    --groups-file "${T_RUN_GROUPS_TSV[$_ri]}" \
                    --composite-matrix-out "$COMPOSITE_MATRIX_OUT" \
                    --profile-tsv-out "$COMPOSITE_PROFILE_TSV" \
                    --profile-plot-png "$COMPOSITE_PROFILE_PNG" \
                    --profile-plot-pdf "$COMPOSITE_PROFILE_PDF" \
                    --title "$COMPOSITE_TITLE" \
                    --ref-label "$T_CENTER_LABEL" \
                    --units-note "$TORNADO_UNITS_NOTE" \
                    >> "$LOG_FILE" 2>&1
                local COMPOSITE_EXIT=$?
                set -e

                if [[ "$COMPOSITE_EXIT" -ne 0 ]]; then
                    err "Composite averaging failed for '$BED_LABEL' / '${T_RUN_LABELS[$_ri]}'. Log: $LOG_FILE"
                    continue
                fi
                ok "Composite matrix written: $(basename "$COMPOSITE_MATRIX_OUT")"
                ok "Composite profile data: $(basename "$COMPOSITE_PROFILE_TSV")"
                ok "Composite trace (mean +/- SD) PNG: $(basename "$COMPOSITE_PROFILE_PNG")"
                [[ -s "$COMPOSITE_PROFILE_PDF" ]] && ok "Composite trace (mean +/- SD) PDF: $(basename "$COMPOSITE_PROFILE_PDF")"

                local -a PLOT_COMPOSITE=(
                    --matrixFile "$COMPOSITE_MATRIX_OUT"
                    --sortRegions "$T_SORT_REGIONS"
                    --samplesLabel "${_run_labels[@]}"
                    --plotTitle "$COMPOSITE_TITLE"
                    --xAxisLabel "Distance from ${T_CENTER_LABEL}"
                    --refPointLabel "$T_CENTER_LABEL"
                    --heatmapHeight 15
                    --heatmapWidth 3
                    --colorList "$TORNADO_HEATMAP_COLORLIST"
                )
                if [[ -n "$T_SORT_USING" ]]; then
                    PLOT_COMPOSITE+=(--sortUsing "$T_SORT_USING")
                    if [[ -n "$T_SORT_REF_GROUP_INDEX" ]]; then
                        PLOT_COMPOSITE+=(--sortUsingSamples "$T_SORT_REF_GROUP_INDEX")
                    fi
                fi

                label "Plotting composite tornado heatmap: ${_group_join}..."
                set +e
                conda run --no-capture-output -n "$ENV_NAME" \
                    plotHeatmap \
                        "${PLOT_COMPOSITE[@]}" \
                        --outFileName "$COMPOSITE_HEATMAP_PNG" \
                        --outFileSortedRegions "$COMPOSITE_SORTED_BED" \
                        --dpi 200 \
                    >> "$LOG_FILE" 2>&1
                local COMPOSITE_HEATMAP_EXIT=$?
                set -e

                if [[ "$COMPOSITE_HEATMAP_EXIT" -eq 0 ]]; then
                    conda run --no-capture-output -n "$ENV_NAME" \
                        plotHeatmap \
                            "${PLOT_COMPOSITE[@]}" \
                            --outFileName "$COMPOSITE_HEATMAP_PDF" \
                        >> "$LOG_FILE" 2>&1 || true

                    ok "Composite tornado PNG: $(basename "$COMPOSITE_HEATMAP_PNG")"
                    [[ -s "$COMPOSITE_HEATMAP_PDF" ]] && ok "Composite tornado PDF: $(basename "$COMPOSITE_HEATMAP_PDF")"
                    ok "Sorted regions BED: $(basename "$COMPOSITE_SORTED_BED")"
                else
                    err "plotHeatmap failed for composite matrix '${T_RUN_LABELS[$_ri]}'. Log: $LOG_FILE"
                fi
            done

            continue
        fi

        # plotHeatmap separates sort direction (--sortRegions) from the value
        # used for sorting (--sortUsing).  Passing "mean" to --sortRegions is
        # invalid, so build the argument list explicitly.
        local -a PLOT_COMMON=(
            --matrixFile "$MATRIX_OUT"
            --sortRegions "$T_SORT_REGIONS"
            --samplesLabel "${SCORE_SHORT_LABELS[@]}"
            --plotTitle "${RESOLVED_CONTRAST} — ${BED_LABEL//_/ }${TORNADO_UNITS_NOTE}"
            --xAxisLabel "Distance from ${T_CENTER_LABEL}"
            --refPointLabel "$T_CENTER_LABEL"
            --heatmapHeight 15
            --heatmapWidth 3
            --colorList "$TORNADO_HEATMAP_COLORLIST"
        )
        if [[ -n "$T_SORT_USING" ]]; then
            PLOT_COMMON+=(--sortUsing "$T_SORT_USING")
            # Restrict which samples' values determine the shared row order
            # (chosen interactively above) -- indices are 1-based and must
            # match the order samples were actually passed to computeMatrix
            # (SCORE_IDS), not the original bundle-wide sample list, since
            # some samples may have been dropped during BAM/bigWig
            # validation above. Falls back to sorting by all samples
            # (deepTools' own default when --sortUsingSamples is omitted)
            # if none of the chosen reference sample(s) survived validation.
            if [[ ${#T_SORT_REF_SAMPLE_NAMES[@]} -gt 0 ]]; then
                local -a SORT_REF_INDICES=()
                local _ref _j
                for _ref in "${T_SORT_REF_SAMPLE_NAMES[@]}"; do
                    for _j in "${!SCORE_IDS[@]}"; do
                        if [[ "${SCORE_IDS[$_j]}" == "$_ref" ]]; then
                            SORT_REF_INDICES+=("$((_j + 1))")
                            break
                        fi
                    done
                done
                if [[ ${#SORT_REF_INDICES[@]} -gt 0 ]]; then
                    PLOT_COMMON+=(--sortUsingSamples "${SORT_REF_INDICES[@]}")
                else
                    warn "None of the chosen reference sample(s) survived BAM/bigWig validation for '$BED_LABEL' -- sorting by all samples instead."
                fi
            fi
        fi

        label "Plotting tornado heatmap..."
        set +e
        conda run --no-capture-output -n "$ENV_NAME" \
            plotHeatmap \
                "${PLOT_COMMON[@]}" \
                --outFileName "$HEATMAP_PNG" \
                --outFileSortedRegions "$SORTED_BED" \
                --dpi 200 \
            >> "$LOG_FILE" 2>&1
        local HEATMAP_EXIT=$?
        set -e

        if [[ "$HEATMAP_EXIT" -eq 0 ]]; then
            conda run --no-capture-output -n "$ENV_NAME" \
                plotHeatmap \
                    "${PLOT_COMMON[@]}" \
                    --outFileName "$HEATMAP_PDF" \
                >> "$LOG_FILE" 2>&1 || true

            ok "Tornado PNG: $(basename "$HEATMAP_PNG")"
            [[ -s "$HEATMAP_PDF" ]] && ok "Tornado PDF: $(basename "$HEATMAP_PDF")"
            ok "Sorted regions BED: $(basename "$SORTED_BED")"

            # ── Standalone profile (line-graph) plot, larger than the thin
            # strip plotHeatmap draws above each heatmap column. --perGroup
            # is used deliberately: this pipeline always hands computeMatrix
            # exactly one region file per run (one BED per all/up/down
            # subset), so deepTools' *default* profile layout -- one panel
            # per sample, one line per region-group -- would degenerate into
            # single-line panels here (verified against deepTools' own
            # parserCommon.py: with only one region group, that layout has
            # nothing to compare per panel). --perGroup instead draws one
            # panel for that single region group with every sample as its
            # own line -- the actual useful comparison. Also confirmed via
            # deepTools' source that plotProfile does NOT accept
            # --xAxisLabel (heatmap-only); --refPointLabel already labels
            # the x=0 tick, so it is omitted here rather than guessed at.
            local PROFILE_PNG="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_profile_${TS}.png"
            local PROFILE_PDF="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_profile_${TS}.pdf"
            local -a PROFILE_COMMON=(
                --matrixFile "$MATRIX_OUT"
                --perGroup
                --samplesLabel "${SCORE_SHORT_LABELS[@]}"
                --plotTitle "${RESOLVED_CONTRAST} — ${BED_LABEL//_/ }${TORNADO_UNITS_NOTE}"
                --refPointLabel "$T_CENTER_LABEL"
                --plotHeight 12
                --plotWidth 20
            )

            label "Plotting profile line graph..."
            set +e
            conda run --no-capture-output -n "$ENV_NAME" \
                plotProfile \
                    "${PROFILE_COMMON[@]}" \
                    --outFileName "$PROFILE_PNG" \
                    --dpi 200 \
                >> "$LOG_FILE" 2>&1
            local PROFILE_EXIT=$?
            set -e

            if [[ "$PROFILE_EXIT" -eq 0 ]]; then
                conda run --no-capture-output -n "$ENV_NAME" \
                    plotProfile \
                        "${PROFILE_COMMON[@]}" \
                        --outFileName "$PROFILE_PDF" \
                    >> "$LOG_FILE" 2>&1 || true

                ok "Profile line graph PNG: $(basename "$PROFILE_PNG")"
                [[ -s "$PROFILE_PDF" ]] && ok "Profile line graph PDF: $(basename "$PROFILE_PDF")"
            else
                warn "plotProfile failed (tornado heatmap above is unaffected). Log: $LOG_FILE"
            fi
        else
            err "plotHeatmap failed. Log: $LOG_FILE"
        fi

    done

    unset T_BEDS T_BED_LABELS T_BAMS T_SAMPLE_IDS T_GROUPS
    unset VALID_BAMS VALID_IDS VALID_BIGWIGS SCORE_IDS
    unset VALID_GROUPS SCORE_GROUPS SCORE_SHORT_LABELS
    unset T_SORT_REF_SAMPLE_NAMES GROUP_REP_COUNTER
    unset T_COMPOSITE_GROUP_NAMES T_COMPOSITE_GROUP_MEMBERS
    unset T_COMPOSITE_GROUP_INDICES T_COMPOSITE_GROUP_LABELS
}

run_tornado_individual() { run_tornado "individual"; }
run_tornado_composite()  { run_tornado "composite"; }

# ─────────────────────────────────────────────────────────────
# Main menu loop
# ─────────────────────────────────────────────────────────────

header "Step 4 · Explorer Menu"
echo -e "  ${DIM}Explorer output root: $EXPLORE_OUT${RESET}"
echo -e "  ${DIM}Each analysis run is saved in its own numbered, timestamped subfolder.${RESET}"
echo -e "  ${DIM}Run as many analyses, in any order, as you like.${RESET}"

while true; do
    blank
    echo -e "  ${BOLD}Choose an analysis:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  PCA plot"
    if [[ "$HAS_VST" != "TRUE" ]]; then
        echo -e "        ${DIM}(unavailable — no VST matrix in bundle)${RESET}"
    fi
    echo -e "    ${CYAN}2${RESET}.  Re-plot contrast — volcano + MA at custom thresholds"
    echo -e "    ${CYAN}3${RESET}.  Peak accessibility boxplot — raw + normalized side by side"
    echo -e "    ${CYAN}4${RESET}.  Sample correlation heatmap"
    if [[ "$HAS_VST" != "TRUE" ]]; then
        echo -e "        ${DIM}(unavailable — no VST matrix in bundle)${RESET}"
    fi
    echo -e "    ${CYAN}5${RESET}.  Motif summary — top known TF motifs for a contrast"
    if [[ "$HAS_HOMER" != "TRUE" ]]; then
        echo -e "        ${DIM}(unavailable — HOMER path not found)${RESET}"
    fi
    echo -e "    ${CYAN}6${RESET}.  Annotation distribution plots — feature composition bar, TSS/TES, gene body, width, chromosome"
    echo -e "    ${CYAN}7${RESET}.  GO/KEGG summary — re-draw enrichment dot plot"
    if [[ "$HAS_GO" != "TRUE" && "$HAS_KEGG" != "TRUE" ]]; then
        echo -e "        ${DIM}(unavailable — no saved GO or KEGG results)${RESET}"
    elif [[ "$HAS_KEGG" != "TRUE" ]]; then
        echo -e "        ${DIM}(GO available; KEGG unavailable)${RESET}"
    fi
    echo -e "    ${CYAN}8${RESET}.  Tornado plot — deepTools heatmap around differential peaks"
    echo -e "        ${DIM}(center: peak center / start / end via deepTools TSS/TES mode · window · sort · grouping · Up/Down split)${RESET}"
    echo -e "    ${CYAN}9${RESET}.  Composite tornado plot — mean per replicate group, +/- SD shading on the trace"
    echo -e "        ${DIM}(same options as #8, but you define replicate groups and each panel/line is a group mean)${RESET}"
    echo -e "    ${CYAN}q${RESET}.  Quit"
    blank

    read -p "  Choice [1/2/3/4/5/6/7/8/9/q]: " MENU_CHOICE

    case "$MENU_CHOICE" in
        1) run_pca_analysis ;;
        2) run_replot_analysis ;;
        3) run_peak_boxplot_analysis ;;
        4) run_correlation_heatmap ;;
        5) run_motif_summary ;;
        6) run_annotation_plots ;;
        7) run_go_kegg_replot ;;
        8) run_tornado_individual ;;
        9) run_tornado_composite ;;
        q|Q) break ;;
        *) err "Enter 1, 2, 3, 4, 5, 6, 7, 8, 9, or q." ;;
    esac
done

blank
header "Explorer Session Complete"
echo -e "  ${BOLD}Organized analysis folders saved under:${RESET}  $EXPLORE_OUT"
blank
echo -e "  ${GREEN}${BOLD}Done.${RESET}"
blank
