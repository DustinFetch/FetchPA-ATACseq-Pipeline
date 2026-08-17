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
ENV_NAME="pepatac"

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

    local TAG="replot_${RUN_ID}_$$"
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local SAFE_LABEL="${RESOLVED_CONTRAST//[^A-Za-z0-9_]/_}"
    local ANALYSIS_OUT_DIR="$EXPLORE_OUT/02_Contrast_Replots/${SAFE_LABEL}_${TS}"
    mkdir -p "$ANALYSIS_OUT_DIR"
    local R_SCRIPT="$ANALYSIS_OUT_DIR/.run_${TAG}.R"
    local LOG_FILE="$ANALYSIS_OUT_DIR/${TAG}.log"
    local OUT_VOL_PDF="$ANALYSIS_OUT_DIR/${SAFE_LABEL}_volcano_fdr${NEW_FDR}_lfc${NEW_FC}_${TS}.pdf"
    local OUT_VOL_PNG="${OUT_VOL_PDF%.pdf}.png"
    local OUT_MA_PDF="${OUT_VOL_PDF/_volcano_/_MA_}"
    local OUT_MA_PNG="${OUT_MA_PDF%.pdf}.png"
    local OUT_CSV="${OUT_VOL_PDF%_volcano_*}_replot_${TS}.csv"

    local RESOLVED_CONTRAST_R OUT_CSV_R OUT_VOL_PDF_R OUT_VOL_PNG_R OUT_MA_PDF_R OUT_MA_PNG_R
    RESOLVED_CONTRAST_R="$(r_string_literal "$RESOLVED_CONTRAST")"
    OUT_CSV_R="$(r_string_literal "$OUT_CSV")"
    OUT_VOL_PDF_R="$(r_string_literal "$OUT_VOL_PDF")"
    OUT_VOL_PNG_R="$(r_string_literal "$OUT_VOL_PNG")"
    OUT_MA_PDF_R="$(r_string_literal "$OUT_MA_PDF")"
    OUT_MA_PNG_R="$(r_string_literal "$OUT_MA_PNG")"

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

res <- b\$all_results[[label]]
if (is.null(res)) stop("Contrast '", label, "' not found in bundle.")

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
label_peaks <- res_plot[res_plot\$Sig, ]
label_peaks\$peak_id <- paste0(label_peaks\$Chr, ":", label_peaks\$Start, "-", label_peaks\$End)

subtitle_vol <- if (fc_cutoff > 0) {
    sprintf("%d up, %d down  (FDR < %.4g  |  |log2FC| >= %.3g)",
            n_up, n_down, fdr_cutoff, fc_cutoff)
} else {
    sprintf("%d up, %d down  (FDR < %.4g)", n_up, n_down, fdr_cutoff)
}

p_vol <- ggplot(res_plot, aes(x=log2FoldChange, y=neg_log10_padj, color=Direction)) +
    geom_point(alpha=0.5, size=1.2) +
    geom_text_repel(data=label_peaks, aes(label=peak_id),
                    size=2.2, max.overlaps=15, show.legend=FALSE) +
    geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed",
               color="grey50", linewidth=0.5) +
    { if (fc_cutoff > 0)
        geom_vline(xintercept=c(-fc_cutoff, fc_cutoff), linetype="dashed",
                   color="grey50", linewidth=0.5)
      else list() } +
    scale_color_manual(values=c(Up="firebrick3", Down="steelblue3", NS="grey70")) +
    labs(title=paste0("Volcano: ", label),
         subtitle=subtitle_vol,
         x="log2 Fold Change", y="-log10(adjusted p-value)") +
    theme_bw(base_size=13)

ggsave(${OUT_VOL_PDF_R}, p_vol, width=8, height=6)
ggsave(${OUT_VOL_PNG_R}, p_vol, width=8, height=6, dpi=150)
cat("Volcano saved:", ${OUT_VOL_PDF_R}, "\n")

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
             y="log2 Fold Change") +
        theme_bw(base_size=13)
    ggsave(${OUT_MA_PDF_R}, p_ma, width=8, height=6)
    ggsave(${OUT_MA_PNG_R}, p_ma, width=8, height=6, dpi=150)
    cat("MA plot saved:", ${OUT_MA_PDF_R}, "\n")
} else {
    cat("  [NOTE] MA plot skipped: no 'Conc' column in results.\n")
}
RSCRIPT_EOF

    label "Running re-plot..."
    if run_r_script "$R_SCRIPT" "$LOG_FILE"; then
        ok "Volcano PDF: $OUT_VOL_PDF"
        ok "Volcano PNG: $OUT_VOL_PNG"
        [[ -f "$OUT_MA_PDF" ]] && ok "MA PDF: $OUT_MA_PDF"
        [[ -f "$OUT_MA_PNG" ]] && ok "MA PNG: $OUT_MA_PNG"
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
    echo -e "    ${CYAN}1${RESET}.  Genomic feature composition (pie)"
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
        p_feat <- ggplot(feature_df, aes(x="", y=Count, fill=Feature)) +
            geom_col(width=1, color="white", linewidth=0.25) +
            coord_polar(theta="y") +
            geom_text(aes(label=sprintf("%s\n%.1f%%", Feature, Percent)),
                      position=position_stack(vjust=0.5), size=3) +
            labs(title=paste0(label, " - Genomic Feature Annotation"),
                 subtitle=paste0("ChIPseeker / ", txdb_pkg), fill="Feature") +
            theme_void(base_size=11) +
            theme(plot.title=element_text(face="bold"))
        save_both(p_feat, file.path(out_dir, "genomic_features_pie"), 9, 7)
        cat("  Saved: genomic_features_pie (pdf/png) + genomic_features_counts.csv\n")
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
    header "Tornado Plot"

    # ── deepTools check ──────────────────────────────────────────
    # computeMatrix consumes the PEPATAC smoothShift bigWig tracks directly.
    # bamCoverage is not needed — PEPATAC already produced normalised bigWigs.
    if ! conda run --no-capture-output -n "$ENV_NAME" bash -c \
            "command -v computeMatrix && command -v plotHeatmap" >/dev/null 2>&1; then
        err "deepTools (computeMatrix / plotHeatmap) not found in the pepatac environment."
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

    # ── Centering ────────────────────────────────────────────────
    blank
    echo -e "  ${BOLD}Center on:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  Peak centers"
    echo -e "    ${CYAN}2${RESET}.  Region starts  ${DIM}(deepTools TSS mode)${RESET}"
    echo -e "    ${CYAN}3${RESET}.  Region ends    ${DIM}(deepTools TES mode)${RESET}"
    blank
    local T_CENTER T_CENTER_LABEL
    while true; do
        read -p "  Choice [1/2/3, default 1]: " T_C
        T_C="${T_C:-1}"
        case "$T_C" in
            # deepTools accepts exactly TSS, TES, or center.  "midpoint" is
            # not a valid --referencePoint value.
            1) T_CENTER="center"; T_CENTER_LABEL="peak center";  break ;;
            2) T_CENTER="TSS";    T_CENTER_LABEL="region start"; break ;;
            3) T_CENTER="TES";    T_CENTER_LABEL="region end";   break ;;
            *) err "Enter 1, 2, or 3." ;;
        esac
    done
    ok "Centering: $T_CENTER_LABEL"

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

    # ── Sample grouping ──────────────────────────────────────────
    blank
    echo -e "  ${BOLD}Sample ordering:${RESET}"
    echo -e "    ${CYAN}1${RESET}.  Original sample order"
    echo -e "    ${CYAN}2${RESET}.  Order by condition group (case then control)"
    blank
    local T_SPLIT=false
    while true; do
        read -p "  Choice [1/2, default 2]: " T_G
        T_G="${T_G:-2}"
        case "$T_G" in
            1) T_SPLIT=false; break ;;
            2) T_SPLIT=true;  break ;;
            *) err "Enter 1 or 2." ;;
        esac
    done

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
    local TS
    TS="$(date '+%Y%m%d_%H%M%S')"
    local TAG="tornado_${RUN_ID}_$$"
    local TORNADO_OUT="$EXPLORE_OUT/08_Tornado_Plots/${SAFE_LABEL}_${TS}"
    mkdir -p "$TORNADO_OUT"
    local PREP_R="$TORNADO_OUT/.${TAG}_prep.R"
    local PREP_OUT="$TORNADO_OUT/.${TAG}_prep.tsv"
    local LOG_FILE="$TORNADO_OUT/${TAG}.log"

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
beds  <- character(0)
blabs <- character(0)
if (subset == "all") {
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
    declare -a VALID_BAMS=() VALID_IDS=()
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
    done

    if [[ ${#VALID_BAMS[@]} -eq 0 ]]; then
        err "No valid BAM files found — aborting tornado."
        return
    fi

    # computeMatrix requires bigWig score files.  Use the smoothShift bigWigs
    # produced by PEPATAC during sample processing — these are ATAC-seq optimised
    # (Tn5 shift-corrected, read-count normalised) and live in the same aligned_*
    # folder as each BAM, named <SampleID>_smooth_shift.bw.
    declare -a VALID_BIGWIGS=() SCORE_IDS=()

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
    done

    if [[ ${#VALID_BIGWIGS[@]} -eq 0 ]]; then
        err "No smoothShift bigWig tracks found — aborting tornado."
        return
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
                ORDERED_BED="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_ordered_by_peak_score_${TS}.bed"
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
                ORDERED_BED="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_ordered_genomically_${TS}.bed"
                sort -k1,1V -k2,2n -k3,3n "$BED_FILE" > "$ORDERED_BED"
                BED_FOR_MATRIX="$ORDERED_BED"
                ;;
            keep|mean)
                ;;
        esac

        local N_PEAKS
        N_PEAKS=$(wc -l < "$BED_FOR_MATRIX")
        local MATRIX_OUT="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_matrix.gz"
        local HEATMAP_PNG="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_tornado_${TS}.png"
        local HEATMAP_PDF="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_tornado_${TS}.pdf"
        local SORTED_BED="$TORNADO_OUT/${SAFE_LABEL}_${BED_LABEL}_sorted_regions_${TS}.bed"

        blank
        label "Computing matrix: $BED_LABEL ($N_PEAKS peaks, ${#VALID_BIGWIGS[@]} bigWigs)..."
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
            >> "$LOG_FILE" 2>&1
        local MATRIX_EXIT=$?
        set -e

        if [[ "$MATRIX_EXIT" -ne 0 ]]; then
            err "computeMatrix failed for '$BED_LABEL'. Log: $LOG_FILE"
            continue
        fi
        ok "Matrix written: $(basename "$MATRIX_OUT")"

        # plotHeatmap separates sort direction (--sortRegions) from the value
        # used for sorting (--sortUsing).  Passing "mean" to --sortRegions is
        # invalid, so build the argument list explicitly.
        local -a PLOT_COMMON=(
            --matrixFile "$MATRIX_OUT"
            --sortRegions "$T_SORT_REGIONS"
            --samplesLabel "${SCORE_IDS[@]}"
            --plotTitle "${RESOLVED_CONTRAST} — ${BED_LABEL//_/ }"
            --xAxisLabel "Distance from ${T_CENTER_LABEL}"
            --refPointLabel "$T_CENTER_LABEL"
            --heatmapHeight 15
            --heatmapWidth 3
        )
        if [[ -n "$T_SORT_USING" ]]; then
            PLOT_COMMON+=(--sortUsing "$T_SORT_USING")
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
        else
            err "plotHeatmap failed. Log: $LOG_FILE"
        fi

    done

    unset T_BEDS T_BED_LABELS T_BAMS T_SAMPLE_IDS T_GROUPS
    unset VALID_BAMS VALID_IDS VALID_BIGWIGS SCORE_IDS
}

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
    echo -e "    ${CYAN}6${RESET}.  Annotation distribution plots — feature pie, TSS/TES, gene body, width, chromosome"
    echo -e "    ${CYAN}7${RESET}.  GO/KEGG summary — re-draw enrichment dot plot"
    if [[ "$HAS_GO" != "TRUE" && "$HAS_KEGG" != "TRUE" ]]; then
        echo -e "        ${DIM}(unavailable — no saved GO or KEGG results)${RESET}"
    elif [[ "$HAS_KEGG" != "TRUE" ]]; then
        echo -e "        ${DIM}(GO available; KEGG unavailable)${RESET}"
    fi
    echo -e "    ${CYAN}8${RESET}.  Tornado plot — deepTools heatmap around differential peaks"
    echo -e "        ${DIM}(center: peak center / start / end via deepTools TSS/TES mode · window · sort · grouping · Up/Down split)${RESET}"
    echo -e "    ${CYAN}q${RESET}.  Quit"
    blank

    read -p "  Choice [1/2/3/4/5/6/7/8/q]: " MENU_CHOICE

    case "$MENU_CHOICE" in
        1) run_pca_analysis ;;
        2) run_replot_analysis ;;
        3) run_peak_boxplot_analysis ;;
        4) run_correlation_heatmap ;;
        5) run_motif_summary ;;
        6) run_annotation_plots ;;
        7) run_go_kegg_replot ;;
        8) run_tornado ;;
        q|Q) break ;;
        *) err "Enter 1, 2, 3, 4, 5, 6, 7, 8, or q." ;;
    esac
done

blank
header "Explorer Session Complete"
echo -e "  ${BOLD}Organized analysis folders saved under:${RESET}  $EXPLORE_OUT"
blank
echo -e "  ${GREEN}${BOLD}Done.${RESET}"
blank
