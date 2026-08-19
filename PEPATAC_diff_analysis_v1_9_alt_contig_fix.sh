#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# PEPATAC DIFFERENTIAL ANALYSIS
# Interactive differential accessibility analysis for
# samples processed by run.sh (PEPATAC local runner).
#
# Flow:
#   1. Locate a PEPATAC run output folder.
#   2. Discover samples from the autodetected R sample sheet.
#   2b. Explicitly include/exclude samples before DiffBind sees them.
#   3. Assign group labels to included samples interactively.
#   4. Define contrasts (group A vs group B).
#   5. Select output folder.
#   6. Show a final run summary and confirm.
#   7. Check required R packages (stops if any are missing; set
#      PEPATAC_DIFF_AUTOINSTALL=TRUE to install them instead).
#   8. Build DiffBind only from included samples, then run DESeq2.
#   9. Per-contrast: diff peaks CSV, ChIPseeker annotation +
#      feature/TSS/TES stats, HOMER motif enrichment, GO/KEGG
#      enrichment. This script computes only — no plots are
#      generated here; use PEPATAC_explore.sh for all graphs
#      (PCA, volcano/MA, annotation distributions, GO/KEGG
#      dot plots, motif summaries, tornado plots).
#  10. Write a final cross-contrast summary.
#  11. Write explorer_bundle.rds for use with PEPATAC_explore.sh.
# ============================================================

DIFF_VERSION="1.9-alt-contig-fix-mle-map-fix"
SCRIPT_VERSION="$DIFF_VERSION"
REFERENCE_SNAPSHOT_SCHEMA_SUPPORTED=4
RUN_STARTED="$(date '+%Y-%m-%d %H:%M:%S')"
RUN_ID="$(date '+%Y%m%d_%H%M%S')"

ENV_NAME="pepatac"
PEPATAC_DIR="$HOME/pepatac"
REFGENIE_CONFIG="$HOME/refgenie/refgenie.yaml"

# ─────────────────────────────────────────────────────────────
# Terminal colors / formatting  (identical to run.sh)
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

command_exists() { command -v "$1" >/dev/null 2>&1; }

# ─────────────────────────────────────────────────────────────
# Path normalisation  (shared with run.sh / explore.sh / install.sh)
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

# ─────────────────────────────────────────────────────────────
# CSV parsing + safe R-literal embedding
# ─────────────────────────────────────────────────────────────
#
# This script generates an R script by interpolating Bash values (sample
# IDs, BAM/peak paths, group and contrast labels) directly into R source
# text. Two related risks live here:
#   1. Reading the autodetected sample sheet with `awk -F'","'` breaks the
#      moment any field contains a literal comma or an escaped quote --
#      it's a naive split on the 3-character separator, not real CSV.
#   2. Embedding a Bash value into an R string literal via `printf
#      '"%s",'` with no escaping means a value containing a `"` or `\`
#      can break out of the R string and get parsed as R syntax --
#      an actual code-injection boundary, not just a cosmetic parsing bug.
# Both are fixed below rather than patched around.

# csv_fields LINE
# RFC4180-aware CSV row splitter: correctly handles fields with embedded
# commas and doubled-quote escapes ("" -> "), unlike a naive split on the
# literal '","' separator. Emits one field per line (quotes already
# stripped) so the caller can read them with `mapfile`.
csv_fields() {
    awk '
    {
        line = $0
        n = length(line)
        field = ""
        inq = 0
        for (i = 1; i <= n; i++) {
            c = substr(line, i, 1)
            if (inq) {
                if (c == "\"") {
                    if (substr(line, i + 1, 1) == "\"") { field = field "\""; i++ }
                    else { inq = 0 }
                } else field = field c
            } else {
                if (c == "\"") inq = 1
                else if (c == ",") { print field; field = "" }
                else field = field c
            }
        }
        print field
    }' <<< "$1"
}

# url_encode STRING
# Percent-encodes STRING to [A-Za-z0-9._~-] plus %XX escapes. Used to make
# any value that will be embedded into the generated R script's source
# text safe by construction: the encoded output can never contain a `"`,
# `\`, or any other character with meaning in R syntax, so there is
# nothing for a value to inject with even if it contains one. The R side
# decodes with utils::URLdecode() (base R, no extra package) before use --
# see .pepatac_url_decode() in the generated R script.
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

# csv_quote VALUE
# RFC4180 field escaping: doubles any embedded double-quote so VALUE can
# be safely wrapped in "..." without corrupting the CSV row. Used when
# writing the DiffBind sample sheet, which DiffBind reads back with its
# own CSV parser -- a different injection surface than the R-literal
# embedding fixed elsewhere in this script (r_vector_literal /
# r_string_literal): this file is only ever read back as data, so an
# unescaped quote here shifts columns or breaks a row rather than
# executing anything, but it's still real corruption on a value
# containing a literal `"` (e.g. a path with a quote in it).
csv_quote() {
    printf '%s' "${1//\"/\"\"}"
}

in_env() {
    conda run --no-capture-output -n "$ENV_NAME" "$@"
}

in_env_clean() {
    conda run --no-capture-output -n "$ENV_NAME" \
        env -u R_ARCH -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$@"
}

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
  Differential Accessibility Analysis
BANNER
echo -e "${RESET}"
echo -e "  ${DIM}PEPATAC Differential Runner  ·  DiffBind + DESeq2${RESET}"
echo -e "  ${DIM}Select samples → assign groups → define contrasts → take off${RESET}"
echo ""

# ─────────────────────────────────────────────────────────────
# Source conda
# ─────────────────────────────────────────────────────────────

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [[ -f "$CONDA_SH" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_SH"
else
    die "Conda not found at $CONDA_SH. Run install.sh first."
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    die "Conda environment '$ENV_NAME' not found. Run install.sh first."
fi

# ─────────────────────────────────────────────────────────────
# STEP 1 — Locate PEPATAC run output folder
# ─────────────────────────────────────────────────────────────

header "Step 1 · PEPATAC Run Folder"

echo -e "  Enter the path to a completed PEPATAC run output folder."
echo -e "  ${DIM}This folder should contain samples_for_R_autodetected.csv${RESET}"
echo -e "  ${DIM}which is written automatically by run.sh.${RESET}"
blank

while true; do
    IFS= read -r -e -p "  Run output folder: " RUN_DIR
    normalize_path_var RUN_DIR

    if [[ -z "$RUN_DIR" ]]; then
        err "Please enter a path."
        continue
    fi

    if [[ ! -d "$RUN_DIR" ]]; then
        err "Directory not found: $RUN_DIR"
        continue
    fi

    AUTODETECTED_SHEET="$RUN_DIR/samples_for_R_autodetected.csv"
    if [[ ! -f "$AUTODETECTED_SHEET" ]]; then
        warn "samples_for_R_autodetected.csv not found in $RUN_DIR"
        echo -e "  ${DIM}Expected: $AUTODETECTED_SHEET${RESET}"
        read -p "  Use this folder anyway (you will enter sample paths manually)? [y/N]: " FORCE
        if [[ "${FORCE,,}" == "y" ]]; then
            AUTODETECTED_SHEET=""
            break
        fi
        continue
    fi

    ok "Run folder found: $RUN_DIR"
    ok "Autodetected sample sheet: $AUTODETECTED_SHEET"
    break
done

# ─────────────────────────────────────────────────────────────
# STEP 2 — Discover samples
# ─────────────────────────────────────────────────────────────

header "Step 2 · Sample Discovery"

declare -a SAMPLE_IDS
declare -a SAMPLE_BAMS
declare -a SAMPLE_PEAKS
declare -a SAMPLE_STATUS_ARR
declare -a SAMPLE_PEAK_CALLERS

# infer_peak_caller PEAK_FILE
# Mirrors the runner's own peak_caller_for() logic. Used as a fallback when
# a sample's actual peak format wasn't available from an upstream CSV --
# either it was entered manually (no CSV at all) or it came from an older
# sheet written before the runner tracked this column. Never guess "narrow"
# just because that's the common case; a wrong guess here is exactly the
# bug this exists to prevent.
infer_peak_caller() {
    local peak_file="$1"
    case "$(basename "${peak_file:-}")" in
        *.narrowPeak) echo "narrow" ;;
        "") echo "narrow" ;;   # nothing to check -- keep the prior default
        *)  echo "bed" ;;
    esac
}

# Samples removed before group assignment are tracked separately so the
# run summary and output records can prove exactly what DiffBind did not see.
declare -a EXCLUDED_SAMPLE_IDS=()
declare -a EXCLUDED_SAMPLE_REASONS=()

if [[ -n "$AUTODETECTED_SHEET" ]]; then
    # An older sheet (written before the runner tracked peak format) won't
    # have this column at all. Check the header rather than assuming column
    # 6 is always PeakCaller -- reading the wrong column silently would be
    # exactly the kind of bug this is meant to fix, not avoid.
    HEADER_LINE=$(head -n 1 "$AUTODETECTED_SHEET")
    HAS_PEAKCALLER_COL=false
    [[ "$HEADER_LINE" == *"PeakCaller"* ]] && HAS_PEAKCALLER_COL=true

    # Parse the CSV using a dedicated file descriptor (fd 3) so that the
    # read loop never competes with terminal stdin for the awk subshells.
    #
    # csv_fields() does a real RFC4180 split (respects quoted commas and
    # doubled-quote escapes) instead of naively splitting on the literal
    # '","' separator, which broke the moment any field contained a comma.
    #
    # Status's column index also now depends on HAS_PEAKCALLER_COL, matching
    # caller's existing header-based guard. The previous code read Status
    # from a hardcoded column 7 unconditionally -- correct when PeakCaller
    # is present (7 columns), but wrong for an older 6-column sheet without
    # it, where Status is actually column 6. That mismatch was silent: a
    # legacy sheet would have gotten an empty Status for every sample.
    while IFS= read -r line <&3; do
        [[ "$line" =~ ^[[:space:]]*$ ]] && continue
        [[ "$line" =~ ^\"?SampleID ]] && continue

        mapfile -t _csv_f < <(csv_fields "$line")
        sid="${_csv_f[0]:-}"
        bam="${_csv_f[3]:-}"
        peaks="${_csv_f[4]:-}"
        if $HAS_PEAKCALLER_COL; then
            caller="${_csv_f[5]:-}"
            status="${_csv_f[6]:-}"
        else
            caller=""
            status="${_csv_f[5]:-}"
        fi

        [[ -z "$sid" ]] && continue

        SAMPLE_IDS+=("$sid")
        SAMPLE_BAMS+=("$bam")
        SAMPLE_PEAKS+=("$peaks")
        SAMPLE_STATUS_ARR+=("$status")
        SAMPLE_PEAK_CALLERS+=("$caller")
    done 3< "$AUTODETECTED_SHEET"

    if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
        die "No samples found in $AUTODETECTED_SHEET"
    fi

    DISCOVERED_SAMPLE_COUNT=${#SAMPLE_IDS[@]}

    # Duplicate SampleID detection happens once, after Step 2b, once every
    # legitimate way to resolve a duplicate (status-based filtering, and
    # manual exclusion) has had a chance to run -- see that check for why.

    echo -e "  Found ${BOLD}${#SAMPLE_IDS[@]}${RESET} sample(s) from the autodetected sheet:\n"

    SKIPPED_COUNT=0
    for i in "${!SAMPLE_IDS[@]}"; do
        sid="${SAMPLE_IDS[$i]}"
        bam="${SAMPLE_BAMS[$i]}"
        peak="${SAMPLE_PEAKS[$i]}"
        status="${SAMPLE_STATUS_ARR[$i]}"

        bam_ok=false
        peak_ok=false
        status_ok=false
        [[ -f "$bam" ]]  && bam_ok=true
        [[ -f "$peak" ]] && peak_ok=true
        [[ "${status^^}" == "PASS" ]] && status_ok=true

        if $status_ok && $bam_ok && $peak_ok; then
            printf "    ${GREEN}✔${RESET}  %-36s  Status ✔  BAM ✔  Peak ✔\n" "$sid"
        else
            printf "    ${YELLOW}⚠${RESET}  %-36s  Status %s [%s]  BAM %s  Peak %s\n" \
                "$sid" \
                "$( $status_ok && echo '✔' || echo '✘' )" \
                "${status:-UNKNOWN}" \
                "$( $bam_ok    && echo '✔' || echo '✘' )" \
                "$( $peak_ok   && echo '✔' || echo '✘' )"
            (( SKIPPED_COUNT++ )) || true
        fi
    done

    blank

    if [[ $SKIPPED_COUNT -gt 0 ]]; then
        warn "$SKIPPED_COUNT sample(s) failed one or more checks (Status≠PASS, missing BAM, or missing peak file)."
        echo -e "  ${DIM}Samples with Status≠PASS may have incomplete or partial outputs.${RESET}"
        echo -e "  ${DIM}Including them risks introducing low-quality data into the differential analysis.${RESET}"
        blank
        read -p "  Exclude flagged samples and continue with passing samples only? [Y/n]: " EXCL
        if [[ "${EXCL,,}" != "n" ]]; then
            VALID_IDS=()
            VALID_BAMS=()
            VALID_PEAKS=()
            VALID_STATUS=()
            VALID_CALLERS=()
            for i in "${!SAMPLE_IDS[@]}"; do
                s_status="${SAMPLE_STATUS_ARR[$i]:-}"
                if [[ "${s_status^^}" == "PASS" ]] \
                   && [[ -f "${SAMPLE_BAMS[$i]}" ]] \
                   && [[ -f "${SAMPLE_PEAKS[$i]}" ]]; then
                    VALID_IDS+=("${SAMPLE_IDS[$i]}")
                    VALID_BAMS+=("${SAMPLE_BAMS[$i]}")
                    VALID_PEAKS+=("${SAMPLE_PEAKS[$i]}")
                    VALID_STATUS+=("${s_status}")
                    VALID_CALLERS+=("${SAMPLE_PEAK_CALLERS[$i]:-}")
                else
                    bam_state="$( [[ -f "${SAMPLE_BAMS[$i]}" ]] && echo present || echo missing )"
                    peak_state="$( [[ -f "${SAMPLE_PEAKS[$i]}" ]] && echo present || echo missing )"
                    EXCLUDED_SAMPLE_IDS+=("${SAMPLE_IDS[$i]}")
                    EXCLUDED_SAMPLE_REASONS+=(
                        "automatic pre-analysis exclusion: status=${s_status:-UNKNOWN}; BAM=${bam_state}; peak=${peak_state}"
                    )
                fi
            done
            SAMPLE_IDS=("${VALID_IDS[@]}")
            SAMPLE_BAMS=("${VALID_BAMS[@]}")
            SAMPLE_PEAKS=("${VALID_PEAKS[@]}")
            SAMPLE_STATUS_ARR=("${VALID_STATUS[@]}")
            SAMPLE_PEAK_CALLERS=("${VALID_CALLERS[@]}")
            if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
                die "No valid samples remain after filtering. Check PEPATAC status and file paths."
            fi
            ok "Continuing with ${#SAMPLE_IDS[@]} valid sample(s)."
        else
            warn "Proceeding with all samples including flagged ones — interpret results cautiously."
        fi
    fi

else
    # Manual sample entry fallback.
    echo -e "  No autodetected sheet found. Enter samples manually."
    echo -e "  ${DIM}Type 'done' as the sample ID when finished.${RESET}"
    blank

    while true; do
        read -p "  Sample ID (or 'done'): " sid
        [[ "$sid" == "done" ]] && break
        [[ -z "$sid" ]] && continue

        IFS= read -r -e -p "    BAM file path:  " bam
        normalize_path_var bam
        IFS= read -r -e -p "    Peak file path: " peak
        normalize_path_var peak

        if [[ ! -f "$bam" ]];  then warn "BAM not found: $bam";  fi
        if [[ ! -f "$peak" ]]; then warn "Peak not found: $peak"; fi

        this_caller="$(infer_peak_caller "$peak")"
        if [[ "$this_caller" == "bed" ]]; then
            warn "Peak file doesn't look like narrowPeak format (no .narrowPeak extension): $peak"
            echo -e "  ${DIM}Recorded as PeakCaller=bed -- DiffBind will be told to parse this as a${RESET}"
            echo -e "  ${DIM}plain BED file, not narrowPeak, so its columns aren't misread.${RESET}"
        fi

        SAMPLE_IDS+=("$sid")
        SAMPLE_BAMS+=("$bam")
        SAMPLE_PEAKS+=("$peak")
        SAMPLE_STATUS_ARR+=("MANUAL")
        SAMPLE_PEAK_CALLERS+=("$this_caller")
    done

    [[ ${#SAMPLE_IDS[@]} -eq 0 ]] && die "No samples entered."
    DISCOVERED_SAMPLE_COUNT=${#SAMPLE_IDS[@]}
    ok "Entered ${#SAMPLE_IDS[@]} sample(s) manually."
fi

# A sample can reach here with no recorded peak-caller value if it came
# from an older sheet written before the runner tracked this column. Infer
# it the same way the manual-entry path does, rather than leaving it blank
# (which would otherwise reach the DiffBind sheet as an empty PeakCaller
# field further down).
for i in "${!SAMPLE_IDS[@]}"; do
    if [[ -z "${SAMPLE_PEAK_CALLERS[$i]:-}" ]]; then
        SAMPLE_PEAK_CALLERS[$i]="$(infer_peak_caller "${SAMPLE_PEAKS[$i]}")"
    fi
done

# ─────────────────────────────────────────────────────────────
# STEP 2b — Explicit sample inclusion/exclusion
# ─────────────────────────────────────────────────────────────

header "Step 2b · Sample Selection"

echo -e "  Choose which samples to include in ${BOLD}this run${RESET}."
echo -e "  ${DIM}DiffBind builds a single shared object from every sample it loads —${RESET}"
echo -e "  ${DIM}the consensus peak set, read counts, normalization, and PCA all reflect${RESET}"
echo -e "  ${DIM}every sample present. Any sample not part of the comparisons you are${RESET}"
echo -e "  ${DIM}running now will contaminate those shared structures.${RESET}"
echo -e "  ${DIM}This script is scoped to one set of comparisons per run. If you want${RESET}"
echo -e "  ${DIM}a different comparison, re-run and exclude accordingly.${RESET}"
echo -e "  ${DIM}Excluded samples are removed before group assignment and before the${RESET}"
echo -e "  ${DIM}DiffBind sample sheet is written — they have zero influence on consensus${RESET}"
echo -e "  ${DIM}peaks, MIN_OVERLAP, counting, normalization, PCA, differential testing,${RESET}"
echo -e "  ${DIM}GO/KEGG, motif backgrounds, or the explorer bundle.${RESET}"
blank

echo -e "  ${BOLD}Currently available samples:${RESET}"
for i in "${!SAMPLE_IDS[@]}"; do
    printf "    ${CYAN}%3d${RESET}  %s\n" "$((i+1))" "${SAMPLE_IDS[$i]}"
done
blank
echo -e "  ${DIM}Enter sample numbers separated by spaces or commas.${RESET}"
echo -e "  ${DIM}Ranges are allowed, for example: 4 8   or   4,8   or   4-6${RESET}"
echo -e "  ${DIM}Press Enter to include every sample shown above.${RESET}"
blank

while true; do
    read -p "  Samples to exclude from all downstream analysis [default: none]: " EXCLUDE_INPUT
    EXCLUDE_INPUT="${EXCLUDE_INPUT//,/ }"

    declare -A _EXCLUDE_IDX=()
    selection_error=""

    if [[ -n "${EXCLUDE_INPUT//[[:space:]]/}" ]] && [[ "${EXCLUDE_INPUT,,}" != "none" ]]; then
        for token in $EXCLUDE_INPUT; do
            if [[ "$token" =~ ^[0-9]+$ ]]; then
                n="$token"
                if (( n < 1 || n > ${#SAMPLE_IDS[@]} )); then
                    selection_error="sample number '$token' is outside 1-${#SAMPLE_IDS[@]}"
                    break
                fi
                _EXCLUDE_IDX["$((n-1))"]=1
            elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
                range_start="${BASH_REMATCH[1]}"
                range_end="${BASH_REMATCH[2]}"
                if (( range_start < 1 || range_end > ${#SAMPLE_IDS[@]} || range_start > range_end )); then
                    selection_error="invalid range '$token' (valid range is 1-${#SAMPLE_IDS[@]})"
                    break
                fi
                for (( n=range_start; n<=range_end; n++ )); do
                    _EXCLUDE_IDX["$((n-1))"]=1
                done
            else
                selection_error="could not understand '$token'"
                break
            fi
        done
    fi

    if [[ -n "$selection_error" ]]; then
        err "$selection_error"
        unset _EXCLUDE_IDX
        continue
    fi

    if (( ${#_EXCLUDE_IDX[@]} == 0 )); then
        ok "All ${#SAMPLE_IDS[@]} currently available samples will be included."
        unset _EXCLUDE_IDX
        break
    fi

    if (( ${#SAMPLE_IDS[@]} - ${#_EXCLUDE_IDX[@]} < 2 )); then
        err "At least two samples must remain after exclusion."
        unset _EXCLUDE_IDX
        continue
    fi

    blank
    echo -e "  ${BOLD}${YELLOW}Samples selected for complete exclusion:${RESET}"
    for i in "${!SAMPLE_IDS[@]}"; do
        if [[ -n "${_EXCLUDE_IDX[$i]+x}" ]]; then
            echo -e "    ${YELLOW}✘${RESET}  ${SAMPLE_IDS[$i]}"
        fi
    done
    blank
    read -p "  Confirm these samples should have zero downstream influence? [y/N]: " CONFIRM_EXCLUSION
    if [[ "${CONFIRM_EXCLUSION,,}" != "y" ]]; then
        warn "Exclusion selection discarded. Choose again."
        unset _EXCLUDE_IDX
        continue
    fi

    INCLUDED_IDS=()
    INCLUDED_BAMS=()
    INCLUDED_PEAKS=()
    INCLUDED_STATUS=()
    INCLUDED_CALLERS=()

    for i in "${!SAMPLE_IDS[@]}"; do
        if [[ -n "${_EXCLUDE_IDX[$i]+x}" ]]; then
            EXCLUDED_SAMPLE_IDS+=("${SAMPLE_IDS[$i]}")
            EXCLUDED_SAMPLE_REASONS+=("user excluded before group assignment and DiffBind sample-sheet creation")
        else
            INCLUDED_IDS+=("${SAMPLE_IDS[$i]}")
            INCLUDED_BAMS+=("${SAMPLE_BAMS[$i]}")
            INCLUDED_PEAKS+=("${SAMPLE_PEAKS[$i]}")
            INCLUDED_STATUS+=("${SAMPLE_STATUS_ARR[$i]:-UNKNOWN}")
            INCLUDED_CALLERS+=("${SAMPLE_PEAK_CALLERS[$i]:-}")
        fi
    done

    SAMPLE_IDS=("${INCLUDED_IDS[@]}")
    SAMPLE_BAMS=("${INCLUDED_BAMS[@]}")
    SAMPLE_PEAKS=("${INCLUDED_PEAKS[@]}")
    SAMPLE_STATUS_ARR=("${INCLUDED_STATUS[@]}")
    SAMPLE_PEAK_CALLERS=("${INCLUDED_CALLERS[@]}")

    unset INCLUDED_IDS INCLUDED_BAMS INCLUDED_PEAKS INCLUDED_STATUS INCLUDED_CALLERS _EXCLUDE_IDX

    ok "Sample selection applied before DiffBind: ${#SAMPLE_IDS[@]} included."
    break
done

blank
echo -e "  ${BOLD}Included in every downstream step (${#SAMPLE_IDS[@]}):${RESET}"
for sid in "${SAMPLE_IDS[@]}"; do
    echo -e "    ${GREEN}✔${RESET}  $sid"
done

if (( ${#EXCLUDED_SAMPLE_IDS[@]} > 0 )); then
    blank
    echo -e "  ${BOLD}Excluded before DiffBind (${#EXCLUDED_SAMPLE_IDS[@]}):${RESET}"
    for i in "${!EXCLUDED_SAMPLE_IDS[@]}"; do
        echo -e "    ${RED}✘${RESET}  ${EXCLUDED_SAMPLE_IDS[$i]}"
        echo -e "       ${DIM}${EXCLUDED_SAMPLE_REASONS[$i]}${RESET}"
    done
fi

blank
ok "Only the included samples will proceed to group assignment."

# Duplicate SampleID check. Placed here deliberately -- after status-based
# filtering and any manual exclusion above, rather than right after
# parsing -- because both of those are legitimate ways to resolve a
# spurious duplicate (excluding one of the two entries by number). Checking
# earlier would block that path. Checking here catches only duplicates
# that are still present in what's actually about to go to DiffBind.
#
# This is a hard rejection with no "continue anyway" option, because the R
# script builds the DiffBind sample sheet and unconditionally stop()s on
# any exact duplicate SampleID regardless of what's chosen here -- so
# offering a choice at this point was always going to be overruled later.
declare -A _SEEN_IDS
DUPLICATE_IDS_FOUND=()
for sid in "${SAMPLE_IDS[@]}"; do
    if [[ -n "${_SEEN_IDS[$sid]+x}" ]]; then
        DUPLICATE_IDS_FOUND+=("$sid")
    fi
    _SEEN_IDS["$sid"]=1
done
unset _SEEN_IDS

if [[ ${#DUPLICATE_IDS_FOUND[@]} -gt 0 ]]; then
    blank
    err "Duplicate SampleIDs found among the samples about to be analyzed:"
    for dup in "${DUPLICATE_IDS_FOUND[@]}"; do
        err "  - $dup"
    done
    blank
    echo -e "  ${DIM}Exact duplicate SampleIDs are never allowed here -- the R script that${RESET}"
    echo -e "  ${DIM}builds the DiffBind sample sheet hard-stops on this regardless, so it's${RESET}"
    echo -e "  ${DIM}rejected now instead of wasting a full run getting there.${RESET}"
    echo -e "  ${DIM}Biological replicates need unique SampleIDs with a shared Condition --${RESET}"
    echo -e "  ${DIM}not duplicate IDs. Two samples both literally named 'WT_S1', for example,${RESET}"
    echo -e "  ${DIM}are not the same as two distinct replicate names like 'WT_S1' and 'WT_S2'${RESET}"
    echo -e "  ${DIM}assigned to the same Condition -- that second case is fine and expected.${RESET}"
    echo -e "  ${DIM}If this came from the autodetected sheet, check the PEPATAC run for a${RESET}"
    echo -e "  ${DIM}lane-merging or naming collision. Rename the conflicting sample(s), or${RESET}"
    echo -e "  ${DIM}exclude one of them in Step 2b above, then re-run this script.${RESET}"
    die "Aborted: duplicate SampleIDs must be resolved before continuing."
fi

# ─────────────────────────────────────────────────────────────
# STEP 3 — Assign group labels
# ─────────────────────────────────────────────────────────────

header "Step 3 · Group Assignment"

echo -e "  Assign a group label to each included sample."
echo -e "  ${DIM}Use the same label for samples that belong to the same condition.${RESET}"
echo -e "  ${DIM}Example: 'control', 'treated', 'KO', 'WT', 'drug_high', etc.${RESET}"
blank

while true; do
    declare -a SAMPLE_GROUPS
    declare -A GROUP_MEMBERS      # group_name -> space-separated indices
    declare -A GROUP_RAW_INPUTS   # group_name -> newline-separated distinct raw (pre-sanitize) inputs seen

    for i in "${!SAMPLE_IDS[@]}"; do
        sid="${SAMPLE_IDS[$i]}"
        while true; do
            echo -e "  Group for ${BOLD}${sid}${RESET}:"
            read -p "    > " grp_raw
            grp="${grp_raw//[^A-Za-z0-9_]/_}"   # replace any non-alphanumeric (except _) with underscore
            if [[ -z "$grp" ]]; then
                err "Group label cannot be empty."
                continue
            fi
            break
        done
        SAMPLE_GROUPS[$i]="$grp"
        GROUP_MEMBERS["$grp"]="${GROUP_MEMBERS[$grp]:-} $i"
        # Track every distinct raw string that sanitized to this group, so a
        # collision like 'WT-24h' / 'WT 24h' both becoming 'WT_24h' can be
        # surfaced below instead of silently merging two intended-different
        # groups with no trace of why.
        if [[ "${GROUP_RAW_INPUTS[$grp]:-}" != *$'\n'"${grp_raw}"$'\n'* ]]; then
            GROUP_RAW_INPUTS["$grp"]="${GROUP_RAW_INPUTS[$grp]:-}"$'\n'"${grp_raw}"$'\n'
        fi
    done

    blank
    echo -e "  ${BOLD}Group summary:${RESET}"
    # Collect unique groups in order of first appearance
    declare -a _SUMMARY_GROUPS
    for i in "${!SAMPLE_IDS[@]}"; do
        grp="${SAMPLE_GROUPS[$i]}"
        already=false
        for ug in "${_SUMMARY_GROUPS[@]:-}"; do
            [[ "$ug" == "$grp" ]] && already=true && break
        done
        $already || _SUMMARY_GROUPS+=("$grp")
    done
    # Print each group, its members, and -- if more than one distinct raw
    # input sanitized down to this same group name -- an explicit heads-up,
    # so a typo-driven merge gets caught here rather than discovered later.
    declare -a COLLIDED_GROUPS=()
    for grp in "${_SUMMARY_GROUPS[@]}"; do
        echo -e "  ${CYAN}${grp}${RESET}"
        for i in "${!SAMPLE_IDS[@]}"; do
            [[ "${SAMPLE_GROUPS[$i]}" == "$grp" ]] && echo -e "      ${DIM}${SAMPLE_IDS[$i]}${RESET}"
        done
        raw_variants_str="${GROUP_RAW_INPUTS[$grp]:-}"
        mapfile -t raw_variants < <(printf '%s' "$raw_variants_str" | sed '/^$/d')
        if [[ ${#raw_variants[@]} -gt 1 ]]; then
            warn "Multiple different inputs sanitized to the same group '$grp':"
            for rv in "${raw_variants[@]}"; do
                echo -e "      ${DIM}typed \"${rv}\"${RESET}"
            done
            COLLIDED_GROUPS+=("$grp")
        fi
    done
    unset _SUMMARY_GROUPS

    # A sanitization collision is a genuinely ambiguous situation -- it can
    # be intentional (typo-tolerant merging) or a real mistake (two
    # different conditions silently combined). Either way, it can't be
    # waved through by the general yes/no below: it needs its own explicit,
    # unambiguous confirmation, or group assignment restarts from scratch.
    if [[ ${#COLLIDED_GROUPS[@]} -gt 0 ]]; then
        blank
        warn "${#COLLIDED_GROUPS[@]} group(s) above were reached by more than one differently-typed input."
        read -p "  Type MERGE to confirm this is intentional, or press Enter to redo group assignment: " MERGE_CONFIRM
        if [[ "${MERGE_CONFIRM^^}" != "MERGE" ]]; then
            blank
            warn "Redoing group assignment."
            unset SAMPLE_GROUPS GROUP_MEMBERS GROUP_RAW_INPUTS COLLIDED_GROUPS
            blank
            continue
        fi
        ok "Confirmed: collision(s) above are intentional."
    fi
    unset COLLIDED_GROUPS

    blank
    read -p "  Are these group assignments correct? [Y/n]: " CONFIRM_GROUPS
    if [[ "${CONFIRM_GROUPS,,}" == "n" ]]; then
        blank
        warn "Redoing group assignment."
        unset SAMPLE_GROUPS GROUP_MEMBERS GROUP_RAW_INPUTS
        blank
        continue
    fi
    break
done

# Collect unique group names.
declare -a UNIQUE_GROUPS
for grp in "${SAMPLE_GROUPS[@]}"; do
    already=false
    for ug in "${UNIQUE_GROUPS[@]:-}"; do
        [[ "$ug" == "$grp" ]] && already=true && break
    done
    $already || UNIQUE_GROUPS+=("$grp")
done

# Warn if any group has only 1 sample (no replicates).
for grp in "${UNIQUE_GROUPS[@]}"; do
    count=0
    for sg in "${SAMPLE_GROUPS[@]}"; do
        [[ "$sg" == "$grp" ]] && (( count++ )) || true
    done
    if [[ $count -eq 1 ]]; then
        warn "Group '$grp' has only 1 sample (no replicates)."
        echo -e "  ${DIM}DESeq2 requires >=2 samples per group for dispersion estimation.${RESET}"
        echo -e "  ${DIM}Loading into DiffBind will still succeed, but ANY contrast that uses${RESET}"
        echo -e "  ${DIM}group '$grp' will be skipped by the R script (a [SKIP] message prints,${RESET}"
        echo -e "  ${DIM}but no DESeq2 test runs and no output is written for that contrast)${RESET}"
        echo -e "  ${DIM}unless it gains a 2nd sample.${RESET}"
    fi
done

# ─────────────────────────────────────────────────────────────
# STEP 4 — Define contrasts
# ─────────────────────────────────────────────────────────────

header "Step 4 · Define Contrasts"

echo -e "  Define one or more contrasts (pairwise comparisons)."
echo -e "  ${DIM}Each contrast compares two groups: 'case' vs 'control'.${RESET}"
echo -e "  ${DIM}Positive fold change = more accessible in the CASE group.${RESET}"
blank
echo -e "  Available groups:"
for grp in "${UNIQUE_GROUPS[@]}"; do
    echo -e "    ${CYAN}•${RESET}  $grp"
done
blank

declare -a CONTRAST_CASES
declare -a CONTRAST_CONTROLS
declare -a CONTRAST_LABELS

CONTRAST_NUM=0

while true; do
    blank
    echo -e "  ${BOLD}Contrast $((CONTRAST_NUM + 1))${RESET}"
    echo -e "  ${DIM}Enter 'done' for the case group when finished adding contrasts.${RESET}"
    blank

    # Case group.
    while true; do
        read -p "  Case group (or 'done' to finish): " CASE_GRP
        [[ "$CASE_GRP" == "done" ]] && break 2
        CASE_GRP="${CASE_GRP//[^A-Za-z0-9_]/_}"
        found=false
        for ug in "${UNIQUE_GROUPS[@]}"; do
            [[ "$ug" == "$CASE_GRP" ]] && found=true && break
        done
        $found && break
        err "Group '$CASE_GRP' not found. Choose from the list above."
    done

    # Control group.
    while true; do
        read -p "  Control group: " CTRL_GRP
        CTRL_GRP="${CTRL_GRP//[^A-Za-z0-9_]/_}"
        if [[ "$CTRL_GRP" == "$CASE_GRP" ]]; then
            err "Case and control cannot be the same group."
            continue
        fi
        found=false
        for ug in "${UNIQUE_GROUPS[@]}"; do
            [[ "$ug" == "$CTRL_GRP" ]] && found=true && break
        done
        $found && break
        err "Group '$CTRL_GRP' not found. Choose from the list above."
    done

    # Optional label. Must be unique after sanitization: this label becomes
    # the per-contrast output directory name ($DIFF_OUT/$label/), so two
    # different contrasts that sanitize to the same label (e.g. 'WT-24h_vs_KO'
    # and 'WT 24h_vs_KO' both becoming 'WT_24h_vs_KO') would silently
    # overwrite one contrast's results with the other's -- no error, no
    # warning, just a missing set of outputs discovered later. Hard-reject
    # a collision and re-prompt, the same way case/control are validated
    # above, rather than offering a choice that would just get overwritten.
    while true; do
        DEFAULT_LABEL="${CASE_GRP}_vs_${CTRL_GRP}"
        read -p "  Label for this contrast [default: ${DEFAULT_LABEL}]: " CLABEL
        CLABEL="${CLABEL:-$DEFAULT_LABEL}"
        CLABEL="${CLABEL//[^A-Za-z0-9_]/_}"
        label_collision=false
        for existing_label in "${CONTRAST_LABELS[@]:-}"; do
            [[ "$existing_label" == "$CLABEL" ]] && label_collision=true && break
        done
        if $label_collision; then
            err "Label '$CLABEL' collides with a contrast already defined (contrast labels become output folder names -- reusing one would silently overwrite that contrast's results). Enter a different label."
            continue
        fi
        break
    done

    CONTRAST_CASES+=("$CASE_GRP")
    CONTRAST_CONTROLS+=("$CTRL_GRP")
    CONTRAST_LABELS+=("$CLABEL")
    (( CONTRAST_NUM++ )) || true

    ok "Added: ${CASE_GRP} vs ${CTRL_GRP}  (label: ${CLABEL})"
done

[[ $CONTRAST_NUM -eq 0 ]] && die "No contrasts defined."

blank
echo -e "  ${BOLD}Defined contrasts:${RESET}"
for i in "${!CONTRAST_LABELS[@]}"; do
    printf "    ${CYAN}%2d${RESET}.  %-30s  (%s vs %s)\n" \
        "$((i+1))" "${CONTRAST_LABELS[$i]}" \
        "${CONTRAST_CASES[$i]}" "${CONTRAST_CONTROLS[$i]}"
done

# ─────────────────────────────────────────────────────────────
# STEP 5 — Genome and validated annotation profile
# ─────────────────────────────────────────────────────────────

header "Step 5 · Reference Genome"

SUPPORTED_GENOMES=("mm10" "hg38" "rn7" "dm6" "danRer11")
CUSTOM_GENOME_REGISTRY_ROOT="$HOME/pepatac_custom_genomes"
IS_CUSTOM_GENOME=false
CUSTOM_REFERENCE_SOURCE_FILE=""
CUSTOM_PROFILE_ASSEMBLY=""
CUSTOM_FASTA_PATH=""
CUSTOM_FASTA_SHA256=""
CUSTOM_CHROM_SIZES_PATH=""
CUSTOM_CHROM_SIZES_SHA256=""
CUSTOM_TXDB_SQLITE_PATH=""
CUSTOM_TXDB_SQLITE_SHA256=""
CUSTOM_ORGDB_SQLITE_PATH=""
CUSTOM_ORGDB_SQLITE_SHA256=""
CUSTOM_TSS_BED_PATH=""
CUSTOM_TSS_BED_SHA256=""
CUSTOM_FEATURE_BED_PATH=""
CUSTOM_FEATURE_BED_SHA256=""
CUSTOM_PROFILE_BUILD_CONF_PATH=""
CUSTOM_PROFILE_BUILD_CONF_SHA256=""
TXDB_SOURCE_VERSION=""
ORGDB_SOURCE_VERSION=""

custom_genome_registry_file() {
    printf '%s/%s/genome.conf\n' "$CUSTOM_GENOME_REGISTRY_ROOT" "$1"
}

ucsc_profile_fasta_url() { printf 'https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.fa.gz\n' "$1" "$1"; }
ucsc_profile_chrom_sizes_url() { printf 'https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes\n' "$1" "$1"; }

sha256_of() {
    local f="$1"
    [[ -f "$f" ]] || { printf '\n'; return 0; }
    command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required for reference-integrity checks."
    sha256sum "$f" | awk '{print $1}'
}

verify_registered_sha256() {
    local f="$1" expected="$2" label_text="$3"
    [[ -n "$expected" ]] || die "No registered SHA-256 is available for $label_text: $f"
    [[ -f "$f" ]] || die "$label_text not found: $f"
    local actual
    actual="$(sha256_of "$f")"
    [[ "$actual" == "$expected" ]] \
        || die "$label_text hash does not match the run snapshot."$'\n'"       File: $f"$'\n'"       Snapshot: $expected"$'\n'"       Actual:   $actual"$'\n'"       Differential analysis will not use altered annotation/reference bytes."
}

resolve_custom_genome_registry_path() {
    local genome="$1" snapshot_file="${RUN_DIR:-}/reference_snapshot/genome.conf" snapshot_name=""
    if [[ -n "${RUN_DIR:-}" && -f "$snapshot_file" ]]; then
        snapshot_name="$(awk -F= '$1=="CUSTOM_GENOME_NAME"{print substr($0,index($0,"=")+1); exit}' "$snapshot_file" 2>/dev/null || true)"
        if [[ "$snapshot_name" == "$genome" ]]; then
            printf '%s\n' "$snapshot_file"
            return 0
        fi
    fi
    custom_genome_registry_file "$genome"
}

load_custom_genome_registry() {
    local genome="$1"
    local reg_file="${2:-$(custom_genome_registry_file "$genome")}" line key val
    [[ -f "$reg_file" ]] || return 1

    local -A rg=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        key="${line%%=*}"; val="${line#*=}"
        [[ "$key" =~ ^[A-Za-z0-9_]+$ ]] && rg["$key"]="$val"
    done < "$reg_file"

    local snapshot_schema="${rg[REFERENCE_SNAPSHOT_SCHEMA]:-0}"
    local registry_schema="${rg[CUSTOM_REFERENCE_REGISTRY_SCHEMA]:-0}"
    if [[ "$reg_file" == */reference_snapshot/genome.conf ]]; then
        [[ "$snapshot_schema" =~ ^[0-9]+$ ]] || die "Reference snapshot has an invalid schema: $snapshot_schema"
        case "$snapshot_schema" in
            3|4) ;;
            *)
                if (( snapshot_schema > REFERENCE_SNAPSHOT_SCHEMA_SUPPORTED )); then
                    die "Reference snapshot schema $snapshot_schema is newer than this differential script understands (maximum $REFERENCE_SNAPSHOT_SCHEMA_SUPPORTED)."
                fi
                die "Reference snapshot schema ${snapshot_schema:-0} predates validated UCSC/TxDb/OrgDb profiles. Re-run sample processing with runner v1.25 or later."
                ;;
        esac
    else
        [[ "$registry_schema" =~ ^[0-9]+$ ]] || die "Custom profile registry has an invalid schema: $registry_schema"
        (( registry_schema >= 3 )) \
            || die "Custom profile registry schema ${registry_schema:-0} predates validated UCSC/TxDb/OrgDb profiles. Re-register this assembly with runner v1.25 or later."
    fi

    CUSTOM_REG_GENOME_NAME="${rg[CUSTOM_GENOME_NAME]:-$genome}"
    CUSTOM_REG_UCSC_ASSEMBLY="${rg[CUSTOM_UCSC_ASSEMBLY]:-$genome}"
    CUSTOM_REG_FASTA_SOURCE="${rg[CUSTOM_FASTA_SOURCE]:-}"
    CUSTOM_REG_CHROM_SOURCE="${rg[CUSTOM_CHROM_SIZES_SOURCE]:-}"
    CUSTOM_REG_FASTA="${rg[CUSTOM_FASTA_CACHE]:-}"
    CUSTOM_REG_FASTA_SHA256="${rg[CUSTOM_FASTA_SHA256]:-}"
    CUSTOM_REG_CHROM_SIZES="${rg[CUSTOM_CHROM_SIZES]:-}"
    CUSTOM_REG_CHROM_SIZES_SHA256="${rg[CUSTOM_CHROM_SIZES_SHA256]:-}"
    CUSTOM_REG_TXDB_PKG="${rg[CUSTOM_TXDB_PACKAGE]:-}"
    CUSTOM_REG_TXDB_VERSION="${rg[CUSTOM_TXDB_PACKAGE_VERSION]:-}"
    CUSTOM_REG_ORGDB_PKG="${rg[CUSTOM_ORGDB_PACKAGE]:-}"
    CUSTOM_REG_ORGDB_VERSION="${rg[CUSTOM_ORGDB_PACKAGE_VERSION]:-}"
    CUSTOM_REG_TXDB_SQLITE="${rg[CUSTOM_TXDB_SQLITE]:-}"
    CUSTOM_REG_TXDB_SQLITE_SHA256="${rg[CUSTOM_TXDB_SQLITE_SHA256]:-}"
    CUSTOM_REG_ORGDB_SQLITE="${rg[CUSTOM_ORGDB_SQLITE]:-}"
    CUSTOM_REG_ORGDB_SQLITE_SHA256="${rg[CUSTOM_ORGDB_SQLITE_SHA256]:-}"
    CUSTOM_REG_TSS_BED="${rg[CUSTOM_TSS_BED]:-}"
    CUSTOM_REG_TSS_BED_SHA256="${rg[CUSTOM_TSS_BED_SHA256]:-}"
    CUSTOM_REG_FEATURE_BED="${rg[CUSTOM_FEATURE_BED]:-}"
    CUSTOM_REG_FEATURE_BED_SHA256="${rg[CUSTOM_FEATURE_BED_SHA256]:-}"
    CUSTOM_REG_PROFILE_BUILD_CONF="${rg[CUSTOM_PROFILE_BUILD_CONF]:-}"
    CUSTOM_REG_PROFILE_BUILD_CONF_SHA256="${rg[CUSTOM_PROFILE_BUILD_CONF_SHA256]:-}"
    CUSTOM_REG_KEGG_ORG="${rg[CUSTOM_KEGG_ORG]:-}"
    CUSTOM_REG_SOURCE_FILE="$reg_file"
    return 0
}

# Try to auto-detect genome/profile name from the run manifest.
DETECTED_GENOME=""
RUN_MANIFEST="$RUN_DIR/run_manifest.txt"
if [[ -f "$RUN_MANIFEST" ]]; then
    DETECTED_GENOME=$(grep -m1 "^Genome:" "$RUN_MANIFEST" 2>/dev/null | awk '{print $2}' || true)
fi

if [[ -n "$DETECTED_GENOME" ]]; then
    echo -e "  Detected genome from run manifest: ${BOLD}${DETECTED_GENOME}${RESET}"
    read -p "  Use this genome? [Y/n]: " USE_DETECTED
    if [[ "${USE_DETECTED,,}" != "n" ]]; then
        GENOME="$DETECTED_GENOME"
        ok "Using genome: $GENOME"
    else
        DETECTED_GENOME=""
    fi
fi

if [[ -z "$DETECTED_GENOME" ]]; then
    echo -e "  Built-in genomes:"
    for i in "${!SUPPORTED_GENOMES[@]}"; do
        printf "    ${CYAN}%2d${RESET}.  %s\n" "$((i+1))" "${SUPPORTED_GENOMES[$i]}"
    done
    echo -e "  ${DIM}Or type a validated UCSC/TxDb/OrgDb profile name registered by run.sh.${RESET}"
    blank

    while true; do
        read -p "  Enter genome/profile name or number: " GENOME_INPUT
        if [[ "$GENOME_INPUT" =~ ^[0-9]+$ ]]; then
            IDX=$((GENOME_INPUT - 1))
            if [[ $IDX -ge 0 && $IDX -lt ${#SUPPORTED_GENOMES[@]} ]]; then
                GENOME="${SUPPORTED_GENOMES[$IDX]}"
            else
                err "Invalid number."
                continue
            fi
        else
            GENOME="$GENOME_INPUT"
        fi

        valid=false
        for sg in "${SUPPORTED_GENOMES[@]}"; do
            [[ "$sg" == "$GENOME" ]] && valid=true && break
        done
        if $valid; then
            ok "Genome: $GENOME"
            break
        elif load_custom_genome_registry "$GENOME" "$(resolve_custom_genome_registry_path "$GENOME")"; then
            ok "Genome: $GENOME (validated user-defined profile)"
            break
        else
            err "Unsupported genome/profile: $GENOME."
            err "No validated run snapshot or profile registry was found."
            echo -e "  ${DIM}Register the UCSC assembly with matching TxDb and OrgDb through run.sh first.${RESET}"
        fi
    done
fi

# Built-ins use installed, audited annotation packages by default (tier 3
# of the resolution order below). If this run has its own frozen,
# fingerprinted TxDb/OrgDb (schema-4 snapshot written by run.sh v1.26+),
# that's preferred instead (tier 1) -- see try_load_run_annotation_snapshot.
# User-defined profiles always use their own frozen TxDb/OrgDb SQLite
# files (tier 1/2, unchanged from before).
#
# Resolution order, matching the order actually implemented below:
#   1. Schema-4 generic frozen annotation snapshot (built-in or custom)
#   2. Schema-3 custom frozen-profile snapshot (legacy custom runs)
#   3. Legacy built-in installed TxDb/OrgDb fallback (no snapshot at all)
#   4. Otherwise hard-stop (custom genomes only -- see the case statement;
#      a custom profile with no usable registry/snapshot is a die(), not a
#      fallback, since there is no "installed package" to fall back to)
TXDB_PKG=""
ORGDB_PKG=""
KEGG_ORG=""

# try_load_run_annotation_snapshot GENOME
# Tier 1. Looks for this run's own reference_snapshot/genome.conf, and
# uses its frozen TxDb/OrgDb only if: it identifies this exact genome, its
# schema is 4 (schema 3 predates built-in annotation snapshots entirely,
# so it never has these fields), and its recorded ANNOTATION_STATUS is
# "validated" (not "disabled_by_override" -- see run.sh's
# PEPATAC_ALLOW_MISSING_ANNOTATION_QC). Hashes are verified immediately
# via the existing die-on-mismatch verify_registered_sha256 -- a hash
# mismatch on a supposedly-immutable frozen file is a real integrity
# problem, not something to silently fall back past.
try_load_run_annotation_snapshot() {
    local genome="$1" snapshot_file="${RUN_DIR:-}/reference_snapshot/genome.conf"
    [[ -n "${RUN_DIR:-}" && -f "$snapshot_file" ]] || return 1

    local -A ss=() line key val
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        key="${line%%=*}"; val="${line#*=}"
        [[ "$key" =~ ^[A-Za-z0-9_]+$ ]] && ss["$key"]="$val"
    done < "$snapshot_file"

    [[ "${ss[CUSTOM_GENOME_NAME]:-}" == "$genome" ]] || return 1

    local schema="${ss[REFERENCE_SNAPSHOT_SCHEMA]:-0}"
    [[ "$schema" =~ ^[0-9]+$ ]] || die "Reference snapshot has an invalid schema value: $schema"
    # Schema < 4 predates built-in annotation snapshots entirely -- a
    # genuine legacy case, fall back to tier 3. Schema > 4 is newer than
    # this script understands -- that IS malformed/unsupported, not legacy,
    # so it dies rather than silently falling back past data it can't read.
    (( schema < 4 )) && return 1
    (( schema == 4 )) || die "Reference snapshot schema $schema is newer than this differential script understands for built-in annotation (maximum 4)."

    local status="${ss[ANNOTATION_STATUS]:-}"
    case "$status" in
        disabled_by_override)
            warn "This run's snapshot recorded annotation QC as disabled for '$genome'."
            warn "Reason: ${ss[ANNOTATION_QC_FAILURE_REASON]:-<not recorded>}"
            warn "Falling back to the currently-installed TxDb/OrgDb packages (tier 3, legacy behavior)."
            warn "This run lacks frozen annotation provenance -- results may differ from a run with validated QC."
            return 1
            ;;
        validated)
            ;;
        *)
            # A schema-4 snapshot with a missing or unrecognized status is
            # not a legitimate "no snapshot" case -- it's a corrupted or
            # incomplete one. Silently falling back to installed packages
            # here would hide exactly the kind of drift this snapshot
            # exists to catch, so this hard-stops instead.
            die "Reference snapshot for '$genome' has schema 4 but an unrecognized or missing ANNOTATION_STATUS ('${status:-<empty>}'). This indicates a corrupted or incomplete snapshot -- refusing to silently fall back to installed packages."
            ;;
    esac

    ANNOTATION_SNAPSHOT_TXDB_PKG="${ss[ANNOTATION_TXDB_PACKAGE]:-}"
    ANNOTATION_SNAPSHOT_TXDB_VERSION="${ss[ANNOTATION_TXDB_PACKAGE_VERSION]:-}"
    ANNOTATION_SNAPSHOT_TXDB_SQLITE="${ss[ANNOTATION_TXDB_SQLITE]:-}"
    ANNOTATION_SNAPSHOT_TXDB_SQLITE_SHA256="${ss[ANNOTATION_TXDB_SQLITE_SHA256]:-}"
    ANNOTATION_SNAPSHOT_ORGDB_PKG="${ss[ANNOTATION_ORGDB_PACKAGE]:-}"
    ANNOTATION_SNAPSHOT_ORGDB_VERSION="${ss[ANNOTATION_ORGDB_PACKAGE_VERSION]:-}"
    ANNOTATION_SNAPSHOT_ORGDB_SQLITE="${ss[ANNOTATION_ORGDB_SQLITE]:-}"
    ANNOTATION_SNAPSHOT_ORGDB_SQLITE_SHA256="${ss[ANNOTATION_ORGDB_SQLITE_SHA256]:-}"
    # Chrom sizes is what lets this script independently repeat the
    # TxDb-vs-assembly-length cross-check for a built-in genome too, the
    # same way it already can for a custom profile via CUSTOM_CHROM_SIZES.
    # Optional field (blank on a snapshot written before this was added) --
    # not required for tier 1 to succeed, just for the extra cross-check.
    ANNOTATION_SNAPSHOT_CHROM_SIZES="${ss[ANNOTATION_CHROM_SIZES]:-}"
    ANNOTATION_SNAPSHOT_CHROM_SIZES_SHA256="${ss[ANNOTATION_CHROM_SIZES_SHA256]:-}"

    if [[ -z "$ANNOTATION_SNAPSHOT_TXDB_SQLITE" || -z "$ANNOTATION_SNAPSHOT_ORGDB_SQLITE" ]]; then
        die "Reference snapshot for '$genome' claims ANNOTATION_STATUS=validated but is missing required TxDb/OrgDb SQLite paths. This indicates a corrupted or incomplete snapshot."
    fi
    verify_registered_sha256 "$ANNOTATION_SNAPSHOT_TXDB_SQLITE" "$ANNOTATION_SNAPSHOT_TXDB_SQLITE_SHA256" "Run TxDb SQLite"
    verify_registered_sha256 "$ANNOTATION_SNAPSHOT_ORGDB_SQLITE" "$ANNOTATION_SNAPSHOT_ORGDB_SQLITE_SHA256" "Run OrgDb SQLite"
    [[ -n "$ANNOTATION_SNAPSHOT_CHROM_SIZES" ]] && \
        verify_registered_sha256 "$ANNOTATION_SNAPSHOT_CHROM_SIZES" "$ANNOTATION_SNAPSHOT_CHROM_SIZES_SHA256" "Run annotation chrom.sizes"
    return 0
}

case "$GENOME" in
    hg38)     TXDB_PKG="TxDb.Hsapiens.UCSC.hg38.knownGene"; ORGDB_PKG="org.Hs.eg.db"; KEGG_ORG="hsa" ;;
    mm10)     TXDB_PKG="TxDb.Mmusculus.UCSC.mm10.knownGene"; ORGDB_PKG="org.Mm.eg.db"; KEGG_ORG="mmu" ;;
    rn7)      TXDB_PKG="TxDb.Rnorvegicus.UCSC.rn7.refGene"; ORGDB_PKG="org.Rn.eg.db"; KEGG_ORG="rno" ;;
    dm6)      TXDB_PKG="TxDb.Dmelanogaster.UCSC.dm6.ensGene"; ORGDB_PKG="org.Dm.eg.db"; KEGG_ORG="dme" ;;
    danRer11) TXDB_PKG="TxDb.Drerio.UCSC.danRer11.refGene"; ORGDB_PKG="org.Dr.eg.db"; KEGG_ORG="dre" ;;
    *)
        IS_CUSTOM_GENOME=true
        CUSTOM_REFERENCE_SOURCE_FILE="$(resolve_custom_genome_registry_path "$GENOME")"
        load_custom_genome_registry "$GENOME" "$CUSTOM_REFERENCE_SOURCE_FILE" \
            || die "Validated profile for '$GENOME' disappeared since selection."

        if [[ "$CUSTOM_REFERENCE_SOURCE_FILE" == "$RUN_DIR/reference_snapshot/genome.conf" ]]; then
            ok "Using this run's immutable validated profile snapshot: $CUSTOM_REFERENCE_SOURCE_FILE"
        else
            if [[ "${PEPATAC_ALLOW_LIVE_CUSTOM_REGISTRY:-0}" != "1" ]]; then
                die "No immutable validated-profile snapshot was found in $RUN_DIR/reference_snapshot/genome.conf."$'\n'"       Refusing to use the live registry because it may describe a newer profile version."$'\n'"       Re-run sample processing with runner v1.25+, or set PEPATAC_ALLOW_LIVE_CUSTOM_REGISTRY=1 only for deliberate recovery."
            fi
            warn "Legacy recovery override enabled: using the current live profile registry."
        fi

        [[ "$CUSTOM_REG_GENOME_NAME" == "$GENOME" ]] \
            || die "Profile-name mismatch: run genome is '$GENOME' but snapshot identifies '$CUSTOM_REG_GENOME_NAME'."
        CUSTOM_PROFILE_ASSEMBLY="$CUSTOM_REG_UCSC_ASSEMBLY"
        [[ "$CUSTOM_REG_FASTA_SOURCE" == "$(ucsc_profile_fasta_url "$CUSTOM_PROFILE_ASSEMBLY")" && \
           "$CUSTOM_REG_CHROM_SOURCE" == "$(ucsc_profile_chrom_sizes_url "$CUSTOM_PROFILE_ASSEMBLY")" ]] \
            || die "Validated profile snapshot is not tied to the official UCSC sequence sources for '$CUSTOM_PROFILE_ASSEMBLY'."
        [[ "$CUSTOM_REG_TXDB_PKG" == *".UCSC.${CUSTOM_PROFILE_ASSEMBLY}."* ]] \
            || die "Validated profile TxDb package '$CUSTOM_REG_TXDB_PKG' does not identify UCSC assembly '$CUSTOM_PROFILE_ASSEMBLY'."
        CUSTOM_FASTA_PATH="$CUSTOM_REG_FASTA"; CUSTOM_FASTA_SHA256="$CUSTOM_REG_FASTA_SHA256"
        CUSTOM_CHROM_SIZES_PATH="$CUSTOM_REG_CHROM_SIZES"; CUSTOM_CHROM_SIZES_SHA256="$CUSTOM_REG_CHROM_SIZES_SHA256"
        CUSTOM_TXDB_SQLITE_PATH="$CUSTOM_REG_TXDB_SQLITE"; CUSTOM_TXDB_SQLITE_SHA256="$CUSTOM_REG_TXDB_SQLITE_SHA256"
        CUSTOM_ORGDB_SQLITE_PATH="$CUSTOM_REG_ORGDB_SQLITE"; CUSTOM_ORGDB_SQLITE_SHA256="$CUSTOM_REG_ORGDB_SQLITE_SHA256"
        CUSTOM_TSS_BED_PATH="$CUSTOM_REG_TSS_BED"; CUSTOM_TSS_BED_SHA256="$CUSTOM_REG_TSS_BED_SHA256"
        CUSTOM_FEATURE_BED_PATH="$CUSTOM_REG_FEATURE_BED"; CUSTOM_FEATURE_BED_SHA256="$CUSTOM_REG_FEATURE_BED_SHA256"
        CUSTOM_PROFILE_BUILD_CONF_PATH="$CUSTOM_REG_PROFILE_BUILD_CONF"; CUSTOM_PROFILE_BUILD_CONF_SHA256="$CUSTOM_REG_PROFILE_BUILD_CONF_SHA256"
        TXDB_PKG="$CUSTOM_REG_TXDB_PKG"; TXDB_SOURCE_VERSION="$CUSTOM_REG_TXDB_VERSION"
        ORGDB_PKG="$CUSTOM_REG_ORGDB_PKG"; ORGDB_SOURCE_VERSION="$CUSTOM_REG_ORGDB_VERSION"
        KEGG_ORG="$CUSTOM_REG_KEGG_ORG"

        [[ -n "$CUSTOM_PROFILE_ASSEMBLY" && -n "$TXDB_PKG" && -n "$TXDB_SOURCE_VERSION" && -n "$ORGDB_PKG" && -n "$ORGDB_SOURCE_VERSION" ]] \
            || die "Validated profile snapshot is incomplete: UCSC assembly plus TxDb/OrgDb package names and versions are required."
        verify_registered_sha256 "$CUSTOM_FASTA_PATH" "$CUSTOM_FASTA_SHA256" "Profile FASTA"
        verify_registered_sha256 "$CUSTOM_CHROM_SIZES_PATH" "$CUSTOM_CHROM_SIZES_SHA256" "Profile chromosome sizes"
        verify_registered_sha256 "$CUSTOM_TXDB_SQLITE_PATH" "$CUSTOM_TXDB_SQLITE_SHA256" "Frozen TxDb SQLite"
        verify_registered_sha256 "$CUSTOM_ORGDB_SQLITE_PATH" "$CUSTOM_ORGDB_SQLITE_SHA256" "Frozen OrgDb SQLite"
        verify_registered_sha256 "$CUSTOM_TSS_BED_PATH" "$CUSTOM_TSS_BED_SHA256" "Profile TSS BED"
        verify_registered_sha256 "$CUSTOM_FEATURE_BED_PATH" "$CUSTOM_FEATURE_BED_SHA256" "Profile feature BED"
        verify_registered_sha256 "$CUSTOM_PROFILE_BUILD_CONF_PATH" "$CUSTOM_PROFILE_BUILD_CONF_SHA256" "Profile build manifest"
        ok "All immutable profile hashes match this run's snapshot."
        ;;
esac

# Tier 1 for built-ins: prefer this run's own frozen, fingerprinted
# TxDb/OrgDb over whatever happens to be installed right now, if this
# run's snapshot has one. Custom genomes already resolved their frozen
# TxDb/OrgDb inside the case statement above (tiers 1/2 there, unchanged);
# this only applies to the five built-ins, which previously had no
# snapshot-based path at all (always tier 3, unconditionally).
ANNOTATION_SNAPSHOT_TXDB_PKG=""; ANNOTATION_SNAPSHOT_TXDB_VERSION=""
ANNOTATION_SNAPSHOT_TXDB_SQLITE=""; ANNOTATION_SNAPSHOT_TXDB_SQLITE_SHA256=""
ANNOTATION_SNAPSHOT_ORGDB_PKG=""; ANNOTATION_SNAPSHOT_ORGDB_VERSION=""
ANNOTATION_SNAPSHOT_ORGDB_SQLITE=""; ANNOTATION_SNAPSHOT_ORGDB_SQLITE_SHA256=""
BUILTIN_ANNOTATION_FROM_SNAPSHOT=false
if ! $IS_CUSTOM_GENOME; then
    if try_load_run_annotation_snapshot "$GENOME"; then
        TXDB_PKG="$ANNOTATION_SNAPSHOT_TXDB_PKG"; TXDB_SOURCE_VERSION="$ANNOTATION_SNAPSHOT_TXDB_VERSION"
        ORGDB_PKG="$ANNOTATION_SNAPSHOT_ORGDB_PKG"; ORGDB_SOURCE_VERSION="$ANNOTATION_SNAPSHOT_ORGDB_VERSION"
        BUILTIN_ANNOTATION_FROM_SNAPSHOT=true
        ok "Using this run's frozen validated TxDb/OrgDb for $GENOME (matches its upstream PEPATAC QC)."
    else
        ok "No usable frozen annotation snapshot for $GENOME -- loading the currently-installed TxDb/OrgDb packages."
    fi
fi

# Single source of truth for what actually gets loaded in R below,
# regardless of which tier supplied it -- collapses what used to be a
# hard custom-vs-built-in branch into "do we have a frozen SQLite for
# this run, or not." See the R body's loadDb()/library() branch.
FROZEN_TXDB_SQLITE=""; FROZEN_TXDB_SQLITE_SHA256=""
FROZEN_ORGDB_SQLITE=""; FROZEN_ORGDB_SQLITE_SHA256=""
FROZEN_CHROM_SIZES=""; FROZEN_CHROM_SIZES_SHA256=""
USE_FROZEN_ANNOTATION=false
if $IS_CUSTOM_GENOME; then
    FROZEN_TXDB_SQLITE="$CUSTOM_TXDB_SQLITE_PATH"; FROZEN_TXDB_SQLITE_SHA256="$CUSTOM_TXDB_SQLITE_SHA256"
    FROZEN_ORGDB_SQLITE="$CUSTOM_ORGDB_SQLITE_PATH"; FROZEN_ORGDB_SQLITE_SHA256="$CUSTOM_ORGDB_SQLITE_SHA256"
    FROZEN_CHROM_SIZES="$CUSTOM_CHROM_SIZES_PATH"; FROZEN_CHROM_SIZES_SHA256="$CUSTOM_CHROM_SIZES_SHA256"
    USE_FROZEN_ANNOTATION=true
elif $BUILTIN_ANNOTATION_FROM_SNAPSHOT; then
    FROZEN_TXDB_SQLITE="$ANNOTATION_SNAPSHOT_TXDB_SQLITE"; FROZEN_TXDB_SQLITE_SHA256="$ANNOTATION_SNAPSHOT_TXDB_SQLITE_SHA256"
    FROZEN_ORGDB_SQLITE="$ANNOTATION_SNAPSHOT_ORGDB_SQLITE"; FROZEN_ORGDB_SQLITE_SHA256="$ANNOTATION_SNAPSHOT_ORGDB_SQLITE_SHA256"
    # May be blank for a snapshot written before ANNOTATION_CHROM_SIZES
    # existed -- the R-side cross-check below already handles a blank
    # CUSTOM_CHROM_SIZES by simply not running (same as it always has for
    # a built-in genome with no cross-check data available).
    FROZEN_CHROM_SIZES="$ANNOTATION_SNAPSHOT_CHROM_SIZES"; FROZEN_CHROM_SIZES_SHA256="$ANNOTATION_SNAPSHOT_CHROM_SIZES_SHA256"
    USE_FROZEN_ANNOTATION=true
fi

if $IS_CUSTOM_GENOME; then
    ok "UCSC assembly:        $CUSTOM_PROFILE_ASSEMBLY"
    ok "Coordinate annotation: $TXDB_PKG ${TXDB_SOURCE_VERSION:-<version not recorded>} (frozen TxDb)"
    ok "Gene metadata / GO:  $ORGDB_PKG ${ORGDB_SOURCE_VERSION:-<version not recorded>} (frozen OrgDb)"
    [[ -n "$KEGG_ORG" ]] && ok "KEGG organism code: $KEGG_ORG" \
        || warn "No KEGG organism code registered -- KEGG enrichment will be skipped."
elif $BUILTIN_ANNOTATION_FROM_SNAPSHOT; then
    ok "Coordinate annotation: $TXDB_PKG ${TXDB_SOURCE_VERSION:-<version not recorded>} (frozen TxDb, this run's snapshot)"
    ok "Gene metadata / GO:  $ORGDB_PKG ${ORGDB_SOURCE_VERSION:-<version not recorded>} (frozen OrgDb, this run's snapshot)"
else
    ok "Coordinate annotation: $TXDB_PKG (installed package)"
    ok "Gene metadata / GO:  $ORGDB_PKG (installed package)"
fi

# Frozen SQLite objects (custom, or built-in with a snapshot) don't need
# the specific TxDb/OrgDb data package installed at all -- AnnotationDbi::
# loadDb() only needs the generic TxDb/OrgDb S4 classes, which are already
# a hard Bioconductor dependency. Only the "no frozen SQLite available"
# case needs the specific package present.
TXDB_INSTALL_PKG="$TXDB_PKG"
ORGDB_INSTALL_PKG="$ORGDB_PKG"
if $USE_FROZEN_ANNOTATION; then
    TXDB_INSTALL_PKG=""
    ORGDB_INSTALL_PKG=""
fi

# ─────────────────────────────────────────────────────────────
# STEP 5b — Optional HOMER motif analysis
# ─────────────────────────────────────────────────────────────

header "Step 5b · Optional HOMER Motifs"

HOMER_CONFIGURE="$(conda run -n "$ENV_NAME" find "$HOME/miniconda3/envs/$ENV_NAME" \
    -name configureHomer.pl 2>/dev/null | head -1)"
HOMER_FIND_MOTIFS="$(conda run -n "$ENV_NAME" find "$HOME/miniconda3/envs/$ENV_NAME" \
    -name findMotifsGenome.pl 2>/dev/null | head -1)"

# The conda-installed HOMER version, recorded in provenance below so a
# future HOMER release (motif databases and algorithm both change over
# time) doesn't leave published motif results with no record of what
# actually produced them.
HOMER_VERSION="$(conda list -n "$ENV_NAME" 2>/dev/null | awk '$1=="homer"{print $2; exit}')"
HOMER_VERSION="${HOMER_VERSION:-unknown}"
HOMER_GENOME_INFO="not applicable (motif analysis disabled)"

RUN_MOTIF_ANALYSIS=false
HOMER_GENOME_READY=false
MOTIF_MODE="known"
# What gets passed as findMotifsGenome.pl's <genome> argument. For the five
# built-ins this is the installed HOMER genome name ($GENOME). HOMER also
# accepts a raw FASTA path in that same argument position, which is what a
# custom genome uses instead -- there's no HOMER-catalog name for it to
# install, so the install step below is skipped entirely for custom genomes.
HOMER_GENOME_ARG="$GENOME"

echo -e "  ${DIM}Peak annotation uses ChIPseeker regardless of this choice.${RESET}"
echo -e "  ${DIM}HOMER is used only for optional motif enrichment.${RESET}"
blank
read -p "  Run HOMER motif analysis on significant peaks? [Y/n]: " RUN_MOTIF_INPUT

if [[ "${RUN_MOTIF_INPUT,,}" == "n" ]]; then
    ok "HOMER motif analysis disabled."
elif [[ -z "$HOMER_CONFIGURE" || -z "$HOMER_FIND_MOTIFS" ]]; then
    warn "HOMER tools were not found in the pepatac environment."
    echo -e "  ${DIM}Peak annotation and GO/KEGG will still run normally.${RESET}"
elif $IS_CUSTOM_GENOME; then
    # A user-defined profile does not rely on a HOMER catalog entry.
    # findMotifsGenome.pl takes a FASTA path
    # directly in place of a genome name -- use the same FASTA run.sh built.
    if [[ -n "$CUSTOM_FASTA_PATH" && -f "$CUSTOM_FASTA_PATH" ]]; then
        HOMER_GENOME_ARG="$CUSTOM_FASTA_PATH"
        HOMER_GENOME_READY=true
        HOMER_GENOME_INFO="custom genome (direct FASTA, no HOMER genome package): $CUSTOM_FASTA_PATH"
        ok "Using local FASTA directly for HOMER (no install needed): $CUSTOM_FASTA_PATH"
    else
        warn "Registered FASTA for '$GENOME' not found: ${CUSTOM_FASTA_PATH:-<none recorded>}"
        warn "HOMER motif analysis will be skipped for this run."
    fi
else
    HOMER_GENOME_DIR="$(find "$HOME/miniconda3/envs/$ENV_NAME/share/homer" \
        -type d -name "$GENOME" 2>/dev/null | head -1)"

    if [[ -d "$HOMER_GENOME_DIR" ]]; then
        ok "HOMER motif genome '$GENOME' is installed: $HOMER_GENOME_DIR"
        HOMER_GENOME_READY=true
        HOMER_GENOME_INFO="$GENOME (HOMER genome package: $HOMER_GENOME_DIR, last modified $(date -r "$HOMER_GENOME_DIR" '+%Y-%m-%d' 2>/dev/null || echo unknown))"
    else
        warn "HOMER motif genome '$GENOME' is not installed."
        echo -e "  ${DIM}Only motif enrichment needs this extra genome package (~1-3 GB).${RESET}"
        blank
        read -p "  Download the HOMER motif genome '$GENOME' now? [Y/n]: " DL_HOMER
        if [[ "${DL_HOMER,,}" != "n" ]]; then
            if in_env perl "$HOMER_CONFIGURE" -install "$GENOME"; then
                ok "HOMER motif genome '$GENOME' installed successfully."
                HOMER_GENOME_READY=true
                # HOMER_GENOME_DIR was captured before install, when the
                # package didn't exist yet -- still empty here. Re-resolve
                # it now that configureHomer.pl has actually created it, or
                # the provenance record below would save a blank path.
                HOMER_GENOME_DIR="$(find "$HOME/miniconda3/envs/$ENV_NAME/share/homer" \
                    -type d -name "$GENOME" 2>/dev/null | head -1)"
                HOMER_GENOME_INFO="$GENOME (HOMER genome package: $HOMER_GENOME_DIR, installed $(date '+%Y-%m-%d'))"
            else
                warn "HOMER genome download failed; motif analysis will be skipped."
            fi
        else
            warn "HOMER motif analysis skipped for this run."
        fi
    fi
fi

# Mode selection applies uniformly whichever branch above set
# HOMER_GENOME_READY -- custom (FASTA-direct) or built-in (installed name).
if $HOMER_GENOME_READY; then
    RUN_MOTIF_ANALYSIS=true
    blank
    echo -e "  ${BOLD}Which motif analysis?${RESET}"
    echo -e "    ${CYAN}1${RESET}.  Known motifs only  ${DIM}(fast)${RESET}"
    echo -e "    ${CYAN}2${RESET}.  De novo only       ${DIM}(slow)${RESET}"
    echo -e "    ${CYAN}3${RESET}.  Both               ${DIM}(slowest)${RESET}"
    blank
    while true; do
        read -p "  Choice [1-3, default: 1]: " MOTIF_CHOICE
        MOTIF_CHOICE="${MOTIF_CHOICE:-1}"
        case "$MOTIF_CHOICE" in
            1) MOTIF_MODE="known";  break ;;
            2) MOTIF_MODE="denovo"; break ;;
            3) MOTIF_MODE="both";   break ;;
            *) err "Enter 1, 2, or 3." ;;
        esac
    done
    ok "HOMER motif mode: $MOTIF_MODE"
fi

# ─────────────────────────────────────────────────────────────
# STEP 6 — CPU Threads for DiffBind
# ─────────────────────────────────────────────────────────────
# NOTE: When appending diff_analysis.sh to run.sh, replace this
# entire step with:  DIFF_THREADS="$THREADS"
# so DiffBind automatically inherits the run-wide thread count.
# ─────────────────────────────────────────────────────────────

header "Step 6 · Analysis Parameters"

echo -e "  Configure the key thresholds for differential accessibility analysis."
blank

# ── MIN_OVERLAP ──────────────────────────────────────────────
echo -e "  ${BOLD}MIN_OVERLAP${RESET} — minimum number of samples a peak must appear in"
echo -e "  to enter the consensus peak set."
echo -e "  ${DIM}• 2 = peak must be called in at least 2 samples (DiffBind default; recommended).${RESET}"
echo -e "  ${DIM}• Higher values → more stringent consensus, may discard condition-specific peaks.${RESET}"
echo -e "  ${DIM}• For n=6 samples (2 replicates × 3 conditions), use 2 to retain replicated peaks.${RESET}"
blank

while true; do
    read -p "  MIN_OVERLAP [default: 2]: " MIN_OVERLAP_INPUT
    if [[ -z "$MIN_OVERLAP_INPUT" ]]; then
        MIN_OVERLAP_VAL=2
        break
    elif [[ "$MIN_OVERLAP_INPUT" =~ ^[0-9]+$ ]] && [[ "$MIN_OVERLAP_INPUT" -ge 1 ]]; then
        MIN_OVERLAP_VAL="$MIN_OVERLAP_INPUT"
        break
    else
        err "Please enter a positive integer."
    fi
done
ok "MIN_OVERLAP = $MIN_OVERLAP_VAL"
blank

# ── FDR cutoff ───────────────────────────────────────────────
echo -e "  ${BOLD}FDR cutoff${RESET} — adjusted p-value threshold for calling a peak significant."
echo -e "  ${DIM}Standard: 0.05 (5% false discovery rate).${RESET}"
blank

while true; do
    read -p "  FDR cutoff [default: 0.05]: " FDR_INPUT
    if [[ -z "$FDR_INPUT" ]]; then
        FDR_CUTOFF_VAL="0.05"
        break
    elif [[ "$FDR_INPUT" =~ ^0?\.[0-9]+$ ]] || [[ "$FDR_INPUT" =~ ^1(\.0+)?$ ]]; then
        FDR_CUTOFF_VAL="$FDR_INPUT"
        break
    else
        err "Please enter a number between 0 and 1 (e.g. 0.05)."
    fi
done
ok "FDR cutoff = $FDR_CUTOFF_VAL"
blank

# ── FC cutoff ────────────────────────────────────────────────
echo -e "  ${BOLD}Fold-change cutoff${RESET} — absolute log2 fold change required for significance."
echo -e "  ${DIM}Applied jointly with FDR: a peak must meet BOTH thresholds.${RESET}"
echo -e "  ${DIM}• 0.585 ≈ 1.5-fold change (recommended for ATAC-seq).${RESET}"
echo -e "  ${DIM}• 1.0   = 2-fold change (more stringent).${RESET}"
echo -e "  ${DIM}• 0     = no fold-change filter (FDR only).${RESET}"
blank

while true; do
    read -p "  |log2FC| cutoff [default: 0.585]: " FC_INPUT
    if [[ -z "$FC_INPUT" ]]; then
        FC_CUTOFF_VAL="0.585"
        break
    elif [[ "$FC_INPUT" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        FC_CUTOFF_VAL="$FC_INPUT"
        break
    else
        err "Please enter a non-negative number (e.g. 0.585)."
    fi
done
ok "FC cutoff = $FC_CUTOFF_VAL (|log2FC| ≥ $FC_CUTOFF_VAL)"
blank

# DiffBind's dba.count() drops any consensus site whose activity (RPKM-based
# by default since DiffBind 3.0) doesn't clear this value on at least one
# sample. Default 1 matches DiffBind's own default and is left untouched for
# normal runs. Low-depth runs (e.g. the smoke-test harness) can hit
# "No sites have activity greater than filter value" purely from shallow
# coverage; PEPATAC_DIFFBIND_FILTER lets an automated caller override this
# without adding an interactive prompt for everyday use.
DIFFBIND_FILTER_VAL="${PEPATAC_DIFFBIND_FILTER:-1}"

# dba.normalize()'s `background` argument controls whether DESeq2-native RLE
# normalization is computed from the consensus peakset alone (FALSE) or from
# larger genome-wide background bins (TRUE). DiffBind's own current guidance
# favors background=TRUE for ATAC-seq-style comparisons: peak-only RLE is
# more vulnerable to composition bias when a large fraction of the genome
# moves in the same direction under treatment, which is common in ATAC
# diff-accessibility (broad opening or closing). TRUE is the default here;
# PEPATAC_DIFFBIND_BACKGROUND=FALSE overrides it for anyone with a specific
# reason to normalize on peaks only (e.g. matching a prior analysis that
# used background=FALSE, or very few peaks are expected to move).
DIFFBIND_BACKGROUND_VAL="${PEPATAC_DIFFBIND_BACKGROUND:-TRUE}"
case "${DIFFBIND_BACKGROUND_VAL,,}" in
    true|t|1|yes)  DIFFBIND_BACKGROUND_VAL="TRUE" ;;
    false|f|0|no)  DIFFBIND_BACKGROUND_VAL="FALSE" ;;
    *) die "PEPATAC_DIFFBIND_BACKGROUND must be TRUE or FALSE (got: $DIFFBIND_BACKGROUND_VAL)" ;;
esac

# dba.count()'s `summits` argument controls whether/how consensus peaks get
# re-centered and trimmed to a fixed width around each peak's point of
# highest read pileup. As of DiffBind 3.0, its OWN default is summits=200
# (a 401bp window) unless summits=FALSE is passed explicitly.
#
# Pinned to 200 here rather than left unset, so this pipeline -- not
# whatever DiffBind's default happens to be on a given install -- owns the
# published value. Leaving this blank would mean a future DiffBind release
# changing its own default silently changes this pipeline's peak geometry
# too, while the provenance table below keeps describing the OLD default
# text. Pinning to today's actual DiffBind default changes nothing about
# current behavior; it just freezes it, the same way the blacklist/Miniconda
# downloads elsewhere in this pipeline are pinned to a known-good snapshot
# rather than "whatever is current."
#
# Some ATAC-seq guidance recommends summits=75 (151bp windows) for
# narrower, less background-diluted regions; PEPATAC_DIFFBIND_SUMMITS=FALSE
# instead preserves each sample's original variable-width, merged peak
# boundaries. Override with either as needed.
DIFFBIND_SUMMITS_VAL="${PEPATAC_DIFFBIND_SUMMITS:-200}"

# GO/KEGG enrichment needs a background/"universe" gene set to test against.
# The preferred universe is genes near this experiment's own accessible
# (consensus) peaks -- testing against the whole genome/database instead
# can inflate significance for terms whose genes were never reachable by
# chromatin accessibility in this cell type to begin with. That fallback
# only happens when fewer than 5 background genes map, which is uncommon
# but not impossible (assembly mismatch, a sparse peak set, etc.).
# Default: skip GO/KEGG entirely rather than silently produce results
# against a scientifically different universe than intended. Override with
# PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE=TRUE to proceed anyway (e.g.
# for a quick look, understanding the caveat) -- every output row is
# tagged with which universe was actually used either way (UniverseType
# column), so results can never be silently ambiguous about this later.
ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL="${PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE:-FALSE}"
case "${ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL,,}" in
    true|t|1|yes)  ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL="TRUE" ;;
    false|f|0|no)  ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL="FALSE" ;;
    *) die "PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE must be TRUE or FALSE (got: $ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL)" ;;
esac

header "Step 7 · CPU Threads"

AVAIL_CORES=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
DEFAULT_DIFF_THREADS=$(( AVAIL_CORES > 4 ? AVAIL_CORES - 2 : AVAIL_CORES ))

echo -e "  Available CPU cores: ${BOLD}${AVAIL_CORES}${RESET}"
echo -e "  ${DIM}Threads are used by DiffBind during read counting and DESeq2 analysis.${RESET}"
echo -e "  ${DIM}Recommended: $DEFAULT_DIFF_THREADS (leaving 2 cores for system)${RESET}"
blank

while true; do
    read -p "  Threads to use [default: ${DEFAULT_DIFF_THREADS}]: " DIFF_THREADS_INPUT
    if [[ -z "$DIFF_THREADS_INPUT" ]]; then
        DIFF_THREADS="$DEFAULT_DIFF_THREADS"
        break
    elif [[ "$DIFF_THREADS_INPUT" =~ ^[0-9]+$ ]] && [[ "$DIFF_THREADS_INPUT" -ge 1 ]]; then
        DIFF_THREADS="$DIFF_THREADS_INPUT"
        break
    else
        err "Please enter a positive integer."
    fi
done

ok "Using $DIFF_THREADS threads for DiffBind."

# ─────────────────────────────────────────────────────────────
# STEP 8 — Output folder
# ─────────────────────────────────────────────────────────────

header "Step 8 · Output Folder"

DEFAULT_DIFF_OUT="$RUN_DIR/diff_analysis_${RUN_ID}"
echo -e "  Where should differential analysis results be written?"
echo -e "  ${DIM}Default: $DEFAULT_DIFF_OUT${RESET}"
blank

IFS= read -r -e -p "  Output folder [press Enter for default]: " DIFF_OUT
normalize_path_var DIFF_OUT
[[ -z "$DIFF_OUT" ]] && DIFF_OUT="$DEFAULT_DIFF_OUT"

ok "Output directory: $DIFF_OUT"
echo -e "  ${DIM}It will be created before analysis begins.${RESET}"

# ─────────────────────────────────────────────────────────────
# STEP 9 — Run summary + confirmation
# ─────────────────────────────────────────────────────────────

header "Step 9 · Run Summary"

echo -e "  ${BOLD}Run folder:${RESET}      $RUN_DIR"
echo -e "  ${BOLD}Output folder:${RESET}   $DIFF_OUT"
echo -e "  ${BOLD}Genome:${RESET}          $GENOME$( $IS_CUSTOM_GENOME && echo " (custom)" || echo "" )"
echo -e "  ${BOLD}Coordinate annotation:${RESET} $( $USE_FROZEN_ANNOTATION && echo "$TXDB_PKG ${TXDB_SOURCE_VERSION:-} (frozen validated TxDb)" || echo "$TXDB_PKG (installed package)" )"
echo -e "  ${BOLD}Gene metadata / GO:${RESET} $ORGDB_PKG${ORGDB_SOURCE_VERSION:+ $ORGDB_SOURCE_VERSION}$( $USE_FROZEN_ANNOTATION && echo " (frozen validated OrgDb)" || echo " (installed package)" )"
echo -e "  ${BOLD}KEGG organism code:${RESET} ${KEGG_ORG:-<none -- KEGG enrichment skipped>}"
echo -e "  ${BOLD}Threads:${RESET}         $DIFF_THREADS"
echo -e "  ${BOLD}Discovered samples:${RESET} $DISCOVERED_SAMPLE_COUNT"
echo -e "  ${BOLD}Included samples:${RESET}   ${#SAMPLE_IDS[@]}"
echo -e "  ${BOLD}Excluded samples:${RESET}   ${#EXCLUDED_SAMPLE_IDS[@]}"
echo -e "  ${BOLD}Groups:${RESET}          ${#UNIQUE_GROUPS[@]}  (${UNIQUE_GROUPS[*]})"
echo -e "  ${BOLD}Contrasts:${RESET}       $CONTRAST_NUM"
echo -e "  ${BOLD}MIN_OVERLAP:${RESET}     $MIN_OVERLAP_VAL"
echo -e "  ${BOLD}FDR cutoff:${RESET}      $FDR_CUTOFF_VAL"
echo -e "  ${BOLD}|log2FC| cutoff:${RESET} $FC_CUTOFF_VAL"
echo -e "  ${BOLD}Normalization:${RESET}   DBA_NORM_RLE, background=$DIFFBIND_BACKGROUND_VAL  (override via PEPATAC_DIFFBIND_BACKGROUND)"
echo -e "  ${BOLD}Peak re-centering:${RESET} ${DIFFBIND_SUMMITS_VAL:-DiffBind default (401bp window)}  (override via PEPATAC_DIFFBIND_SUMMITS)"
echo -e "  ${BOLD}Enrichment universe fallback:${RESET} ${ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL}  (if TRUE, proceeds on the OrgDb default when <5 background genes map, instead of skipping GO/KEGG; override via PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE)"
echo -e "  ${DIM}Note: each contrast is normalized/tested independently (no other loaded${RESET}"
echo -e "  ${DIM}group influences its statistics), but every contrast shares one consensus${RESET}"
echo -e "  ${DIM}peak set built once from all included samples -- a deliberate tradeoff${RESET}"
echo -e "  ${DIM}for consistent, comparable peak sets across comparisons, worth stating in${RESET}"
echo -e "  ${DIM}any write-up of these results.${RESET}"
blank

for i in "${!CONTRAST_LABELS[@]}"; do
    echo -e "    ${CYAN}$((i+1)).${RESET}  ${CONTRAST_LABELS[$i]}  →  ${CONTRAST_CASES[$i]} vs ${CONTRAST_CONTROLS[$i]}"
done
blank

echo -e "  ${BOLD}Samples included in DiffBind:${RESET}"
for sid in "${SAMPLE_IDS[@]}"; do
    echo -e "    ${GREEN}✔${RESET}  $sid"
done
if (( ${#EXCLUDED_SAMPLE_IDS[@]} > 0 )); then
    blank
    echo -e "  ${BOLD}Samples excluded before DiffBind:${RESET}"
    for i in "${!EXCLUDED_SAMPLE_IDS[@]}"; do
        echo -e "    ${RED}✘${RESET}  ${EXCLUDED_SAMPLE_IDS[$i]}"
        echo -e "       ${DIM}${EXCLUDED_SAMPLE_REASONS[$i]}${RESET}"
    done
fi
blank

echo -e "  ${DIM}Per contrast, the runner will produce:${RESET}"
echo -e "  ${DIM}  • DiffBind consensus peak count matrix${RESET}"
echo -e "  ${DIM}  • DESeq2 differential accessibility results (CSV, all peaks + significant)${RESET}"
echo -e "  ${DIM}  • PCA coordinates across all samples (TSV)${RESET}"
echo -e "  ${DIM}  • ChIPseeker annotated peaks + feature/TSS/TES/width stats (TSV)${RESET}"
echo -e "  ${DIM}  • GO/KEGG enrichment results (CSV: full tested set + significant-only subset)${RESET}"
echo -e "  ${DIM}  • HOMER motifs + GO/KEGG run 3x per contrast: combined, Up-only, Down-only${RESET}"
echo -e "  ${DIM}  • Cross-contrast summary table${RESET}"
echo -e "  ${DIM}This script computes only — no plots are generated.${RESET}"
echo -e "  ${DIM}Run PEPATAC_explore.sh afterward for PCA, volcano/MA, annotation,${RESET}"
echo -e "  ${DIM}GO/KEGG, motif, and tornado plots, regenerated on demand.${RESET}"
blank

read -p "  Start differential analysis? [Y/n]: " FINAL_CONFIRM
[[ "${FINAL_CONFIRM,,}" == "n" ]] && die "Aborted by user."

# ─────────────────────────────────────────────────────────────
# STEP 10 — Create output directory + write DiffBind sample sheet
# ─────────────────────────────────────────────────────────────

header "Step 10 · Preparing Output"

# print_peak_format_alert
# Prints a maximally visible banner naming every sample whose peak file
# wasn't real narrowPeak output (PeakCaller=bed). Called once here, as
# part of the decision to proceed, and again in the final run summary so
# it can't be missed by anyone who only reads the end of the log. Both
# calls read from the same BED_FALLBACK_* arrays, so the warning can't
# drift between the two places it's shown.
print_peak_format_alert() {
    local border="════════════════════════════════════════════════════════"
    blank
    echo -e "${BOLD}${RED}${border}${RESET}"
    echo -e "${BOLD}${RED}  ⚠  PEAK FORMAT WARNING -- NOT ALL SAMPLES USED narrowPeak  ⚠${RESET}"
    echo -e "${BOLD}${RED}${border}${RESET}"
    blank
    echo -e "  ${BOLD}${#BED_FALLBACK_SAMPLES[@]} of ${#SAMPLE_IDS[@]}${RESET} sample(s) did not have a real narrowPeak file."
    echo -e "  PEPATAC's peak calling for these fell back to a summits/generic BED file"
    echo -e "  instead -- this usually means peak calling produced no peaks, or the"
    echo -e "  narrowPeak output was missing or corrupted for that sample."
    blank
    for i in "${!BED_FALLBACK_SAMPLES[@]}"; do
        echo -e "    ${RED}✘${RESET}  ${BOLD}${BED_FALLBACK_SAMPLES[$i]}${RESET}"
        echo -e "       ${DIM}${BED_FALLBACK_PEAKS[$i]}${RESET}"
    done
    blank
    echo -e "  These are correctly being sent to DiffBind as PeakCaller=bed, not narrow,"
    echo -e "  so the file's columns won't be misread as narrowPeak's p-value/score/summit"
    echo -e "  fields. But a BED fallback means weaker or different peak evidence for these"
    echo -e "  samples than the rest of the cohort -- something upstream already didn't go"
    echo -e "  as expected for them. Check their PEPATAC peak-calling logs before trusting"
    echo -e "  results that depend heavily on these samples."
    echo -e "${BOLD}${RED}${border}${RESET}"
    blank
}

mkdir -p "$DIFF_OUT"
DIFFBIND_SHEET="$DIFF_OUT/diffbind_samplesheet.csv"
INCLUDED_SAMPLES_FILE="$DIFF_OUT/included_samples.tsv"
EXCLUDED_SAMPLES_FILE="$DIFF_OUT/excluded_samples.tsv"
DIFF_LOG="$DIFF_OUT/diff_analysis_${RUN_ID}.log"
DIFF_R_SCRIPT="$DIFF_OUT/run_diffbind.R"

# Write explicit sample-selection records before creating the DiffBind sheet.
{
    printf 'SampleID\tBAM\tPeaks\tPeakCaller\tStatus\n'
    for i in "${!SAMPLE_IDS[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "${SAMPLE_IDS[$i]}" "${SAMPLE_BAMS[$i]}" "${SAMPLE_PEAKS[$i]}" \
            "${SAMPLE_PEAK_CALLERS[$i]:-narrow}" "${SAMPLE_STATUS_ARR[$i]:-UNKNOWN}"
    done
} > "$INCLUDED_SAMPLES_FILE"

{
    printf 'SampleID\tReason\n'
    for i in "${!EXCLUDED_SAMPLE_IDS[@]}"; do
        printf '%s\t%s\n' "${EXCLUDED_SAMPLE_IDS[$i]}" "${EXCLUDED_SAMPLE_REASONS[$i]}"
    done
} > "$EXCLUDED_SAMPLES_FILE"

# Which samples are about to be told to DiffBind as PeakCaller=bed instead
# of narrow. Computed from the final, already-filtered SAMPLE_IDS so this
# only reflects what's actually going into the analysis.
BED_FALLBACK_SAMPLES=()
BED_FALLBACK_PEAKS=()
for i in "${!SAMPLE_IDS[@]}"; do
    if [[ "${SAMPLE_PEAK_CALLERS[$i]:-narrow}" == "bed" ]]; then
        BED_FALLBACK_SAMPLES+=("${SAMPLE_IDS[$i]}")
        BED_FALLBACK_PEAKS+=("${SAMPLE_PEAKS[$i]}")
    fi
done

if [[ ${#BED_FALLBACK_SAMPLES[@]} -gt 0 ]]; then
    print_peak_format_alert
    read -p "  Proceed with these sample(s) included as-is? [y/N]: " BED_FALLBACK_CONFIRM
    if [[ "${BED_FALLBACK_CONFIRM,,}" != "y" ]]; then
        die "Aborted. Investigate the sample(s) listed above, or exclude them in Step 2b, then re-run."
    fi
    warn "Continuing with ${#BED_FALLBACK_SAMPLES[@]} bed-fallback sample(s) -- this will be shown again in the final summary."
fi

# Write DiffBind sample sheet from the already-filtered arrays only.
{
    echo "SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller"
    declare -A GROUP_REPLICATE_COUNT
    for i in "${!SAMPLE_IDS[@]}"; do
        grp="${SAMPLE_GROUPS[$i]}"
        GROUP_REPLICATE_COUNT["$grp"]=$(( ${GROUP_REPLICATE_COUNT["$grp"]:-0} + 1 ))
        rep="${GROUP_REPLICATE_COUNT[$grp]}"
        bam="${SAMPLE_BAMS[$i]}"
        peak="${SAMPLE_PEAKS[$i]}"
        caller="${SAMPLE_PEAK_CALLERS[$i]:-narrow}"
        echo "\"$(csv_quote "${SAMPLE_IDS[$i]}")\",\"$(csv_quote "$grp")\",${rep},\"$(csv_quote "$bam")\",\"$(csv_quote "$peak")\",\"$(csv_quote "$caller")\""
    done
} > "$DIFFBIND_SHEET"

EXPECTED_SHEET_ROWS=${#SAMPLE_IDS[@]}
ACTUAL_SHEET_ROWS=$(( $(wc -l < "$DIFFBIND_SHEET") - 1 ))
if [[ "$ACTUAL_SHEET_ROWS" -ne "$EXPECTED_SHEET_ROWS" ]]; then
    die "DiffBind sample-sheet verification failed: expected $EXPECTED_SHEET_ROWS included samples, found $ACTUAL_SHEET_ROWS rows."
fi

ok "Included sample record: $INCLUDED_SAMPLES_FILE"
ok "Excluded sample record: $EXCLUDED_SAMPLES_FILE"
ok "DiffBind sample sheet: $DIFFBIND_SHEET"
ok "Verified: DiffBind sheet contains exactly $ACTUAL_SHEET_ROWS included sample(s)."

# ─────────────────────────────────────────────────────────────
# STEP 11 — Check / install required R packages
# ─────────────────────────────────────────────────────────────

header "Step 11 · R Package Check"

# Default: a missing package stops the run rather than silently installing
# whatever CRAN/Bioconductor currently ships. The installer is meant to be
# the authority on the tested environment (frozen annotation, pinned CLI
# tools, etc.) -- quietly patching the R side at analysis time would let a
# published run use different package versions than what was validated,
# with no record of it happening. Set PEPATAC_DIFF_AUTOINSTALL=TRUE to
# restore the old convenience behavior for exploratory/iterative use.
DIFF_AUTOINSTALL_VAL="${PEPATAC_DIFF_AUTOINSTALL:-FALSE}"
case "${DIFF_AUTOINSTALL_VAL,,}" in
    true|t|1|yes)  DIFF_AUTOINSTALL_VAL="TRUE" ;;
    false|f|0|no)  DIFF_AUTOINSTALL_VAL="FALSE" ;;
    *) die "PEPATAC_DIFF_AUTOINSTALL must be TRUE or FALSE (got: $DIFF_AUTOINSTALL_VAL)" ;;
esac

in_env_clean Rscript --vanilla - "$DIFF_OUT" "$TXDB_INSTALL_PKG" "$ORGDB_INSTALL_PKG" "$DIFF_AUTOINSTALL_VAL" << 'RCHECK'
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
selected_txdb <- args[2]
selected_orgdb <- args[3]
autoinstall <- toupper(args[4]) == "TRUE"

# Full dependency list for DiffBind + DESeq2 + RLE normalization.
# edgeR and limma are required by DiffBind's DBA_NORM_RLE normalization.
# The Bioconductor packages must be installed via BiocManager, not CRAN.
#
# selected_txdb/selected_orgdb are populated only for the five built-in
# profiles. User-defined profiles load immutable TxDb/OrgDb SQLite snapshots,
# so their source packages do not need to remain installed for later analysis.
# Filter(nzchar, ...) keeps the package check valid in both cases.
bioc_pkgs <- c(
    "DiffBind",           # core differential binding
    "DESeq2",             # statistical engine
    "edgeR",              # required by DiffBind RLE normalization
    "limma",              # required by DiffBind
    "BiocParallel",       # parallel processing
    "GenomicRanges",      # genomic interval handling
    "Rsamtools",          # BAM file reading
    "SummarizedExperiment",
    "AnnotationDbi",
    "GenomicFeatures",
    "GenomeInfoDb",
    "IRanges",
    "ChIPseeker",         # assembly-aware peak annotation
    "BiocManager",
    "clusterProfiler",    # GO and KEGG enrichment analysis
    "enrichplot",         # clusterProfiler visualizations
    "GO.db",
    Filter(nzchar, c(selected_txdb, selected_orgdb))
)
cran_pkgs <- c(
    "ggplot2",
    "ggrepel",
    "dplyr",
    "tidyr",
    "RSQLite"
)
all_needed <- c(bioc_pkgs, cran_pkgs)

missing <- all_needed[!vapply(all_needed, requireNamespace, logical(1), quietly=TRUE)]

if (length(missing) > 0) {
    if (!autoinstall) {
        stop("Required R package(s) missing: ", paste(missing, collapse=", "),
             "\nRe-run PEPATAC_install.sh to bring the environment back to its tested state.",
             "\n(Or set PEPATAC_DIFF_AUTOINSTALL=TRUE to install missing packages from",
             " current CRAN/Bioconductor instead -- convenient, but the resulting package",
             " versions are then unvalidated and won't match a pinned/tested environment.)")
    }
    cat("  Installing missing R packages:", paste(missing, collapse=", "), "\n")
    options(repos=c(CRAN="https://cloud.r-project.org"))
    cran_miss <- intersect(missing, cran_pkgs)
    bioc_miss  <- intersect(missing, bioc_pkgs)
    if (length(cran_miss) > 0) {
        cat("  Installing from CRAN:", paste(cran_miss, collapse=", "), "\n")
        install.packages(cran_miss, quiet=TRUE)
    }
    if (length(bioc_miss) > 0) {
        cat("  Installing from Bioconductor:", paste(bioc_miss, collapse=", "), "\n")
        BiocManager::install(bioc_miss, ask=FALSE, update=FALSE)
    }
    # Verify everything installed successfully.
    still_missing <- all_needed[!vapply(all_needed, requireNamespace, logical(1), quietly=TRUE)]
    if (length(still_missing) > 0) {
        stop("Failed to install: ", paste(still_missing, collapse=", "),
             "\nTry running install.sh again or install manually.")
    }
    cat("  All packages installed successfully.\n")
} else {
    cat("  All required R packages are present.\n")
}

# Print versions of installed analysis packages. Source TxDb/OrgDb package
# versions for user-defined profiles come from the immutable snapshot and are
# printed by the generated analysis script.
cat("\n  Key package versions:\n")
version_pkgs <- c("DiffBind", "DESeq2", "ChIPseeker", "GenomicFeatures",
                   "clusterProfiler", "enrichplot",
                   Filter(nzchar, c(selected_txdb, selected_orgdb)))
for (pkg in version_pkgs) {
    tryCatch(
        cat(sprintf("    %-22s %s\n", pkg,
                    as.character(packageVersion(pkg)))),
        error = function(e) cat(sprintf("    %-22s (version unavailable)\n", pkg))
    )
}
RCHECK

ok "R packages ready."

# ─────────────────────────────────────────────────────────────
# STEP 12 — Write the R analysis script
# ─────────────────────────────────────────────────────────────

header "Step 12 · Writing R Analysis Script"

# Serialize contrast arrays, sample metadata, and free-text paths into a
# format R can read SAFELY. Every value below either came from interactive
# free-text input (group/contrast labels) or is a filesystem path (BAMs,
# output dir) -- naive `printf '"%s",'` embedding of these directly into R
# source text is an injection boundary (a value containing `"` or `\` can
# break out of the R string literal). r_vector_literal/r_string_literal
# percent-encode every value so the embedded text can only ever contain
# [A-Za-z0-9._~-] and %XX escapes -- nothing that can be interpreted as R
# syntax -- then the R side decodes with utils::URLdecode() before use.
# See the definitions of these two functions near the top of this script.
CONTRASTS_CASES_R="$(r_vector_literal "${CONTRAST_CASES[@]}")"
CONTRASTS_CTRLS_R="$(r_vector_literal "${CONTRAST_CONTROLS[@]}")"
CONTRASTS_LABELS_R="$(r_vector_literal "${CONTRAST_LABELS[@]}")"

# Serialize sample IDs, BAM paths, and group assignments for the R script.
# These are stored in the explorer bundle so tornado plots can find BAMs
# bundle so the tornado explorer menu can find BAMs without re-running.
SAMPLE_IDS_R="$(r_vector_literal "${SAMPLE_IDS[@]}")"
SAMPLE_BAMS_R="$(r_vector_literal "${SAMPLE_BAMS[@]}")"
SAMPLE_GROUPS_R="$(r_vector_literal "${SAMPLE_GROUPS[@]}")"

# DIFF_OUT and DIFFBIND_SHEET are free-text (user-typed / normalize_path'd)
# filesystem paths embedded directly into R source -- same injection
# boundary as above, same fix.
DIFF_OUT_R="$(r_string_literal "$DIFF_OUT")"
DIFFBIND_SHEET_R="$(r_string_literal "$DIFFBIND_SHEET")"

# User-defined profile metadata and immutable SQLite paths are user-provided
# or filesystem-derived values, so encode them before embedding in R source.
GENOME_R="$(r_string_literal "$GENOME")"
TXDB_PKG_R="$(r_string_literal "$TXDB_PKG")"
ORGDB_PKG_R="$(r_string_literal "$ORGDB_PKG")"
KEGG_ORG_R="$(r_string_literal "$KEGG_ORG")"
PROFILE_ASSEMBLY_R="$(r_string_literal "${CUSTOM_PROFILE_ASSEMBLY:-$GENOME}")"
TXDB_SOURCE_VERSION_R="$(r_string_literal "$TXDB_SOURCE_VERSION")"
ORGDB_SOURCE_VERSION_R="$(r_string_literal "$ORGDB_SOURCE_VERSION")"
# CUSTOM_TXDB_SQLITE/CUSTOM_ORGDB_SQLITE (R names, unchanged, to minimize
# churn in the R body below) are fed from the unified FROZEN_* bash
# variables now, not the custom-only *_PATH ones -- this is what lets a
# built-in genome with its own frozen snapshot reach the same
# AnnotationDbi::loadDb() path a custom profile already used. See
# USE_FROZEN_ANNOTATION_R just below for the branch condition itself.
CUSTOM_TXDB_SQLITE_R="$(r_string_literal "$FROZEN_TXDB_SQLITE")"
CUSTOM_ORGDB_SQLITE_R="$(r_string_literal "$FROZEN_ORGDB_SQLITE")"
CUSTOM_CHROM_SIZES_R="$(r_string_literal "$FROZEN_CHROM_SIZES")"
CUSTOM_FASTA_SHA256_R="$(r_string_literal "$CUSTOM_FASTA_SHA256")"
CUSTOM_CHROM_SIZES_SHA256_R="$(r_string_literal "$FROZEN_CHROM_SIZES_SHA256")"
CUSTOM_TXDB_SQLITE_SHA256_R="$(r_string_literal "$FROZEN_TXDB_SQLITE_SHA256")"
CUSTOM_ORGDB_SQLITE_SHA256_R="$(r_string_literal "$FROZEN_ORGDB_SQLITE_SHA256")"
USE_FROZEN_ANNOTATION_R=$( $USE_FROZEN_ANNOTATION && echo "TRUE" || echo "FALSE" )
HOMER_CONFIGURE_R="$(r_string_literal "$HOMER_CONFIGURE")"
HOMER_FIND_MOTIFS_R="$(r_string_literal "$HOMER_FIND_MOTIFS")"
HOMER_GENOME_ARG_R="$(r_string_literal "$HOMER_GENOME_ARG")"
HOMER_VERSION_R="$(r_string_literal "$HOMER_VERSION")"
HOMER_GENOME_INFO_R="$(r_string_literal "$HOMER_GENOME_INFO")"

cat > "$DIFF_R_SCRIPT" << RSCRIPT_EOF
# ============================================================
# PEPATAC Differential Accessibility Analysis
# Generated by diff_analysis.sh ${SCRIPT_VERSION}
# Run ID: ${RUN_ID}
# Started: ${RUN_STARTED}
# ============================================================

suppressPackageStartupMessages({
    library(DiffBind)
    library(DESeq2)
    library(dplyr)
    library(tidyr)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(AnnotationDbi)
})

# ── Safe-embedding decode helpers ─────────────────────────────
# Counterparts to Bash's url_encode/r_vector_literal/r_string_literal.
# base R only (utils::URLdecode) -- no extra package dependency.
.pepatac_url_decode <- function(x) utils::URLdecode(x)
.pepatac_decode_vec <- function(x) vapply(x, .pepatac_url_decode, character(1), USE.NAMES = FALSE)

# ── Parameters ──────────────────────────────────────────────
DIFF_OUT        <- ${DIFF_OUT_R}
DIFFBIND_SHEET  <- ${DIFFBIND_SHEET_R}
GENOME          <- ${GENOME_R}
TXDB_PKG        <- ${TXDB_PKG_R}
ORGDB           <- ${ORGDB_PKG_R}
KEGG_ORG        <- ${KEGG_ORG_R}
PROFILE_ASSEMBLY <- ${PROFILE_ASSEMBLY_R}
TXDB_SOURCE_VERSION <- ${TXDB_SOURCE_VERSION_R}
ORGDB_SOURCE_VERSION <- ${ORGDB_SOURCE_VERSION_R}
IS_CUSTOM_GENOME <- $( $IS_CUSTOM_GENOME && echo "TRUE" || echo "FALSE" )
# TRUE whenever a frozen SQLite is actually available to load -- true for
# every custom profile, and for a built-in genome only when this run has
# its own validated annotation snapshot (see USE_FROZEN_ANNOTATION_R in
# bash above). IS_CUSTOM_GENOME stays reserved for things that are
# genuinely about the custom-profile system specifically (KEGG defaults,
# the profile's own FASTA/chrom.sizes provenance, the profile-vs-TxDb
# seqlevels cross-check below) rather than "is a frozen SQLite in use."
USE_FROZEN_ANNOTATION <- ${USE_FROZEN_ANNOTATION_R}
CUSTOM_TXDB_SQLITE <- ${CUSTOM_TXDB_SQLITE_R}
CUSTOM_ORGDB_SQLITE <- ${CUSTOM_ORGDB_SQLITE_R}
CUSTOM_CHROM_SIZES <- ${CUSTOM_CHROM_SIZES_R}
CUSTOM_FASTA_SHA256 <- ${CUSTOM_FASTA_SHA256_R}
CUSTOM_CHROM_SIZES_SHA256 <- ${CUSTOM_CHROM_SIZES_SHA256_R}
CUSTOM_TXDB_SQLITE_SHA256 <- ${CUSTOM_TXDB_SQLITE_SHA256_R}
CUSTOM_ORGDB_SQLITE_SHA256 <- ${CUSTOM_ORGDB_SQLITE_SHA256_R}
HAS_ORGDB       <- nzchar(ORGDB)
HAS_KEGG        <- nzchar(KEGG_ORG)
CONTRAST_CASES  <- ${CONTRASTS_CASES_R}
CONTRAST_CTRLS  <- ${CONTRASTS_CTRLS_R}
CONTRAST_LABELS <- ${CONTRASTS_LABELS_R}
DIFF_THREADS    <- ${DIFF_THREADS}    # set interactively; see Step 6 note in diff_analysis.sh
EXPECTED_SAMPLE_COUNT <- ${#SAMPLE_IDS[@]}

FDR_CUTOFF  <- ${FDR_CUTOFF_VAL}   # set interactively: adjusted p-value threshold
FC_CUTOFF   <- ${FC_CUTOFF_VAL}   # set interactively: |log2FC| threshold (applied jointly with FDR)
MIN_OVERLAP <- ${MIN_OVERLAP_VAL}  # set interactively: min samples a peak must appear in for consensus
DIFFBIND_FILTER <- ${DIFFBIND_FILTER_VAL}  # min site activity for dba.count(); override via PEPATAC_DIFFBIND_FILTER
DIFFBIND_BACKGROUND <- ${DIFFBIND_BACKGROUND_VAL}  # dba.normalize() background=; override via PEPATAC_DIFFBIND_BACKGROUND
DIFFBIND_SUMMITS <- "${DIFFBIND_SUMMITS_VAL}"  # dba.count() summits=; "" = DiffBind's own default; override via PEPATAC_DIFFBIND_SUMMITS
ALLOW_DEFAULT_ENRICHMENT_UNIVERSE <- ${ALLOW_DEFAULT_ENRICHMENT_UNIVERSE_VAL}  # override via PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE

RUN_MOTIF_ANALYSIS <- $( [[ "$RUN_MOTIF_ANALYSIS" == "true" ]] && echo "TRUE" || echo "FALSE" )
MOTIF_MODE         <- "${MOTIF_MODE}"          # "known", "denovo", or "both" (set interactively)
HOMER_CONFIGURE    <- ${HOMER_CONFIGURE_R}     # path to configureHomer.pl in conda env
HOMER_FIND_MOTIFS  <- ${HOMER_FIND_MOTIFS_R}   # path to findMotifsGenome.pl in conda env
HOMER_GENOME_ARG   <- ${HOMER_GENOME_ARG_R}    # installed HOMER genome name, or a raw FASTA path for a custom genome
HOMER_VERSION      <- ${HOMER_VERSION_R}       # conda-installed HOMER package version
HOMER_GENOME_INFO  <- ${HOMER_GENOME_INFO_R}   # HOMER genome package identity/path used for motif finding

# Sample metadata — stored in the
# explorer bundle so tornado plots can find BAMs interactively.
SAMPLE_IDS    <- ${SAMPLE_IDS_R}
SAMPLE_BAMS   <- ${SAMPLE_BAMS_R}
SAMPLE_GROUPS <- ${SAMPLE_GROUPS_R}
# Named BAM vector: names = SampleID, values = BAM path.
names(SAMPLE_BAMS)   <- SAMPLE_IDS
names(SAMPLE_GROUPS) <- SAMPLE_IDS

# ── Coordinate and gene annotation objects ───────────────────
# Frozen SQLite (custom profile, or a built-in genome with its own run
# snapshot) loads via AnnotationDbi::loadDb(). Otherwise -- a built-in
# genome with no snapshot for this run -- loads the currently-installed
# package, same as before this fix existed.
if (USE_FROZEN_ANNOTATION) {
    cat(sprintf("  Loading frozen validated TxDb: %s\n", CUSTOM_TXDB_SQLITE))
    TXDB <- suppressWarnings(suppressMessages(AnnotationDbi::loadDb(CUSTOM_TXDB_SQLITE)))
    if (!inherits(TXDB, "TxDb")) stop("Frozen TxDb SQLite did not load as a TxDb object.")
    ORGDB_OBJECT <- suppressWarnings(suppressMessages(AnnotationDbi::loadDb(CUSTOM_ORGDB_SQLITE)))
    if (!inherits(ORGDB_OBJECT, "OrgDb")) stop("Frozen OrgDb SQLite did not load as an OrgDb object.")
    TXDB_PKG_LABEL <- sprintf("%s %s (frozen SQLite)", TXDB_PKG, TXDB_SOURCE_VERSION)
    ORGDB_LABEL <- sprintf("%s %s (frozen SQLite)", ORGDB, ORGDB_SOURCE_VERSION)
} else {
    suppressPackageStartupMessages({
        library(TXDB_PKG, character.only=TRUE)
        library(ORGDB, character.only=TRUE)
    })
    TXDB <- getExportedValue(TXDB_PKG, TXDB_PKG)
    ORGDB_OBJECT <- getExportedValue(ORGDB, ORGDB)
    TXDB_PKG_LABEL <- TXDB_PKG
    ORGDB_LABEL <- ORGDB
}
if (!HAS_KEGG) {
    cat("  No KEGG organism code registered for this genome -- KEGG enrichment will be skipped.\n")
}

# ── Helper: safe directory creation ─────────────────────────
make_dir <- function(d) { dir.create(d, showWarnings=FALSE, recursive=TRUE); d }

# ── Helper: write TSV diagnostics ───────────────────────────
write_tsv <- function(x, path) {
    write.table(x, file=path, sep="\t", quote=FALSE, row.names=FALSE, na="")
}

# ── Annotation provenance and TxDb ↔ OrgDb compatibility ─────
annotation_provenance_dir <- make_dir(file.path(DIFF_OUT, "diagnostics", "annotation_provenance"))

TXDB_METADATA <- tryCatch(AnnotationDbi::metadata(TXDB), error=function(e) data.frame())
if (nrow(TXDB_METADATA) > 0) {
    write_tsv(TXDB_METADATA, file.path(annotation_provenance_dir, "txdb_metadata.tsv"))
}

txdb_genome_rows <- which(tolower(as.character(TXDB_METADATA\$name)) == "genome")
TXDB_REPORTED_GENOME <- if (length(txdb_genome_rows) > 0) {
    as.character(TXDB_METADATA\$value[txdb_genome_rows[1]])
} else {
    NA_character_
}
if (USE_FROZEN_ANNOTATION && (is.na(TXDB_REPORTED_GENOME) || !nzchar(TXDB_REPORTED_GENOME))) {
    stop("Frozen TxDb does not report its UCSC assembly in metadata.")
}
if (!is.na(TXDB_REPORTED_GENOME) && nzchar(TXDB_REPORTED_GENOME) &&
    !identical(tolower(TXDB_REPORTED_GENOME), tolower(PROFILE_ASSEMBLY))) {
    stop(sprintf("TxDb assembly mismatch: profile assembly is '%s' but %s reports '%s'.",
                 PROFILE_ASSEMBLY, TXDB_PKG, TXDB_REPORTED_GENOME))
}
# CUSTOM_CHROM_SIZES now holds the frozen chrom.sizes for either source
# (custom profile, or a built-in genome using a run snapshot -- see
# FROZEN_CHROM_SIZES in bash above), so this independently repeats the
# exact seqlevel/length check for both, not just custom. It's still
# skipped for a built-in genome with no snapshot at all (tier 3, legacy
# install-based loading) or a snapshot written before this field existed
# -- CUSTOM_CHROM_SIZES is simply blank in either case, same as before.
# run.sh's build_profile_annotation_assets()/validate_frozen_annotation_
# profile() already performed the equivalent check once at build time;
# this is an independent re-check at the point of use.
if (nzchar(CUSTOM_CHROM_SIZES)) {
    chrom <- read.delim(CUSTOM_CHROM_SIZES, header=FALSE, stringsAsFactors=FALSE,
                        col.names=c("seqname", "length"))
    if (!nrow(chrom) || anyDuplicated(chrom\$seqname)) {
        stop("Snapshot chrom.sizes is empty or contains duplicate sequence names.")
    }
    si <- GenomeInfoDb::seqinfo(TXDB)
    txseq <- as.character(GenomeInfoDb::seqlevels(si))
    txlen <- GenomeInfoDb::seqlengths(si)
    # AnnotationDbi::loadDb() (used above to load this frozen TxDb) always
    # returns the FULL original sequence set -- any active-seqlevels
    # restriction applied when the TxDb was originally built (e.g. to
    # exclude alt/patch contigs not present in the alignment reference) is
    # a session-level view filter, not something saveDb()/loadDb() persist
    # into the SQLite file. So a reloaded TxDb reporting sequences beyond
    # chrom.sizes is expected and normal here, not a sign of corruption --
    # only a length disagreement on a sequence present in BOTH is.
    overlap <- intersect(txseq, chrom\$seqname)
    if (!length(overlap)) {
        stop("Frozen TxDb and snapshot chrom.sizes share NO sequence names at all. This is a genuine mismatch, not just extra alt/patch contigs.")
    }
    idx <- match(overlap, chrom\$seqname)
    tidx <- match(overlap, txseq)
    if (any(is.na(txlen[tidx]))) {
        stop("Frozen TxDb is missing one or more chromosome lengths for sequences shared with chrom.sizes; exact assembly validation is not possible.")
    }
    if (any(txlen[tidx] != chrom\$length[idx])) {
        bad <- which(txlen[tidx] != chrom\$length[idx])[1]
        stop(sprintf("Frozen TxDb length mismatch for %s: TxDb=%s, chrom.sizes=%s.",
                     overlap[bad], txlen[tidx][bad], chrom\$length[idx[bad]]))
    }
}

normalize_annotation_ids <- function(ids, keytype) {
    ids <- as.character(ids)
    if (keytype %in% c("ENSEMBL", "REFSEQ")) {
        ids <- sub("\\\\.[0-9]+$", "", ids)
    }
    ids
}

choose_orgdb_keytype <- function(txdb, orgdb, txdb_pkg, orgdb_label, diagnostics_dir) {
    tx_ids <- unique(as.character(AnnotationDbi::keys(txdb, keytype="GENEID")))
    tx_ids <- tx_ids[!is.na(tx_ids) & nzchar(tx_ids)]
    if (length(tx_ids) == 0) stop("The selected TxDb contains no GENEID keys.")

    available <- AnnotationDbi::keytypes(orgdb)
    preferred <- if (grepl("knownGene", txdb_pkg, fixed=TRUE)) {
        c("ENTREZID", "ENSEMBL", "REFSEQ", "SYMBOL", "FLYBASE", "ZFIN")
    } else if (grepl("ensGene", txdb_pkg, fixed=TRUE)) {
        c("ENSEMBL", "FLYBASE", "ENTREZID", "SYMBOL", "REFSEQ", "ZFIN")
    } else {
        c("ENTREZID", "REFSEQ", "ENSEMBL", "SYMBOL", "FLYBASE", "ZFIN")
    }
    candidates <- intersect(preferred, available)
    if (length(candidates) == 0) {
        stop("No compatible identifier keytypes are shared with ", orgdb_label, ".")
    }

    rows <- lapply(candidates, function(kt) {
        org_keys <- tryCatch(as.character(AnnotationDbi::keys(orgdb, keytype=kt)),
                             error=function(e) character(0))
        normalized <- normalize_annotation_ids(tx_ids, kt)
        data.frame(
            Keytype=kt,
            TxDb_gene_IDs=length(tx_ids),
            Exact_matches=sum(tx_ids %in% org_keys),
            Normalized_matches=sum(normalized %in% org_keys),
            Mapping_rate=max(mean(tx_ids %in% org_keys), mean(normalized %in% org_keys)),
            stringsAsFactors=FALSE
        )
    })
    diag <- do.call(rbind, rows)
    diag <- diag[order(-diag\$Mapping_rate, match(diag\$Keytype, preferred)), , drop=FALSE]
    write_tsv(diag, file.path(diagnostics_dir, "txdb_orgdb_keytype_compatibility.tsv"))

    best <- diag[1, , drop=FALSE]
    cat(sprintf("  Annotation ID join: %s GENEID → %s %s (%.1f%% overlap)\n",
                txdb_pkg, orgdb_label, best\$Keytype, 100 * best\$Mapping_rate))
    if (best\$Mapping_rate < 0.70) {
        stop(sprintf(
            paste0("TxDb/OrgDb compatibility is too low (%.1f%%; minimum 70%%). ",
                   "This usually means the organism or identifier type is wrong. See %s."),
            100 * best\$Mapping_rate,
            file.path(diagnostics_dir, "txdb_orgdb_keytype_compatibility.tsv")
        ))
    }
    if (best\$Mapping_rate < 0.90) {
        warning(sprintf("TxDb/OrgDb compatibility is %.1f%%; annotation will continue, but unmapped genes are reported.",
                        100 * best\$Mapping_rate))
    }
    as.character(best\$Keytype)
}

TXDB_GENE_KEYTYPE <- if (HAS_ORGDB) {
    choose_orgdb_keytype(TXDB, ORGDB_OBJECT, TXDB_PKG_LABEL, ORGDB_LABEL, annotation_provenance_dir)
} else {
    NA_character_
}

safe_pkg_version <- function(pkg) {
    if (!nzchar(pkg)) return(NA_character_)
    tryCatch(as.character(packageVersion(pkg)), error=function(e) NA_character_)
}
TXDB_RECORDED_VERSION <- if (USE_FROZEN_ANNOTATION) TXDB_SOURCE_VERSION else safe_pkg_version(TXDB_PKG)
ORGDB_RECORDED_VERSION <- if (USE_FROZEN_ANNOTATION) ORGDB_SOURCE_VERSION else safe_pkg_version(ORGDB)

annotation_source <- data.frame(
    Genome=GENOME,
    UCSC_assembly=PROFILE_ASSEMBLY,
    TxDb_package=TXDB_PKG,
    TxDb_version=TXDB_RECORDED_VERSION,
    TxDb_storage=if (USE_FROZEN_ANNOTATION) CUSTOM_TXDB_SQLITE else "installed package",
    TxDb_reported_genome=TXDB_REPORTED_GENOME,
    OrgDb_package=ORGDB,
    OrgDb_version=ORGDB_RECORDED_VERSION,
    OrgDb_storage=if (USE_FROZEN_ANNOTATION) CUSTOM_ORGDB_SQLITE else "installed package",
    # Reference_FASTA_SHA256/Chrom_sizes_SHA256 stay IS_CUSTOM_GENOME-only,
    # not USE_FROZEN_ANNOTATION: these describe the GENOME's own sequence
    # assets, which for a built-in genome (even one using a frozen
    # annotation snapshot) live in Refgenie or the local-build cache --
    # untouched by, and out of scope for, this annotation-parity fix.
    Reference_FASTA_SHA256=if (IS_CUSTOM_GENOME) CUSTOM_FASTA_SHA256 else NA_character_,
    Chrom_sizes_SHA256=if (USE_FROZEN_ANNOTATION) CUSTOM_CHROM_SIZES_SHA256 else NA_character_,
    TxDb_SQLite_SHA256=if (USE_FROZEN_ANNOTATION) CUSTOM_TXDB_SQLITE_SHA256 else NA_character_,
    OrgDb_SQLite_SHA256=if (USE_FROZEN_ANNOTATION) CUSTOM_ORGDB_SQLITE_SHA256 else NA_character_,
    TxDb_gene_keytype=TXDB_GENE_KEYTYPE,
    Coordinate_authority=if (USE_FROZEN_ANNOTATION) "Frozen validated TxDb SQLite" else "Installed TxDb package",
    Gene_metadata_authority=if (USE_FROZEN_ANNOTATION) "Frozen validated OrgDb SQLite" else "Installed OrgDb package",
    KEGG_organism=if (HAS_KEGG) KEGG_ORG else "none (KEGG enrichment skipped)",
    Motif_tool=if (RUN_MOTIF_ANALYSIS) sprintf("HOMER findMotifsGenome.pl (HOMER v%s)", HOMER_VERSION) else "disabled",
    Motif_genome_package=if (RUN_MOTIF_ANALYSIS) HOMER_GENOME_INFO else NA_character_,
    stringsAsFactors=FALSE
)
write_tsv(annotation_source,
          file.path(annotation_provenance_dir, "annotation_sources.tsv"))
cat(sprintf("  Coordinate annotation: %s\n", TXDB_PKG_LABEL))
cat(sprintf("  Gene metadata / GO:  %s\n", ORGDB_LABEL))

GENE_RANGES <- suppressWarnings(GenomicFeatures::genes(TXDB))

first_mapped_value <- function(ids, column) {
    ids <- as.character(ids)
    result <- rep(NA_character_, length(ids))
    if (is.null(ORGDB_OBJECT)) {
        # Defensive fallback for legacy bundles. Validated user profiles and
        # all built-ins require an OrgDb, so this branch should not normally run.
        return(result)
    }
    valid <- !is.na(ids) & nzchar(ids) & ids != "." & toupper(ids) != "NA"
    if (!any(valid) || !(column %in% AnnotationDbi::columns(ORGDB_OBJECT))) {
        return(result)
    }

    mapping_ids <- normalize_annotation_ids(ids[valid], TXDB_GENE_KEYTYPE)
    unique_ids <- unique(mapping_ids[!is.na(mapping_ids) & nzchar(mapping_ids)])
    out <- tryCatch(
        AnnotationDbi::mapIds(ORGDB_OBJECT,
            keys=unique_ids, column=column, keytype=TXDB_GENE_KEYTYPE,
            multiVals="first"),
        error=function(e) {
            warning(sprintf("OrgDb mapping to %s failed: %s", column, conditionMessage(e)))
            setNames(rep(NA_character_, length(unique_ids)), unique_ids)
        }
    )
    result[valid] <- unname(as.character(out[mapping_ids]))
    result
}

standardize_chipseeker_feature <- function(annotation) {
    x <- as.character(annotation)
    out <- rep("Other", length(x))
    out[grepl("^Promoter", x, ignore.case=TRUE)] <- "Promoter"
    out[grepl("5' UTR|5UTR", x, ignore.case=TRUE)] <- "5' UTR"
    out[grepl("3' UTR|3UTR", x, ignore.case=TRUE)] <- "3' UTR"
    out[grepl("Exon", x, ignore.case=TRUE)] <- "Exon"
    out[grepl("Intron", x, ignore.case=TRUE)] <- "Intron"
    out[grepl("Downstream", x, ignore.case=TRUE)] <- "Downstream"
    out[grepl("Intergenic", x, ignore.case=TRUE)] <- "Distal Intergenic"
    out[is.na(x) | !nzchar(x)] <- "Unannotated"
    out
}

annotate_peak_set <- function(peaks_df, output_tsv, set_label) {
    required <- c("PeakID", "Chr", "Start", "End")
    missing_cols <- setdiff(required, colnames(peaks_df))
    if (length(missing_cols) > 0) {
        stop("Peak annotation input is missing: ", paste(missing_cols, collapse=", "))
    }

    peaks_df\$PeakID <- as.character(peaks_df\$PeakID)
    peaks_df\$Chr <- as.character(peaks_df\$Chr)
    peaks_df\$Start <- as.integer(peaks_df\$Start)
    peaks_df\$End <- as.integer(peaks_df\$End)
    if (anyDuplicated(peaks_df\$PeakID)) stop("PeakID values must be unique for annotation.")

    valid_coords <- !is.na(peaks_df\$Start) & !is.na(peaks_df\$End) &
                    peaks_df\$Start >= 1 & peaks_df\$End >= peaks_df\$Start
    valid_seq <- peaks_df\$Chr %in% GenomeInfoDb::seqlevels(TXDB)
    valid <- valid_coords & valid_seq

    seq_match_rate <- mean(valid_seq)
    cat(sprintf("  %s: %.1f%% of peaks use chromosome names present in %s\n",
                set_label, 100 * seq_match_rate, TXDB_PKG_LABEL))
    if (seq_match_rate < 0.70) {
        stop(sprintf("Only %.1f%% of peak chromosome names match the TxDb; likely assembly/style mismatch.",
                     100 * seq_match_rate))
    }
    if (any(!valid_seq)) {
        bad_seq <- sort(unique(peaks_df\$Chr[!valid_seq]))
        writeLines(bad_seq,
            file.path(annotation_provenance_dir,
                      paste0(gsub("[^A-Za-z0-9]+", "_", set_label),
                             "_unmatched_seqlevels.txt")))
        warning(sprintf("%s: %d peaks on unmatched TxDb sequences will remain unannotated.",
                        set_label, sum(!valid_seq)))
    }

    result <- data.frame(
        PeakID=peaks_df\$PeakID,
        Chr=peaks_df\$Chr,
        Start=peaks_df\$Start,
        End=peaks_df\$End,
        Width=peaks_df\$End - peaks_df\$Start + 1L,
        AnnotationDetail=rep(NA_character_, nrow(peaks_df)),
        Feature=rep("Unannotated", nrow(peaks_df)),
        GeneID=rep(NA_character_, nrow(peaks_df)),
        GeneIDKeytype=rep(TXDB_GENE_KEYTYPE, nrow(peaks_df)),
        SYMBOL=rep(NA_character_, nrow(peaks_df)),
        ENTREZID=rep(NA_character_, nrow(peaks_df)),
        ENSEMBL=rep(NA_character_, nrow(peaks_df)),
        GENENAME=rep(NA_character_, nrow(peaks_df)),
        TranscriptID=rep(NA_character_, nrow(peaks_df)),
        DistanceToTSS=rep(NA_real_, nrow(peaks_df)),
        DistanceToTES=rep(NA_real_, nrow(peaks_df)),
        GeneBodyPosition=rep(NA_real_, nrow(peaks_df)),
        GeneChr=rep(NA_character_, nrow(peaks_df)),
        GeneStart=rep(NA_integer_, nrow(peaks_df)),
        GeneEnd=rep(NA_integer_, nrow(peaks_df)),
        GeneStrand=rep(NA_character_, nrow(peaks_df)),
        stringsAsFactors=FALSE
    )

    if (any(valid)) {
        gr <- GenomicRanges::GRanges(
            seqnames=peaks_df\$Chr[valid],
            ranges=IRanges::IRanges(start=peaks_df\$Start[valid], end=peaks_df\$End[valid])
        )
        names(gr) <- peaks_df\$PeakID[valid]
        peak_anno <- ChIPseeker::annotatePeak(
            gr, TxDb=TXDB, tssRegion=c(-3000, 3000), verbose=FALSE
        )
        ann <- as.data.frame(peak_anno)
        required_ann_cols <- c("seqnames", "start", "end", "annotation", "geneId")
        missing_ann_cols <- setdiff(required_ann_cols, colnames(ann))
        if (length(missing_ann_cols) > 0) {
            stop("ChIPseeker output is missing required columns: ",
                 paste(missing_ann_cols, collapse=", "))
        }
        if (nrow(ann) != sum(valid)) {
            dropped <- sum(valid) - nrow(ann)
            dropped_pct <- 100 * dropped / sum(valid)
            # ChIPseeker does not guarantee 1:1 output rows on large peak sets --
            # a small number of peaks can be silently dropped internally (this is
            # a known ChIPseeker behavior, not specific to this pipeline). Losing
            # a handful of peaks out of hundreds of thousands should not fail the
            # entire run; those peaks are left annotated as "Unannotated" below
            # instead. A large drop is a different story and still stops the run,
            # since that would point at something actually wrong (wrong TxDb,
            # coordinate system mismatch, etc.) rather than this normal edge case.
            if (dropped_pct > 2) {
                stop(sprintf(
                    "ChIPseeker returned %d rows for %d valid peaks (%.2f%% dropped) -- too large to treat as a benign mismatch.",
                    nrow(ann), sum(valid), dropped_pct))
            }
            warning(sprintf(
                "%s: ChIPseeker returned %d rows for %d valid peaks (%d peaks / %.3f%% will be marked Unannotated).",
                set_label, nrow(ann), sum(valid), dropped, dropped_pct))
        }

        # ChIPseeker may sort its GRanges internally. Restore the original peak
        # order explicitly before attaching differential statistics.
        ann_ids <- paste0(as.character(ann\$seqnames), ":", ann\$start, "-", ann\$end)
        expected_ids <- peaks_df\$PeakID[valid]
        reorder_idx <- match(expected_ids, ann_ids)
        if (anyNA(reorder_idx)) {
            unmatched <- sum(is.na(reorder_idx))
            unmatched_pct <- 100 * unmatched / length(reorder_idx)
            if (unmatched_pct > 2) {
                stop(sprintf(
                    "Could not match %d/%d (%.2f%%) ChIPseeker output rows back to input peak IDs -- too large to treat as a benign mismatch.",
                    unmatched, length(reorder_idx), unmatched_pct))
            }
            warning(sprintf(
                "%s: %d peak(s) (%.3f%%) could not be matched back to a ChIPseeker output row and will be marked Unannotated.",
                set_label, unmatched, unmatched_pct))
        }
        ann <- ann[reorder_idx, , drop=FALSE]

        idx <- which(valid)
        gene_ids <- as.character(ann\$geneId)
        result\$AnnotationDetail[idx] <- as.character(ann\$annotation)
        result\$Feature[idx] <- standardize_chipseeker_feature(ann\$annotation)
        result\$GeneID[idx] <- gene_ids
        if ("transcriptId" %in% colnames(ann)) result\$TranscriptID[idx] <- as.character(ann\$transcriptId)
        if ("distanceToTSS" %in% colnames(ann)) result\$DistanceToTSS[idx] <- as.numeric(ann\$distanceToTSS)

        result\$SYMBOL[idx] <- first_mapped_value(gene_ids, "SYMBOL")
        result\$ENTREZID[idx] <- first_mapped_value(gene_ids, "ENTREZID")
        result\$ENSEMBL[idx] <- first_mapped_value(gene_ids, "ENSEMBL")
        result\$GENENAME[idx] <- first_mapped_value(gene_ids, "GENENAME")

        gene_idx <- match(gene_ids, names(GENE_RANGES))
        have_gene <- !is.na(gene_idx)
        if (any(have_gene)) {
            target <- idx[have_gene]
            gg <- GENE_RANGES[gene_idx[have_gene]]
            # gg can contain the same gene multiple times (many peaks can
            # share a nearest gene), so its inherited names (gene IDs) are
            # not unique. as.data.frame(gg) would use those names as
            # row.names and fail with "duplicate row.names" on any peak
            # set large enough for genes to have >1 nearby peak (i.e.
            # essentially always). Pull fields directly instead.
            gstart <- as.integer(BiocGenerics::start(gg))
            gend <- as.integer(BiocGenerics::end(gg))
            gstrand <- as.character(BiocGenerics::strand(gg))
            mid <- floor((result\$Start[target] + result\$End[target]) / 2)
            glen <- pmax(1, gend - gstart)

            result\$GeneChr[target] <- as.character(GenomeInfoDb::seqnames(gg))
            result\$GeneStart[target] <- gstart
            result\$GeneEnd[target] <- gend
            result\$GeneStrand[target] <- gstrand
            result\$DistanceToTES[target] <- ifelse(
                gstrand == "+", mid - gend,
                ifelse(gstrand == "-", gstart - mid, NA_real_)
            )
            pos <- ifelse(
                gstrand == "+", (mid - gstart) / glen,
                ifelse(gstrand == "-", (gend - mid) / glen, NA_real_)
            )
            pos[pos < 0 | pos > 1] <- NA_real_
            result\$GeneBodyPosition[target] <- pos
        }

        assigned <- unique(gene_ids[!is.na(gene_ids) & nzchar(gene_ids) & gene_ids != "." & toupper(gene_ids) != "NA"])
        if (length(assigned) > 0) {
            if (HAS_ORGDB) {
                org_keys <- AnnotationDbi::keys(ORGDB_OBJECT, keytype=TXDB_GENE_KEYTYPE)
                mapping_ids <- normalize_annotation_ids(assigned, TXDB_GENE_KEYTYPE)
                assigned_rate <- mean(mapping_ids %in% org_keys)
                cat(sprintf("  %s: TxDb → OrgDb mapping for assigned genes: %.1f%% (%d unique genes)\n",
                            set_label, 100 * assigned_rate, length(assigned)))
                if (assigned_rate < 0.70) {
                    stop(sprintf("%s annotation gene mapping fell below 70%% (%.1f%%).",
                                 set_label, 100 * assigned_rate))
                }
            } else {
                # No OrgDb to check gene IDs against -- TxDb-based coordinate
                # annotation (nearest gene, TSS/TES distance, gene body
                # position) is entirely unaffected by this and already
                # written to \`result\` above; this block only ever validated
                # the TxDb↔OrgDb ID join, which doesn't apply here.
                cat(sprintf("  %s: %d gene(s) assigned by coordinate (no OrgDb registered -- ID-mapping check skipped)\n",
                            set_label, length(assigned)))
            }
        }
    }

    extra_cols <- setdiff(colnames(peaks_df), c("PeakID", "Chr", "Start", "End"))
    if (length(extra_cols) > 0) result <- cbind(result, peaks_df[, extra_cols, drop=FALSE])
    write_tsv(result, output_tsv)
    cat(sprintf("  %s: ChIPseeker annotation saved → %s (%d peaks)\n",
                set_label, basename(output_tsv), nrow(result)))
    result
}

# ── Helper: normalize names for safe matching ────────────────
# This lets us match sample names even if R/DiffBind converts characters
# like "-" into "." in column names.
normalize_name <- function(x) {
    tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

# ── Helper: identify true sample-count columns for PCA ───────
find_pca_sample_columns <- function(count_mat, samples_df, pca_dir) {
    sample_ids <- as.character(samples_df\$SampleID)
    table_cols <- colnames(count_mat)

    selected_cols <- rep(NA_character_, length(sample_ids))
    match_method  <- rep(NA_character_, length(sample_ids))

    # 1) Exact match: best and preferred.
    exact_idx <- match(sample_ids, table_cols)
    exact_ok  <- !is.na(exact_idx)
    selected_cols[exact_ok] <- table_cols[exact_idx[exact_ok]]
    match_method[exact_ok]  <- "exact"

    # 2) make.names() match: common when column names contain hyphens.
    remaining <- which(is.na(selected_cols))
    if (length(remaining) > 0) {
        safe_table  <- make.names(table_cols, unique=FALSE)
        safe_sample <- make.names(sample_ids, unique=FALSE)
        safe_idx <- match(safe_sample[remaining], safe_table)
        safe_ok <- !is.na(safe_idx)

        selected_cols[remaining[safe_ok]] <- table_cols[safe_idx[safe_ok]]
        match_method[remaining[safe_ok]]  <- "make.names"
    }

    # 3) punctuation-insensitive match: last controlled fallback.
    remaining <- which(is.na(selected_cols))
    if (length(remaining) > 0) {
        norm_table  <- normalize_name(table_cols)
        norm_sample <- normalize_name(sample_ids)

        for (ii in remaining) {
            hits <- which(norm_table == norm_sample[ii])
            if (length(hits) == 1) {
                selected_cols[ii] <- table_cols[hits]
                match_method[ii]  <- "punctuation_insensitive"
            }
        }
    }

    map_df <- data.frame(
        SampleID = sample_ids,
        PCAColumn = selected_cols,
        MatchMethod = match_method,
        stringsAsFactors = FALSE
    )

    # Diagnostics: every column in the DiffBind peak table.
    col_diag <- data.frame(
        Column = table_cols,
        Class = vapply(count_mat, function(x) paste(class(x), collapse=","), character(1)),
        IsNumeric = vapply(count_mat, is.numeric, logical(1)),
        SelectedForPCA = table_cols %in% selected_cols,
        stringsAsFactors = FALSE
    )

    write_tsv(map_df, file.path(pca_dir, "pca_sample_column_map.tsv"))
    write_tsv(col_diag, file.path(pca_dir, "diffbind_peak_table_column_classes.tsv"))

    if (anyNA(selected_cols)) {
        missing <- map_df\$SampleID[is.na(map_df\$PCAColumn)]
        cat("  [WARN] PCA skipped: could not confidently find sample-count columns for:\n")
        cat(paste0("    - ", missing, collapse="\n"), "\n")
        cat("  Diagnostics written:\n")
        cat("    ", file.path(pca_dir, "pca_sample_column_map.tsv"), "\n", sep="")
        cat("    ", file.path(pca_dir, "diffbind_peak_table_column_classes.tsv"), "\n", sep="")
        return(NULL)
    }

    if (any(duplicated(selected_cols))) {
        cat("  [WARN] PCA skipped: multiple samples mapped to the same count column.\n")
        cat("  Diagnostics written:\n")
        cat("    ", file.path(pca_dir, "pca_sample_column_map.tsv"), "\n", sep="")
        return(NULL)
    }

    non_numeric <- selected_cols[!vapply(count_mat[, selected_cols, drop=FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
        cat("  [WARN] PCA skipped: matched sample columns are not numeric:\n")
        cat(paste0("    - ", non_numeric, collapse="\n"), "\n")
        cat("  Diagnostics written:\n")
        cat("    ", file.path(pca_dir, "diffbind_peak_table_column_classes.tsv"), "\n", sep="")
        return(NULL)
    }

    selected_cols
}

# ── Step A: Load sample sheet + build DiffBind object ───────
cat("\n[1/6] Loading sample sheet and counting reads...\n")

# Hard guard: the generated sheet is the complete boundary of the analysis.
# Any sample absent from this sheet cannot enter consensus construction,
# normalization, PCA, contrasts, enrichment backgrounds, or the explorer bundle.
sample_sheet_check <- read.csv(DIFFBIND_SHEET, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(sample_sheet_check) != EXPECTED_SAMPLE_COUNT) {
    stop(sprintf(
        "Sample-isolation check failed before DiffBind: expected %d included rows, found %d.",
        EXPECTED_SAMPLE_COUNT, nrow(sample_sheet_check)
    ))
}
if (anyDuplicated(sample_sheet_check\$SampleID)) {
    stop("Sample-isolation check failed: duplicate SampleID values in DiffBind sample sheet.")
}

dba_obj <- dba(sampleSheet=DIFFBIND_SHEET)

n_loaded_samples <- nrow(dba_obj\$samples)
cat("  Samples loaded:", n_loaded_samples, "\n")

loaded_ids <- as.character(dba_obj\$samples\$SampleID)
sheet_ids  <- as.character(sample_sheet_check\$SampleID)
if (n_loaded_samples != EXPECTED_SAMPLE_COUNT || !setequal(loaded_ids, sheet_ids)) {
    stop(
        "Sample-isolation check failed after DiffBind load.\n",
        "Sheet IDs: ", paste(sheet_ids, collapse=", "), "\n",
        "Loaded IDs: ", paste(loaded_ids, collapse=", ")
    )
}
cat("  Sample isolation verified: DiffBind loaded only the explicitly included samples.\n")

# Save the sample metadata exactly as DiffBind interpreted it.
diagnostics_dir <- make_dir(file.path(DIFF_OUT, "diagnostics"))
write_tsv(as.data.frame(dba_obj\$samples), file.path(diagnostics_dir, "diffbind_samples_loaded.tsv"))

cat(sprintf("  MIN_OVERLAP: %d (set interactively; peak must appear in at least this many samples)\n", MIN_OVERLAP))
cat(sprintf("  DIFFBIND_FILTER: %s (min site activity kept by dba.count(); override via PEPATAC_DIFFBIND_FILTER)\n", DIFFBIND_FILTER))

# dba.count() does its own internal parallel counting via the DBA object's
# own "cores" config option (set below)
# it does NOT use the BiocParallel backend registered later for dba.analyze().
# Left unset, DiffBind defaults to BiocParallel::multicoreWorkers() (every
# logical core on the machine), ignoring the thread count chosen
# interactively above. That can spawn far more parallel BAM-reading workers
# than there is RAM for — especially under WSL2, where the VM often has
# much less memory available than the host — and counting silently
# thrashes instead of speeding up. Set it explicitly so the counting step
# actually uses what was chosen.
dba_obj\$config\$RunParallel <- DIFF_THREADS > 1
dba_obj\$config\$cores       <- DIFF_THREADS
cat(sprintf("  DiffBind counting cores set to: %d\n", DIFF_THREADS))

# Build consensus peak set. MIN_OVERLAP controls how stringent
# the consensus is: a peak must appear in at least this many samples.
#
# IMPORTANT, worth being explicit about: this consensus peak set is built
# ONCE, here, across every included sample -- before Step D subsets down
# to each contrast's own two groups for normalization and dispersion
# estimation. That means the *statistics* for a given contrast are fully
# isolated to its own two groups (no other loaded group can influence a
# comparison it isn't part of -- see Step D), but the *feature universe
# being tested* is shared across every contrast: a peak only present in a
# group that isn't part of a given comparison can still occupy a "slot" in
# that comparison's tested peak set, and the same shared peak set is what
# feeds the accessible-peak background used for GO/KEGG enrichment later.
# This is a deliberate, common design choice (every contrast shares one
# consistent, comparable peak universe) rather than an oversight -- the
# alternative, recounting consensus peaks separately per contrast, would
# be slower and would make tested peak sets differ between comparisons,
# which has its own downsides. It's stated here explicitly so it's a
# documented tradeoff rather than an implicit one.
#
# summits=: see the PEPATAC_DIFFBIND_SUMMITS comment in diff_analysis.sh.
# Pinned to 200 (401bp window) by default -- this pipeline's own default,
# not a deferral to whatever DiffBind currently defaults to. Only reaches
# the DiffBind-default branch below if a user explicitly sets
# PEPATAC_DIFFBIND_SUMMITS="" to opt back into that.
count_args <- list(dba_obj, minOverlap=MIN_OVERLAP, score=DBA_SCORE_READS,
                   filter=DIFFBIND_FILTER, bParallel=DIFF_THREADS > 1)
if (nzchar(DIFFBIND_SUMMITS)) {
    summits_val <- if (toupper(DIFFBIND_SUMMITS) %in% c("TRUE", "FALSE")) {
        as.logical(DIFFBIND_SUMMITS)
    } else {
        suppressWarnings(as.numeric(DIFFBIND_SUMMITS))
    }
    if (is.na(summits_val)) {
        stop(sprintf("PEPATAC_DIFFBIND_SUMMITS must be TRUE, FALSE, or a number (got: %s)",
                     DIFFBIND_SUMMITS))
    }
    count_args\$summits <- summits_val
    summits_recorded <- sprintf("%s (overridden via PEPATAC_DIFFBIND_SUMMITS)", DIFFBIND_SUMMITS)
    cat(sprintf("  Peak re-centering (summits): %s\n", summits_recorded))
} else {
    summits_recorded <- "DiffBind's own default, not pinned (explicit PEPATAC_DIFFBIND_SUMMITS=\"\" override)"
    cat(sprintf("  Peak re-centering (summits): %s\n", summits_recorded))
}
dba_obj <- do.call(dba.count, count_args)

# Peak re-centering (summits=) is otherwise only visible in the console
# log while the run is happening -- an upstream DiffBind default change
# would then silently change feature geometry with no record of it after
# the fact. Persisted here alongside the other dba.count() parameters that
# together determine the counted feature set.
write_tsv(
    data.frame(
        Summits=summits_recorded,
        MinOverlap=MIN_OVERLAP,
        Score_metric="DBA_SCORE_READS",
        Filter_threshold=DIFFBIND_FILTER,
        Background_normalization_bins=DIFFBIND_BACKGROUND,
        stringsAsFactors=FALSE
    ),
    file.path(diagnostics_dir, "diffbind_count_parameters.tsv")
)

# Retrieve the counted consensus peak table once for diagnostics/PCA.
# IMPORTANT: this table is not assumed to be purely numeric.
count_mat <- dba.peakset(dba_obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
cat("  Consensus peaks:", nrow(count_mat), "\n")
cat("  Count matrix built.\n")

# ── Save full raw count matrix (TSV.GZ + RDS) ────────────────
# The run summary promises this file. It is useful for downstream
# work (clustering, custom plots) without re-running DiffBind.
count_mat_path_rds <- file.path(diagnostics_dir, "consensus_peak_raw_counts.rds")
count_mat_path_tsv <- file.path(diagnostics_dir, "consensus_peak_raw_counts.tsv.gz")
saveRDS(count_mat, count_mat_path_rds)
tryCatch({
    con <- gzfile(count_mat_path_tsv, "w")
    write.table(count_mat, con, sep="\t", quote=FALSE, row.names=FALSE)
    close(con)
    cat(sprintf("  Raw count matrix saved: %d peaks × %d columns\n",
                nrow(count_mat), ncol(count_mat)))
}, error = function(e) {
    cat("  [WARN] Could not write gzipped count matrix:", conditionMessage(e), "\n")
    cat("  RDS version still available at:", count_mat_path_rds, "\n")
})

# Save a small preview so the PCA input is inspectable without opening huge files.
write_tsv(head(count_mat, 25), file.path(diagnostics_dir, "diffbind_peak_table_preview_first25.tsv"))

# ── Write and annotate the consensus-peak background ──────────
# The same TxDb/OrgDb annotation function is used here and for each
# significant peak set. This keeps GO foreground and universe definitions
# on one coordinate system and one identifier mapping.
bg_bed_path <- file.path(diagnostics_dir, "all_consensus_peaks.bed")
bg_ann_path <- file.path(diagnostics_dir, "all_consensus_peaks_annotated.tsv")

coord_cols_present <- all(c("Chr","Start","End") %in% names(count_mat))
if (!coord_cols_present) {
    lc <- tolower(names(count_mat))
    if (all(c("chr","start","end") %in% lc)) {
        names(count_mat)[lc == "chr"]   <- "Chr"
        names(count_mat)[lc == "start"] <- "Start"
        names(count_mat)[lc == "end"]   <- "End"
        coord_cols_present <- TRUE
    }
}

if (coord_cols_present) {
    bg_peaks <- data.frame(
        PeakID=paste0(count_mat\$Chr, ":", count_mat\$Start, "-", count_mat\$End),
        Chr=as.character(count_mat\$Chr),
        Start=as.integer(count_mat\$Start),
        End=as.integer(count_mat\$End),
        stringsAsFactors=FALSE
    )

    write.table(
        # DiffBind/GRanges coordinates are 1-based; BED starts are 0-based.
        data.frame(bg_peaks\$Chr, pmax(0L, bg_peaks\$Start - 1L), bg_peaks\$End,
                   bg_peaks\$PeakID, 0, "."),
        bg_bed_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
    )
    cat(sprintf("  Background BED written: %d peaks → %s\n",
                nrow(bg_peaks), basename(bg_bed_path)))

    background_annotation <- annotate_peak_set(
        bg_peaks, bg_ann_path, "Consensus-background peaks"
    )
    rm(background_annotation, bg_peaks)
    invisible(gc())
} else {
    stop(sprintf(
        "Could not identify Chr/Start/End columns in the DiffBind count matrix. Available: %s",
        paste(names(count_mat), collapse=", ")
    ))
}

# ── Checkpoint: save counted DiffBind object ─────────────────
# Read-counting from BAMs is the slowest step (30–60 min for 6 samples).
# Saving here means a downstream crash doesn't require a full recount.
counted_rds <- file.path(diagnostics_dir, "diffbind_counted.rds")
saveRDS(dba_obj, counted_rds)
cat(sprintf("  Checkpoint saved: %s\n", basename(counted_rds)))

# ── Step B: Whole-experiment normalization overview (informational) ──
cat("\n[2/6] Computing a whole-experiment normalization overview...\n")
# IMPORTANT: this does NOT normalize dba_obj itself, and is NOT what any
# individual contrast's DESeq2 model actually uses below. Each contrast in
# Step D is normalized independently, on a subset containing only its own
# two groups, using DiffBind's modern design-based contrast mode
# (dba.contrast(design="~Condition", ...)) rather than the legacy
# group1=/group2= form. That distinction matters here specifically:
# DiffBind's own documentation states that when dba.contrast() is set up
# with design=FALSE (which is what supplying group1=/group2= without
# design= does), dba.normalize()'s DESeq2 default reverts to DBA_NORM_LIB
# "for backwards compatibility" -- i.e. contrast mode and normalization
# are coupled in DiffBind's legacy path in a way they are not in the
# modern path. A single normalization computed once, globally, and then
# handed to a legacy-mode contrast risked silently not being the
# normalization actually used for the real per-contrast statistics, even
# though it was logged and saved as if it were. Normalizing per-contrast,
# on the modern design path, removes that ambiguity rather than requiring
# it to be proven safe.
#
# What follows here is deliberately just a whole-experiment sanity check
# (one badly-behaved library is often visible here before any contrast
# runs) -- not a value fed into any contrast's actual statistics. For the
# factors a specific contrast actually used, see that contrast's own
# <label>/diagnostics/diffbind_normalization_factors.tsv, written in Step D.
overview_obj <- tryCatch(
    dba.normalize(dba_obj, method=DBA_DESEQ2, normalize=DBA_NORM_RLE,
                 background=DIFFBIND_BACKGROUND),
    error = function(e) {
        cat("  [WARN] Whole-experiment normalization overview failed:", conditionMessage(e), "\n")
        NULL
    }
)

if (!is.null(overview_obj)) {
    tryCatch({
        norm_info <- dba.normalize(overview_obj, method=DBA_DESEQ2, bRetrieve=TRUE)
        if (!is.null(norm_info\$norm.factors)) {
            nf_df <- data.frame(
                Sample     = overview_obj\$samples\$SampleID,
                # Full precision -- this file is also the fallback source of
                # truth for the explorer if a contrast-specific file isn't
                # available. Rounding here would only ever lose precision
                # for no benefit; a display-only rounded column can be added
                # by any downstream consumer that wants one.
                NormFactor = as.numeric(norm_info\$norm.factors),
                stringsAsFactors = FALSE
            )
            cat("  Whole-experiment normalization factors (informational only -- see note above):\n")
            print(nf_df, row.names=FALSE)
            nf_range <- max(nf_df\$NormFactor) / min(nf_df\$NormFactor)
            if (nf_range > 5) {
                cat(sprintf("  [WARN] Normalization factor range is %.1fx — inspect FRiP and library sizes.\n", nf_range))
            }
            write_tsv(nf_df, file.path(diagnostics_dir, "diffbind_normalization_factors_whole_experiment.tsv"))
        }
    }, error = function(e) {
        cat("  [WARN] Could not retrieve whole-experiment normalization factors:", conditionMessage(e), "\n")
    })
}
cat("  (per-contrast factors are computed independently in Step D below)\n")

# ── Step C: PCA across all samples (done once) ───────────────
cat("\n[3/6] PCA across all samples...\n")
pca_dir <- make_dir(file.path(DIFF_OUT, "pca"))

# Use DiffBind's own sample table — authoritative, version-independent.
samples_df <- as.data.frame(dba_obj\$samples)
n_samples  <- nrow(samples_df)

if (n_samples >= 3) {
    # PCA is intentionally built as an ATAC-seq analogue of your RNA-seq DESeq2 PCA:
    #   RNA-seq: genes × samples raw counts  -> DESeq2 transform -> PCA
    #   ATAC:    peaks × samples raw counts  -> DESeq2 transform -> PCA
    #
    # Differential testing below still uses DiffBind/DESeq2 on the raw peak counts.
    # This PCA uses a DESeq2 variance-stabilizing transformation only for visualization/QC.
    # It does NOT use CPM and it does NOT feed transformed values into differential testing.
    #
    # Crucially, we select ONLY columns that map to actual SampleID values.
    # We do NOT use "everything except coordinates", because DiffBind peak tables
    # can contain non-count metadata columns.
    sample_cols <- find_pca_sample_columns(count_mat, samples_df, pca_dir)

    if (!is.null(sample_cols)) {
        raw_counts_df <- count_mat[, sample_cols, drop=FALSE]
        raw_counts <- as.matrix(raw_counts_df)
        storage.mode(raw_counts) <- "numeric"
        colnames(raw_counts) <- as.character(samples_df\$SampleID)

        # Build peak IDs only from coordinate columns that actually exist.
        have_coords <- all(c("Chr","Start","End") %in% names(count_mat))
        peak_ids <- if (have_coords) {
            paste0(count_mat\$Chr, ":", count_mat\$Start, "-", count_mat\$End)
        } else {
            paste0("peak_", seq_len(nrow(raw_counts)))
        }
        rownames(raw_counts) <- make.unique(peak_ids)

        # Clean the matrix before handing it to DESeq2.
        # DESeq2 expects non-negative integer counts.
        raw_counts[!is.finite(raw_counts)] <- NA_real_
        if (anyNA(raw_counts)) {
            cat("  [WARN] PCA count matrix contained NA/Inf values; replacing them with 0 before DESeq2 VST.\n")
            raw_counts[is.na(raw_counts)] <- 0
        }
        if (any(raw_counts < 0, na.rm=TRUE)) {
            cat("  [WARN] PCA count matrix contained negative values; clipping them to 0 before DESeq2 VST.\n")
            raw_counts[raw_counts < 0] <- 0
        }

        # Remove peaks with no counts in all samples.
        keep_nonzero <- rowSums(raw_counts, na.rm=TRUE) > 0
        raw_counts <- raw_counts[keep_nonzero, , drop=FALSE]

        # Round only for the DESeq2 PCA transform. The original DiffBind object is unchanged.
        deseq_counts <- round(raw_counts)
        storage.mode(deseq_counts) <- "integer"

        sample_meta <- data.frame(
            Sample = as.character(samples_df\$SampleID),
            Condition = factor(as.character(samples_df\$Condition)),
            Replicate = if ("Replicate" %in% names(samples_df)) {
                factor(as.character(samples_df\$Replicate))
            } else {
                factor(as.character(seq_len(n_samples)))
            },
            PCAColumn = sample_cols,
            RawLibraryCountForPCA = as.numeric(colSums(raw_counts, na.rm=TRUE)),
            stringsAsFactors = FALSE
        )
        rownames(sample_meta) <- sample_meta\$Sample
        write_tsv(sample_meta, file.path(pca_dir, "pca_sample_metadata.tsv"))

        matrix_summary_pre <- data.frame(
            Metric = c(
                "samples",
                "consensus_peaks_before_filtering",
                "peaks_after_nonzero_filter",
                "min_raw_library_count_for_pca",
                "max_raw_library_count_for_pca",
                "pca_transform"
            ),
            Value = c(
                ncol(deseq_counts),
                nrow(count_mat),
                nrow(deseq_counts),
                min(sample_meta\$RawLibraryCountForPCA),
                max(sample_meta\$RawLibraryCountForPCA),
                "DESeq2 variance-stabilizing transformation; size factors estimated with type='poscounts'"
            )
        )
        write_tsv(matrix_summary_pre, file.path(pca_dir, "pca_matrix_summary_pre_vst.tsv"))

        if (nrow(deseq_counts) < 2) {
            cat("  [WARN] PCA skipped: fewer than 2 nonzero consensus peaks after filtering.\n")
        } else if (ncol(deseq_counts) < 3) {
            cat(sprintf("  [WARN] PCA skipped: %d sample(s) — PCA requires at least 3.\n", ncol(deseq_counts)))
        } else {
            # Use design=~1 and blind=TRUE because this is QC PCA, not differential testing.
            # type='poscounts' is safer for sparse ATAC peak matrices than the default
            # median-ratio size factor method, which can fail when many peaks contain zeros.
            dds_pca <- DESeqDataSetFromMatrix(
                countData = deseq_counts,
                colData   = sample_meta,
                design    = ~ 1
            )

            vst_mat <- NULL
            vst_method <- NA_character_
            vst_error <- NULL

            dds_pca <- tryCatch(
                estimateSizeFactors(dds_pca, type="poscounts"),
                error = function(e) {
                    vst_error <<- paste("estimateSizeFactors(type='poscounts') failed:", conditionMessage(e))
                    NULL
                }
            )

            if (is.null(dds_pca)) {
                cat("  [WARN] PCA skipped:", vst_error, "\n")
                write_tsv(data.frame(Error=vst_error), file.path(pca_dir, "pca_deseq2_vst_error.tsv"))
            } else {
                # Try the fast DESeq2::vst() first. It can fail on very small feature sets,
                # so fall back to varianceStabilizingTransformation(), which is slower but
                # more tolerant for small test matrices.
                vsd <- tryCatch({
                    vst_method <<- "DESeq2::vst(blind=TRUE)"
                    vst(dds_pca, blind=TRUE)
                }, error = function(e1) {
                    cat("  [WARN] DESeq2::vst() failed; falling back to varianceStabilizingTransformation().\n")
                    cat("         vst() message:", conditionMessage(e1), "\n")
                    tryCatch({
                        vst_method <<- "DESeq2::varianceStabilizingTransformation(blind=TRUE)"
                        varianceStabilizingTransformation(dds_pca, blind=TRUE)
                    }, error = function(e2) {
                        vst_error <<- paste(
                            "Both DESeq2 VST methods failed.",
                            "vst():", conditionMessage(e1),
                            "varianceStabilizingTransformation():", conditionMessage(e2)
                        )
                        NULL
                    })
                })

                if (is.null(vsd)) {
                    cat("  [WARN] PCA skipped:", vst_error, "\n")
                    write_tsv(data.frame(Error=vst_error), file.path(pca_dir, "pca_deseq2_vst_error.tsv"))
                } else {
                    vst_mat <- assay(vsd)

                    # Remove zero-variance peaks after VST. We use scale.=FALSE, matching
                    # DESeq2 plotPCA behavior on transformed values.
                    peak_var <- apply(vst_mat, 1, var, na.rm=TRUE)
                    keep_var <- is.finite(peak_var) & peak_var > 0
                    vst_mat <- vst_mat[keep_var, , drop=FALSE]

                    matrix_summary <- data.frame(
                        Metric = c(
                            "samples",
                            "consensus_peaks_before_filtering",
                            "peaks_after_nonzero_filter",
                            "peaks_after_vst_variance_filter",
                            "min_raw_library_count_for_pca",
                            "max_raw_library_count_for_pca",
                            "vst_method"
                        ),
                        Value = c(
                            ncol(deseq_counts),
                            nrow(count_mat),
                            nrow(deseq_counts),
                            nrow(vst_mat),
                            min(sample_meta\$RawLibraryCountForPCA),
                            max(sample_meta\$RawLibraryCountForPCA),
                            vst_method
                        )
                    )
                    write_tsv(matrix_summary, file.path(pca_dir, "pca_matrix_summary.tsv"))

                    size_factor_df <- data.frame(
                        Sample = names(sizeFactors(dds_pca)),
                        DESeq2SizeFactor = as.numeric(sizeFactors(dds_pca)),
                        stringsAsFactors = FALSE
                    )
                    write_tsv(size_factor_df, file.path(pca_dir, "pca_deseq2_size_factors.tsv"))

                    # Store the full transformed matrix as RDS for reproducibility without
                    # creating a massive text file. Also write a small preview for inspection.
                    saveRDS(vst_mat, file.path(pca_dir, "pca_deseq2_vst_matrix.rds"))
                    preview_n <- min(1000, nrow(vst_mat))
                    preview_df <- data.frame(PeakID=rownames(vst_mat)[seq_len(preview_n)],
                                             vst_mat[seq_len(preview_n), , drop=FALSE],
                                             check.names=FALSE)
                    write_tsv(preview_df, file.path(pca_dir, "pca_deseq2_vst_matrix_preview_first1000.tsv"))

                    if (nrow(vst_mat) < 2) {
                        cat("  [WARN] PCA skipped: fewer than 2 variable peaks after DESeq2 VST.\n")
                    } else {
                        pca_res <- prcomp(t(vst_mat), center=TRUE, scale.=FALSE)

                        if (ncol(pca_res\$x) < 2) {
                            cat("  [WARN] PCA skipped: fewer than 2 principal components available.\n")
                        } else {
                            pca_df <- as.data.frame(pca_res\$x[, 1:2, drop=FALSE])
                            pca_df\$Sample <- rownames(pca_df)
                            pca_df <- left_join(pca_df, sample_meta, by="Sample")

                            var_exp <- round(100 * pca_res\$sdev^2 / sum(pca_res\$sdev^2), 1)

                            write_tsv(pca_df, file.path(pca_dir, "pca_all_samples_coordinates.tsv"))
                            write_tsv(
                                data.frame(PC=paste0("PC", seq_along(var_exp)), PercentVariance=var_exp),
                                file.path(pca_dir, "pca_percent_variance.tsv")
                            )

                            cat("  PCA coordinates computed from DESeq2 VST-transformed peak counts.\n")
                            cat("  PCA diagnostics written to:", pca_dir, "\n")
                            cat("  Use the explorer (menu item 1) to plot PCA interactively.\n")
                        }
                    }
                }
            }
        }
    }
} else {
    cat(sprintf("  Skipping PCA: %d sample(s) — PCA requires at least 3.\n", n_samples))
}

# ── Step D: Per-contrast differential analysis ───────────────
cat("\n[4/6] Running per-contrast DESeq2 analysis...\n")

all_results <- list()

for (ci in seq_along(CONTRAST_LABELS)) {
    case_grp  <- CONTRAST_CASES[ci]
    ctrl_grp  <- CONTRAST_CTRLS[ci]
    label     <- CONTRAST_LABELS[ci]

    cat(sprintf("\n  Contrast %d/%d: %s (%s vs %s)\n",
        ci, length(CONTRAST_LABELS), label, case_grp, ctrl_grp))

    # Guard: DESeq2 requires ≥2 samples per group to estimate dispersion.
    conditions <- dba_obj\$samples\$Condition
    n_case <- sum(conditions == case_grp)
    n_ctrl <- sum(conditions == ctrl_grp)
    if (n_case < 2 || n_ctrl < 2) {
        cat(sprintf("  [SKIP] %s: DESeq2 requires ≥2 samples per group (%s=%d, %s=%d). Skipping.\n",
            label, case_grp, n_case, ctrl_grp, n_ctrl))
        next
    }

    # Check that both groups exist in the loaded sample sheet.
    if (!case_grp %in% conditions || !ctrl_grp %in% conditions) {
        cat(sprintf("  [SKIP] Groups '%s' or '%s' not found in loaded samples.\n",
            case_grp, ctrl_grp))
        next
    }

    # Subset to a fresh DBA object containing only this contrast's two
    # groups. Counts are preserved from the single global dba.count() call
    # earlier (dba() with a sample-level mask keeps counts; it's peak-level
    # masks that would force a recount) -- only the sample set is
    # restricted here, so no other loaded group can influence dispersion
    # estimation or normalization for this specific comparison.
    #
    # Built directly from the Condition column rather than DiffBind's own
    # named masks (dba_obj\$masks[[case_grp]]). DiffBind auto-generates a
    # mask for every unique value across several attribute categories at
    # once -- Tissue, Factor, Condition, Treatment, Caller, Replicate --
    # so a group name that happens to also match, say, a Replicate number
    # used elsewhere in the sheet could resolve to the wrong mask. Reading
    # samples\$Condition directly removes that indirection entirely: the
    # mask is guaranteed to reflect exactly the Condition column, which is
    # the column this script actually assigns group labels into.
    dba_sample_conditions <- as.character(dba_obj\$samples\$Condition)
    contrast_mask <- dba_sample_conditions %in% c(case_grp, ctrl_grp)
    dba_sub <- tryCatch(
        dba(dba_obj, mask=contrast_mask),
        error = function(e) {
            cat(sprintf("  [SKIP] %s: could not subset to this contrast's samples (%s).\n",
                        label, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(dba_sub)) next

    # Modern design-based contrast on the two-group subset. Deliberately
    # NOT group1=/group2= (which puts dba.contrast() into DiffBind's
    # pre-3.0 "no-design" backward-compatibility mode -- see the note at
    # the top of Step B for why that specifically matters here: DiffBind's
    # own documentation ties normalization defaults to which contrast mode
    # is in use).
    dba_sub <- tryCatch(
        dba.contrast(dba_sub, design="~Condition",
                    contrast=c("Condition", case_grp, ctrl_grp)),
        error = function(e) {
            cat(sprintf("  [SKIP] %s: dba.contrast(design=) failed (%s).\n",
                        label, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(dba_sub)) next

    # Normalize THIS SUBSET explicitly. These are the factors that will
    # actually be used for this contrast's DESeq2 model -- not inherited
    # from the whole-experiment overview computed in Step B, which never
    # touched dba_obj and was informational only.
    dba_sub <- dba.normalize(dba_sub, method=DBA_DESEQ2,
                             normalize=DBA_NORM_RLE,
                             background=DIFFBIND_BACKGROUND)

    # bBlacklist=FALSE, bGreylist=FALSE: without these, dba.analyze()'s own
    # defaults (DBA\$config, normally TRUE) can silently call dba.blacklist()
    # itself before analyzing -- DiffBind auto-detects the genome and
    # applies its OWN blacklist if one is available for it, independent of
    # whatever blacklist decision was already made upstream. Two problems
    # follow from that: (1) blacklist handling is already an explicit
    # runner-level choice (PEPATAC_run_flexible_paths.sh's USE_BLACKLIST),
    # already baked into the peaks/BAMs before this script ever sees them --
    # a second, independent blacklist application inside DiffBind can
    # silently override that choice regardless of what the runner decided;
    # (2) if it removes consensus intervals, the matrix dba.analyze() ends
    # up testing differs from the one that was just normalized above, which
    # is exactly the kind of change that would make normalization factors
    # saved before this call describe a matrix that's no longer the one
    # actually analyzed. Disabling DiffBind's hidden layer keeps the
    # blacklist decision entirely where it was already made -- upstream, in
    # the runner -- rather than adding a second, redundant decision point
    # here that could disagree with it.
    BiocParallel::register(BiocParallel::MulticoreParam(workers=DIFF_THREADS))
    dba_sub <- dba.analyze(dba_sub, method=DBA_DESEQ2,
                           bBlacklist=FALSE, bGreylist=FALSE,
                           bParallel=DIFF_THREADS > 1)

    # Also retrieve the UNSHRUNKEN MLE log2FoldChange, alongside DiffBind's
    # own Fold (which dba.report() below may return apeglm/ashr-shrunk --
    # see DiffBind's own docs). Verified on real data: plotting the shrunken
    # Fold against this same contrast's Wald p-value can make the effect-
    # size/significance relationship look far tighter than it really is,
    # because shrinkage strength is itself a function of the same standard
    # error that drives the p-value -- two related-but-distinct quantities
    # that can visually read as one. The MLE estimate is the actual tested
    # coefficient (stat = MLE_log2FC / lfcSE, exactly) and is the more
    # conventional volcano-plot x-axis. Both are kept here, not one chosen
    # for the user: MLE is noisier but direct/unregularized; shrunken Fold
    # is more stable for ranking/visualization but can overstate how
    # deterministic the effect-size/significance relationship looks. The
    # explorer's re-plot menu lets users pick between them (see
    # run_replot_analysis() in PEPATAC_explore...sh). Retrieved here, once,
    # right after the real analysis -- re-running this from the explorer
    # later would mean reloading the BAM-counting checkpoint from scratch.
    # Non-fatal on failure: mle_cols stays NULL, the explorer just falls
    # back to shrunken Fold automatically (see the merge below).
    mle_cols <- NULL
    tryCatch({
        mle_dds <- dba.analyze(dba_sub, bRetrieveAnalysis=DBA_DESEQ2)
        if (!is.null(mle_dds)) {
            mle_res_df <- as.data.frame(
                DESeq2::results(mle_dds, contrast=c("Condition", case_grp, ctrl_grp))
            )
            peaks_coords <- dba.peakset(dba_sub, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
            # dba.peakset()'s Chr/Start/End capitalization varies by DiffBind
            # version/context -- don't assume, look it up (same defensive
            # pattern the explorer's peak-boxplot code already uses).
            find_coord_col <- function(df, candidates) {
                hit <- intersect(candidates, names(df))
                if (length(hit) == 0) NA_character_ else hit[1]
            }
            pc_chr   <- find_coord_col(peaks_coords, c("Chr", "chr", "CHR", "seqnames"))
            pc_start <- find_coord_col(peaks_coords, c("Start", "start", "START"))
            pc_end   <- find_coord_col(peaks_coords, c("End", "end", "END"))
            if (is.na(pc_chr) || is.na(pc_start) || is.na(pc_end)) {
                cat(sprintf("  [WARN] %s: could not retrieve unshrunken MLE log2FC -- dba.peakset() coordinate columns not recognized (found: %s). Only the DiffBind Fold will be available for plotting.\n",
                            label, paste(names(peaks_coords), collapse=", ")))
            } else {
                # DiffBind's design-based DESeq2 path (pv.DEinitDESeq2, used
                # whenever dba.contrast() was called with design=/contrast=
                # as this pipeline does) can drop peaks below its own
                # internal count filter before ever building the
                # DESeqDataSet -- separately from dba.count()'s own
                # DIFFBIND_FILTER, and separately from dba.report()'s own
                # per-contrast filtering below. When that happens, it does
                # NOT renumber the surviving rows 1..N: it sets
                # rownames(counts) <- which(keep), i.e. the ORIGINAL row
                # index into the peakset that was passed in, preserved
                # exactly (verified against DiffBind's actual current
                # source, R/analyze_deseq2.R, and confirmed
                # bRetrieveAnalysis=DBA_DESEQ2 returns that object with zero
                # further modification -- DBA.R's return(DBA\$DESeq2\$DEdata)).
                # So rather than assuming every consensus peak survived into
                # the DESeq2 fit (the assumption that broke on the K4meTest
                # dataset -- 54,461 peaks vs 54,437 DESeq2 rows), read those
                # rownames back as indices into peaks_coords directly, and
                # validate every part of that assumption explicitly rather
                # than trusting it implicitly:
                #   1. rownames(mle_res_df) must be IDENTICAL to
                #      rownames(mle_dds) -- DESeq2::results() is documented
                #      to preserve row order/names from its input, but this
                #      whole fix exists to stop trusting implicit row-order
                #      assumptions, so it is asserted here instead.
                #   2. every rowname must match ^[1-9][0-9]*\$ literally
                #      (a plain positive integer, no leading zero, no
                #      whitespace/decimal) BEFORE it is ever passed to
                #      as.integer() -- so a rowname format some other
                #      DiffBind version uses can't be silently
                #      misinterpreted as an index.
                #   3. the resulting indices must be unique and in range.
                # If a future DiffBind version changes this internal
                # behavior, any one of these checks failing means this
                # fails loudly (a [WARN], MLE simply unavailable) rather
                # than silently mismatching peaks to the wrong MLE values.
                mle_row_ids      <- rownames(mle_dds)
                mle_res_row_ids  <- rownames(mle_res_df)
                row_ids_match    <- identical(mle_res_row_ids, mle_row_ids)
                row_id_format_ok <- length(mle_row_ids) > 0 &&
                                    all(grepl("^[1-9][0-9]*\$", mle_row_ids))
                mle_idx <- if (row_id_format_ok) suppressWarnings(as.integer(mle_row_ids)) else integer(0)

                idx_valid <- row_ids_match &&
                             row_id_format_ok &&
                             length(mle_idx) == nrow(mle_res_df) &&
                             length(mle_idx) > 0 &&
                             !anyNA(mle_idx) &&
                             !anyDuplicated(mle_idx) &&
                             all(mle_idx >= 1 & mle_idx <= nrow(peaks_coords))

                if (idx_valid) {
                    mle_coords <- peaks_coords[mle_idx, , drop=FALSE]
                    mle_cols <- data.frame(
                        Chr = mle_coords[[pc_chr]],
                        Start = mle_coords[[pc_start]],
                        End = mle_coords[[pc_end]],
                        MLE_log2FoldChange = mle_res_df\$log2FoldChange,
                        MLE_lfcSE = mle_res_df\$lfcSE,
                        MLE_stat = mle_res_df\$stat,
                        stringsAsFactors = FALSE
                    )

                    # Coordinates must uniquely identify each mapped peak,
                    # or the coordinate-based attachment further below
                    # (match(), specifically chosen over merge() so a
                    # duplicate key can never silently multiply rows) would
                    # have an ambiguous target. Refuse rather than guess.
                    mle_key <- paste(mle_cols\$Chr, mle_cols\$Start, mle_cols\$End, sep="\t")
                    if (anyDuplicated(mle_key)) {
                        cat(sprintf("  [WARN] %s: could not retrieve unshrunken MLE log2FC -- mapped MLE coordinates are not unique (%d duplicate coordinate key(s)); refusing ambiguous coordinate attachment. Only the DiffBind Fold will be available for plotting.\n",
                                    label, sum(duplicated(mle_key))))
                        mle_cols <- NULL
                    } else {
                        n_filtered <- nrow(peaks_coords) - nrow(mle_res_df)
                        cat(sprintf("  Unshrunken MLE log2FC also retrieved for %s -- available to the explorer as an alternate x-axis.\n",
                                    label))
                        cat(sprintf("    Consensus peaks:               %d\n", nrow(peaks_coords)))
                        cat(sprintf("    Peaks represented in DESeq2:   %d\n", nrow(mle_res_df)))
                        cat(sprintf("    Filtered before DESeq2 model:  %d\n", n_filtered))
                        cat(sprintf("    MLE coordinate mapping:        VERIFIED (rownames(mle_dds) used as original-peak indices, validated format/row-order/unique/in-range/complete)\n"))

                        mle_diag_dir <- make_dir(file.path(DIFF_OUT, label, "diagnostics"))
                        map_summary <- data.frame(
                            Metric = c("Consensus_peaks", "Peaks_represented_in_DESeq2",
                                       "Filtered_before_DESeq2_model", "Mapping_verified"),
                            Value = c(nrow(peaks_coords), nrow(mle_res_df), n_filtered, "TRUE"),
                            stringsAsFactors = FALSE
                        )
                        write_tsv(map_summary, file.path(mle_diag_dir, "diffbind_mle_mapping_summary.tsv"))

                        if (n_filtered > 0) {
                            omitted_idx <- setdiff(seq_len(nrow(peaks_coords)), mle_idx)
                            omitted_df <- peaks_coords[omitted_idx, c(pc_chr, pc_start, pc_end), drop=FALSE]
                            names(omitted_df) <- c("Chr", "Start", "End")
                            write_tsv(omitted_df, file.path(mle_diag_dir, "mle_peaks_filtered_before_deseq2.tsv"))
                            # Coordinates only -- this shows WHICH peaks were
                            # dropped before the DESeq2 fit, not why (that
                            # would require also saving the filter
                            # score/threshold DiffBind's internal filterFun
                            # actually applied).
                            cat(sprintf("    Filtered peak coordinates saved (which peaks were dropped, not why): mle_peaks_filtered_before_deseq2.tsv\n"))
                        }
                    }
                } else {
                    cat(sprintf("  [WARN] %s: could not retrieve unshrunken MLE log2FC -- rownames(mle_dds) did not validate as usable original-peak indices (row order matches results()=%s, format all ^[1-9][0-9]*=%s, n=%d vs %d DESeq2 results, has NA=%s, has duplicates=%s, all in-range=%s). Only the DiffBind Fold will be available for plotting.\n",
                                label, row_ids_match, row_id_format_ok, length(mle_idx), nrow(mle_res_df),
                                if (length(mle_idx) > 0) anyNA(mle_idx) else NA,
                                if (length(mle_idx) > 0) anyDuplicated(mle_idx) > 0 else NA,
                                if (length(mle_idx) > 0) all(mle_idx >= 1 & mle_idx <= nrow(peaks_coords)) else NA))
                }
            }
        } else {
            # dba.analyze(bRetrieveAnalysis=DBA_DESEQ2) returning NULL is not
            # an R error, so the tryCatch below never fires for this case --
            # log it explicitly rather than leave it silent.
            cat(sprintf("  [WARN] %s: could not retrieve unshrunken MLE log2FC -- dba.analyze(bRetrieveAnalysis=DBA_DESEQ2) returned NULL (no error thrown). Only the DiffBind Fold will be available for plotting.\n",
                        label))
        }
    }, error = function(e) {
        cat(sprintf("  [WARN] %s: could not retrieve unshrunken MLE log2FC (%s). Only the DiffBind Fold will be available for plotting.\n",
                    label, conditionMessage(e)))
    })

    # Save this contrast's actual normalization factors AFTER analyzing,
    # not before -- retrieving them post-dba.analyze() (rather than
    # immediately after the dba.normalize() call above) proves they
    # describe the object that was actually analyzed, not a pre-analysis
    # snapshot that could have been invalidated by a hidden step in
    # between. With bBlacklist/bGreylist explicitly disabled above, no such
    # hidden step exists anymore -- but retrieving after dba.analyze()
    # rather than before costs nothing and removes any doubt about it,
    # here or in any future change to this sequence.
    #
    # So the explorer -- and anyone auditing results later -- can find the
    # exact per-sample factors this specific comparison used. These can
    # legitimately differ from another contrast's factors for a sample
    # that appears in both, since each contrast is normalized
    # independently on its own two-group subset; that's expected, not a
    # bug, under the modern per-contrast model.
    contrast_diag_dir <- make_dir(file.path(DIFF_OUT, label, "diagnostics"))
    tryCatch({
        norm_info_c <- dba.normalize(dba_sub, method=DBA_DESEQ2, bRetrieve=TRUE)
        if (!is.null(norm_info_c\$norm.factors)) {
            nf_df_c <- data.frame(
                Sample     = dba_sub\$samples\$SampleID,
                NormFactor = as.numeric(norm_info_c\$norm.factors),
                stringsAsFactors = FALSE
            )
            write_tsv(nf_df_c, file.path(contrast_diag_dir, "diffbind_normalization_factors.tsv"))
            cat(sprintf("  Normalization factors (this contrast's own, RLE, background=%s):\n",
                        DIFFBIND_BACKGROUND))
            print(nf_df_c, row.names=FALSE)
        }
    }, error = function(e) {
        cat(sprintf("  [WARN] %s: could not save this contrast's normalization factors: %s\n",
                    label, conditionMessage(e)))
    })

    # Extract results at the chosen FDR cutoff (th=1 pulls everything).
    # precision=0: dba.report()'s default (precision=2:3 for DataType=
    # DBA_DATA_FRAME) rounds Fold to 2 decimal places and signif()s p-value/
    # FDR to 3 significant figures before returning them. FC_CUTOFF/FDR_CUTOFF
    # below are compared directly against these values, so without
    # precision=0 a peak genuinely near either threshold can be classified
    # differently than the full-precision DESeq2 result would give -- e.g. a
    # true log2FC of 0.5849 rounds to 0.58 and fails ">= 0.585" even though
    # the real value was within 0.0001 of the cutoff. precision=0 stores full
    # precision throughout (bundle, CSVs, plots), matching what's actually
    # tested against the thresholds.
    res_dba <- dba.report(dba_sub, method=DBA_DESEQ2,
                          th=1, bUsePval=FALSE, precision=0,
                          DataType=DBA_DATA_FRAME)

    if (is.null(res_dba) || nrow(res_dba) == 0) {
        cat("  [WARN] No results returned for this contrast.\n")
        next
    }

    # Standardize column names across DiffBind versions.
    # Avoid regex strings here because this R script is generated from Bash.
    # Direct name replacement is safer and avoids escape-sequence errors.
    nms <- names(res_dba)
    nms[nms == "Fold"] <- "log2FoldChange"
    nms[nms %in% c("p.value", "pvalue", "p-value")] <- "pvalue"
    nms[nms == "FDR"] <- "padj"
    names(res_dba) <- nms

    # Attach the unshrunken MLE columns retrieved above (if that succeeded)
    # by genomic coordinate. match() is used instead of merge() specifically
    # so a duplicate coordinate key can never silently multiply res_dba's
    # rows -- the MLE mapping block above already refuses to build mle_cols
    # with duplicate keys, but this re-checks independently rather than
    # relying on that guarantee holding across an edit to either block.
    # match() also guarantees res_dba's own row order/count is unchanged,
    # which merge() does not strictly guarantee; asserted explicitly below
    # rather than assumed.
    if (!is.null(mle_cols)) {
        res_key <- paste(res_dba\$Chr, res_dba\$Start, res_dba\$End, sep="\t")
        mle_key <- paste(mle_cols\$Chr, mle_cols\$Start, mle_cols\$End, sep="\t")
        if (anyDuplicated(mle_key)) {
            cat(sprintf("  [WARN] %s: MLE coordinate keys are not unique at attachment time -- skipping MLE attachment to the report.\n", label))
        } else {
            n_before_mle_attach <- nrow(res_dba)
            mle_match <- match(res_key, mle_key)
            res_dba\$MLE_log2FoldChange <- mle_cols\$MLE_log2FoldChange[mle_match]
            res_dba\$MLE_lfcSE          <- mle_cols\$MLE_lfcSE[mle_match]
            res_dba\$MLE_stat           <- mle_cols\$MLE_stat[mle_match]
            stopifnot(nrow(res_dba) == n_before_mle_attach)
            n_mle_matched <- sum(!is.na(mle_match))
            cat(sprintf("  MLE columns attached to DiffBind report: %d/%d rows matched by coordinate.\n",
                        n_mle_matched, nrow(res_dba)))
            if (n_mle_matched != nrow(res_dba)) {
                cat(sprintf("  [WARN] %s: %d DiffBind report rows had no MLE coordinate match; those MLE fields are NA.\n",
                            label, nrow(res_dba) - n_mle_matched))
            }
        }
    }

    res_dba\$Contrast  <- label
    # Significance requires BOTH FDR and fold-change thresholds (set interactively).
    # FC_CUTOFF == 0 means FDR-only filtering.
    # NOTE: log2FoldChange here is DiffBind's own (possibly shrunken) Fold
    # column, NOT MLE_log2FoldChange attached above -- the MLE columns are
    # an alternate plotting/diagnostic quantity only and do not change
    # which peaks this pipeline calls significant.
    if (FC_CUTOFF > 0) {
        res_dba\$Sig <- !is.na(res_dba\$padj) &
                       res_dba\$padj < FDR_CUTOFF &
                       abs(res_dba\$log2FoldChange) >= FC_CUTOFF
    } else {
        res_dba\$Sig <- !is.na(res_dba\$padj) & res_dba\$padj < FDR_CUTOFF
    }
    res_dba\$Direction <- ifelse(res_dba\$Sig & res_dba\$log2FoldChange > 0,
                                "Up", ifelse(res_dba\$Sig & res_dba\$log2FoldChange < 0,
                                "Down", "NS"))

    all_results[[label]] <- res_dba

    # ── Contrast output directory ──────────────────────────
    c_dir <- make_dir(file.path(DIFF_OUT, label))

    # ── CSV: all peaks ─────────────────────────────────────
    csv_path <- file.path(c_dir, paste0(label, "_all_peaks.csv"))
    write.csv(res_dba, csv_path, row.names=FALSE)
    cat(sprintf("  Saved: %s\n", basename(csv_path)))

    # ── CSV: significant peaks only ────────────────────────
    sig_res <- res_dba[res_dba\$Sig, ]
    sig_csv <- file.path(c_dir, paste0(label, "_significant_peaks.csv"))
    write.csv(sig_res, sig_csv, row.names=FALSE)
    sig_label <- if (FC_CUTOFF > 0) {
        sprintf("FDR<%.3f & |log2FC|>=%.3f", FDR_CUTOFF, FC_CUTOFF)
    } else {
        sprintf("FDR<%.3f", FDR_CUTOFF)
    }
    cat(sprintf("  Significant peaks (%s): %d  (up: %d, down: %d)\n",
        sig_label,
        sum(res_dba\$Sig),
        sum(res_dba\$Direction == "Up"),
        sum(res_dba\$Direction == "Down")))
}

# ── Step E: ChIPseeker peak annotation ───────────────────────
cat("\n[5a/6] Annotating significant peaks with ChIPseeker...\n")

for (ci in seq_along(CONTRAST_LABELS)) {
    label   <- CONTRAST_LABELS[ci]
    res_dba <- all_results[[label]]
    if (is.null(res_dba)) next

    sig_res <- res_dba[res_dba\$Sig, , drop=FALSE]
    if (nrow(sig_res) == 0) {
        cat(sprintf("  %s: no significant peaks to annotate.\n", label))
        next
    }

    c_dir    <- file.path(DIFF_OUT, label)
    bed_path <- file.path(c_dir, paste0(label, "_sig_peaks.bed"))
    ann_path <- file.path(c_dir, paste0(label, "_annotated_peaks.tsv"))

    peak_ids <- paste0(sig_res\$Chr, ":", sig_res\$Start, "-", sig_res\$End)
    extra <- sig_res[, setdiff(colnames(sig_res), c("Chr", "Start", "End")), drop=FALSE]
    peak_input <- cbind(
        data.frame(
            PeakID=peak_ids,
            Chr=as.character(sig_res\$Chr),
            Start=as.integer(sig_res\$Start),
            End=as.integer(sig_res\$End),
            stringsAsFactors=FALSE
        ),
        extra
    )

    write.table(
        # DiffBind/GRanges coordinates are 1-based; BED starts are 0-based.
        data.frame(peak_input\$Chr, pmax(0L, peak_input\$Start - 1L), peak_input\$End,
                   peak_input\$PeakID, 0, "."),
        bed_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
    )

    # Directional BEDs (Up-only / Down-only), alongside the combined one
    # above. Opening and closing chromatin often reflect different TF
    # programs; pooling both directions into one HOMER motif search can
    # dilute or cancel signal that a directional search would catch. Used
    # by the directional HOMER split in Step E2 below.
    for (dir_tag in c("Up", "Down")) {
        dir_idx <- which(peak_input\$Direction == dir_tag)
        if (length(dir_idx) == 0) next
        dir_bed_path <- file.path(c_dir, paste0(label, "_", tolower(dir_tag), "_peaks.bed"))
        write.table(
            data.frame(peak_input\$Chr[dir_idx], pmax(0L, peak_input\$Start[dir_idx] - 1L),
                       peak_input\$End[dir_idx], peak_input\$PeakID[dir_idx], 0, "."),
            dir_bed_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
        )
    }

    annotate_peak_set(peak_input, ann_path, label)
}

# ── Step E2: HOMER motif analysis ────────────────────────────
cat("\n[5b/6] HOMER motif analysis...\n")

# Runs findMotifsGenome.pl for a single BED file and reports/parses the
# results into motif_dir. Factored out so the combined + directional
# (Up-only / Down-only) runs below share one implementation instead of
# tripling this block.
run_homer_for_bed <- function(bed_path, motif_dir, label, run_note,
                               find_motifs_pl, bg_arg) {
    cat(sprintf("\n  Motifs: %s\n", run_note))

    if (!file.exists(bed_path)) {
        cat(sprintf("  %s: peaks BED not found — skipping.\n", run_note))
        return(invisible(NULL))
    }

    n_sig <- as.integer(system(
        sprintf("wc -l < %s", shQuote(bed_path)), intern=TRUE))
    if (is.na(n_sig) || n_sig < 5) {
        cat(sprintf("  %s: too few peaks (%d) for motif analysis.\n",
                    run_note, max(0L, n_sig, na.rm=TRUE)))
        return(invisible(NULL))
    }

    dir.create(motif_dir, showWarnings=FALSE, recursive=TRUE)

    # Build findMotifsGenome.pl command.
    # -size given: use exact peak coordinates rather than re-centering.
    # -p: threads.
    # MOTIF_MODE (chosen interactively) decides which HOMER analyses run:
    #   known  -> -nomotif  (known-motif enrichment only; fast)
    #   denovo -> -noknown  (de novo discovery only; slow)
    #   both   -> (no flag; HOMER runs both; slowest)
    motif_flag <- switch(MOTIF_MODE,
        known  = "-nomotif",
        denovo = "-noknown",
        both   = "",
        "-nomotif"   # safe fallback: known only
    )
    mode_label <- switch(MOTIF_MODE,
        known  = "known motifs only",
        denovo = "de novo only",
        both   = "known + de novo",
        "known motifs only"
    )
    motif_cmd <- sprintf(
        "perl %s %s %s %s -size given -bg %s -p %d %s > %s 2>&1",
        shQuote(find_motifs_pl),
        shQuote(bed_path),
        shQuote(HOMER_GENOME_ARG),
        shQuote(motif_dir),
        bg_arg,
        max(1L, DIFF_THREADS),
        motif_flag,
        shQuote(file.path(motif_dir, "homer_motif_stdout.txt"))
    )

    cat(sprintf("  %s: running findMotifsGenome.pl (%d peaks, %s)...\n",
                run_note, n_sig, mode_label))
    motif_ret <- system(motif_cmd)

    if (motif_ret == 0) {
        # Output files depend on which analyses were requested.
        known_html  <- file.path(motif_dir, "knownResults.html")
        known_txt   <- file.path(motif_dir, "knownResults.txt")
        denovo_html <- file.path(motif_dir, "homerResults.html")

        # ── Known-motif enrichment (present unless de novo-only) ──
        if (MOTIF_MODE %in% c("known", "both")) {
            if (file.exists(known_txt) && file.size(known_txt) > 0) {
                # Parse top 10 known motifs for the log.
                known_df <- tryCatch(
                    read.delim(known_txt, header=TRUE, stringsAsFactors=FALSE,
                               check.names=FALSE, quote=""),
                    error = function(e) NULL
                )
                if (!is.null(known_df) && nrow(known_df) > 0) {
                    # Column names vary by HOMER version — find p-value col.
                    pval_col <- grep("P-value|pvalue|p.value|Log P-value",
                                     colnames(known_df),
                                     ignore.case=TRUE, value=TRUE)[1]
                    name_col <- grep("Motif Name|motif.name|name",
                                     colnames(known_df),
                                     ignore.case=TRUE, value=TRUE)[1]
                    pct_col  <- grep("% of Target|Target%|target.percent",
                                     colnames(known_df),
                                     ignore.case=TRUE, value=TRUE)[1]

                    cat(sprintf("  %s: top known motifs:\n", run_note))
                    top_n <- min(10, nrow(known_df))
                    for (ri in seq_len(top_n)) {
                        mname <- if (!is.na(name_col)) known_df[ri, name_col] else "unknown"
                        mpval <- if (!is.na(pval_col)) known_df[ri, pval_col] else "?"
                        mpct  <- if (!is.na(pct_col))  known_df[ri, pct_col]  else "?"
                        cat(sprintf("    %2d. %-40s  p=%s  target%%=%s\n",
                                    ri, mname, mpval, mpct))
                    }

                    # Save clean top-motifs summary CSV.
                    write.csv(known_df[seq_len(min(50, nrow(known_df))), ],
                        file.path(motif_dir, paste0(label, "_known_motifs_top50.csv")),
                        row.names=FALSE)
                }
                if (file.exists(known_html)) {
                    cat(sprintf("  %s: known motifs (browser): %s\n", run_note, known_html))
                }
            } else {
                cat(sprintf("  %s: known-motif enrichment produced no output.\n", run_note))
                cat(sprintf("     Check: %s\n",
                            file.path(motif_dir, "homer_motif_stdout.txt")))
            }
        }

        # ── De novo discovery (present unless known-only) ─────────
        if (MOTIF_MODE %in% c("denovo", "both")) {
            if (file.exists(denovo_html)) {
                cat(sprintf("  %s: de novo motifs (browser): %s\n", run_note, denovo_html))
            } else {
                cat(sprintf("  %s: de novo discovery produced no homerResults.html.\n", run_note))
                cat(sprintf("     Check: %s\n",
                            file.path(motif_dir, "homer_motif_stdout.txt")))
            }
        }

        cat(sprintf("  %s: motif results → %s/\n", run_note, basename(motif_dir)))
    } else {
        cat(sprintf("  %s: findMotifsGenome.pl failed (exit %d).\n",
                    run_note, motif_ret))
        cat(sprintf("     Check: %s\n",
                    file.path(motif_dir, "homer_motif_stdout.txt")))
    }
    invisible(NULL)
}

if (!RUN_MOTIF_ANALYSIS) {
    cat("  Motif analysis skipped (disabled at startup).\n")
} else {

    # Use the findMotifsGenome.pl path found directly by Bash at startup.
    find_motifs_pl <- HOMER_FIND_MOTIFS

    if (!nzchar(find_motifs_pl) || !file.exists(find_motifs_pl)) {
        cat(sprintf("  [WARN] findMotifsGenome.pl not found — skipping motifs.\n"))
        cat(sprintf("         Searched inside: %s\n",
                    dirname(dirname(HOMER_CONFIGURE))))
    } else {

        # Background: all consensus peaks (correct ATAC-seq background).
        # Same background used for every run below (combined + directional).
        bg_bed <- file.path(DIFF_OUT, "diagnostics", "all_consensus_peaks.bed")
        bg_arg <- if (file.exists(bg_bed)) shQuote(bg_bed) else "automatic"

        for (ci in seq_along(CONTRAST_LABELS)) {
            label <- CONTRAST_LABELS[ci]
            c_dir <- file.path(DIFF_OUT, label)

            # Combined (both directions pooled) -- kept as the default
            # overview output, same behavior as before this change.
            run_homer_for_bed(
                file.path(c_dir, paste0(label, "_sig_peaks.bed")),
                file.path(c_dir, "motifs"),
                label, label, find_motifs_pl, bg_arg
            )

            # Directional (Up-only / Down-only). Opening and closing
            # chromatin often reflect different TF programs; pooling both
            # into a single motif search can dilute or cancel either
            # signal. Skipped automatically (via run_homer_for_bed's own
            # file-exists / n>=5 checks) when a contrast has too few
            # peaks in one direction, or none.
            for (dir_tag in c("Up", "Down")) {
                dir_bed <- file.path(c_dir, paste0(label, "_", tolower(dir_tag), "_peaks.bed"))
                run_homer_for_bed(
                    dir_bed,
                    file.path(c_dir, paste0("motifs_", tolower(dir_tag))),
                    label, paste0(label, " (", dir_tag, "-only)"),
                    find_motifs_pl, bg_arg
                )
            }
        }
    }
}

# ── Step E3: ChIPseeker annotation summary tables ─────────────
cat("\n[5c/6] Summarizing ChIPseeker annotation distributions...\n")

feature_order <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron",
                   "Downstream", "Distal Intergenic", "Other", "Unannotated")

for (ci in seq_along(CONTRAST_LABELS)) {
    label <- CONTRAST_LABELS[ci]
    c_dir <- file.path(DIFF_OUT, label)
    ann_path <- file.path(c_dir, paste0(label, "_annotated_peaks.tsv"))

    if (!file.exists(ann_path)) {
        cat(sprintf("  %s: annotation table not found — skipping annotation summary.\n", label))
        next
    }

    ann <- tryCatch(
        read.delim(ann_path, header=TRUE, stringsAsFactors=FALSE,
                   quote="", comment.char="", check.names=FALSE),
        error=function(e) {
            cat(sprintf("  %s: could not read annotation table: %s\n",
                        label, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(ann) || nrow(ann) == 0) next

    summary_dir <- file.path(c_dir, "annotation_summary")
    dir.create(summary_dir, showWarnings=FALSE, recursive=TRUE)

    # Genomic feature composition — the only aggregate worth precomputing here.
    # Everything else (TSS/TES distance, gene-body position, peak width,
    # chromosome distribution) is a direct per-peak column already sitting in
    # ann_path, so the explorer re-derives those histograms from that table
    # on demand rather than duplicating the computation here.
    feat <- factor(ann\$Feature, levels=feature_order)
    feature_df <- as.data.frame(table(feat), stringsAsFactors=FALSE)
    colnames(feature_df) <- c("Feature", "Count")
    feature_df <- feature_df[feature_df\$Count > 0, , drop=FALSE]
    feature_df\$Percent <- 100 * feature_df\$Count / sum(feature_df\$Count)
    write.csv(feature_df, file.path(summary_dir, "genomic_features_counts.csv"), row.names=FALSE)
    cat(sprintf("  %s: genomic feature counts saved (%d features).\n", label, nrow(feature_df)))
}

# ── Step F: GO and KEGG enrichment analysis ───────────────────
cat("\n[5d/6] Running GO and KEGG enrichment analysis...\n")

# direction = NULL reads all genes in the annotation table (combined,
# both directions pooled). direction = "Up"/"Down" filters to peaks with
# that Direction first. Every significant-peak annotation table carries a
# Direction column (see extra_cols in annotate_peak_set / Step E), so a
# directional request finding no such column is treated as zero genes for
# that direction rather than silently falling back to the combined set.
read_annotation_genes <- function(annotation_tsv, direction = NULL) {
    if (!file.exists(annotation_tsv) || file.size(annotation_tsv) == 0) {
        return(list(entrez=character(0), symbols=character(0)))
    }
    ann <- tryCatch(
        read.delim(annotation_tsv, header=TRUE, stringsAsFactors=FALSE,
                   quote="", comment.char="", check.names=FALSE),
        error=function(e) {
            warning(sprintf("Could not read annotation table %s: %s",
                            annotation_tsv, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(ann)) return(list(entrez=character(0), symbols=character(0)))

    if (!is.null(direction)) {
        if ("Direction" %in% colnames(ann)) {
            ann <- ann[ann\$Direction == direction, , drop=FALSE]
        } else {
            return(list(entrez=character(0), symbols=character(0)))
        }
    }

    entrez <- if ("ENTREZID" %in% colnames(ann)) as.character(ann\$ENTREZID) else character(0)
    symbols <- if ("SYMBOL" %in% colnames(ann)) as.character(ann\$SYMBOL) else character(0)
    entrez <- sort(unique(entrez[!is.na(entrez) & nzchar(entrez) & entrez != "."]))
    symbols <- sort(unique(symbols[!is.na(symbols) & nzchar(symbols) & symbols != "."]))
    list(entrez=entrez, symbols=symbols)
}

# run_note: what gets printed in log messages (e.g. "label (Up-only)").
# file_tag: appended to output filenames (e.g. "up" -> "<label>_up_BP.csv");
# "" (default, combined run) keeps the original unsuffixed filenames.
run_go_kegg <- function(gene_entrez, bg_entrez, gene_symbols, label, c_dir,
                        run_note = label, file_tag = "") {
    gene_entrez <- sort(unique(as.character(gene_entrez)))
    gene_entrez <- gene_entrez[!is.na(gene_entrez) & nzchar(gene_entrez)]
    if (!is.null(bg_entrez)) {
        bg_entrez <- sort(unique(as.character(bg_entrez)))
        bg_entrez <- bg_entrez[!is.na(bg_entrez) & nzchar(bg_entrez)]
    }

    if (length(gene_entrez) < 5) {
        cat(sprintf("  %s: too few mapped Entrez genes (%d) for enrichment — skipping.\n",
                    run_note, length(gene_entrez)))
        return(invisible(NULL))
    }

    # Recorded as its own column in every output CSV below -- not just
    # mentioned in the console log -- so which universe a given enrichment
    # result rests on is machine-readable from the file itself.
    universe_type <- if (is.null(bg_entrez)) {
        "OrgDb_whole_database_default (PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE override)"
    } else {
        "Accessible_peaks (ATAC consensus background)"
    }

    cat(sprintf("  %s: %d foreground Entrez genes; background=%s\n",
                run_note, length(gene_entrez),
                if (is.null(bg_entrez)) "OrgDb default" else length(bg_entrez)))

    go_dir <- file.path(c_dir, "GO")
    dir.create(go_dir, showWarnings=FALSE, recursive=TRUE)

    file_label <- if (nzchar(file_tag)) paste0(label, "_", file_tag) else label

    enrichr_path <- file.path(c_dir, paste0(file_label, "_genes_for_enrichr.txt"))
    writeLines(sort(unique(gene_symbols)), enrichr_path)

    if (HAS_ORGDB) {
        for (ont in c("BP", "MF", "CC")) {
            ont_name <- c(BP="Biological_Process",
                          MF="Molecular_Function",
                          CC="Cellular_Component")[[ont]]
            tryCatch({
                # pvalueCutoff/qvalueCutoff = 1: return every term clusterProfiler
                # actually tested, not just the ones clearing a threshold.
                # clusterProfiler computes p-values for all terms with any gene
                # overlap regardless of these cutoffs, so loosening them to 1
                # costs nothing at runtime -- it just stops silently discarding
                # the non-passing rows before they ever reach a CSV. The output
                # was previously labeled "all terms tested" while actually only
                # containing the pre-filtered subset; this makes that literally
                # true. A separate _significant.csv (FDR<0.05) is still written
                # as a convenience subset.
                go_res <- clusterProfiler::enrichGO(
                    gene=gene_entrez,
                    universe=bg_entrez,
                    OrgDb=ORGDB_OBJECT,
                    keyType="ENTREZID",
                    ont=ont,
                    pAdjustMethod="BH",
                    pvalueCutoff=1,
                    qvalueCutoff=1,
                    readable=TRUE
                )
                res_df <- if (!is.null(go_res)) as.data.frame(go_res) else NULL
                if (!is.null(res_df) && nrow(res_df) > 0) {
                    res_df\$UniverseType <- universe_type
                    write.csv(res_df,
                        file.path(go_dir, paste0(file_label, "_", ont_name, ".csv")),
                        row.names=FALSE)
                    sig_df <- res_df[!is.na(res_df\$p.adjust) & res_df\$p.adjust < 0.05, , drop=FALSE]
                    write.csv(sig_df,
                        file.path(go_dir, paste0(file_label, "_", ont_name, "_significant.csv")),
                        row.names=FALSE)
                    cat(sprintf("  %s: GO %s — %d terms tested (full table), %d significant (FDR<0.05, see _significant.csv)\n",
                                run_note, ont, nrow(res_df), nrow(sig_df)))
                } else {
                    cat(sprintf("  %s: GO %s — enrichGO returned no result.\n", run_note, ont))
                }
            }, error=function(e) {
                cat(sprintf("  %s: GO %s failed: %s\n", run_note, ont, conditionMessage(e)))
            })
        }
    } else {
        cat(sprintf("  %s: GO skipped -- no OrgDb registered for this genome.\n", run_note))
    }

    if (HAS_KEGG) {
        tryCatch({
            kegg_res <- clusterProfiler::enrichKEGG(
                gene=gene_entrez,
                universe=bg_entrez,
                organism=KEGG_ORG,
                pAdjustMethod="BH",
                pvalueCutoff=1,
                qvalueCutoff=1
            )
            kegg_df <- if (!is.null(kegg_res)) as.data.frame(kegg_res) else NULL
            if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
                kegg_df\$UniverseType <- universe_type
                write.csv(kegg_df,
                    file.path(go_dir, paste0(file_label, "_KEGG_pathways.csv")),
                    row.names=FALSE)
                sig_kegg <- kegg_df[!is.na(kegg_df\$p.adjust) & kegg_df\$p.adjust < 0.05, , drop=FALSE]
                write.csv(sig_kegg,
                    file.path(go_dir, paste0(file_label, "_KEGG_pathways_significant.csv")),
                    row.names=FALSE)
                cat(sprintf("  %s: KEGG — %d pathways tested (full table), %d significant (FDR<0.05, see _significant.csv)\n",
                            run_note, nrow(kegg_df), nrow(sig_kegg)))
            } else {
                cat(sprintf("  %s: KEGG — enrichKEGG returned no result.\n", run_note))
            }
        }, error=function(e) {
            cat(sprintf("  %s: KEGG failed: %s\n", run_note, conditionMessage(e)))
        })
    } else {
        cat(sprintf("  %s: KEGG skipped -- no KEGG organism code registered for this genome.\n", run_note))
    }
}

bg_ann_path <- file.path(DIFF_OUT, "diagnostics", "all_consensus_peaks_annotated.tsv")
bg_gene_info <- read_annotation_genes(bg_ann_path)
bg_entrez <- if (length(bg_gene_info\$entrez) >= 5) bg_gene_info\$entrez else NULL
# Referenced again just before the "complete" banner so this doesn't only
# show up once, early, where a long run's output could bury it.
enrichment_universe_is_default <- is.null(bg_entrez)
# Default: skip GO/KEGG entirely when the universe would silently swap to
# the OrgDb whole-database default -- overridable per PEPATAC_ALLOW_
# DEFAULT_ENRICHMENT_UNIVERSE. When the accessible-gene background is
# available, this is always TRUE regardless of the override setting.
run_go_kegg_enabled <- !enrichment_universe_is_default || ALLOW_DEFAULT_ENRICHMENT_UNIVERSE

print_universe_fallback_alert <- function() {
    cat("\n")
    if (run_go_kegg_enabled) {
        cat("  ============================================================\n")
        cat("  GO/KEGG ENRICHMENT UNIVERSE FELL BACK TO THE ORGDB DEFAULT\n")
        cat("  (proceeding anyway: PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE=TRUE)\n")
        cat("  ============================================================\n")
    } else {
        cat("  ============================================================\n")
        cat("  GO/KEGG ENRICHMENT SKIPPED -- NO ACCESSIBLE-GENE BACKGROUND\n")
        cat("  ============================================================\n")
    }
    cat("  Fewer than 5 Entrez genes mapped from the consensus-peak\n")
    cat("  background annotation. Testing against the OrgDb whole-genome/\n")
    cat("  whole-database default universe instead of an ATAC-accessible-\n")
    cat("  genes-only universe computes p-values against every gene in the\n")
    cat("  database, not just genes near an accessible peak in this\n")
    cat("  experiment -- which can inflate significance for terms whose\n")
    cat("  genes aren't chromatin-accessible here at all.\n")
    if (run_go_kegg_enabled) {
        cat("  Every output row below is tagged UniverseType=\n")
        cat("  \"OrgDb_whole_database_default\" so this can't be missed later.\n")
    } else {
        cat("  GO/KEGG enrichment is being skipped for every contrast below\n")
        cat("  as a result (this does not affect DiffBind/DESeq2 differential\n")
        cat("  accessibility results, HOMER motifs, or anything else).\n")
        cat("  Set PEPATAC_ALLOW_DEFAULT_ENRICHMENT_UNIVERSE=TRUE to proceed\n")
        cat("  anyway and accept that caveat.\n")
    }
    cat("  If this is unexpected, check why so few genes mapped:\n")
    cat(sprintf("  %s\n", bg_ann_path))
    cat("  ============================================================\n\n")
}

if (!is.null(bg_entrez)) {
    cat(sprintf("  Accessible-gene background: %d Entrez genes from consensus peaks\n",
                length(bg_entrez)))
    writeLines(bg_gene_info\$symbols,
        file.path(DIFF_OUT, "diagnostics", "background_genes_for_enrichr.txt"))
} else {
    print_universe_fallback_alert()
}

if (!run_go_kegg_enabled) {
    cat("\n  Skipping GO/KEGG enrichment for all contrasts (see alert above).\n")
} else {
    for (ci in seq_along(CONTRAST_LABELS)) {
        label <- CONTRAST_LABELS[ci]
        c_dir <- file.path(DIFF_OUT, label)
        ann_path <- file.path(c_dir, paste0(label, "_annotated_peaks.tsv"))
        cat(sprintf("\n  GO/KEGG: %s\n", label))

        fg <- read_annotation_genes(ann_path)
        if (length(fg\$entrez) == 0) {
            cat(sprintf("  %s: no mapped Entrez genes in ChIPseeker annotation — skipping.\n", label))
            next
        }

        # Combined (both directions pooled) -- kept as the default overview
        # output, same behavior as before this change.
        run_go_kegg(fg\$entrez, bg_entrez, fg\$symbols, label, c_dir, run_note=label)

        # Directional (Up-only / Down-only). Opening and closing chromatin
        # often reflect different regulatory programs; pooling both directions
        # into one enrichment test can dilute or cancel either signal.
        for (dir_tag in c("Up", "Down")) {
            fg_dir <- read_annotation_genes(ann_path, direction=dir_tag)
            if (length(fg_dir\$entrez) == 0) {
                cat(sprintf("  %s (%s-only): no mapped Entrez genes — skipping.\n", label, dir_tag))
                next
            }
            run_go_kegg(fg_dir\$entrez, bg_entrez, fg_dir\$symbols, label, c_dir,
                        run_note=paste0(label, " (", dir_tag, "-only)"),
                        file_tag=tolower(dir_tag))
        }
    }
}

# ── Step G: Cross-contrast summary ───────────────────────────
cat("\n[6/6] Writing cross-contrast summary...\n")

summary_rows <- lapply(CONTRAST_LABELS, function(lbl) {
    r <- all_results[[lbl]]
    if (is.null(r)) return(NULL)
    data.frame(
        Contrast       = lbl,
        Total_peaks    = nrow(r),
        Sig_peaks      = sum(r\$Sig, na.rm=TRUE),
        Up_peaks       = sum(r\$Direction == "Up",   na.rm=TRUE),
        Down_peaks     = sum(r\$Direction == "Down",  na.rm=TRUE),
        FDR_cutoff     = FDR_CUTOFF,
        log2FC_cutoff  = FC_CUTOFF,
        Normalization  = sprintf("DBA_NORM_RLE(background=%s)", DIFFBIND_BACKGROUND),
        stringsAsFactors = FALSE
    )
})
summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

# Guard: if every contrast was skipped (most commonly because a group had
# fewer than 2 samples -- see the [SKIP] messages in Step D above), don't
# fall through into writing a corrupt/empty summary and printing a
# "complete" banner over zero actual results.
if (is.null(summary_df) || nrow(summary_df) == 0) {
    stop(paste(
        "No contrast produced results -- every contrast was skipped.",
        "This is most commonly because a group had fewer than 2 samples",
        "(DESeq2 requires >=2 replicates per group for dispersion",
        "estimation). Check the [SKIP] messages above, fix group",
        "assignments/replicates, and re-run."
    ))
}

summary_csv <- file.path(DIFF_OUT, "contrast_summary.csv")
write.csv(summary_df, summary_csv, row.names=FALSE)

cat("\n  Cross-contrast summary:\n\n")
print(summary_df, row.names=FALSE)
cat(sprintf("\n  Summary saved: %s\n", summary_csv))

# Shown once already (right where the fallback happened) and repeated here
# on purpose -- someone skimming straight to the end of a long run's
# output should not miss that GO/KEGG enrichment rested on the OrgDb
# default universe rather than an ATAC-accessible-genes-only background.
if (enrichment_universe_is_default) {
    print_universe_fallback_alert()
}

cat("\n══════════════════════════════════════════════════\n")
cat("  Differential analysis complete.\n")
cat(sprintf("  Results: %s\n", DIFF_OUT))
cat("══════════════════════════════════════════════════\n\n")

# ── Explorer bundle ───────────────────────────────────────────
# Saves everything the explorer needs to disk in one place so the
# interactive explorer script (PEPATAC_explore.sh) can load it
# without re-running DiffBind or re-reading BAMs.
#
# Contents (lean — no DiffBind object, no raw count matrix):
#   vst_matrix    — peaks × samples VST matrix (already in memory from PCA)
#   all_results   — named list of per-contrast result data frames (th=1, all peaks)
#   col_data      — sample metadata (SampleID, Condition, Replicate)
#   contrast_*    — contrast definition vectors
#   fdr/fc cutoff — original run thresholds (explorer can override for re-plots)
#   homer_*       — absolute paths to HOMER tools (validated at run time)
#   genome / txdb / orgdb — annotation provenance and GO/KEGG re-plots
#
# vst_matrix will be NULL when the global PCA step was skipped (< 3 samples
# or VST failed). The explorer detects this and disables PCA / correlation
# heatmap gracefully, but all other menu items still work.

cat("\nSaving explorer bundle...\n")

tryCatch({
    # Collect sample metadata from DiffBind's own table — authoritative.
    explorer_col_data <- tryCatch(
        as.data.frame(dba_obj\$samples)[, intersect(
            c("SampleID", "Condition", "Replicate"),
            colnames(as.data.frame(dba_obj\$samples))
        ), drop = FALSE],
        error = function(e) {
            warning("Could not extract sample metadata from dba_obj: ", conditionMessage(e))
            NULL
        }
    )

    # vst_mat may or may not exist depending on whether the PCA step succeeded.
    # Check for it safely — it lives in the pca/ subdirectory RDS if the global
    # variable was cleaned up, but it should still be in scope at this point.
    vst_for_bundle <- if (exists("vst_mat") && !is.null(vst_mat)) {
        vst_mat
    } else {
        # Try loading from the PCA RDS written earlier in this session.
        pca_rds <- file.path(DIFF_OUT, "pca", "pca_deseq2_vst_matrix.rds")
        if (file.exists(pca_rds)) {
            tryCatch(readRDS(pca_rds),
                     error = function(e) { warning("Could not reload VST from PCA RDS: ", conditionMessage(e)); NULL })
        } else {
            NULL
        }
    }

    go_result_files <- list.files(
        DIFF_OUT, pattern="_GO_(BP|MF|CC)[.]csv$",
        recursive=TRUE, full.names=TRUE
    )
    kegg_result_files <- list.files(
        DIFF_OUT, pattern="_KEGG_pathways[.]csv$",
        recursive=TRUE, full.names=TRUE
    )
    has_go_results <- length(go_result_files) > 0
    has_kegg_results <- length(kegg_result_files) > 0

    # Only advertise contrasts that actually completed. CONTRAST_LABELS/
    # CASES/CTRLS above are every contrast the user REQUESTED, but Step D's
    # dba.analyze() loop uses \`next\` to skip past ones that fail its
    # <2-replicates-per-group guard (or a subset/contrast/analyze error),
    # leaving no entry in all_results[[label]] for them. Advertising the
    # full requested list here would let the Explorer offer a skipped
    # contrast and then crash on b\$all_results[[label]] finding nothing.
    # Filtering all three vectors together keeps their index alignment intact.
    completed_mask <- CONTRAST_LABELS %in% names(all_results)
    bundle_cases   <- CONTRAST_CASES[completed_mask]
    bundle_ctrls   <- CONTRAST_CTRLS[completed_mask]
    bundle_labels  <- CONTRAST_LABELS[completed_mask]

    explorer_bundle <- list(
        schema_version   = 3L,
        generated        = as.character(Sys.time()),
        run_id           = "${RUN_ID}",
        script_version   = "${SCRIPT_VERSION}",
        diff_out         = DIFF_OUT,
        genome           = GENOME,
        # User-defined profiles retain the source package labels while loading
        # immutable frozen database objects. Availability of enrichment menus
        # is controlled separately by has_go/has_kegg result flags.
        txdb             = TXDB_PKG_LABEL,
        orgdb            = ORGDB_LABEL,
        has_go           = has_go_results,
        has_kegg         = has_kegg_results,
        txdb_gene_keytype = TXDB_GENE_KEYTYPE,
        annotation_sources = annotation_source,
        contrast_cases   = bundle_cases,
        contrast_ctrls   = bundle_ctrls,
        contrast_labels  = bundle_labels,
        fdr_cutoff       = FDR_CUTOFF,
        fc_cutoff        = FC_CUTOFF,
        col_data         = explorer_col_data,
        vst_matrix       = vst_for_bundle,
        all_results      = all_results,
        homer_find_motifs = HOMER_FIND_MOTIFS,
        homer_configure  = HOMER_CONFIGURE,
        # Sample BAM paths and group assignments — used by the explorer
        # tornado menu (item 7) to run deepTools without re-running DiffBind.
        sample_ids       = SAMPLE_IDS,
        sample_bams      = SAMPLE_BAMS,
        sample_groups    = SAMPLE_GROUPS
    )

    bundle_path <- file.path(DIFF_OUT, "explorer_bundle.rds")
    saveRDS(explorer_bundle, bundle_path)

    vst_note <- if (!is.null(vst_for_bundle)) {
        sprintf("%d peaks × %d samples", nrow(vst_for_bundle), ncol(vst_for_bundle))
    } else {
        "NULL (PCA was skipped — explorer PCA/correlation will be unavailable)"
    }
    cat(sprintf("  Explorer bundle saved: %s\n", basename(bundle_path)))
    cat(sprintf("  VST matrix: %s\n", vst_note))
    cat(sprintf("  Contrasts in bundle: %d\n", length(all_results)))
    cat(sprintf("  GO results available: %s\n", has_go_results))
    cat(sprintf("  KEGG results available: %s\n", has_kegg_results))
    cat(sprintf("  Tip: run PEPATAC_explore.sh and point it at:\n  %s\n", DIFF_OUT))

}, error = function(e) {
    cat("  [WARN] Explorer bundle could not be saved:", conditionMessage(e), "\n")
    cat("  The differential analysis results are still complete.\n")
    cat("  Re-run the analysis to generate the explorer bundle.\n")
})

# ── Session info for reproducibility ─────────────────────────
tryCatch({
    si_path <- file.path(DIFF_OUT, "R_session_info.txt")
    sink(si_path)
    cat("PEPATAC Differential Analysis — R Session Info\n")
    cat(sprintf("Run ID: ${RUN_ID}\n"))
    cat(sprintf("Script version: ${SCRIPT_VERSION}\n\n"))
    print(sessionInfo())
    sink()
    cat(sprintf("  Session info saved: %s\n", basename(si_path)))
}, error = function(e) {
    if (sink.number() > 0) sink()
    cat("  [WARN] Could not save session info:", conditionMessage(e), "\n")
})
RSCRIPT_EOF

ok "R script written: $DIFF_R_SCRIPT"

# ─────────────────────────────────────────────────────────────
# STEP 13 — Run the analysis
# ─────────────────────────────────────────────────────────────

header "Step 13 · Running Analysis"

echo -e "  ${DIM}Log: $DIFF_LOG${RESET}"
blank

set +e
in_env_clean Rscript --vanilla "$DIFF_R_SCRIPT" 2>&1 | tee "$DIFF_LOG"
R_EXIT="${PIPESTATUS[0]}"
set -e

blank
if [[ "$R_EXIT" -eq 0 ]]; then
    ok "Analysis completed successfully."
else
    err "R script exited with code $R_EXIT. Check log: $DIFF_LOG"
    exit "$R_EXIT"
fi

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────

header "Run Complete"

echo -e "  ${BOLD}Results written to:${RESET}   $DIFF_OUT"
echo -e "  ${BOLD}Log:${RESET}                  $DIFF_LOG"
echo -e "  ${BOLD}R script:${RESET}             $DIFF_R_SCRIPT"
echo -e "  ${BOLD}Included samples:${RESET}      $INCLUDED_SAMPLES_FILE"
echo -e "  ${BOLD}Excluded samples:${RESET}      $EXCLUDED_SAMPLES_FILE"
echo -e "  ${BOLD}DiffBind sheet:${RESET}       $DIFFBIND_SHEET"
echo -e "  ${BOLD}DiffBind sample count:${RESET} ${#SAMPLE_IDS[@]} included; ${#EXCLUDED_SAMPLE_IDS[@]} excluded before loading"
blank
echo -e "  ${BOLD}Per-contrast folders:${RESET}"
for label in "${CONTRAST_LABELS[@]}"; do
    echo -e "    ${CYAN}•${RESET}  $DIFF_OUT/$label/"
    echo -e "        ${DIM}├── ${label}_significant_peaks.csv${RESET}"
    echo -e "        ${DIM}├── ${label}_up_peaks.bed / _down_peaks.bed  ← directional BEDs${RESET}"
    echo -e "        ${DIM}├── ${label}_genes_for_enrichr.txt (+ _up / _down)  ← paste into enrichr.com${RESET}"
    echo -e "        ${DIM}├── ${label}_annotated_peaks.tsv${RESET}"
    echo -e "        ${DIM}├── GO/  ← clusterProfiler GO + KEGG, combined + Up-only + Down-only,${RESET}"
    echo -e "        ${DIM}│         full tested set (*.csv) + significant subset (*_significant.csv)${RESET}"
    echo -e "        ${DIM}├── annotation_summary/genomic_features_counts.csv${RESET}"
    if $RUN_MOTIF_ANALYSIS; then
        case "$MOTIF_MODE" in
            known)  echo -e "        ${DIM}└── motifs/, motifs_up/, motifs_down/  → knownResults.html  ← known TF motifs, combined + directional${RESET}" ;;
            denovo) echo -e "        ${DIM}└── motifs/, motifs_up/, motifs_down/  → homerResults.html  ← de novo motifs, combined + directional${RESET}" ;;
            both)   echo -e "        ${DIM}└── motifs/, motifs_up/, motifs_down/  → knownResults.html + homerResults.html${RESET}" ;;
        esac
    fi
done
blank
echo -e "  ${BOLD}Cross-contrast summary:${RESET}  $DIFF_OUT/contrast_summary.csv"
echo -e "  ${BOLD}PCA (all samples):${RESET}       $DIFF_OUT/pca/"
echo -e "  ${BOLD}Diagnostics + checkpoints:${RESET} $DIFF_OUT/diagnostics/"
echo -e "    ${DIM}• diffbind_counted.rds       — post-count checkpoint (reload to skip BAM recounting)${RESET}"
echo -e "    ${DIM}• consensus_peak_raw_counts.rds / .tsv.gz${RESET}"
echo -e "    ${DIM}• diffbind_normalization_factors_whole_experiment.tsv  — informational only,${RESET}"
echo -e "    ${DIM}  NOT what any individual contrast used (see below)${RESET}"
echo -e "    ${DIM}• annotation_provenance/annotation_sources.tsv${RESET}"
echo -e "    ${DIM}• annotation_provenance/txdb_orgdb_keytype_compatibility.tsv${RESET}"
echo -e "  ${BOLD}Per-contrast normalization:${RESET} <contrast>/diagnostics/diffbind_normalization_factors.tsv"
echo -e "    ${DIM}(the real factors that contrast's own DESeq2 model used — each contrast is${RESET}"
echo -e "    ${DIM}normalized independently, so these can differ for a sample shared across contrasts)${RESET}"
echo -e "  ${BOLD}Session info:${RESET}            $DIFF_OUT/R_session_info.txt"
echo -e "  ${BOLD}Explorer bundle:${RESET}         $DIFF_OUT/explorer_bundle.rds"
blank
echo -e "  ${DIM}This script only computed results — no plots were generated.${RESET}"
echo -e "  ${DIM}Next step: run PEPATAC_explore.sh and point it at this output folder.${RESET}"
echo -e "  ${DIM}The explorer gives you interactive PCA, re-plots at custom thresholds,${RESET}"
echo -e "  ${DIM}peak boxplots, correlation heatmaps, motif summaries, GO/KEGG views,${RESET}"
echo -e "  ${DIM}annotation distribution plots (feature pie, TSS/TES, gene body, width,${RESET}"
echo -e "  ${DIM}chromosome), and tornado plots — any settings, regenerated as many${RESET}"
echo -e "  ${DIM}times as you like without re-running DiffBind.${RESET}"
blank

# Shown once already (as a confirmation gate before the DiffBind sheet was
# written) and repeated here on purpose -- someone skimming straight to
# the end of a long run's output should not be able to miss that some of
# these results rest on a peak-calling fallback for one or more samples.
if [[ ${#BED_FALLBACK_SAMPLES[@]} -gt 0 ]]; then
    print_peak_format_alert
fi

echo -e "  ${GREEN}${BOLD}Differential analysis complete.${RESET}"
blank
