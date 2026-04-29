set -euo pipefail

# --- Configuration Loading ---
# This script expects 'config.yaml' to be in the project root directory.
CONFIG_FILE="../../config.yaml"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found." >&2
    echo "Please copy 'config.yaml.template' to 'config.yaml' and configure it." >&2
    exit 1
fi

# Load parameters from config.yaml using yq
PROJECT_BASE_DIR=$(yq '.PROJECT_BASE_DIR' "$CONFIG_FILE")
BASE_DATA_STORAGE_DIR=$(yq '.BASE_DATA_STORAGE_DIR' "$CONFIG_FILE")
TRIMMED_OUTPUT_SUBDIR=$(yq '.DIR_SETTINGS.TRIMMED_OUTPUT_SUBDIR' "$CONFIG_FILE")

# Check if yq command was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to parse config.yaml. Make sure yq is installed and config.yaml is valid." >&2
    exit 1
fi

tissue="Hepatic"

replicates=("1" "2" "3")

# compound="AZA" #
# compid="Azathioprine"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
#

# dose=("The")
# timeTHE=("000")
# replicates=("3")
# 
# compound="DIC" #
# compid="Diclofenac"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# dose=("The" "Tox")
# timeTHE=("002" "008" "024" "072" "168" "240" "336")
# timeTOX=("002" "008" "024" "072" "168" "240")
# replicates=("1" "2" "3")
# 
# compound="DIC" #
# compid="Diclofenac"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# 
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
# 
# 
# compound="ISO" #
# compid="Isoniazid"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# compound="MTX" #
# compid="Methotrexate"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 

# dose=("Tox")
# timeTOX=("000")
# replicates=("1" "2" "3")
# 
# compound="CYC" #
# compid="Cyclosporin"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# compound="ISO" #
# compid="Isoniazid"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# dose=("The" "Tox")
# timeTHE=("008" "024" "072" "168" "240" "336")
# timeTOX=("002" "008" "024" "072" "168" "240")
# 
# compound="PHE" #
# compid="Phenytoin"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"

# dose=("Tox")
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
# timeTOX=("168" "240")
# 
# compound="RIF" #
# compid="Rifampicin"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"

# dose=("The")
# compound="VPA" #
# compid="Valproic_Acid"
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
# timeTOX=("000" "002" "008" "024" "072" "168" "240")
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd "$PROJECT_BASE_DIR"
# 
# 
timeCTRL=("002" "008" "024" "072" "168" "240" "336")
trimmed_dir="$TRIMMED_OUTPUT_SUBDIR"

compound="ConDMSO" #
compid="Con_0.1_DMSO"
. ./scripts/ciri2/ciri2v6_multiple_ctrl_wout_bwa.sh
cd "$PROJECT_BASE_DIR"

# dose=("Tox")
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
# timeTOX=("002" "008" "024" "072" "168" "240")
# 
# 
# compound="APA" #
# compid="Acetaminophen"
# . ./scripts/ciri2/ciri2v6.1_APA.sh

