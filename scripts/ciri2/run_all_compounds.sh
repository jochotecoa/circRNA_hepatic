tissue="Hepatic"
dose=("The" "Tox")
timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
timeTOX=("002" "008" "024" "072" "168" "240")
replicates=("1" "2" "3")


# compound="AZA" #
# compid="Azathioprine"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="CYC" #
compid="Cyclosporin"
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="DIC" #
compid="Diclofenac"
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="ISO" #
compid="Isoniazid"
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="MTX" #
compid="Methotrexate"
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="PHE" #
compid="Phenytoin"
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="RIF" #
compid="Rifampicin"
. ./scripts/ciri2/ciri2v6.1_multiple.sh


compound="VPA" #
compid="Valproic_Acid"
timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
timeTOX=("000" "002" "008" "024" "072" "168" "240")
. ./scripts/ciri2/ciri2v6.1_multiple.sh

compound="ConDMSO" #
compid="Con_0.1_DMSO"
timeCTRL=("002" "008" "024" "072" "168" "240" "336")
. ./scripts/ciri2/ciri2v6.1_multiple.sh
# 
# compound="APA" #
# compid="Acetaminophen"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh

