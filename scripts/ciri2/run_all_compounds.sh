tissue="Hepatic"

# compound="AZA" #
# compid="Azathioprine"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd /share/script/hecatos/juantxo/circRNA_hepatic
#
# compound="CYC" #
# compid="Cyclosporin"
# . ./scripts/ciri2/ciri2v6.1_multiple.sh
# cd /share/script/hecatos/juantxo/circRNA_hepatic


dose=("The")
timeTHE=("000")
replicates=("3")

compound="DIC" #
compid="Diclofenac"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

dose=("The" "Tox")
timeTHE=("002" "008" "024" "072" "168" "240" "336")
timeTOX=("002" "008" "024" "072" "168" "240")
replicates=("1" "2" "3")

compound="DIC" #
compid="Diclofenac"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic


timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")


compound="ISO" #
compid="Isoniazid"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="MTX" #
compid="Methotrexate"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="PHE" #
compid="Phenytoin"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="RIF" #
compid="Rifampicin"
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="VPA" #
compid="Valproic_Acid"
timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
timeTOX=("000" "002" "008" "024" "072" "168" "240")
. ./scripts/ciri2/ciri2v6.1_multiple.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="ConDMSO" #
compid="Con_0.1_DMSO"
timeCTRL=("002" "008" "024" "072" "168" "240" "336")
. ./scripts/ciri2/ciri2v6.1_multiple_ctrl.sh
cd /share/script/hecatos/juantxo/circRNA_hepatic

compound="APA" #
compid="Acetaminophen"
timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
timeTOX=("002" "008" "024" "072" "168" "240")
. ./scripts/ciri2/ciri2v6.1_APA.sh

