#!/bin/sh

DB=uniref100
method=crnpred

Usage () {
echo "Usage: $0 [-m method] [-d DB] input_seq"
echo "  input_seq: protein sequence in the FASTA format."
echo "  method: 'crnpred2k', 'crnpred5k' or 'linear' (default: crnpred2k)"
echo "  DB: BLAST sequence database (default: uniref100)"
exit 2
}

#Usage
args=`getopt m:d:h $*`
set -- $args

for i in $*; do
  case "$i" in
    -m) method=$2; shift 2;;
    -d) DB=$2; shift 2;;
    -h) Usage ;;
    --) shift; break;;
  esac
done
inseq=${1}

if [ "hoge${inseq}" = "hoge" ]; then
  Usage
fi

prof=`basename ${inseq}`.d.prof
fout=`basename ${inseq}`.d.out

blastpgp -d $DB -h 0.0005 -j 3 -i $inseq -Q $prof > /dev/null 2>&1 

echo "# Input sequence" > $fout
cat $inseq >> $fout
echo '' >> $fout
echo '' >> $fout

echo "*** SS, CN, RWCO predictions..."
if [ "$method" = "crnpred2k" ]; then
  echo "   prediction by CRN2000 (takes some time...)"
  ${CRNPRED_DIR}/bin/xpredm2000 $prof >> $fout
elif [ "$method" = "crnpred5k" ]; then
  echo "   prediction by CRN5000 (takes some time...)"
  ${CRNPRED_DIR}/bin/xpredm5000 $prof >> $fout
else
  echo "   prediction by linear method"
  ${CRNPRED_DIR}/bin/lpredm $prof >> $fout
fi
