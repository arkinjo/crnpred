#!/usr/bin/env zsh

ndim=${1:-"2000"}
nd=${2:-"20"}
\rm -f wm.tmp wi.tmp wss.tmp wcn.tmp wrwco.tmp

ls Mat/wmat${ndim}_${nd}_1.0* | sed 's:^:w2000/:' > wm.tmp
ls Mat/win${ndim}_21.0.01* | sed 's:^:w2000/:' > wi.tmp
for typ in ss cn rwco ; do
    ls w_${typ}/xwoutall* | sed 's:^:w2000/:' > w${typ}.tmp
done

paste wm.tmp wi.tmp wss.tmp wcn.tmp wrwco.tmp

\rm -f wm.tmp wi.tmp wss.tmp wcn.tmp wrwco.tmp
