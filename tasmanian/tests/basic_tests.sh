#!/bin/bash

Dir=$(dirname "${BASH_SOURCE[0]}")
[ $(diff ${Dir}/test_test.table.csv <(cat ${Dir}/G_T.sorted.sam |\
                                    python ${Dir}/intersections_test.py -b ${Dir}/bedfile.bedGraph |\
                                    python ${Dir}/tasmanian_test.py -r ${Dir}/small_region.fa 2>/dev/null) |\
wc -l) -eq 0 ]\
\
	&& echo "super" || exit 1
