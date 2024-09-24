#!/bin/bash

Dir=$(dirname "${BASH_SOURCE[0]}")

# Generate new table
cat ${Dir}/G_T.sorted.sam |\
python ${Dir}/intersections_test.py -b ${Dir}/bedfile.bedGraph |\
python ${Dir}/tasmanian_test.py -r ${Dir}/small_region.fa 2>/dev/null > ${Dir}/new_test.table.csv

# Compare new table with existing one
diff_output=$(diff ${Dir}/test_test.table.csv ${Dir}/new_test.table.csv)
diff_lines=$(echo -n "$diff_output" | wc -l)

if [ $diff_lines -eq 0 ]; then
    echo "passed"
else
    echo "failed"
    exit 1
fi
