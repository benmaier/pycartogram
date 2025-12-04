#!/bin/bash
# Fix duplicate LC_RPATH in cCartogram

CCARTOGRAM_PATH="/Users/bfmaier/miniconda3/lib/python3.12/site-packages/cCartogram.cpython-312-darwin.so"

echo "Current rpaths:"
otool -l "$CCARTOGRAM_PATH" | grep -A2 LC_RPATH

echo ""
echo "Removing duplicate rpath..."
install_name_tool -delete_rpath '/Users/bfmaier/miniconda3/lib' "$CCARTOGRAM_PATH"

echo ""
echo "Rpaths after fix:"
otool -l "$CCARTOGRAM_PATH" | grep -A2 LC_RPATH
