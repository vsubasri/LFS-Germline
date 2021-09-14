#!/usr/bin/env bash
for f in final_results/*.telseq; do
    tail -n +3 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
