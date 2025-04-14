#!/bin/bash

DIR=external/neko/src
OPT=external/neko/flinter_rc.yml

# Get the score printed between the arrows
for file in $(find $DIR -name "*.f90"); do

    score=$(flint score -r $OPT $file 2>/dev/null |
        grep -oP '(?<=\>\|)[^\|\<]+(?=\|\<)')
    report=$(flint lint -r $OPT $file 2>/dev/null)
    if [ -z "$report" ]; then
        continue
    else
        echo "$report"
    fi

    while true; do
        echo "$score: $file"
        if [ "$score" == "10.00" ]; then
            break
        else
            read -p "How to proceede: " proceede
        fi
        case $proceede in
        [nN])
            echo "Moving to the next file"
            break
            ;;
        [pP])
            flint lint -r $OPT $file 2>/dev/null
            ;;
        [uU])
            score=$(flint score -r $OPT $file 2>/dev/null |
                grep -oP '(?<=\>\|)[^\|\<]+(?=\|\<)')
            ;;
        [qQ])
            echo "Quitting the script"
            exit 0
            ;;
        *)
            echo "Please specify if we should [u]pdate the score,"
            echo "[p]rint the detailed flint report, go to [n]ext file"
            echo "or [q]uit the script"
            ;;
        esac
    done

done
