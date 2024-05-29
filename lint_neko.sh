#!/bin/bash

DIR=external/neko/src
OPT=external/neko/flinter_rc.yml

# Get the score printed between the arrows
for file in $(find $DIR -name "*.f90"); do

    score=$(flint score -r $OPT $file 2>/dev/null |
        grep -oP '(?<=\>\|)[^\|\<]+(?=\|\<)')

    if [ -z "$score" ]; then
        echo "Error: $file"
        cat /tmp/flint_error
        continue
    elif [ -f /tmp/flint_error ]; then
        rm /tmp/flint_error
    fi
    if [ "$score" != "10.00" ]; then
        report=$(flint lint -r $OPT $file 2>/dev/null)
        if [ -z "$report" ]; then
            flint stats -r $OPT $file 2>/dev/null
        else
            echo "$report"
        fi
    fi
    echo "$score: $file"

    while true; do
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
        [uU])
            report=$(flint lint -r $OPT $file 2>/dev/null)
            if [ -z "$report" ]; then
                flint stats -r $OPT $file 2>/dev/null
            else
                echo "$report"
            fi
            score=$(flint score -r $OPT $file 2>/dev/null |
                grep -oP '(?<=\>\|)[^\|\<]+(?=\|\<)')
            echo "$score: $file"
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
