#!/bin/bash

MAIN_DIR=$(realpath -m "$(dirname "$0")/../")
OPT=$MAIN_DIR/flinter_rc.yml

# Check if the flint command is available
if ! command -v flint &>/dev/null; then
    echo "Flint could not be found. Please install it first."
    exit 1
fi

# List the modified files from develop branch
FILES=$(git diff --name-only develop -- $MAIN_DIR | grep -E '\.f90$')
if [ -z "$FILES" ]; then
    echo "No modified files found."
    exit 0
fi

# Get the score printed between the arrows
for file in $FILES; do
    file=$(realpath -m "$file")

    score=$(flint score -r $OPT $file 2>/dev/null |
        grep -oP '(?<=\>\|)[^\|\<]+(?=\|\<)')

    if [ -z "$score" ]; then
        continue
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
