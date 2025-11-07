# ============================================================================ #
# Robust argument parsing
# This function is used to parse the command line arguments and set the
# appropriate variables. It uses the `getopt` command to handle both short
# and long options.

function robust-getopt {
    # Define the options
    local NAME=$1
    local OPT=$2
    local OPTIONS=$3
    local PARSED
    shift 3

    local OS_TYPE=$(uname)

    # Use GNU getopt if available, fallback to BSD getopt on macOS
    if [ "$OS_TYPE" == "Linux" ]; then
        if command -v getopt >/dev/null 2>&1; then
            PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")

            if [ $? -ne 0 ]; then
                echo "Error: Failed to parse options using GNU getopt." >&2
                exit 1
            fi
        else
            echo "Error: GNU getopt is not installed. Please install it and try again." >&2
            exit 1
        fi
    elif [ "$OS_TYPE" == "Darwin" ]; then
        if command -v gnugetopt >/dev/null 2>&1; then
            PARSED=$(gnugetopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")

            if [ $? -ne 0 ]; then
                echo "Error: Failed to parse options using GNU getopt on macOS." >&2
                exit 1
            fi
        # Check if the user has installed GNU getopt via Homebrew
        elif command -v gnu-getopt >/dev/null 2>&1; then
            PARSED=$(gnu-getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")

            if [ $? -ne 0 ]; then
                echo "Error: Failed to parse options using GNU getopt on macOS." >&2
                exit 1
            fi
        # Fallback to BSD getopt if GNU getopt is not available
        # Note: BSD getopt does not support long options
        elif command -v getopt >/dev/null 2>&1; then
            echo "Warning: macOS uses BSD getopt, long options are not supported." >&2
            PARSED=$(getopt $OPT "$@")

            if [ $? -ne 0 ]; then
                echo "Error: Failed to parse options using BSD getopt." >&2
                exit 1
            fi
        else
            echo "Error: Neither GNU getopt nor BSD getopt is available. Please ensure one is installed and try again." >&2
            exit 1
        fi
    else
        echo "Error: Unsupported operating system: $OS_TYPE" >&2
        exit 1
    fi

    # Check if the parsed options are empty
    if [ -z "$PARSED" ]; then
        echo "Error: No options provided." >&2
        exit 1
    fi

    # Return the parsed options
    echo "$PARSED"
    return $?
}
