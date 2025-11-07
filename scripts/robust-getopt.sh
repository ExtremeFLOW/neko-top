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

    local GETOPT=""
    local GNU_GETOPT=""

    # Check if the OS is Linux or Darwin (macOS)
    if [ "$OS_TYPE" == "Linux" ]; then
        if command -v getopt >/dev/null 2>&1; then
            GETOPT=$(command -v getopt)
        else
            echo "Error: getopt is not installed. Please install it and try again." >&2
            exit 1
        fi
    elif [ "$OS_TYPE" == "Darwin" ]; then

        if command -v gnugetopt >/dev/null 2>&1; then
            GETOPT=$(command -v gnugetopt)
        elif command -v gnu-getopt >/dev/null 2>&1; then
            GETOPT=$(command -v gnu-getopt)
        elif command -v getopt >/dev/null 2>&1; then
            GETOPT=$(command -v getopt)
        else
            echo "Error: Neither GNU getopt nor BSD getopt is available. Please ensure one is installed and try again." >&2
            exit 1
        fi
    else
        echo "Error: Unsupported operating system: $OS_TYPE" >&2
    fi

    $GETOPT -T >/dev/null
    if [ $? -eq 4 ]; then
        PARSED=$($GETOPT --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
    else
        PARSED=$($GETOPT $OPT "$@")
    fi

    if [ $? -ne 0 ]; then
        echo "Error: Failed to parse options." >&2
        exit 1
    fi

    # Return the parsed options
    echo "$PARSED"
    return $?
}
