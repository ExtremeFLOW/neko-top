# ============================================================================ #
# Robust argument parsing
# This function is used to parse the command line arguments and set the
# appropriate variables. It uses the `getopt` command to handle both short
# and long options.

function robust-getopt {
    # Define the options
    local NAME=$1    # The name of the script
    local OPT=$2     # Short options, e.g., "h,t,c,q,d:"
    local OPTIONS=$3 # Long options, e.g., "help,test,procs:"
    local PARSED     # The parsed options to be returned
    shift 3

    local GETOPT=""
    local OS_TYPE=$(uname)

    # Check if the OS is Linux or Darwin (macOS)
    if [ "$OS_TYPE" == "Linux" ]; then
        if command -v getopt >/dev/null 2>&1; then
            GETOPT=$(command -v getopt)
        else
            echo "Error: getopt is not installed." >&2
            echo "Please install it and try again." >&2
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
            echo "Error: Neither GNU getopt nor BSD getopt is available." >&2
            echo "Please ensure one is installed and try again." >&2
            exit 1
        fi
    else
        echo "Error: Unsupported operating system: $OS_TYPE" >&2
    fi

    $GETOPT -T >/dev/null
    if [ $? -eq 4 ]; then
        PARSED=$($GETOPT --options=$OPT --longoptions=$OPTIONS --name $NAME -- "$@")
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
