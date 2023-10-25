# Simple setup to compile the current example.
set(DRIVER_TYPE "user")
set(EXTRA_SOURCES
    ${NEKO_DIR}/examples/turb_channel/turb_channel.f90
)

build_example(${DRIVER_TYPE})
