################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= ./SAfish
# Cuda source files (compiled with cudacc)
CUFILES		:= SAfish.cu

################################################################################
# Rules and targets

# location of cuda sdk
CUDA_SDK	= /opt/NVIDIA_GPU_Computing_SDK

# override variables that common.mk needs in order to find files in the SDK
ROOTDIR		= ${CUDA_SDK}/C/common
BINDIR		= ./bin

include ${CUDA_SDK}/C/common/common.mk
