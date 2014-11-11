################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= ./fish_gpu

# Cuda source files (compiled with cudacc)
CUFILES		:= fish.cu

################################################################################
# Rules and targets

# location of cuda sdk
CUDA_SDK	= /opt/NVIDIA_GPU_Computing_SDK

# override variables that common.mk needs in order to find files in the SDK
ROOTDIR		= ${CUDA_SDK}/C/common
BINDIR		= ./bin

include ${CUDA_SDK}/C/common/common.mk
