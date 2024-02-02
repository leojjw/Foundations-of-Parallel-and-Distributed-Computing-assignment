compiler  		=g++

NVCC 			=nvcc
NVCC_FLAGS 		=-std=c++11 -O3 -Wno-deprecated-gpu-targets

LIBRARIES 		=-L/${CUDA_DIR}/lib64 -lcudart -lm

ifeq ($(mode),debug)
	cflags    	+=-std=c++17 -g3 -ggdb  -fsanitize=address -fno-omit-frame-pointer -march=native  -fopenmp
	linkflags 	+= -fsanitize=address -fopenmp -L/${CUDA_DIR}/lib64 -lcudart -lm
else
	cflags    	+=-std=c++17 -O3 -march=native  -fopenmp -fopt-info-vec -fopt-info-inline -fopenmp
	linkflags 	+=-fopenmp -L/${CUDA_DIR}/lib64 -lcudart -lm
endif

SrcDir	  		=./src
ObjDir    		=./obj
HeadRoot		=./include
HeadDir         +=$(HeadRoot) $(foreach dir,$(HeadRoot),$(wildcard $(dir)/*))
source    		=$(foreach dir,$(SrcDir),$(wildcard $(dir)/*.cpp))
head      		=$(foreach dir,$(HeadDir),$(wildcard $(dir)/*.hpp))
object    		=$(patsubst %.cpp,$(ObjDir)/%.o,$(notdir $(source)))

target 	  		=HeatDifu
NO_COLOR		=\033[0m
OK_COLOR		=\033[32;01m


$(target):$(object) ./obj/gpu_cal.o $(head)
	$(compiler) -o $(target) $(object) ./obj/gpu_cal.o $(linkflags) $(lib)
	@printf "$(OK_COLOR)Compiling Is Successful!\nExecutable File: $(target) $(NO_COLOR)\n"

$(ObjDir)/%.o:$(SrcDir)/%.cpp $(head)
	$(compiler) -c $(cflags) $< -o $@ -I $(HeadRoot)

./obj/gpu_cal.o:$(HeadRoot)/Vector/gpu_cal.cu
	$(NVCC) $(NVCC_FLAGS) -c $^ -o $@

.PHONY:run 
run:$(target)
	@printf "$(OK_COLOR)$(target) is executing $(NO_COLOR)\n"
	./$(target)

.PHONY:clean	 
clean:
	-rm $(object) ./obj/gpu_cal.o $(target)

.PHONY:clean_all	 
clean_all:
	-rm $(object) ./obj/gpu_cal.o $(target) ./output/*

.PHONY:plot
plot:
	python3 Plot.py
