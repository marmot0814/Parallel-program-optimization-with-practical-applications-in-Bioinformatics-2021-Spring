TARGET          := main.out

SRC_DIR         := src/
INC_DIR         := include/
OBJ_DIR         := objs/

CC              := g++-10
CCFLAGS         := -I$(INC_DIR) -std=c++17 -Ofast -g3 -I/usr/local/include/libsimdpp-2.1 -fconcepts-ts -DSIMDPP_ARCH_X86_AVX2 -mavx2

SRC_FILES       := $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJ_FILES       := $(patsubst $(SRC_DIR)%.cpp,$(OBJ_DIR)%.o,$(SRC_FILES))

all: $(shell mkdir -p $(OBJ_DIR)) $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CC) $(CCFLAGS) -o $@ $^

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	$(shell mkdir -p $(dir $@))
	$(CC) $(CCFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

