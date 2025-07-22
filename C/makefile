# Variables
CC = icx
CFLAGS = -Wno-unused-command-line-argument -lm -llapacke -lblas -qopenmp -O3 -march=native -fp-model precise
SRC_DIR = src
BUILD_DIR = build
APP_DIR = app
MY_PROGRAM = HNCMSA.out

# Source and object files
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES))

# Targets
all: $(MY_PROGRAM)

$(MY_PROGRAM): $(OBJECTS) $(BUILD_DIR)/main.o
	$(CC) $^ $(CFLAGS) -o $@

$(BUILD_DIR)/main.o: $(APP_DIR)/main.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(MY_PROGRAM)

#.PHONY: all clean