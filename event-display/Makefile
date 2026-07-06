ROOT_CONFIG ?= root-config
CXX ?= g++
AR ?= ar

BIN_DIR := bin
BUILD_DIR := build
LIB_DIR := lib
SRC_DIR := src
INC_DIR := include
APP_DIR := app
WEB_DIR := web

CXXFLAGS += -g -O2 -Wall -Wextra -Wno-sign-compare -std=c++17
CPPFLAGS += -I. -I$(INC_DIR) -I$(SRC_DIR) $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs) -lGeom -lEG
LDFLAGS += $(shell $(ROOT_CONFIG) --ldflags)

CORE_SRCS := $(wildcard $(SRC_DIR)/*.cc)
CORE_OBJS := $(CORE_SRCS:$(SRC_DIR)/%.cc=$(BUILD_DIR)/src/%.o)
CORE_DEPS := $(CORE_OBJS:.o=.d)

APP_SRCS := $(wildcard $(APP_DIR)/*.cc)
APP_OBJS := $(APP_SRCS:$(APP_DIR)/%.cc=$(BUILD_DIR)/app/%.o)
APP_DEPS := $(APP_OBJS:.o=.d)

LIB_ANALYZER := $(LIB_DIR)/libAnalyzer.a
EXPORT_BIN := $(BIN_DIR)/tpc_event_export

.PHONY: all lib clean distclean convert geometry serve

all: $(EXPORT_BIN)

lib: $(LIB_ANALYZER)

$(EXPORT_BIN): $(APP_OBJS) $(LIB_ANALYZER)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(APP_OBJS) -L$(LIB_DIR) -lAnalyzer $(ROOT_LIBS)

$(LIB_ANALYZER): $(CORE_OBJS)
	@mkdir -p $(LIB_DIR)
	$(AR) rcs $@ $^

$(BUILD_DIR)/src/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/app/%.o: $(APP_DIR)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -MP -c $< -o $@

convert: $(EXPORT_BIN)
	$(EXPORT_BIN) -i $(INPUT) -o $(or $(OUTPUT),$(WEB_DIR)/display.root) -e $(or $(EVENT),0) -n $(or $(COUNT),1) $(if $(TREE),-t $(TREE),)

geometry: $(EXPORT_BIN)
	$(EXPORT_BIN) -o $(WEB_DIR)/display.root

serve:
	python3 app/tpc_event_server.py --host 127.0.0.1 --port 8080

clean:
	rm -rf $(BUILD_DIR) $(LIB_ANALYZER)

distclean: clean
	rm -rf $(BIN_DIR) $(WEB_DIR)/display.root

-include $(CORE_DEPS) $(APP_DEPS)
