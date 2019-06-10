TEST_BIN_DIR = bin

TEST_OUT_DIR = test_output

$(TEST_BIN_DIR)/test_rational: $(TEST_BIN_DIR) test_rational.cpp include/ext/rational.h include/ext/rational.tcc
	g++ -std=gnu++17 -I./include -I../include -I../polynomial/include -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_rational test_rational.cpp -lquadmath

clean:
	rm -rf $(TEST_BIN_DIR)/*

test: $(TEST_OUT_DIR) $(TEST_BIN_DIR)/test_rational
	./$(TEST_BIN_DIR)/test_rational > $(TEST_OUT_DIR)/test_rational.txt
	fi

$(TEST_BIN_DIR):
	if test ! -d $(TEST_BIN_DIR); then \
	  mkdir $(TEST_BIN_DIR); \
	fi

$(TEST_OUT_DIR):
	if test ! -d $(TEST_OUT_DIR); then \
	  mkdir $(TEST_OUT_DIR); \
	fi
