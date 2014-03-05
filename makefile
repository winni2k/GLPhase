all:
	$(MAKE) -C src clean
	$(MAKE) -C src

test:
	$(MAKE) -C src test

oxford:
	$(MAKE) -C src oxford

.PHONY: all test oxford
