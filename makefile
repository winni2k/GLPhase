all:
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean

debug:
	$(MAKE) -C src debug

test:
	$(MAKE) -C src test
	t/regression_tests.pl

oxford:
	$(MAKE) -C src oxford

.PHONY: all test oxford
