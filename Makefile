all:
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean

debug:
	$(MAKE) -C src debug

test:
	$(MAKE) -C src test

disttest:
	$(MAKE) -C src disttest

oxford:
	$(MAKE) -C src oxford

.PHONY: all test oxford
