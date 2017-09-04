.PHONY: lint init
.DEFAULT_GOAL := lint

init:
	conda install numpy scipy matplotlib flake8 pytest

lint:
	flake8 .
