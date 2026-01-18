#!/usr/bin/env make

SHELL := bash
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.SECONDARY:
.SUFFIXES:
.SECONDEXPANSION:
#.ONESHELL:
.DELETE_ON_ERROR:
.DEFAULT: all
.PHONY: all clean dirs

SHELL := bash
GTDB := /home/shared/db/gtdb/r226/gtdb_genomes_r226
TSV := $(GTDB)/GTDB.tsv

MIN := 100
MAX := 9999999
KMER := 9
STEP := 5

GROUPS := $(file < groups.txt)
#$(info GROUPS: $(GROUPS))

#GROUP_DIRS := $(GROUPS:%=$(GTDB)/%)
#$(info GROUPP_DIRS: $(GROUP_DIRS))
#vpath %.faa.gz $(GROUP_DIRS)

all : $(GROUPS:%=%/kmers)

#%/kmers : $$(file < $(GTDB)/%/proteomes.txt)
#%/kmers : $$(wildcard $(GTDB)/%/*.faa.gz)

#.SECONDEXPANSION:
#$(DIR)/%/kmers : $(GTDB)/%/proteomes.txt $$(addprefix $(GTDB)/$$*/, $$(file < $(GTDB)/$$*/proteomes.txt))
#	@echo $* = $^
#	@wc -l < $<

#$(DIR)/%/kmers : $(GTDB)/%/proteomes.txt $$(addprefix $(GTDB)/$$*/, $$(file < $(DIR)/$$*/proteomes.txt))
#$(DIR)/%/kmers : $(GTDB)/%/proteomes.txt $$(patsubst %.faa.gz,$(DIR)/$$*/%.kmer, $$(file < $(GTDB)/$$*/proteomes.txt))
.SECONDEXPANSION:
#%/kmers : $(GTDB)/%/proteomes.txt  #$$(file < $(GTDB)/$$*/proteomes.txt)
#%/kmers : $$(patsubst %.faa.gz,./$$*/%.kmer,$$(strip $$(file < $(GTDB)/$$*/proteomes.txt)))
%/kmers : $$(patsubst %.faa.gz,%.kmer,$$(addprefix $$*/,$$(strip $$(file < $(GTDB)/%/proteomes.txt))))
	@echo "TARGET=$@"
	@echo "STEM=$*"
	@echo "DEPS=$^"
	
%.kmer : $(GTDB)/%.faa.gz
	seqkit sliding \
	  -w 0 -W $(KMER) -s $(STEP) $< \
	| tantan -p -x X \
	| seqkit seq -u --seq \
	| grep -v '[*XUOW]' \
	| csvtk uniq -H \
	> $(@)

%/nsamples : $(TSV)
	grep -m 1 -F $(*) $(TSV) \
	| cut -f 4 \
	> $(@)
    
groups.txt : $(TSV)
	cat $< \
	| grep -F 'Bacteria/' \
	| tsvtk filter2 -H -f '$$4 >= $(MIN)' \
	| tsvtk filter2 -H -f '$$4 <= $(MAX)' \
	| cut -f 2 \
	| sed -E 's,/[^/]+$$,,' \
	| sort \
	| uniq \
	> $@
	xargs -a $@ mkdir -p 

clean :
	rm -f groups.txt
	#rm -fr $(DIR)

