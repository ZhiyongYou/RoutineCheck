inc_dir = include/
src_dir = src/
obj_dir = obj/
CCOMPILER  = g++
DEBUG = -O2 #-ggdb  #-g -gstabs -gstabs+ -ggdb -glevel

object = $(obj_dir)astro.o \
		$(obj_dir)WFCTASLC.o \
		$(obj_dir)LedEvtSel.o \
		$(obj_dir)WFCTARec.o \
		$(obj_dir)WFCTAMap.o \
		$(obj_dir)CollectIntensity.o \
		$(obj_dir)EventDict.o \
		$(obj_dir)WFCTAEvent.o

all:main draw

draw: $(obj_dir)draw.a
	$(CCOMPILER) -o draw $(obj_dir)draw.a `root-config --cflags --libs` -lMinuit -lMatrix

$(obj_dir)draw.a: $(obj_dir)draw.o $(object)
	ar -r $(obj_dir)draw.a $(obj_dir)draw.o $(object)

$(obj_dir)draw.o: draw.cc Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

main: $(obj_dir)main.a
	$(CCOMPILER) -o main $(obj_dir)main.a `root-config --cflags --libs` -lMinuit -lMatrix

$(obj_dir)main.a: $(obj_dir)main.o $(object)
	ar -r $(obj_dir)main.a $(obj_dir)main.o $(object)

$(obj_dir)main.o: main.cc Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WFCTASLC.o: $(src_dir)WFCTASLC.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)LedEvtSel.o: $(src_dir)LedEvtSel.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WFCTARec.o: $(src_dir)WFCTARec.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WFCTAMap.o: $(src_dir)WFCTAMap.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)astro.o: $(src_dir)astro.c Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)CollectIntensity.o: $(src_dir)CollectIntensity.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)EventDict.o: $(src_dir)EventDict.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WFCTAEvent.o: $(src_dir)WFCTAEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

.PHONY : clean
 clean :
	rm main draw EventDict_rdict.pcm obj/*
