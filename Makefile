
exe: PSO.cpp Entete.h
	g++ $^ -o $@

run:
	./exe 10 0.1 0.1 0.1 100

run2: exe
	./exe 100 0.8 1.62 1.62 10000 CutTest.txt
